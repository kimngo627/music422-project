import torch
import torch.nn as nn
import torch.nn.functional as F
import torchaudio
import math
import cdpam

import sys


class EncoderBlock(nn.Module):
    def __init__(self, in_channels, out_channels, kernel_size=3, stride=2, padding=1):
        super(EncoderBlock, self).__init__()
        self.conv = nn.Conv1d(in_channels, out_channels, kernel_size, stride, padding)
        self.bn = nn.BatchNorm1d(out_channels)
        self.relu = nn.ReLU()
        
    def forward(self, x):
        x = self.conv(x)
        x = self.bn(x)
        x = self.relu(x)
        return x


class DecoderBlock(nn.Module):
    def __init__(self, in_channels, out_channels, kernel_size=3, stride=2, padding=1, output_padding=1):
        super(DecoderBlock, self).__init__()
        self.conv_transpose = nn.ConvTranspose1d(in_channels, out_channels, kernel_size, stride, padding, output_padding)
        self.bn = nn.BatchNorm1d(out_channels)
        self.relu = nn.ReLU()
        
    def forward(self, x):
        x = self.conv_transpose(x)
        x = self.bn(x)
        x = self.relu(x)
        return x


class VectorQuantizer(nn.Module):
    """
    Vector Quantization layer for explicit bitrate control
    Based on VQ-VAE (van den Oord et al., 2017)
    """
    def __init__(self, embedding_dim, num_embeddings, commitment_cost=0.25):
        super(VectorQuantizer, self).__init__()
        self.embedding_dim = embedding_dim
        self.num_embeddings = num_embeddings
        self.commitment_cost = commitment_cost
        
        # Create codebook
        self.codebook = nn.Embedding(num_embeddings, embedding_dim)
        self.codebook.weight.data.uniform_(-1/num_embeddings, 1/num_embeddings)
        
    def forward(self, inputs):
        # Flatten inputs
        flat_inputs = inputs.permute(0, 2, 1).contiguous().view(-1, self.embedding_dim)
        
        # Calculate distances
        distances = (torch.sum(flat_inputs**2, dim=1, keepdim=True) 
                    + torch.sum(self.codebook.weight**2, dim=1)
                    - 2 * torch.matmul(flat_inputs, self.codebook.weight.t()))
        
        # Quantization
        encoding_indices = torch.argmin(distances, dim=1).unsqueeze(1)
        encodings = torch.zeros(encoding_indices.shape[0], self.num_embeddings, device=inputs.device)
        encodings.scatter_(1, encoding_indices, 1)
        
        # Quantized vectors
        quantized = torch.matmul(encodings, self.codebook.weight).view(inputs.shape[0], inputs.shape[2], self.embedding_dim)
        
        # Loss
        q_latent_loss = F.mse_loss(quantized.permute(0, 2, 1), inputs)
        commitment_loss = F.mse_loss(quantized.permute(0, 2, 1).detach(), inputs)
        loss = q_latent_loss + self.commitment_cost * commitment_loss
        
        # Straight-through estimator
        quantized = inputs + (quantized.permute(0, 2, 1) - inputs).detach()
        
        return quantized, loss, encoding_indices
    
    def quantize(self, inputs):
        """Used for inference only"""
        flat_inputs = inputs.permute(0, 2, 1).contiguous().view(-1, self.embedding_dim)
        distances = (torch.sum(flat_inputs**2, dim=1, keepdim=True) 
                    + torch.sum(self.codebook.weight**2, dim=1)
                    - 2 * torch.matmul(flat_inputs, self.codebook.weight.t()))
        encoding_indices = torch.argmin(distances, dim=1)
        return encoding_indices
    
    def decode(self, indices):
        """Convert indices back to vectors"""
        quantized = self.codebook(indices)
        return quantized


class TemporalEncoder(nn.Module):
    """RNN for modeling temporal dependencies"""
    def __init__(self, input_dim, hidden_dim, num_layers=1, bidirectional=True):
        super(TemporalEncoder, self).__init__()
        self.gru = nn.GRU(
            input_size=input_dim,
            hidden_size=hidden_dim,
            num_layers=num_layers,
            batch_first=True,
            bidirectional=bidirectional
        )
        
        self.hidden_dim = hidden_dim
        self.bidirectional = bidirectional
        self.output_dim = hidden_dim * 2 if bidirectional else hidden_dim
        
    def forward(self, x):
        # Input shape: [batch, channels, time]
        batch, channels, time = x.shape
        x = x.permute(0, 2, 1)  # -> [batch, time, channels]
        
        # Apply GRU
        outputs, _ = self.gru(x)
        
        # Back to original shape
        outputs = outputs.permute(0, 2, 1)  # -> [batch, channels, time]
        
        return outputs


class PerceptualAudioCodec(nn.Module):
    def __init__(self, 
                 target_bitrate=64000,  # bits per second
                 sample_rate=44100,
                 frame_length=1.0,     # seconds
                 latent_dim=64,
                 use_rnn=True):
        super(PerceptualAudioCodec, self).__init__()
        
        # Store parameters
        self.target_bitrate = target_bitrate
        self.sample_rate = sample_rate
        self.frame_length = frame_length
        self.latent_dim = latent_dim
        self.use_rnn = use_rnn
        
        # Calculate sample-level compression
        input_samples = int(sample_rate * frame_length)
        reduction_factor = 32  # Each encoding step reduces by factor of 2 (2^5=32)
        self.latent_samples = input_samples // reduction_factor
        
        # Calculate bits per latent sample to meet target bitrate
        # Formula: target_bitrate = latent_samples * bits_per_latent / frame_length
        bits_per_latent = target_bitrate * frame_length / self.latent_samples
        
        # Calculate required codebook size (2^bits_per_symbol)
        required_codebook_size = 2 ** math.ceil(bits_per_latent)
        self.codebook_size = min(8192, max(64, required_codebook_size))  # Reasonable bounds
        
        # Log configuration
        print(f"Target bitrate: {target_bitrate/1000:.1f} kbps")
        print(f"Latent samples: {self.latent_samples}")
        print(f"Bits per latent: {bits_per_latent:.2f}")
        print(f"Codebook size: {self.codebook_size} (log2: {math.log2(self.codebook_size):.2f} bits)")
        print(f"Estimated actual bitrate: {self.latent_samples * math.log2(self.codebook_size) / frame_length / 1000:.1f} kbps")
        
        # CNN Encoder
        self.encoder = nn.Sequential(
            EncoderBlock(1, 32),
            EncoderBlock(32, 64),
            EncoderBlock(64, 128),
            EncoderBlock(128, 256),
            EncoderBlock(256, latent_dim)
        )
        
        # Optional RNN for temporal modeling
        if use_rnn:
            self.temporal_encoder = TemporalEncoder(
                input_dim=latent_dim,
                hidden_dim=latent_dim // 2,
                num_layers=2,
                bidirectional=True
            )
        
        # Vector Quantizer for bitrate control
        self.quantizer = VectorQuantizer(
            embedding_dim=latent_dim,
            num_embeddings=self.codebook_size
        )
        
        # CNN Decoder
        self.decoder = nn.Sequential(
            DecoderBlock(latent_dim, 256),
            DecoderBlock(256, 128),
            DecoderBlock(128, 64),
            DecoderBlock(64, 32),
            DecoderBlock(32, 1)
        )
        
    def encode(self, x):
        # CNN encoding
        features = self.encoder(x)
        
        # Optional RNN processing
        if self.use_rnn:
            features = self.temporal_encoder(features)
        
        # Vector quantization
        quantized, vq_loss, indices = self.quantizer(features)

        return quantized, vq_loss, indices
    
    def decode(self, x):
        # CNN decoding
        x = self.decoder(x)
        return x
    
    def forward(self, x):
        # Encode

        quantized, vq_loss, _ = self.encode(x)
        
        # Decode
        reconstructed = self.decode(quantized)
        
        return reconstructed, vq_loss
    
    def compress(self, x):
        """Convert audio to compressed representation (indices)"""
        with torch.no_grad():
            features = self.encoder(x)
            if self.use_rnn:
                features = self.temporal_encoder(features)
            indices = self.quantizer.quantize(features)
        return indices
    
    def decompress(self, indices):
        """Convert indices back to audio"""
        with torch.no_grad():
            batch_size = 1
            time_steps = indices.shape[0] // batch_size
            
            # Reshape indices and get embeddings
            indices = indices.view(batch_size, time_steps)
            quantized = self.codebook(indices)
            quantized = quantized.permute(0, 2, 1)  # [batch, embedding_dim, time]
            
            # Decode
            audio = self.decoder(quantized)
        return audio
    
    def get_actual_bitrate(self):
        """Calculate actual bitrate in kbps"""
        bits_per_sample = math.log2(self.codebook_size)
        bits_per_second = self.latent_samples * bits_per_sample / self.frame_length
        return bits_per_second / 1000  # Convert to kbps


# Training function
def train_audio_codec(model, dataloader, epochs=10, learning_rate=0.001):
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode='min', factor=0.5, patience=3, verbose=True
    )
    loss_fn = cdpam.CDPAM()
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = model.to(device)
    
    for epoch in range(epochs):
        model.train()
        total_loss = 0
        perceptual_losses = 0
        vq_losses = 0
        
        for batch_idx, audio in enumerate(dataloader):
            audio = audio.unsqueeze(1).to(device)  # Add channel dimension
            
            # Forward pass
            reconstructed, vq_loss = model(audio)
            vq_loss = vq_loss.sum()
            
            # Calculate perceptual loss
            p_loss = loss_fn.forward(audio.squeeze(1), reconstructed.squeeze(1))
            p_loss = p_loss.sum()

            # Calculate mse loss
            mse_loss = F.mse_loss(audio.squeeze(1), reconstructed.squeeze(1)[:, :audio.size(2)])
            
            # Total loss
            loss = mse_loss + p_loss + 0.01 * vq_loss

            # Backward pass
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            
            # Track losses
            total_loss += loss.item()
            perceptual_losses += p_loss.item()
            vq_losses += vq_loss.item()
            
            if batch_idx % 100 == 0:
                print(f"Epoch {epoch}, Batch {batch_idx}/{len(dataloader)}, "
                      f"Loss: {loss.item():.4f}, "
                      f"Perceptual: {p_loss.item():.4f}, "
                      f"VQ: {vq_loss.item():.4f}")
        
        # End of epoch
        avg_loss = total_loss / len(dataloader)
        avg_perceptual = perceptual_losses / len(dataloader)
        avg_vq = vq_losses / len(dataloader)
        
        print(f"Epoch {epoch} completed. "
              f"Avg Loss: {avg_loss:.4f}, "
              f"Avg Perceptual: {avg_perceptual:.4f}, "
              f"Avg VQ: {avg_vq:.4f}")
        
        # Update learning rate
        scheduler.step(avg_loss)
        
        # Validate bitrate
        print(f"Current bitrate: {model.get_actual_bitrate():.1f} kbps")
    
    return model
