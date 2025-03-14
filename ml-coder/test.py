import torch
import torchaudio
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from tqdm import tqdm
import cdpam


def evaluate_codec(model, test_loader, save_samples=True, num_samples_to_save=5, output_dir='./samples'):
    """
    Evaluate the audio codec on test data
    
    Args:
        model: Trained audio codec model
        test_loader: DataLoader with test audio snippets
        save_samples: Whether to save audio samples
        num_samples_to_save: Number of sample pairs to save
        output_dir: Directory to save samples
        
    Returns:
        dict: Dictionary with evaluation metrics
    """
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = model.to(device)
    model.eval()
    
    # Create output directory if needed
    if save_samples:
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True, parents=True)
    
    # Initialize metrics
    mse_losses = []
    snr_values = []
    compression_ratios = []
    perceptual_losses = []
    
    # Track samples for saving
    samples_to_save = []
    
    with torch.no_grad():
        for batch_idx, audio in enumerate(tqdm(test_loader, desc="Evaluating")):
            # Add channel dimension and move to device
            audio = audio.unsqueeze(1).to(device)
            
            # Run through model
            reconstructed, _ = model(audio)
            
            # Calculate metrics for each sample in batch
            for i in range(audio.size(0)):
                original = audio[i].cpu()
                recon = reconstructed[i].cpu()[:, :original.size(1)]

                # Mean Squared Error
                mse = torch.mean((original - recon) ** 2).item()
                mse_losses.append(mse)
                
                # Signal-to-Noise Ratio (SNR)
                signal_power = torch.mean(original ** 2).item()
                noise_power = mse
                snr = 10 * np.log10(signal_power / noise_power) if noise_power > 0 else 100
                snr_values.append(snr)
                
                # Perceptual loss (if model has perceptual loss function)
                loss_fn = cdpam.CDPAM()
                p_loss = loss_fn.forward(original, recon).item()
                perceptual_losses.append(p_loss)
                
                # Save some samples
                if save_samples and len(samples_to_save) < num_samples_to_save:
                    samples_to_save.append((original.squeeze(0).numpy(), recon.squeeze(0).numpy()))
            
            # Calculate compression ratio
            if hasattr(model, 'get_actual_bitrate'):
                original_bitrate = 16 * model.sample_rate  # 16-bit PCM audio
                compressed_bitrate = model.get_actual_bitrate() * 1000  # convert kbps to bps
                ratio = original_bitrate / compressed_bitrate
                compression_ratios.append(ratio)
            
            # Don't process the entire dataset if we're just looking for samples
            if save_samples and len(samples_to_save) >= num_samples_to_save and batch_idx > 5:
                break
    
    # Save audio samples and visualizations
    if save_samples:
        for i, (original, reconstructed) in enumerate(samples_to_save):
            # Save audio files
            original_path = output_path / f"sample_{i}_original.wav"
            recon_path = output_path / f"sample_{i}_reconstructed.wav"
            
            sample_rate = model.sample_rate if hasattr(model, 'sample_rate') else 44100
            torchaudio.save(original_path, torch.tensor(original).unsqueeze(0), sample_rate)
            torchaudio.save(recon_path, torch.tensor(reconstructed).unsqueeze(0), sample_rate)
            
            # Create and save visualization
            plt.figure(figsize=(12, 8))
            
            # Waveform comparison
            plt.subplot(2, 1, 1)
            plt.plot(original, label='Original')
            plt.plot(reconstructed, label='Reconstructed', alpha=0.8)
            plt.legend()
            plt.title(f'Sample {i} - Waveform Comparison')
            plt.xlabel('Sample')
            plt.ylabel('Amplitude')
            
            # Error
            plt.subplot(2, 1, 2)
            plt.plot(np.abs(original - reconstructed))
            plt.title('Absolute Error')
            plt.xlabel('Sample')
            plt.ylabel('Error Magnitude')
            
            plt.tight_layout()
            plt.savefig(output_path / f"sample_{i}_comparison.png")
            plt.close()
    
    # Calculate average metrics
    avg_mse = np.mean(mse_losses)
    avg_snr = np.mean(snr_values)
    avg_compression = np.mean(compression_ratios) if compression_ratios else None
    avg_perceptual = np.mean(perceptual_losses) if perceptual_losses else None
    
    # Prepare results
    results = {
        'mse': avg_mse,
        'snr_db': avg_snr,
        'compression_ratio': avg_compression,
        'perceptual_loss': avg_perceptual,
        'bitrate_kbps': model.get_actual_bitrate() if hasattr(model, 'get_actual_bitrate') else None
    }
    
    print("\nEvaluation Results:")
    print(f"MSE: {avg_mse:.6f}")
    print(f"SNR: {avg_snr:.2f} dB")
    if avg_compression:
        print(f"Compression Ratio: {avg_compression:.1f}x")
    if avg_perceptual:
        print(f"Perceptual Loss: {avg_perceptual:.6f}")
    if hasattr(model, 'get_actual_bitrate'):
        print(f"Bitrate: {model.get_actual_bitrate():.1f} kbps")
    
    return results


def encode_decode_file(model, file_path, output_path=None):
    """
    Process a single audio file through the codec
    
    Args:
        model: Trained audio codec model
        file_path: Path to input audio file
        output_path: Path to save reconstructed audio (None for auto-generate)
        
    Returns:
        tuple: (original audio, reconstructed audio)
    """
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = model.to(device)
    model.eval()
    
    # Load audio file
    waveform, sample_rate = torchaudio.load(file_path)
    
    # Convert to mono if stereo
    if waveform.shape[0] > 1:
        waveform = torch.mean(waveform, dim=0, keepdim=True)
    
    # Resample if needed
    if hasattr(model, 'sample_rate') and sample_rate != model.sample_rate:
        resampler = torchaudio.transforms.Resample(
            orig_freq=sample_rate,
            new_freq=model.sample_rate
        )
        waveform = resampler(waveform)
        sample_rate = model.sample_rate
    
    # Process in segments if the audio is long
    segment_length = 44100  # 1 second at 44.1kHz
    reconstructed_segments = []
    
    with torch.no_grad():
        for i in range(0, waveform.shape[1], segment_length):
            # Extract segment
            segment = waveform[:, i:i+segment_length]
            
            # Pad if shorter than segment_length
            if segment.shape[1] < segment_length:
                padded = torch.zeros(1, segment_length, device=device)
                padded[:, :segment.shape[1]] = segment
                segment = padded
            
            # Process through model
            segment = segment.to(device)
            reconstructed, _ = model(segment.unsqueeze(0))
            
            # Add to list of reconstructed segments
            reconstructed_segments.append(reconstructed.squeeze(0).cpu())
    
    # Concatenate segments
    reconstructed_waveform = torch.cat(reconstructed_segments, dim=1)
    
    # Trim to original length
    reconstructed_waveform = reconstructed_waveform[:, :waveform.shape[1]]
    
    # Save if output path specified
    if output_path is None:
        input_path = Path(file_path)
        output_path = input_path.parent / f"{input_path.stem}_reconstructed{input_path.suffix}"
    
    torchaudio.save(output_path, reconstructed_waveform, sample_rate)
    
    print(f"Processed {file_path}")
    print(f"Saved reconstructed audio to {output_path}")
    
    return waveform, reconstructed_waveform
