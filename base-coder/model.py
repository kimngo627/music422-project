import torch
import torch.nn as nn

class TransientDetectionModel(nn.Module):
    def __init__(self):
        super(TransientDetectionModel, self).__init__()
        
        # CNN layers for feature extraction
        self.cnn = nn.Sequential(
            # First CNN block
            nn.Conv2d(1, 16, kernel_size=3, stride=1, padding=1),
            nn.BatchNorm2d(16),
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=2, stride=2),
            
            # Second CNN block
            nn.Conv2d(16, 32, kernel_size=3, stride=1, padding=1),
            nn.BatchNorm2d(32),
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=2, stride=2),
            
            # Third CNN block
            nn.Conv2d(32, 64, kernel_size=3, stride=1, padding=1),
            nn.BatchNorm2d(64),
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=2, stride=2),
            
            # Fourth CNN block
            nn.Conv2d(64, 128, kernel_size=3, stride=1, padding=1),
            nn.BatchNorm2d(128),
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=2, stride=2)
        )
        
        # Pre-LSTM layer with corrected input dimension (4096)
        self.pre_lstm = nn.Sequential(
            nn.Linear(4096, 256),  # Changed from 128 to 4096
            nn.LayerNorm(256),
            nn.ReLU()
        )
        
        # LSTM layer
        self.lstm = nn.LSTM(
            input_size=256,
            hidden_size=128,
            num_layers=2,
            batch_first=True,
            dropout=0.3,
            bidirectional=True
        )
        
        # Classifier
        self.classifier = nn.Sequential(
            nn.Linear(256, 64),  # 256 = 128*2 (bidirectional LSTM)
            nn.ReLU(),
            nn.Dropout(0.5),
            nn.Linear(64, 1),
            nn.Sigmoid()
        )
    
    def forward(self, x):
        # Add channel dimension if needed
        if x.dim() == 3:
            x = x.unsqueeze(1)
        
        batch_size = x.size(0)
        
        # Process through CNN
        x = self.cnn(x)
        
        # Reshape for pre-LSTM
        # The flattened dimension must be 4096 to match the checkpoint
        lstm_in = x.permute(0, 3, 1, 2)
        lstm_in = lstm_in.reshape(batch_size, lstm_in.size(1), -1)
        
        # Apply pre-LSTM
        lstm_in = self.pre_lstm(lstm_in)
        
        # Apply LSTM
        lstm_out, _ = self.lstm(lstm_in)
        
        # Take final timestep output
        final_out = lstm_out[:, -1, :]
        
        # Apply classifier
        output = self.classifier(final_out)
        
        return output