import torch
import torch.nn as nn

class TransientDetectionModel(nn.Module):
    def __init__(self):
        super(TransientDetectionModel, self).__init__()
        
        self.cnn = nn.Sequential(
            nn.Conv2d(1, 16, kernel_size=3, stride=1, padding=1),
            nn.BatchNorm2d(16),
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=2, stride=2),  # Reduces spatial dimensions by half
            
            nn.Conv2d(16, 32, kernel_size=3, stride=1, padding=1),
            nn.BatchNorm2d(32),
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=2, stride=2), 
            
            nn.Conv2d(32, 64, kernel_size=3, stride=1, padding=1),
            nn.BatchNorm2d(64),
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=2, stride=2),  
            
            nn.Conv2d(64, 128, kernel_size=3, stride=1, padding=1),
            nn.BatchNorm2d(128),
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=2, stride=2) 
        )
        
        self.lstm = None
        self.pre_lstm = None
        
        self.classifier = nn.Sequential(
            nn.Linear(256, 64),  # 256 = 128*2 (bidirectional LSTM)
            nn.ReLU(),
            nn.Dropout(0.5),
            nn.Linear(64, 1),
            nn.Sigmoid()
        )
    
    def forward(self, x):
        if x.dim() == 3:
            x = x.unsqueeze(1)
            #print(f"Added channel dimension: {x.shape}")
        
        batch_size = x.size(0)
        #print(f"Input shape: {x.shape}")
        
        # Process through CNN
        x = self.cnn(x)
        #print(f"After CNN: {x.shape}")
        
        # Reshape for LSTM: [batch, channels, height, width] -> [batch, width, channels*height]
        lstm_in = x.permute(0, 3, 1, 2)
        lstm_in = lstm_in.reshape(batch_size, lstm_in.size(1), -1)
        #print(f"Reshaped for LSTM: {lstm_in.shape}")
        
        lstm_input_features = lstm_in.size(2)
        
        if self.pre_lstm is None:
            self.pre_lstm = nn.Sequential(
                nn.Linear(lstm_input_features, 256),
                nn.LayerNorm(256),
                nn.ReLU()
            ).to(x.device)
            #print(f"Created pre-LSTM layer with input size {lstm_input_features}")
        
        lstm_in = self.pre_lstm(lstm_in)
        #print(f"After pre-LSTM: {lstm_in.shape}")
        
        if self.lstm is None:
            self.lstm = nn.LSTM(
                input_size=256, 
                hidden_size=128,
                num_layers=2,
                batch_first=True,
                dropout=0.3,
                bidirectional=True
            ).to(x.device)
            #print("Created LSTM with input size 256")
        
        lstm_out, _ = self.lstm(lstm_in)
        #print(f"After LSTM: {lstm_out.shape}")
        
        final_out = lstm_out[:, -1, :]
        
        output = self.classifier(final_out)
        
        return output