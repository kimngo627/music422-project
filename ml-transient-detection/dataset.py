import os
import json
import torch
from torch.utils.data import Dataset
import torchaudio
import numpy as np
import librosa
import matplotlib.pyplot as plt

class TransientDataset(Dataset):
    def __init__(self, audio_dir, labels_dir, context_blocks=0, transform=None):
        self.audio_dir = audio_dir
        self.labels_dir = labels_dir
        self.transform = transform
        self.context_blocks = context_blocks
        self.data = []
        
        self.BLOCK_SIZE = 1024  
        self.HOP_SIZE = 512     
        
        # Spectrogram parameters 
        self.n_fft = 1024      
        self.hop_length = 64  
        self.win_length = 128  
        
        self._load_data()
        
        # Debug: Check first few items
        print("\n==== SPECTROGRAM DEBUG ====")
        for i in [200, 670, 900]:
            spec, label = self.__getitem__(i)
            print(f"Item {i}: shape={spec.shape}, label={label.item()}")
        print("===========================\n")
    
    def _load_data(self):
        """Load audio files and their labels"""
        audio_files = sorted([f for f in os.listdir(self.audio_dir) 
                             if f.endswith(".wav") and not f.startswith("._")])
        
        processed_files_count = 0
        
        for audio_file in audio_files:
            label_file = os.path.join(self.labels_dir, audio_file.replace(".wav", ".json"))
            
            if not os.path.exists(label_file):
                print(f"Warning: No label file found for {audio_file}")
                continue
                
            with open(label_file, 'r') as f:
                label_data = json.load(f)
                
            transient_blocks = label_data.get("transient_blocks", [])
            if not transient_blocks:
                print(f"Note: File {audio_file} has no transient blocks")
                continue
                
            for block in transient_blocks:
                if "transient" in block and block["transient"] is not None:  # Skip unlabeled blocks
                    self.data.append({
                        "audio_file": os.path.join(self.audio_dir, audio_file),
                        "block_index": block["block"],
                        "start_sample": block["start_sample"],
                        "label": block["transient"]
                    })
            
            processed_files_count += 1
        
        print(f"Loaded {len(self.data)} labeled blocks from {processed_files_count} audio files")
    
    def __len__(self):
        return len(self.data)
    
    def __getitem__(self, idx):
        item = self.data[idx]
        
        try:
            waveform, sample_rate = torchaudio.load(item["audio_file"])
            
            # Debug info
            if idx in [200, 670, 900]: 
                print(f"\nItem {idx}:")
                print(f"  Audio file: {item['audio_file']}")
                print(f"  Sample rate: {sample_rate} Hz")
                print(f"  Waveform shape: {waveform.shape}")
                print(f"  Block index: {item['block_index']}")
                print(f"  Start sample: {item['start_sample']}")
            
            start_sample = item["start_sample"]
            end_sample = start_sample + self.BLOCK_SIZE
            
            if end_sample > waveform.shape[1]:
                pad_amount = end_sample - waveform.shape[1]
                waveform = torch.nn.functional.pad(waveform, (0, pad_amount))
            
            audio_block = waveform[0, start_sample:end_sample]
            
            if idx in [200, 670, 900]:
                print(f"  Audio block shape: {audio_block.shape}")
                print(f"  Audio duration: {audio_block.shape[0]/sample_rate:.4f} seconds")
            
            audio_np = audio_block.numpy()
            
            expected_frames = 1 + (self.BLOCK_SIZE - self.win_length) // self.hop_length
            
            if idx in [200, 670, 900]:
                print(f"  STFT parameters: n_fft={self.n_fft}, hop_length={self.hop_length}, win_length={self.win_length}")
                print(f"  Expected time frames: {expected_frames}")
            
            D = librosa.stft(audio_np, 
                            n_fft=self.n_fft, 
                            hop_length=self.hop_length, 
                            win_length=self.win_length)
            
            spec = np.abs(D) ** 2
            log_spec = librosa.power_to_db(spec, ref=np.max)
            
            if idx in [200, 670, 900]:
                print(f"  Spectrogram shape: {log_spec.shape}")
                # Visualize the spectrogram for debugging
                plt.figure(figsize=(10, 4))
                librosa.display.specshow(log_spec, sr=sample_rate, hop_length=self.hop_length, x_axis='time', y_axis='linear')
                plt.colorbar(format='%+2.0f dB')
                plt.title(f'Spectrogram (item {idx}, label={item["label"]})')
                plt.tight_layout()
                plt.show()
            
            features = torch.from_numpy(log_spec).float()
            
            if self.transform:
                features = self.transform(features)
            
            label = torch.tensor(item["label"], dtype=torch.long)
            
            return features, label
            
        except Exception as e:
            print(f"Error loading item {idx}: {e}")
            import traceback
            traceback.print_exc()
            return torch.zeros(self.n_fft//2 + 1, 16), torch.tensor(0, dtype=torch.long)