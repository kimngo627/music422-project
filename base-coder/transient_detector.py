"""
transient_detector.py -- Interface for different transient detection methods
"""

import numpy as np
import torch
from enum import Enum
import librosa

from model import TransientDetectionModel

from odf import ODFTransientDetector

class DetectionMethod(Enum):
    """Enumeration of available transient detection methods"""
    MODEL = "model"                
    ODF = "odf"                 

class BlockType(Enum):
    """Enumeration of block types for block switching"""
    LONG = bin(0)    # 0b00
    SHORT = bin(1)   # 0b01
    START = bin(2)   # 0b10
    STOP = bin(3)    # 0b11

class BlockState:
    """Manages state transitions for block switching based on transient detection"""
    
    def __init__(self):
        self.current_state = BlockType.LONG
        
        self.transitions = {
            # current_state: {transient_detected: next_state, not_detected: next_state}
            BlockType.LONG: {True: BlockType.START, False: BlockType.LONG},
            BlockType.START: {True: BlockType.SHORT, False: BlockType.SHORT},
            BlockType.SHORT: {True: BlockType.SHORT, False: BlockType.STOP},
            BlockType.STOP: {True: BlockType.START, False: BlockType.LONG}
        }
    
    def update(self, transient_detected):
        self.current_state = self.transitions[self.current_state][transient_detected]
        return self.current_state

class TransientDetector:
    """Class for detecting transients in audio data with ml or odf method"""
    
    def __init__(self, method=DetectionMethod.MODEL, model_path=None, device='cpu',
                 threshold_factor=1.5, min_threshold=30.0, sample_rate=48000, block_size=1024):
        """
        Initialize the transient detector
        """
        self.method = method
        self.device = device
        self.block_state = BlockState()
        self.sample_rate = sample_rate
        self.block_size = block_size
        
        if method == DetectionMethod.MODEL:
            if model_path is None:
                raise ValueError("Model path must be provided when using MODEL method")
            self.model = TransientDetectionModel()
            checkpoint = torch.load(model_path, map_location=device)
            if 'model_state_dict' in checkpoint:
                self.model.load_state_dict(checkpoint['model_state_dict'])
            else:
                self.model.load_state_dict(checkpoint)
            self.model.eval()
            self.model.to(device)
            
        elif method == DetectionMethod.ODF:
            self.odf_detector = ODFTransientDetector(
                sample_rate=sample_rate,
                block_size=block_size,
                threshold_factor=threshold_factor,
                min_threshold=min_threshold
            )
    
    def detect(self, audio_block, sample_rate=None):
        """
        Detect if there's a transient in the audio block
        """
        if sample_rate is None:
            sample_rate = self.sample_rate
            
        if self.method == DetectionMethod.MODEL:
            return self._detect_with_model(audio_block, sample_rate)
        elif self.method == DetectionMethod.ODF:
            return self._detect_with_odf(audio_block)
        else:
            raise ValueError(f"Unknown detection method: {self.method}")
    
    def _detect_with_model(self, audio_block, sample_rate):
        """Use trained model for transient detection"""
        D = librosa.stft(audio_block, n_fft=1024, hop_length=64, win_length=128)
        spec = np.abs(D) ** 2
        log_spec = librosa.power_to_db(spec, ref=np.max)
        
        # Convert to tensor and add batch dimension
        features = torch.from_numpy(log_spec).float().unsqueeze(0)  # [1, freq_bins, time_frames]
        features = features.permute(0, 3, 2, 1)

        if features.dim() > 3:
            # If we somehow have more than 3 dimensions, take just the first channel
            features = features[:, 0, :, :]
        #print(features.shape)
        
        with torch.no_grad():
            prediction = self.model(features.to(self.device))
            is_transient = prediction.item() > 0.5
        
        return is_transient
    
    def _detect_with_odf(self, audio_block):
        """Use ODF method for transient detection"""
        return self.odf_detector.detect(audio_block)
    
    def get_block_type(self, next_audio_block, sample_rate=None):
        """
        Detect transients and determine the appropriate block type
        """
        # For ODF method, use its built-in block type management
        if self.method == DetectionMethod.ODF:
            return self.odf_detector.get_block_type(next_audio_block)
        
        # For other methods, use the block state machine
        # Detect transient
        is_transient = self.detect(next_audio_block, sample_rate)
        
        # Update block state based on transient detection
        block_type = self.block_state.update(is_transient)
        
        return block_type
    
if __name__ == "__main__":
    # Test with a sample audio file
    import soundfile as sf
    import matplotlib.pyplot as plt
    import numpy as np
    from blockswitching import SineWindowFunc, TransitionWindow
    
    audio_path = "/Users/kimngo/Downloads/music422/music422-project/base-coder/castanets.wav"
    audio, sr = sf.read(audio_path)
    
    long_block_size = 1024
    short_block_size = 128
    hop_size = long_block_size // 2
    start_block = 115
    end_block = 125
    
    start_sample = start_block * hop_size
    end_sample = (end_block + 2) * hop_size 
    
    if end_sample > len(audio):
        audio = np.pad(audio, (0, end_sample - len(audio)))
    
    segment_audio = audio[start_sample:end_sample]
    segment_length = len(segment_audio)
    
    print(f"Visualizing blocks {start_block}-{end_block}")
    print(f"Sample range: {start_sample}-{end_sample}")
    print(f"Segment length: {segment_length} samples")

    detector = TransientDetector(
        method=DetectionMethod.MODEL,
        model_path="/Users/kimngo/Downloads/music422/music422-project/base-coder/transient_detection_model_final.pth"
    )
    
    block_types = []
    transient_positions = []
    
    for i in range(start_block, end_block + 1):
        block_start = i * hop_size
        block_end = block_start + long_block_size
        next_block_start = block_start + hop_size
        next_block_end = next_block_start + long_block_size
        
        if block_end > len(audio):
            block = np.pad(audio[block_start:], (0, block_end - len(audio)))
        else:
            block = audio[block_start:block_end]
        if next_block_end > len(audio):
            next_block = np.pad(audio[next_block_start:], (0, next_block_end - len(audio)))
        else:
            next_block = audio[next_block_start:next_block_end]
        
        is_transient = detector.detect(next_block, sr)
        if is_transient:
            transient_positions.append(next_block_start - start_sample)  
        
        block_type = detector.get_block_type(next_block, sr)
        block_types.append(block_type)
        print(f"Block {i}: {block_type.name}, Transient: {is_transient}")
    
    plt.figure(figsize=(15, 6))
    
    segment_time = np.arange(segment_length) / sr
    plt.plot(segment_time, segment_audio, 'b', alpha=0.7, label='Audio')
    
    windows = []
    window_positions = []
    
    for i, wtype in enumerate(block_types):
        block_idx = start_block + i
        block_pos = block_idx * hop_size - start_sample
        
        if wtype == BlockType.LONG:
            window_values = SineWindowFunc(long_block_size)
            window_label = "Long Window"
            
        elif wtype == BlockType.SHORT:
            window_values = np.zeros(long_block_size)
            for j in range(long_block_size // short_block_size):
                short_start_idx = j * short_block_size // 2
                short_end_idx = min(short_start_idx + short_block_size, long_block_size)
                if short_end_idx <= short_start_idx:
                    continue
                
                short_window = SineWindowFunc(short_block_size)
                short_length = short_end_idx - short_start_idx
                window_values[short_start_idx:short_end_idx] = short_window[:short_length]
            window_label = "Short Windows"
            
        elif wtype == BlockType.START:
            ones_array = np.ones(long_block_size)
            window_values = TransitionWindow(ones_array, SineWindowFunc, long_block_size, short_block_size, start=True)
            window_label = "Start Transition"
            
        elif wtype == BlockType.STOP:
            ones_array = np.ones(long_block_size)
            window_values = TransitionWindow(ones_array, SineWindowFunc, long_block_size, short_block_size, start=False)
            window_label = "Stop Transition"
            
        windows.append((window_values, window_label))
        window_positions.append(block_pos)
    
    used_labels = set()
    
    for i, ((window_values, label), position) in enumerate(zip(windows, window_positions)):
        window_time = np.arange(position, position + len(window_values)) / sr
        
        plt.plot(window_time, window_values, alpha=0.7)
    
    for pos in transient_positions:
        pos_time = pos / sr
        plt.plot(pos_time, 0.9, 'ro', markersize=8, label='Transient' if 'Transient' not in used_labels else "")
        used_labels.add('Transient')
    
    plt.title("Window Shapes During Block Switching (Blocks {}-{})".format(start_block, end_block))
    plt.xlabel("Time (seconds)")
    plt.ylabel("Amplitude")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()