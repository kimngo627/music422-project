"""
odf.py -- Onset Detection Function (ODF) based transient detector

"""

import numpy as np

class ODFTransientDetector:
    """
    Class for detecting transients in audio data using the Onset Detection Function (ODF)
    approach
    """
    
    def __init__(self, sample_rate=48000, block_size=1024, 
                 threshold_factor=1.5, min_threshold=30.0,
                 history_size=5):
        """
        Initialize the transient detector

        """
        self.sample_rate = sample_rate
        self.block_size = block_size
        self.threshold_factor = threshold_factor
        self.min_threshold = min_threshold
        self.history_size = history_size
        
        self.window = np.hanning(block_size)

        self.bands = [
            (0, 500),        # Low frequency band
            (500, 2000),     # Mid-low frequency band
            (2000, 8000),    # Mid-high frequency band
            (8000, 20000)    # High frequency band
        ]
        
        self.band_indices = []
        for low_freq, high_freq in self.bands:
            low_bin = int(low_freq * block_size / sample_rate)
            high_bin = min(int(high_freq * block_size / sample_rate), block_size // 2)
            self.band_indices.append((low_bin, high_bin))
        
        # Weights for different frequency bands
        self.band_weights = [0, 0, 2.0, 1.5]
        
        # State variables for block-by-block processing
        self.prev_energies = [0] * len(self.bands)
        self.detection_history = []
    
    def detect(self, audio_block):
        """
        Detect if there's a transient in the current audio block
        """
        if len(audio_block) != self.block_size:
            raise ValueError(f"Expected block size {self.block_size}, got {len(audio_block)}")
        
        windowed_data = audio_block * self.window
        fft_data = np.abs(np.fft.fft(windowed_data))
        
        # Calculate band energies
        current_energies = []
        for j, (low_bin, high_bin) in enumerate(self.band_indices):
            # Calculate energy in this band
            band_energy = np.sum(fft_data[low_bin:high_bin]**2)
            band_energy_db = 10 * np.log10(band_energy + 1e-10)
            current_energies.append(band_energy_db)
        
        # Calculate ODF as weighted sum of energy differences across bands
        odf_value = 0
        for j in range(len(self.bands)):
            energy_diff = max(0, current_energies[j] - self.prev_energies[j])
            odf_value += self.band_weights[j] * energy_diff
        
        self.prev_energies = current_energies
        
        is_transient = odf_value > self.min_threshold
        
        return is_transient