import os
import torch
import torchaudio
import numpy as np
from torch.utils.data import Dataset, DataLoader, random_split
from pathlib import Path


class AudioSnippetDataset(Dataset):
    """
    Dataset for loading and preprocessing audio files into fixed-length snippets.
    """
    def __init__(self, 
                 audio_dir, 
                 sample_rate=44100, 
                 snippet_length=1.0,
                 mono=True,
                 normalize=True):
        """
        Args:
            audio_dir (str): Directory containing audio files (.wav)
            sample_rate (int): Target sample rate
            snippet_length (float): Length of snippets in seconds
            mono (bool): Convert to mono if True
            normalize (bool): Normalize audio to range [-1, 1] if True
        """
        self.audio_dir = Path(audio_dir)
        self.sample_rate = sample_rate
        self.snippet_length = snippet_length
        self.mono = mono
        self.normalize = normalize
        
        # Number of samples in each snippet
        self.snippet_samples = int(sample_rate * snippet_length)
        
        # Find all wav files
        self.audio_files = list(self.audio_dir.glob('**/*.wav'))
        
        # Pre-calculate total number of snippets
        self.snippets = []
        for audio_file in self.audio_files:
            metadata = torchaudio.info(audio_file)
            duration = metadata.num_frames / metadata.sample_rate
            num_snippets = max(1, int(duration / snippet_length))
            
            for i in range(num_snippets):
                self.snippets.append((audio_file, i))
                
        print(f"Found {len(self.audio_files)} audio files, {len(self.snippets)} snippets")
    
    def __len__(self):
        return len(self.snippets)
    
    def __getitem__(self, idx):
        """
        Returns a 1-second snippet of audio as a tensor.
        """
        audio_file, snippet_idx = self.snippets[idx]
        
        # Load audio file
        waveform, original_sample_rate = torchaudio.load(audio_file)
        
        # Convert to mono if needed
        if self.mono and waveform.shape[0] > 1:
            waveform = torch.mean(waveform, dim=0, keepdim=True)
        
        # Resample if needed
        if original_sample_rate != self.sample_rate:
            resampler = torchaudio.transforms.Resample(
                orig_freq=original_sample_rate, 
                new_freq=self.sample_rate
            )
            waveform = resampler(waveform)
        
        # Calculate snippet start position
        start_sample = int(snippet_idx * self.snippet_samples)
        
        # Handle case where audio is shorter than snippet length
        if waveform.shape[1] < self.snippet_samples:
            # Pad with zeros
            padded = torch.zeros(1, self.snippet_samples)
            padded[0, :waveform.shape[1]] = waveform[0, :]
            snippet = padded
        else:
            # Extract snippet
            end_sample = min(start_sample + self.snippet_samples, waveform.shape[1])
            snippet = waveform[0, start_sample:end_sample]
            
            # Pad if snippet is shorter than required length
            if snippet.shape[0] < self.snippet_samples:
                padded = torch.zeros(self.snippet_samples)
                padded[:snippet.shape[0]] = snippet
                snippet = padded
        
        # Normalize to [-1, 1]
        if self.normalize:
            if torch.abs(snippet).max() > 0:
                snippet = snippet / torch.abs(snippet).max()
        
        # Remove channel dimension for consistency
        if snippet.dim() == 2:
            snippet = snippet.squeeze(0)
        
        return snippet


def create_train_test_dataloaders(audio_dir, 
                                batch_size=32, 
                                sample_rate=44100, 
                                snippet_length=1.0,
                                test_size=0.2,
                                num_workers=4,
                                random_seed=42):
    """
    Creates train and test DataLoaders for audio snippets.
    
    Args:
        audio_dir (str): Directory containing audio files
        batch_size (int): Batch size
        sample_rate (int): Target sample rate
        snippet_length (float): Length of snippets in seconds
        test_size (float): Proportion of data to use for testing (0.0 to 1.0)
        num_workers (int): Number of worker threads for loading data
        random_seed (int): Seed for reproducible train/test splits
        
    Returns:
        tuple: (train_dataloader, test_dataloader)
    """
    # Create the full dataset
    full_dataset = AudioSnippetDataset(
        audio_dir=audio_dir,
        sample_rate=sample_rate,
        snippet_length=snippet_length
    )
    
    # Calculate split sizes
    dataset_size = len(full_dataset)
    test_length = int(dataset_size * test_size)
    train_length = dataset_size - test_length
    
    # Create the splits
    torch.manual_seed(random_seed)
    train_dataset, test_dataset = random_split(
        full_dataset, 
        [train_length, test_length]
    )
    
    # Reset random seed to avoid affecting other random operations
    torch.manual_seed(torch.initial_seed())
    
    # Create dataloaders
    train_dataloader = DataLoader(
        train_dataset,
        batch_size=batch_size,
        shuffle=True,
        num_workers=num_workers,
        pin_memory=True
    )
    
    test_dataloader = DataLoader(
        test_dataset,
        batch_size=batch_size,
        shuffle=False,  # No need to shuffle test data
        num_workers=num_workers,
        pin_memory=True
    )
    
    print(f"Train set: {len(train_dataset)} snippets, {len(train_dataloader)} batches")
    print(f"Test set: {len(test_dataset)} snippets, {len(test_dataloader)} batches")
    
    return train_dataloader, test_dataloader


def create_train_val_test_dataloaders(audio_dir, 
                                    batch_size=32, 
                                    sample_rate=44100, 
                                    snippet_length=1.0,
                                    val_size=0.1,
                                    test_size=0.1,
                                    num_workers=4,
                                    random_seed=42):
    """
    Creates train, validation, and test DataLoaders for audio snippets.
    
    Args:
        audio_dir (str): Directory containing audio files
        batch_size (int): Batch size
        sample_rate (int): Target sample rate
        snippet_length (float): Length of snippets in seconds
        val_size (float): Proportion of data to use for validation
        test_size (float): Proportion of data to use for testing
        num_workers (int): Number of worker threads for loading data
        random_seed (int): Seed for reproducible splits
        
    Returns:
        tuple: (train_dataloader, val_dataloader, test_dataloader)
    """
    # Create the full dataset
    full_dataset = AudioSnippetDataset(
        audio_dir=audio_dir,
        sample_rate=sample_rate,
        snippet_length=snippet_length
    )
    
    # Calculate split sizes
    dataset_size = len(full_dataset)
    test_length = int(dataset_size * test_size)
    val_length = int(dataset_size * val_size)
    train_length = dataset_size - test_length - val_length
    
    # Create the splits
    torch.manual_seed(random_seed)
    train_dataset, val_dataset, test_dataset = random_split(
        full_dataset, 
        [train_length, val_length, test_length]
    )
    
    # Reset random seed
    torch.manual_seed(torch.initial_seed())
    
    # Create dataloaders
    train_dataloader = DataLoader(
        train_dataset,
        batch_size=batch_size,
        shuffle=True,
        num_workers=num_workers,
        pin_memory=True
    )
    
    val_dataloader = DataLoader(
        val_dataset,
        batch_size=batch_size,
        shuffle=False,
        num_workers=num_workers,
        pin_memory=True
    )
    
    test_dataloader = DataLoader(
        test_dataset,
        batch_size=batch_size,
        shuffle=False,
        num_workers=num_workers,
        pin_memory=True
    )
    
    print(f"Train set: {len(train_dataset)} snippets, {len(train_dataloader)} batches")
    print(f"Validation set: {len(val_dataset)} snippets, {len(val_dataloader)} batches")
    print(f"Test set: {len(test_dataset)} snippets, {len(test_dataloader)} batches")
    
    return train_dataloader, val_dataloader, test_dataloader
