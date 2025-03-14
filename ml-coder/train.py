from arch import train_audio_codec, PerceptualAudioCodec
from audio_process import create_train_test_dataloaders
from test import evaluate_codec
import torch

model = PerceptualAudioCodec(
            target_bitrate=128000,  # 128 kbps
            sample_rate=44100,
            frame_length=1.0,
            latent_dim=256,
            use_rnn=True
        )

train_loader, test_loader = create_train_test_dataloaders(
        audio_dir='../test-items/SQAM_WAV',
        batch_size=4,
        sample_rate=44100,
        snippet_length=1.0,
        test_size=0.2
    )

# trained_model = train_audio_codec(model, train_loader, epochs=50)
# torch.save(trained_model.state_dict(), 'trained_model_latent512_codebook2**20.pth')

model.load_state_dict(torch.load('trained_model_latent256_codebook8192.pth'))

test_results = evaluate_codec(
    model=model,
    test_loader=test_loader,
    save_samples=True,
    num_samples_to_save=5,
    output_dir='./codec_samples_latent512_test'
)

train_results = evaluate_codec(
    model=model,
    test_loader=train_loader,
    save_samples=True,
    num_samples_to_save=5,
    output_dir='./codec_samples_latent512_train'
)

print('LATENT DIM 256')
print(f'test_results: {test_results}')
print(f'train_results: {train_results}')

# model = PerceptualAudioCodec(
#             target_bitrate=128000,  # 128 kbps
#             sample_rate=44100,
#             frame_length=1.0,
#             latent_dim=64,
#             use_rnn=True
#         )

# trained_model = train_audio_codec(model, train_loader, epochs=50)
# torch.save(trained_model.state_dict(), 'trained_model.pth')

# # model.load_state_dict(torch.load('trained_model.pth'))

# test_results = evaluate_codec(
#     model=model,
#     test_loader=test_loader,
#     save_samples=True,
#     num_samples_to_save=5,
#     output_dir='./codec_samples_latent64_test'
# )

# train_results = evaluate_codec(
#     model=model,
#     test_loader=train_loader,
#     save_samples=True,
#     num_samples_to_save=5,
#     output_dir='./codec_samples_latent64_train'
# )
