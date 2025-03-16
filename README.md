# music422-project

For the entropy coding branch, you may generate encodings with the following command arguments:

1. `--input_file`: specify the wav file you'd like to encode
2. `--entropy_coding`: specify the entropy coding style among 
* `"full"`: entropy codes the full .pac file,
* `"block"`: entropy codes per-block into the .pac file, 
* `None`: does not apply any entropy coding (this is the setting by default)
3. `--data_rate`: specify the data rate in kbps; the default is 128
* _NOTE_: This is the data rate for the non-entropy coded pac file. When applying entropy coding you will achieve a lower data rate than specified. We are unable to know in advance what data rate the entropy coder will achieve. Based on emprical compression ratios, the estimate is `150 kbps` for the baseline coder achieves `128 kbps` for the per-block entropy coder (with block size 1024); `110 kbps` for the baseline coder achieves `96 kbps` for the per-block entropy coder (with block size 1024); `185 kbps` for the baseline coder achieves `128 kbps` for the full file entropy coder; `150 kbps` for the baseline coder achieves `96 kbps` for the full file entropy coder.
4. `--block_size`: specify the block size used by the audio coder (and hence, by the entropy coder in the per-block entropy coding setting); the default is 1024

An example command
```
python pacfile.py --input_file castanets.wav --entropy_coding block --data_rate 150 --block_size 1024
```
will output `castanets_150kbps_1024.wav` (which actually would be closer to around 128 kbps since we applied block entropy coding).