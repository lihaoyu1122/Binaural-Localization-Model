clear all
close all
clc

sr = 48000;
noise = 10*transpose(0.01 * randn(1,10*sr));
fft_n = fft(noise);
low_cutoff_freq = 125;
high_cutoff_freq = 500;
fft_n(1:round(low_cutoff_freq*10))=0;
fft_n(round(high_cutoff_freq*10):10*sr+1-round(high_cutoff_freq*10))=0;
fft_n(10*sr+1-round(low_cutoff_freq*10):10*sr)=0;
noise = real(ifft(fft_n));

audiowrite("noise.wav",noise,sr);