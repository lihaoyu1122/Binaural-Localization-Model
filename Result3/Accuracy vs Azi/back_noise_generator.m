clear all
close all
clc

sr = 48000;
for SNR = -18:1:-8
Amp = 10/(19*10^(SNR/20));
noise = Amp*transpose(0.01 * randn(1,10*sr));

rendered_noise_l = zeros(480361,1);
rendered_noise_r = zeros(480361,1);
for theta = ([0:15:135 225:15:345]) / 180*pi
    uiopen('C:\Users\12811\Desktop\BF+HRTF\AWTH_HRTF\HRTF_Data\HRTF_perFrequencyInterpolation.ita',1);
    phi = 90 /180*pi;
    coord = itaCoordinates([1 phi theta],'sph');
    HRTF_interp = HRTF_perFrequencyInterpolation.interp(coord);
    hrir_l_sound = transpose(HRTF_interp.time(:,1));
    hrir_r_sound = transpose(HRTF_interp.time(:,2));
    specific_noise_l = conv(noise,hrir_l_sound);
    specific_noise_r = conv(noise,hrir_r_sound);
    rendered_noise_l = rendered_noise_l + specific_noise_l;
    rendered_noise_r = rendered_noise_r + specific_noise_r;
end

audiowrite(strcat("backnoise_snr",num2str(SNR),".wav"),[rendered_noise_l rendered_noise_r],sr);
end