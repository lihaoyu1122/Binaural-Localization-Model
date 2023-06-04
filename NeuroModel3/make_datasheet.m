%% preparation
clear all
close all
clc

uiopen('C:\Users\12811\Desktop\BF+HRTF\AWTH_HRTF\HRTF_Data\HRTF_perFrequencyInterpolation.ita',1)
sr = HRTF_perFrequencyInterpolation.samplingRate;
time_period = 0.015;%calculate period 15ms
n_period = time_period * sr;
%% compiling database
theta_step = pi/36;
phi_step = pi/36;
n_theta = 2*pi/theta_step;
n_phi = pi/phi_step;
neuro_database = zeros(4,n_period,(n_theta+1)*(n_phi+1));
hrir_l = zeros(1,2*n_period);
hrir_r = zeros(1,2*n_period);
hrtf_l = zeros(1,2*n_period);
hrtf_r = zeros(1,2*n_period);
for i = 0:1:n_theta%theta,0:pi/6:2pi
    for j = 0:1:n_phi%phi,0:pi/6:pi
        theta = theta_step*i;
        phi = phi_step*j;
        coord = itaCoordinates([1 phi theta],'sph');
        HRTF_interp = HRTF_perFrequencyInterpolation.interp(coord);
        hrir_l(1:length(HRTF_perFrequencyInterpolation.time(:,1))) = transpose(HRTF_interp.time(:,1));
        hrir_r(1:length(HRTF_perFrequencyInterpolation.time(:,1))) = transpose(HRTF_interp.time(:,2));
        hrtf_l = fft(hrir_l);
        hrtf_r = fft(hrir_r);
        for k = 1:1:n_period
            freq = (k-1)*sr/(2*n_period);
            ILD = 10*log10(abs(hrtf_r(k))/abs(hrtf_l(k)));
            deltaangle = angle(hrtf_l(k))-angle(hrtf_r(k));
            if(deltaangle>pi)
                deltaangle = deltaangle - 2*pi;
            end
            if(deltaangle<=-pi)
                deltaangle = deltaangle + 2*pi;
            end
            ITD = deltaangle/(2*pi*freq);
            if(freq == 0)
                ITD = 0;
            end
            LSO_r = 2/(1+exp(ILD/5));
            LSO_l = 2/(1+exp(-ILD/5));
            MSO_r = 2^(sin(4*pi*freq*-ITD))*1.5/(1+exp((freq-5000)/500));
            MSO_l = 2^(sin(4*pi*freq*ITD))*1.5/(1+exp((freq-5000)/500));
            neuro_database(:,k,i*(n_phi+1)+j+1)=[LSO_r;MSO_r;MSO_l;LSO_l];
        end
    end
end
%% save file
save('neuro_datasheet','theta_step','phi_step','n_theta','n_phi','neuro_database','time_period','n_period','sr');