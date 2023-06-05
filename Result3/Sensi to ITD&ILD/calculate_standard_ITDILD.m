%% preparation
clear all
close all
clc

uiopen('C:\Users\12811\Desktop\BF+HRTF\AWTH_HRTF\HRTF_Data\HRTF_perFrequencyInterpolation.ita',1)
sr = HRTF_perFrequencyInterpolation.samplingRate;
time_period = 0.015;%calculate period 15ms
n_period = time_period * sr;

%% calc standard ITDILD
theta_step = pi/36;
phi_step = pi/36;
n_theta = 2*pi/theta_step;
n_phi = pi/phi_step;
ITD_database = zeros(1,(n_theta+1)*(n_phi+1));
ILD_database = zeros(1,(n_theta+1)*(n_phi+1));
for i = 0:1:n_theta%theta,0:pi/6:2pi
    for j = 0:1:n_phi%phi,0:pi/6:pi
        theta = theta_step*i;
        phi = phi_step*j;
        coord = itaCoordinates([1 phi theta],'sph');
        HRTF_interp = HRTF_perFrequencyInterpolation.interp(coord);
        hrir_l_sound = zeros(1,2*n_period);%add zeros
        hrir_r_sound = zeros(1,2*n_period);
        hrir_l_sound(1:length(HRTF_perFrequencyInterpolation.time(:,1))) = transpose(HRTF_interp.time(:,1));
        hrir_r_sound(1:length(HRTF_perFrequencyInterpolation.time(:,1))) = transpose(HRTF_interp.time(:,2));
        hrtf_l = fft(hrir_l_sound);
        hrtf_r = fft(hrir_r_sound);
        k_600 = round(600*2*n_period/sr)+1;
        k_700 = round(700*2*n_period/sr)+1;
        k_800 = round(800*2*n_period/sr)+1;
        freq = (k_600-1)*sr/(2*n_period);
        deltaangle = angle(hrtf_l(k_600))-angle(hrtf_r(k_600));
        if(deltaangle>pi)
            deltaangle = deltaangle - 2*pi;
        end
        if(deltaangle<=-pi)
            deltaangle = deltaangle + 2*pi;
        end
        ITD_600 = deltaangle/(2*pi*freq);
        if(freq == 0)
            ITD_600 = 0;
        end
        freq = (k_700-1)*sr/(2*n_period);
        deltaangle = angle(hrtf_l(k_700))-angle(hrtf_r(k_700));
        if(deltaangle>pi)
            deltaangle = deltaangle - 2*pi;
        end
        if(deltaangle<=-pi)
            deltaangle = deltaangle + 2*pi;
        end
        ITD_700 = deltaangle/(2*pi*freq);
        if(freq == 0)
            ITD_700 = 0;
        end
        freq = (k_800-1)*sr/(2*n_period);
        deltaangle = angle(hrtf_l(k_800))-angle(hrtf_r(k_800));
        if(deltaangle>pi)
            deltaangle = deltaangle - 2*pi;
        end
        if(deltaangle<=-pi)
            deltaangle = deltaangle + 2*pi;
        end
        ITD_800 = deltaangle/(2*pi*freq);
        if(freq == 0)
            ITD_800 = 0;
        end
        ITD_database(i*(n_phi+1)+j+1) = -median([ITD_600,ITD_700,ITD_800])*1000000;
        p_l = pwelch(hrir_l_sound,174,87);
        p_r = pwelch(hrir_r_sound,174,87);
        ILD_database(i*(n_phi+1)+j+1) = 10*log(sum(p_r)/sum(p_l));
    end
end
%% Save
save('basis','ITD_database','ILD_database');