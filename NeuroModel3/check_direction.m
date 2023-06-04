%% preparation
clear all
close all
clc

load('neuro_datasheet.mat')
%[audio,Fs] = audioread('C:\Users\12811\Desktop\BF+HRTF\Result\Accuracy vs Bw\1-20.wav');
[audio,Fs] = audioread('handel.wav');

time_period = 0.015;%calculate period 15ms
n_period = time_period * sr;
for theta = (0:5:355) / 180*pi
%% convolute stereo sound
uiopen('C:\Users\12811\Desktop\BF+HRTF\AWTH_HRTF\HRTF_Data\HRTF_perFrequencyInterpolation.ita',1)
% theta = 2*pi-pi/6;
% phi = pi/2;%check case
%theta = 75 /180*pi;
phi = 90 /180*pi;%experiment case
coord = itaCoordinates([1 phi theta],'sph');
HRTF_interp = HRTF_perFrequencyInterpolation.interp(coord);
hrir_l_sound = transpose(HRTF_interp.time(:,1));
hrir_r_sound = transpose(HRTF_interp.time(:,2));
%case AWTH
if(Fs~=sr)
    y_interp = transpose(interp1(1:1:length(audio),audio,(1:1:length(audio)/Fs*sr)/sr*Fs));
    y_interp(isnan(y_interp)) = 0;
else
    y_interp = audio;
end
y_sound_l = conv(y_interp,hrir_l_sound);
y_sound_r = conv(y_interp,hrir_r_sound);
%% find possibility
n_piece = floor(length(y_sound_l)/n_period);
possi_chart = zeros(n_piece,(n_theta+1)*(n_phi+1));
parfor j = 1:1:n_piece
    %extract piece
    piece_l = zeros(1,2*n_period);%add zeros
    piece_r = zeros(1,2*n_period);
    piece_l(1:n_period) = y_sound_l((j-1)*n_period+1:j*n_period);
    piece_r(1:n_period) = y_sound_r((j-1)*n_period+1:j*n_period);
    %calculate spike rate & possibility
    fft_l = fft(piece_l);
    fft_r = fft(piece_r);
    possi_freq = zeros(1,n_period);
    intensity = zeros(1,n_period);
    spike_rate = zeros(4,n_period);
    for k = 1:1:n_period
        freq = (k-1)*sr/(2*n_period);
        ILD = 10*log10(abs(fft_r(k))/abs(fft_l(k)));
        intensity(k) = abs(fft_r(k))*abs(fft_l(k));
        deltaangle = angle(fft_l(k))-angle(fft_r(k));
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
        spike_rate(:,k) = [LSO_r;MSO_r;MSO_l;LSO_l];
    end
    local_possi_chart = zeros(1,(n_theta+1)*(n_phi+1));
    for i = 1:1:(n_theta+1)*(n_phi+1)
        for k = 1:1:n_period
            target_spike_rate = neuro_database(:,k,i);
            possi_freq(k)=transpose(spike_rate(:,k))*target_spike_rate/(norm(spike_rate(:,k))*norm(target_spike_rate));
        end
        local_possi_chart(i) = sum(intensity/sum(intensity).*possi_freq);
    end
    possi_chart(j,:) = local_possi_chart;
end
%% Draw figure
x = 1:1:n_piece;
y = 1:1:(n_theta+1)*(n_phi+1);
[X,Y]=meshgrid(x,y);
mesh(X,Y,transpose(abs(possi_chart)));
title("Direction possibility (sound direction: theta = 11pi/6,phi=pi/2:channel:81)");
ylabel('channel');
xlabel('piece');
zlabel('possibility');
clbar = colorbar;
clbar.Label.String = 'Possibility';
%% Draw sphere
figure(2)
[X,Y,Z]=sphere(72);
X = X(1:2:end,:);
Y = Y(1:2:end,:);
Z = Z(1:2:end,:);
ave_possi = mean(abs(possi_chart));
C = reshape(ave_possi,[37 73]);
surf(X,Y,Z,C);
title("Direction possibility (sound direction: theta = 11pi/6,phi=pi/2,(-1,0,0) corresponds to (0,pi/2),(0,0,-1)to(0,0)); surface represents possibility on downleft point");
xlabel('x');
ylabel('y');
zlabel('z');
clbar = colorbar;
clbar.Label.String = 'Possibility';
%% Draw detailed sphere
figure(3)
[X,Y,Z]=sphere(72);
X = X(1:2:end,:);
Y = Y(1:2:end,:);
Z = Z(1:2:end,:);
ave_possi = mean(abs(possi_chart));
ave_possi_high = (ave_possi - 0.9)*10;
ave_possi_high(ave_possi_high<0)=0;
C = reshape(ave_possi_high,[37 73]);
surf(X,Y,Z,C);
title("Direction possibility (sound direction: theta = 11pi/6,phi=pi/2,(-1,0,0) corresponds to (0,pi/2),(0,0,-1)to(0,0)); surface represents possibility on downleft point");
xlabel('x');
ylabel('y');
zlabel('z');
clbar = colorbar;
clbar.Label.String = 'Possibility';
%% Report
[~,max_channel] = max(ave_possi);
max_theta_i = ceil(max_channel/37)-1;
max_phi_i = max_channel - max_theta_i * 37 - 1;
fprintf('point with highest accuracy:theta:%i degree,phi:%i degree\n',max_theta_i*5,max_phi_i*5);

fileID = fopen('output.txt','a');
fprintf(fileID,'%d %d\n',max_theta_i*5,max_phi_i*5);
fclose(fileID);
end