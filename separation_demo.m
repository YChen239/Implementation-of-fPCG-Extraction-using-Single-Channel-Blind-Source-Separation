% This script is used to load and visualize the real data.
%
% Author: Edgar Lobaton
% Created: 1/22/2019

% Loading audio file
[y0,Fs] = audioread('../Physionet-Real-sufhsdb/f100.wav');

% Selecing a channel since both channels seem to be the same
y = y0(:,1);
t = (1:length(y))/Fs;

pl = audioplayer(y,Fs);
play(pl); pause(10); stop(pl);

%% ANALYSIS

% Produce spectogram
winSz = round(2*Fs);
noverlap = round(0.5*winSz);
nfft = [];
figure(1); h(1) = subplot(2,1,1);
spectrogram(y,winSz,noverlap,nfft,Fs,'yaxis');%Short-time Forier Transform

% Filter based on appropriate band from papers and observed in spectogram
Fc = 500;
[b,a] = butter(6,Fc/(Fs/2));
%figure(2), freqz(b,a,4*2048,Fs);
yp = filter(b,a,y);

% Plotting results
figure(1);  h(2) = subplot(2,1,2);
plot(t,y,t,yp);
legend('Original','Filtered');


[imf,residual,info] = emd(yp,'Interpolation','pchip');

% STFT to each channel
[l , n] = size(imf);
Py = [];
C = [];
for i = 1:n
    Pyn = abs(stft(imf(:,i), 3000)).^2;
    Py = cat(3,Py,Pyn);
    [W,H] = sparse_nmf(Pyn,2); %sparse NMF with k=2 2clusters
    C = [C;istft(sqrt(W(:,1)*H(1,:)), 3000)'];
    C = [C;istft(sqrt(W(:,2)*H(2,:)), 3000)'];
end

%Clustering
Z=abs(C).^2./l;
CZ = Z*Z';
image(CZ,'CDataMapping','scaled')
colorbar

C=real(C);
output = (C(1,:)+C(2,:)+C(3,:)).*10;
%% RESULTS

% Playing new filtered signal
pl = audioplayer(output,Fs);
play(pl); pause(10); stop(pl);
plot(t(1:end-10),output(1,1:end-10),t(1:end-10),yp(1:end-10,1));

audiowrite('f10_output.wav',output,Fs)

subplot(311),plot(t(1:end-10),y(1:end-10,1)),title('Original Signal');
subplot(312),plot(t(1:end-10),yp(1:end-10,1)),title('Lowpass Signal');
subplot(313),plot(t(1:end-10),output(1,1:end-10)),title('Final Signal');

% Look into for segmentation of waveform:
% https://physionet.org/physiotools/hss/