% -----------------------------------------------------------------------------
%  File          : SNR_test.m
%  Author        : Chengcheng Xu and A. Young
%  Created       : 02/20/2017
%  Project       : LCLS II Stripline BPM
% -----------------------------------------------------------------------------
%  Description :
%  Test the SNR of the ADC
% -----------------------------------------------------------------------------
%  Note :
%  2/20/17:  designed for GbE platform.  
%            single stream data format
% ./launch.sh attnsweep_test.py -b1 -s512 -n1 -d /data -Y yaml/current/*_project.yaml/000TopLevel.yaml
out = [];
Nfig=0;

for i = 4:7
% Adding file directory path
Data_directory = '/afs/slac/g/lcls/users/BPM/LCLS_II/Data/20231114-172104/';  % Adding the fold name

dir_name  = datestr(now,'mm_dd_yyyy_HHMMSS');
directory = Data_directory
% directory = ['/afs/slac/g/lcls/users/BPM/LCLS_II/Data/BPM_Res_' dir_nam20211018-103638/e '/'];
%mat_fname = [Data_directory,'BPM_Res_data_' dir_name '.mat'];
%
ctrl = struct; % sets up various parameters that seem to be needed for decoding
ctrl.byte_reverse = 1;
ADC.index=i; %0==chan0,1==chan1, 2=chan2, 3==chan3, bay1==chan1==46


% ch1     ch2     ch3     ch4
% [3]     [2]     [1]     [0]      - Two slot crate
% P7      P6      P5      P4       - SNR 
% stream7 stream6 stream5 stream4  - Attnuation Sweep
% green   blue    yellow  red      - Color of channel
% +x      -x      -y      +y       - Face of AMC


ADC.bay=0;
ch1.fname = [directory, 'Stream', int2str(ADC.index), '_', int2str(ADC.bay), '.bin'];
[ch1.header, ch1.nFrame, ch1.data_raw]  = bin2mat_single_GbE(ch1.fname, ctrl);
% Got the Data
ch1.min = min(ch1.data_raw); 
ch1.max = max(ch1.data_raw); 
ch1.amp = ch1.max-ch1.min;
ch1.mean = mean(ch1.data_raw);

fs = (370)*1e6; % fs in Hertz
Struck_250_fs = 250e6;
ADC_FS_Vpp    = 1.7; % ADC16DX370 fullscale input voltage 
ADC_impedance = 100;
ADC_FS_Vp     = ADC_FS_Vpp/2;
ADC_FS_Vrms   = 0.707*ADC_FS_Vp;
ADC_FS_Pmw    = ((1e3)*(ADC_FS_Vrms^2))/(ADC_impedance);
ADC_FS_dBm    = 10*log10(ADC_FS_Pmw)
bitnum = 16; %16 bit ADC
L=length(ch1.data_raw)

%
% Plot 1
Nfig=Nfig+1;
% figure(Nfig);
subplot(4,3,Nfig);
plot(ch1.data_raw(1:250))
title('Time domain plot dot only');
xlabel('Samples');
ylabel('ADC counts');

%
% Plot 2
Nfig=Nfig+1;
% figure(Nfig);
subplot(4,3,Nfig);
plot(ch1.data_raw(1:250)-ch1.mean)
title('Time domain plot');
xlabel('Samples');
ylabel('ADC counts');
mean_data=ch1.data_raw-ch1.mean;

mean_data_max = max(mean_data)
mean_data_min = min(mean_data)
mean_data_amp = mean_data_max-mean_data_min
dBFS_offset = 20*log10(mean_data_max/(2^(bitnum-1)));

NFFT = 2^nextpow2(L+2); 
fullscale = 2^(bitnum-1);

% Blackman window
cc = 1/0.42;          % correction coefficient
w  = (0.42 - .5*(cos(2*pi*(1:L)'/(L+1)))+(0.16/2)*(cos(2*pi*(1:L)'/(L+1)))); % hardcoded blackman in case you do not signal proc tool box
% w = blackman(L);  % window function (blackman)
% end of Blackman

data_fullscale = mean_data ./ fullscale;
f = fs/2 * linspace(0,1,NFFT/2);

rmsdata_fullscale = data_fullscale - mean(data_fullscale);

Y           = fft(squeeze(rmsdata_fullscale(:,1)) .* w, NFFT) * cc / L;	
amp_spec    = 2 * abs(Y(1:NFFT/2));
amp_spec_dB = 20*log10(amp_spec); % in dBFS

%
% Plot 3
Nfig=Nfig+1;
subplot(4,3, Nfig);
% figure(Nfig);
% plot((f*1e-6), squeeze(amp_spec_dB)); grid on
% plot((f*1e-6), squeeze(amp_spec_dB),'*'); grid on
plot(f*1e-6, amp_spec_dB); grid on % frequency domain only going up to 125MHz
title(['Frequency domain plot, ADC fullscale is ' num2str(ADC_FS_dBm) 'dBm']);
xlabel('Frequency (MHz)');
ylabel('dBFS');

i=1;
fft_gain = 10*log10(L/2);
%fft_gain = 10*log10(f(1)/2); % Calculating the fft gain only going to the 125MHz bandwidth
sig_max = max(amp_spec(10:end));
sig_max_ind=find(amp_spec(:) <= sig_max & (amp_spec(:) >=sig_max));
noise_array=cat(2, amp_spec(1:sig_max_ind-20)', amp_spec(sig_max_ind+20:end)');

noise_sum = sum(noise_array);
noise_sum_dB = 20*log10(noise_sum)
noise_std = std(noise_array);
noise_std_dB = 20*log10(noise_std)

sig_samp = 2; % Take the +/- x sample points from the signal ind
sig_dBm_dist = zeros((sig_samp*2)+1,1)
sig_mw_dist = zeros((sig_samp*2)+1,1)
for i=1:(sig_samp*2)+1
    sig_dBm_dist(i) = ADC_FS_dBm+amp_spec_dB((sig_max_ind-sig_samp)+i-1)
    sig_mw_dist(i)  = dBm2mW(sig_dBm_dist(i));
end
sig_power = mW2dBm(sum(sig_mw_dist))


sig_max_dBFS = ADC_FS_dBm + sig_power;
SNR = sig_max_dBFS - noise_std_dB - fft_gain
out = [out sig_power SNR]
end
