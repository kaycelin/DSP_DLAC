%% 2020-03-25, Add IpwrdBPerRBW_Noise
% %% 2020-03-25, Add Noise DCOffset Cancellation --> Removed
%% 2020-03-26, There is NO Negative Voltage Ampliture of NOISE
%% 2020-03-26, IpwrdB_NoiseFullBand: the total Noise for the FullBand
%% 2020-03-29, ISSUE, flag_NFdB: [complexdata(wf)],'SNR','Noise' ???
%% 2021-03-15, Change Noise from real to complex, 'compx1' is correct ???
%% 2021-03-15, scale_complex = 1/2 ???

function [waveform_awgn, noise] = trx_SNR_g(waveform, SNRdB_FullBand, PwrdB_FullBand_Noise, NFdB, flag_NFdB, weight_Noise, flag_Debug_AddNoise)

% ROW
DIM_FFT = 2;
if size(waveform,1)>size(waveform,2) % COLUMN
    waveform=waveform.'; % switch to ROW
    flag_wf_original = 'COLUMN';
else
    flag_wf_original = 'ROW';
end

if ~exist('fs','var')||~exist('fc','var')||~exist('bw_Carrier','var')
    flag_cal_pwr_method = 'td';
end

Nsamps = length(waveform);

if ~exist('NFdB','var')||isempty(NFdB)
    NFdB = 0;
elseif length(flag_NFdB)==length(waveform) % Ref. waveform
    [EVM, ~] = dsp_evm_timexcorr_inband_g(flag_NFdB, circshift(waveform,0), [], 0); % EVM Full Band, ISSUE ?
    SNRdBin = - [0 + 20*log10(EVM/100)];
    SNRdB_FullBand = SNRdBin-NFdB;
elseif strcmp(flag_NFdB,'SNR')&~isempty(PwrdB_FullBand_Noise)
    SNRdB_FullBand = SNRdB_FullBand-NFdB;
elseif strcmp(flag_NFdB, 'Noise')&~isempty(PwrdB_FullBand_Noise)
    PwrdB_FullBand_Noise = PwrdB_FullBand_Noise+NFdB;    
else
    error('Check Noise input!')
end

%% 2020-03-29, flag_NFdB: [complexdata(wf)],'SNR','Noise'
if exist('SNRdB_FullBand','var')&~isempty(SNRdB_FullBand)
    SNR = 10^(SNRdB_FullBand/10);
    
    switch flag_cal_pwr_method
        case {'td'} % calculate pwr of txwaveform by time domain
            Pwr_waveform = sum(abs(waveform).^2, DIM_FFT)/Nsamps; % time domain
            [Pwr_dB_waveform] = Pwr_Inband_g(fft(waveform, Nsamps, DIM_FFT), [], [], 0.01e6, 'half', 0); % freq domain
            
        case {'FD'} % calculate pwr of txwaveform by freq. domain inband
            [Pwr_dB_waveform] = Pwr_Inband_g(fft(waveform, Nsamps, DIM_FFT), fs, fc, bw_Carrier, 0.01e6, 'half', 0); % freq domain
            [Pwr_dB_waveform,~,Opwr_dB_waveform] = Pwr_Inband_g(fft(waveform, Nsamps, DIM_FFT), fs, fc+[-bw_Carrier bw_Carrier]/2, 0.01e6, 'half', 0); % freq domain
            
            Pwr_waveform = 10^(Pwr_dB_waveform/10);
    end
    Pwr_Noise = Pwr_waveform/SNR;
%     noise = sqrt(Pwr_Noise)* randn(1,Nsamps);
elseif exist('PwrdB_FullBand_Noise','var')&~isempty(PwrdB_FullBand_Noise)
    Pwr_Noise = 10^(PwrdB_FullBand_Noise./10);
%     noise = sqrt(Pwr_Noise/1)* randn(1,Nsamps);
end

%% 2021-03-15, scale_complex = 1/2 ???
scale_complex = 1/2;
% scale_complex = 1;

if ~exist('weight_Noise','var')||isempty(weight_Noise)
%     weight_Noise = 1e-2;
%     weight_Noise = [-1 1] ;
    noiseI = sqrt(Pwr_Noise/scale_complex)*(2*rand(1,Nsamps)-1);
    noiseQ = sqrt(Pwr_Noise/scale_complex)*(2*rand(1,Nsamps)-1);

%     noise = sqrt(Pwr_Noise/1)*(1*rand(1,Nsamps)-0);
% noise=sqrt(Pwr_Noise/1)*ones(1,Nsamps);

% elseif strcmp(weight_Noise,'Debug')
%     noiseI = sqrt(Pwr_Noise/scale_complex)*(2*ones(1,Nsamps)-1);
%     noiseQ = sqrt(Pwr_Noise/scale_complex)*(2*ones(1,Nsamps)-1);
else
    noiseI = sqrt(Pwr_Noise/scale_complex)*(randi([0 1],[1,Nsamps])*2-1).*(2*weight_Noise*rand(1,Nsamps)-weight_Noise+1);
    noiseQ = sqrt(Pwr_Noise/scale_complex)*(randi([0 1],[1,Nsamps])*2-1).*(2*weight_Noise*rand(1,Nsamps)-weight_Noise+1);

%     (randi([0 1],[1,Nsamps])*2-1).*
%     random_index = (2*weight_Noise*rand(1,Nsamps)-weight_Noise);
%     random_index_POS = find(random_index)
%     random_index_NEG = random_index-1;
    
end
%% 2021-03-15, Change Noise from real to complex, 'compx1' is correct ???
if  ~exist('flag_Debug_AddNoise','var')||isempty(flag_Debug_AddNoise)||flag_Debug_AddNoise==0
    flag_Debug_AddNoise = 'compx1';
end
switch flag_Debug_AddNoise
    case {'compx1'}
        noise = noiseI + 1i*noiseQ;
    case {'compx2'}
        noise = 0.5*noiseI + 0.5*1i*noiseQ;
    case {'real'}
        noise = noiseI;
end
% noise = sqrt(Pwr_Noise/1)*randn(1,Nsamps);
% noise = sqrt(Pwr_Noise/1)*(weight_Noise*rand(1,Nsamps)-weight_Noise/2);
% noise = sqrt(Pwr_Noise/1)*(diff(weight_Noise)*rand(1,Nsamps)+weight_Noise(1));

% %% Add Noise DCOffset Cancellation
% NOISE = fft(noise);
% NOISE(1) = 0;
% noise_DCOff = ifft(NOISE);
% waveform_awgn = waveform+noise_DCOff;
waveform_awgn = waveform+noise;

% export and recover flag_wf_original
if strcmp(flag_wf_original,'COLUMN')
    waveform_awgn = waveform_awgn.'; % back to column typ
    noise = noise.';
end
%% check
% % figure(20200326)
% % PLOT_FFT_dB_g2(waveform_awgn, 983040000, Nsamps, ['txwf noise'], 'df', 0);
% [IpwrdB_wf, ~, ~] = Pwr_Inband_g(fft(waveform,[],DIM_FFT), 3.3178e+09, 2e9+[0 160000], 0, 'full', 0);
% [IpwrdB_N, ~, ~] = Pwr_Inband_g(fft(noise,[],DIM_FFT), 3.3178e+09, 2e9+[0 160000], 0, 'full', 0);
% [FpwrdB_N, ~, ~] = Pwr_Inband_g(fft(noise,[],DIM_FFT), [], 2e9+[0 160000], 0, 'full', 0);
% 
% % [IpwrdB_wfN, ~, ~] = Pwr_Inband_g(fft(waveform_awgn), 3.3178e+09, 2e9+[0
% % 160000], 0, 'full', 0);
% [IpwrdB_wfN, ~, ~] = Pwr_Inband_g(fft(waveform_awgn,[],DIM_FFT), 3.3178e+09, 0+[0 3.3178e+09], 0, 'full', 0);
end