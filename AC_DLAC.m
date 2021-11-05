%% PA1, 2021-11-03,
%% PA2, 2021-11-05, issue, why add -6dB ? 
%% PA3, 2021-11-05, AWGN

clear all
close all

fnum = 110301
%% input: signal type
flag_sigType = 'DLAC' % AC/DLAC
flag_fs_extend = 1; % signal fs extend

if strcmpi(flag_sigType, 'DLAC')
    %% input: DL signal
    DLconfig.bw_Channel = '20MHz'
    DLconfig.fs = 122.88e6; % flag_fs_extend
    
    %% output: DL signal
    [sigDL, sigGrid, ConfigDL] = OFDM_Mod_g(DLconfig,'64QAM','off','NR',fnum,[]);
end

if strcmpi(flag_sigType, 'DLAC')||strcmpi(flag_sigType, 'AC')
    %% input: AC signal
    if flag_fs_extend
        fsAC = [23.04e6, 122.88e6]
    else
        fsAC = [23.04e6]*ones(1,2)
    end
    df = 100
    Nsamps = fsAC(2)/df
    Nbr = 4
    NsampsAC = [Nsamps, Nbr]
    dACSubcarrier = 240e3;
    bwChannelAC = {'20MHz',fsAC,dACSubcarrier}
    
    %% output: AC signal
    [sigAC, ConfigAC] = AntCal_genACSource_g100(bwChannelAC, [], [NsampsAC], [], [], fnum+1 ,[]);
end

%% output: carrier configuration
if strcmpi(flag_sigType, 'DLAC')
    if ConfigDL.fs ~= ConfigAC.fs
        error('signal fs should be the same!')
    end
end   
fs = ConfigAC.fs
bwInband = (ConfigAC.bwCarrier/2+0.5e6)*[-1 1]
bwChannel = ConfigAC.bwChannel
Nsamps = size(sigAC,1)
Nbr = size(sigAC,2)
df = fs/Nsamps


%% input: channel filter
FIRc1_Wtype = "kaiser"
FIRc1_Ftype = "LPF"
FIRc1_Order = NaN
FIRc1_fTolerance = -0.1e6
FIRc1_K_AttdB = 60
FIRc1_K_fdelta = 0.5e6
FIRc1_fcutoffL = bwChannel/2
FIRc1_fcutoffH = 0
FIRc1_Export = fs
NCarriers = 1

%% output: channel filter
b_ch = SYM_FIRApp(FIRc1_Wtype,FIRc1_Ftype,FIRc1_Order,FIRc1_K_AttdB,FIRc1_K_fdelta,FIRc1_fTolerance,FIRc1_fcutoffL,FIRc1_fcutoffH,df,bwInband,NCarriers,FIRc1_Export,fnum+2);
b_ch = b_ch{:};

%% output: signal+channel filter
if strcmpi(flag_sigType, 'DLAC')
    sigDLch = conv(sigDL(:,1), b_ch, 'same');
    PLOT_FFT_dB_g(sigDLch, fs, Nsamps, ['sigDLch'], 'df', 'full', 'pwr', [fnum+3, 1,2,1]);
end

for ibr = 1:Nbr
    sigACch(:,ibr) = conv(sigAC(:,ibr), b_ch, 'same');
end
PLOT_FFT_dB_g(sigACch, fs, Nsamps, ['sigACch'], 'df', 'full', 'pwr', [fnum+3, 1,2,2]);

%% Digital gain or output power assignment
if strcmpi(flag_sigType, 'DLAC')
    %% input: DLAC output power
    PodB_DL = -15
    
    %% output: signal DL+gain
    vo_DL = 10.^((PodB_DL)/20);
    vi_DL = ( sqrt(mean(abs(sigDLch).^2)) );
    sigDLchGain = sigDLch.*vo_DL/vi_DL;
    PLOT_FFT_dB_g(sigDLchGain, fs, Nsamps, ['sigDLchGain, output pwrdB:',num2str(PodB_DL)], 'df', 'full', 'pwr', [fnum+4, 1,2,1]); 
    PdBm_DL = 10*log10(mean(abs(sigDLchGain(:,1)).^2))
end

%% input: AC output power ??? 
%% 2021-11-05, issue, why add -6dB ? 
PodB_AC = -15+10*log10(numel(ConfigAC.nSubcarrier)/ConfigDL.NScs)-6 %% 2021-10-27, Equal subcarrier for DLAC

%% output: signal AC+gain
vo_AC = 10.^((PodB_AC)/20);
vi_AC = ( sqrt(mean(abs(sigACch(:,1)).^2)) );
sigACchGain = sigACch.*vo_AC/vi_AC;
PLOT_FFT_dB_g(sigACchGain, fs, Nsamps, ['sigACchGain, output pwrdB:',num2str(PodB_AC)], 'df', 'full', 'pwr', [fnum+4, 1,2,2]);
PdBm_AC = 10*log10(mean(abs(sigACchGain(:,1)).^2))

%% NCO and Carrier combination
if strcmpi(flag_sigType, 'DLAC')
    %% input: fNCO for [DL1 AC DL2]
    fNCO = [-20e6 0 20e6]
    flag_SUMCarrier = 1
    disp_title = 'DL+AC+DL'
    idAC=2
    idDL=[1 3]

    %% output: DLAC signal
    if any(fNCO > fs/2)
        error('check fNCO range!')
    end
    if 0
        for ibr=1:Nbr
            sigDLchGain1(:,ibr) = circshift(sigDLchGain,ibr-1);
        end
    else
        sigDLchGain1 = sigDLchGain;
    end
    sigDLAC_cell = [{sigDLchGain1}; {sigACchGain}; {sigDLchGain1}]; % row
    [~, sigDLACoutput_cell, ~, ~, ~] = SYM_NCOApp(sigDLAC_cell, fNCO, fs, bwInband, flag_SUMCarrier, fnum+5, disp_title, 0);
    sigACout = cell2mat(sigDLACoutput_cell);
else
    idAC=1
    sigACout = sigACchGain;
end

%% Export: DLAC waveform
sigACout;
if 0
    save(['waveform_',flag_sigType,'_3x20MHz_122p88MHz_4Br.mat'],'sigACout')    
end
PLOT_FFT_dB_g(sigACout, fs, Nsamps, ['siganl ',flag_sigType], 'df', 'full', 'pwr', [fnum+50]);

%% input: Branch phase shifter
phsShiftNbrDeg = [0 10 -7 2]

%% output: signal+phaseShifter
sigACphsShift = sigACout.*exp(1i*phsShiftNbrDeg./180*pi);

%% output: Antenna output for Transmitter
sigAC_TXout = sigACphsShift;

%% output: Antenna combination for all branch combination to one branch
sigAC_AntComb= sum(sigAC_TXout,2);

%% 2021-11-05, AWGN
flag_AWGN = 1
if flag_AWGN
    sigAC_AntComb = sum(sigAC_TXout,2);
    
    %% input: AWGN
    SNRdB_FullBand = 50
    [sigAC_AntComb_AWGN, noise] = trx_SNR_g(sigAC_AntComb, SNRdB_FullBand, [], [], 'SNR', 0, 0);
    PLOT_FFT_dB_g(sigAC_AntComb_AWGN, fs, Nsamps, ['waveforme+AWGN SNRdB:',num2str(SNRdB_FullBand)], 'df', 'full', 'pwr', [fnum+70]);

    % export
    sigAC_AntComb = sigAC_AntComb_AWGN;
end

%% output: Antenna combiner input for Receiver
sigAC_RxIn = sigAC_AntComb;

%% Downcoversion by LO

%% NCO shift to Zero and multicarrier divide
if exist('fNCO','var')&&~isempty(fNCO)
    fNCO2Zero = -fNCO;
    flag_NCO2Zero_MultiCarrierDivide = 'DIV'
    idC=1
    [sigAC_RxInNCO2Zero, ~, ~, ~, table_NCO_AC] = SYM_NCOApp(sigAC_RxIn, fix(fNCO2Zero/df)*df, fs, [], flag_NCO2Zero_MultiCarrierDivide, fnum+6, [], []);
else
    sigAC_RxInNCO2Zero = sigAC_RxIn;
end
    
%% DLAC channel filter for demodulation
if strcmpi(flag_sigType, 'DLAC')
    for idC = 1:size(sigAC_RxInNCO2Zero,2)
        sigAC_RxInForDmod(:,idC) = conv(sigAC_RxInNCO2Zero(:,idC), b_ch, 'same');
    end
else
    sigAC_RxInForDmod = sigAC_RxIn;
end
[~,~,plt]=PLOT_FFT_dB_g(sigAC_RxInForDmod, fs, Nsamps, ['sigACout For Demodulation'], 'df', 'full', 'pwr', [fnum+7]);
title('DL1(blue), AC(orange), DL2(yellow)')
%% AC Demodulation
ACDomod_Nshift = 0;
[ACdmod_t0Mean, ACdmod_p0DegMean, ACdmod_phEstDeg, ACdmod_phEstDegDrift, ACdmod_dataCapCor, ACdmod_SNR, ACdmod_dataCapNbrwoPD]...
    = AntCal_phaseDemodulateApp_g100(sigAC_RxInForDmod(:,idAC),ConfigAC,ACDomod_Nshift,fnum+8, []);

%% AC Calculate EQ coefficents and Phase Alignment
%% input: AC EQ number of taps
ACEQ_numFirTaps = 2
ACEQ_NoOfAlignBr = 1

if ConfigAC.Nfft_org ~= ConfigAC.Nfft
    ConfigAC_EQ = ConfigAC;
    ConfigAC_EQ.Nfft = ConfigAC_EQ.Nfft_org; % replace Nfft configuration for EQ solve
    if 0
        ConfigAC_EQ.Nfft = fsAC(1)/ConfigAC_EQ.dfScs
    end
end

%% output: AC signal+phaseAlignment(sigAC_PhsAlign), phaseAlignment EQ filter(ACEQ_taps)
[sigAC_PhsAlign,ACEQ_taps,ACEQ_B] = AntCal_EQphsAlignApp(sigAC,ConfigAC_EQ,[], ACdmod_phEstDeg,ACEQ_numFirTaps,ACEQ_NoOfAlignBr,fnum+8, []);

%% DL(OFDM) demodulation:
if strcmpi(flag_sigType, 'DLAC')
    for idC = idDL
        rxGrid(:,:,idC) = OFDM_Demod_g(ConfigDL, sigAC_RxInForDmod(:,idC)*vi_DL/vo_DL);
        PLOT_Constellation(rxGrid(:,:,idC)/1,[],fnum+8+idC,ConfigDL.MOD,0)
        [evm_DL] = dsp_evm_g(rxGrid(:,:,idC), sigGrid)
    end
end

%% Analysis
flag_analysis = 1
% k=1
if flag_analysis
    x(k) = SNRdB_FullBand
    y1(k,:) = ACdmod_SNR
    y2(k,:) = ACdmod_phEstDegDrift
    y3(k,:) = ACdmod_p0DegMean-phsShiftNbrDeg
    y4(k) = evm_DL
    k=k+1
end
if flag_analysis
    for kk=1:numel(x)
        figure(fnum+90)
        plot(x, y1(kk,:), 'DisplayName', ['SNRdBout:',num2str(round(y1(kk,:),2))] ), hold on, legend
        title(['SNRin vs SNRout vs Branches'])
        xlabel('SNRdBin'), ylabel('SNRdBout')
        
        figure(fnum+91)
        plot(x, y2(kk,:), 'DisplayName', ['phaseDriftdeg:',num2str(round(y2(kk,:),2))] ), hold on, legend
        title(['SNRin vs phaseDriftdeg vs Branches'])
        xlabel('SNRdBin'), ylabel('phaseDriftdeg')

        figure(fnum+92)
        plot(x, y3(kk,:), 'DisplayName', ['mean(phaseDriftdeg) compare to RefPhaseShift:',num2str(round(y3(kk,:),2))] ), hold on, legend
        title(['SNRin vs phaseShifterForEQ vs Branches'])
        xlabel('SNRdBin'), ylabel('mean(phaseDriftdeg) compare to RefPhaseShift')

    end
    
    figure(fnum+93)
    plot(x, y4, 'DisplayName', ['EVM(DL):',num2str(round(y4,2))] ), hold on, legend
    title(['SNRin vs EVM(DL) vs Branches'])
    xlabel('SNRdBin'), ylabel('EVM')
end
