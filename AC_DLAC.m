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






stop=1

    PdBm_DLin = 10*log10(mean(abs(sigDL(:,1)).^2))
    PdBm_DLo = 10*log10(mean(abs(sigAC_RxInForDmod(:,1)).^2))
    PdBm_DLo_n = 10*log10(mean(abs(sigAC_RxInForDmod(:,1)*vi_DL/vo_DL).^2))
    PdBm_DLo_n = 10*log10(mean(abs(sigAC_RxIn(:,1)*1/1).^2))

    sigAC_RxIn


DLconfig.bw_Channel = '20MHz'
DLconfig.fs = 122.88e6;
[sigDL2, ~, Config2] = OFDM_Mod_g(DLconfig,'64QAM','off','NR',fnum,[]);
rxGrid2 = OFDM_Demod_g(ConfigDL, sigDL2);
PLOT_Constellation(rxGrid2,[],103002,[])

PdBm_DLo = 10*log10(mean(abs(sigDL2(:,1)).^2))




ConfigAC2 = ConfigAC
ConfigAC2.Nfft = 512
[ACdmod_t0Mean, ACdmod_p0DegMean, ACdmod_phEstDeg, ACdmod_phEstDegDrift, ACdmod_dataCapCor, ACdmod_SNR, ACdmod_dataCapNbrwoPD]...
    = AntCal_phaseDemodulateApp_g100(sigAC_RxInForDmod,ConfigAC2,ACDomod_Nshift,fnum+7, []);








load('waveform_AC_20MHz_122p88MHz_4Brs.mat')
load('waveform_AC_20MHz_122p88MHz_4Brs_cell.mat')
load('waveform_DLAC_3x20MHz_122p88MHz.mat')
ACconfig = sigAC_cell{2}
sigAC = sigAC_cell{1};
sigOut_DLAC = sigAC;
sigOut_DLAC = signal;


% sigOut_DLAC = sum(sigAC,2);

bwChannel = ACconfig.bwChannel
fs = ACconfig.fs
Nsamps = length(sigAC)
df = Nsamps/fs
bwInband = ACconfig.bwInband
PLOT_FFT_dB_g(sigOut_DLAC, fs, Nsamps, ['sigDLchGain, output pwrdB:',num2str([])], 'df', 'full', 'pwr', [111+4, 1,2,1]);



ACDomod_Nshift = 0;
fnum=101+1;
phsShiftNbrDeg = [0 10 -7 2]

sigAC_phsShift = sigOut_DLAC.*exp(1i*phsShiftNbrDeg./180*pi);

sigAC_sum = sum(sigAC_phsShift,2);


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
b_ch_AC = SYM_FIRApp(FIRc1_Wtype,FIRc1_Ftype,FIRc1_Order,FIRc1_K_AttdB,FIRc1_K_fdelta,FIRc1_fTolerance,FIRc1_fcutoffL,FIRc1_fcutoffH,df,bwInband,NCarriers,FIRc1_Export,fnum+2);
b_ch_AC = b_ch_AC{:};

sigAC_sum_ch(:,1) = conv(sigAC_sum(:,1), b_ch_AC, 'same');
PLOT_FFT_dB_g(sigAC_sum_ch(:,:), fs, Nsamps*1, ['sigAC + channel filter'], 'df', 'full', 'pwr', [fnum+2]);



    %% D7b. AC DeModulation
close all
[ACdmod_t0Mean, ACdmod_p0DegMean, ACdmod_phEstDeg, ACdmod_phEstDegDrift, ACdmod_dataCapCor, ACdmod_SNR, ACdmod_dataCapNbrwoPD]...
    = AntCal_phaseDemodulateApp_g100(sigAC_sum_ch,ACconfig,ACDomod_Nshift,fnum, []);

ConfigAC_2 = ConfigAC
ConfigAC_2.Nfft = 512
[ACdmod_t0Mean, ACdmod_p0DegMean, ACdmod_phEstDeg, ACdmod_phEstDegDrift, ACdmod_dataCapCor, ACdmod_SNR, ACdmod_dataCapNbrwoPD]...
    = AntCal_phaseDemodulateApp_g100(sigAC_sum_ch,ConfigAC_2,ACDomod_Nshift,fnum, []);

%% D7c. AC Calculate EQ coefficents and Phase Alignment
ACEQ_numFirTaps = 3;
ACEQ_NoOfAlignBr = 1;

% delta=fix(122.88e6/23.04e6)
% center = Nsamps/2
% NSubcarrier_half = fix(ACconfig.NSubcarrier/2)
% sigAC_DMC = sigAC(center-NSubcarrier_half:center+NSubcarrier_half);
% PLOT_FFT_dB_g(sigAC_DMC(:,:), fs, Nsamps*1, ['sigAC + channel filter'], 'df', 'full', 'pwr', [fnum+2]);
% 
% 
% 
% 
% PLOT_FFT_dB_g(sigAC(:,:), fs, Nsamps*4, ['sigAC + channel filter'], 'df', 'full', 'pwr', [fnum+2]);
% for ibr=1:size(sigAC,2)
% sigAC_ch(:,ibr) = conv(sigAC(:,ibr), b_ch_AC, 'same');
% end
% PLOT_FFT_dB_g(sigAC_ch(:,:), fs, Nsamps*4, ['sigAC + channel filter'], 'df', 'full', 'pwr', [fnum+2]);
% 
% r=8
% sigAC_DMC = sigAC_ch(1:r:end,:);
% f1 = fftshift(fft(sigAC_ch, [], 2));
% sigAC_DMC = ifft(fftshift(f1(center-NSubcarrier_half:center+NSubcarrier_half,:)), [], 2);
% PLOT_FFT_dB_g(sigAC_DMC(:,:), fs/r, length(sigAC_DMC)*1, ['sigAC + channel filter'], 'df', 'full', 'pwr', [fnum+2]);
% 


%% 
ACconfig.Nfft = 23.04e6/240e3

fnum=fnum+1;
[waveformD7_ACPhsAlign,ACEQ_taps,ACEQ_B] = AntCal_EQphsAlignApp(sigAC,ACconfig,[], ACdmod_phEstDeg,ACEQ_numFirTaps,ACEQ_NoOfAlignBr,fnum, []);
ACPhsShiftNbrDeg_T7_Ref = (ACPhsShiftNbrDeg_T7(ACEQ_NoOfAlignBr)-ACPhsShiftNbrDeg_T7).'.*ones(1,length(ACconfig.nSubcarrier));
for idBR = 1:size(ACPhsShiftNbrDeg_T7_Ref,1)
    plot(ACconfig.nSubcarrier, ACPhsShiftNbrDeg_T7_Ref(idBR,:),'r','DisplayName',['Ref idBR',num2str(idBR),', [',num2str(mean(ACPhsShiftNbrDeg_T7_Ref(idBR,:))),']deg'])
    legend;
end



  %% D7. OFDM DeModulation
        fnum=fnum+1;
        [rxGridD7cell, EVMd7_DemodulationDL] = SYM_DeMODApp(waveformD7cell_DL, waveformD7cell_ref, DLconfig, fs, bwInbandC0, fnum, fnum_dir);
          [rxGridD7cell, EVMd7_DemodulationDL] = SYM_DeMODApp(sigDL_1, sigDL, Config, fs, bwInband, fnum, []);
r=16
sigDL_1 = sigDL(1:r:end);
PLOT_FFT_dB_g(sigDL_1(:,:), fs/r, length(sigDL_1)*1, ['sigAC + channel filter'], 'df', 'full', 'pwr', [fnum+2]);
PLOT_FFT_dB_g(sigDL(:,:), fs, length(sigDL)*1, ['sigAC + channel filter'], 'df', 'full', 'pwr', [fnum+2]);
sigDL_1 = fftshift(fft(sigDL));

Nfft = Config.Nfft
df_Scs = Config.dfScs
NScs = Config.NScs
f = (-Nfft/2:Nfft/2-1)*df_Scs;
f_Scs = [(-NScs/2:-1),(1:NScs/2)]*df_Scs; % exclude DC
ind_Scs = find(ismember(f,f_Scs)==1);
sigDL_2 = zeros(23.04e6/30e3,1);

Nfft2 = 23.04e6/30e3
df_Scs = Config.dfScs
NScs = Config.NScs
f = (-Nfft2/2:Nfft2/2-1)*df_Scs;
f_Scs2 = [(-NScs/2:-1),(1:NScs/2)]*df_Scs; % exclude DC
ind_Scs2 = find(ismember(f,f_Scs2)==1);

sigDL_2(ind_Scs2) = sigDL_1(ind_Scs);
sigDL_3 = ifft(fftshift(sigDL_2));

PLOT_FFT_dB_g(sigDL_3(:,:), fs, length(sigDL_3)*1, ['sigAC + channel filter'], 'df', 'full', 'pwr', [fnum+2]);
%% !
fs=Config.fs
bwInband=[]
fnum=1030
[rxGridD7cell, EVMd7_DemodulationDL] = SYM_DeMODApp(sigDL, [], Config, fs, bwInband, fnum, []);

%% !!
[sigDL, ~, Config] = OFDM_Mod_g('20MHz','64QAM','off','LTE',fnum,[]);
rxGrid = OFDM_Demod_g(Config, sigDL);
PLOT_Constellation(rxGrid,[],103001)

[sigDL1, ~, Config1] = OFDM_Mod_g('20MHz','QPSK','off','NR',fnum,[]);
rxGrid1 = OFDM_Demod_g(Config1, sigDL1);
PLOT_Constellation(rxGrid1,[],103001)

DLconfig.bw_Channel = '20MHz'
DLconfig.fs = 122.88e6;
[sigDL2, ~, Config2] = OFDM_Mod_g(DLconfig,'64QAM','off','NR',fnum,[]);
rxGrid2 = OFDM_Demod_g(Config2, sigDL2);
rxGrid2b = OFDM_Demod_g(Config2, sigDL2);
PLOT_Constellation(rxGrid2,[],103002,[])
PLOT_Constellation(rxGrid2b,[],103002,[],0)

