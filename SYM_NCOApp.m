function [waveform_NCO_idCcell, waveform_NCO_Combcell, bwInbandOut_idC, evmNCOComb, tableIOputNCO_idC] = NCOApp(waveform, fNCO, fs, bwInband, flag_SUMCarrier, fnum, disp_title, fnum_save_dir)
% flag_SUMCarrier = SUMon/SUMoff

% if exist('flag_SUMCarrier','var')&&~isempty(flag_SUMCarrier)&&strcmp(flag_SUMCarrier,'SUM')
%     fno=1;
% else
%     fno=0;
% end
if ~exist('fnum_save_dir')||isempty(fnum_save_dir)||fnum_save_dir==0
    fnum_save_dir = 0;
else
    fnum_save_dir = 1;
end

if strcmp(flag_SUMCarrier,'SUM')||isempty(flag_SUMCarrier)||all(flag_SUMCarrier==1)
    flag_SUMCarrier = 'SUM';
    if iscell(waveform)
        NCarriers = size(waveform,1);
    else
        NCarriers = size(waveform,3);
    end
    if NCarriers == 1
        NCarriers = numel(fNCO);
    end
elseif strcmp(flag_SUMCarrier,'DIV')||flag_SUMCarrier==0 % divide
    flag_SUMCarrier = 'DIV';
    NCarriers = numel(fNCO);
elseif ~exist('flag_SUMCarrier','var')||isempty(flag_SUMCarrier)||(NCarriers==1)
    flag_SUMCarrier = 'NA';
end

% if NCarriers==1
%     flag_SUMCarrier = [];
% end

for idC=1:NCarriers
    if iscell(waveform)
        flag_Format_Cell = 1;
        if size(waveform,1)>1
            waveform_idC = waveform{idC};
        else
            waveform_idC = waveform{:};
        end
    else
        flag_Format_Cell = 0;
        if size(waveform,3)>1
            waveform_idC = waveform(:,:,idC);
        else size(waveform,3)==1
            waveform_idC = waveform;
        end
    end
    Nsamps = length(waveform_idC);
    
    NCarriersOfFS=size(fs,1);
    if NCarriersOfFS~=1
        fs_idC = fs(idC,:);
    else
        fs_idC = fs;
    end
    df = fs_idC/Nsamps;
    
    if isrow(fNCO)
        fNCO = fNCO(:)
    end
    NCarriersOfFNCO=size(fNCO,1);
    if NCarriersOfFNCO~=1
        fNCO_idC = ceil(fNCO(idC,:)/df)*df;
    else
        fNCO_idC = ceil(fNCO/df)*df;
    end
    
    if exist('bwInband','var')&&~isempty(bwInband)
        NCarriersOfIbw=size(bwInband,1);
        if NCarriersOfIbw~=1
            %         bwInband_idC = bwInband(idC,:);
            bwInband_idC = ceil(bwInband(idC,:)/df)*df;
        else
            %         bwInband_idC = bwInband;
            bwInband_idC = ceil(bwInband/df)*df;
        end
        
        bwInbandOut_idC(idC,:) = fNCO_idC+bwInband_idC;
    else
        bwInbandOut_idC = [];
    end
    
    fNCO_MHz(idC,:) = fNCO_idC/1e6;
    disp(['NCO with ',num2str(fNCO_MHz(idC,:)),'MHz ===================================================='])
    
%     Nsamps = length(waveform_idC);
    Nbr = size(waveform_idC,2);
    [~,DIMFFT] = max(size(waveform_idC));
    
    %% Generate Multi-Carriers
    if DIMFFT ==1
        t=([0:Nsamps-1]/fs_idC).'; % COLUMN
    else
        t=([0:Nsamps-1]/fs_idC); % ROW
    end
    waveform_NCO_idC=waveform_idC.*exp(1i*2*pi*fNCO_idC*t);
    waveform_NCO_idCcell{idC,:}=waveform_NCO_idC;
    
    % plot
    if exist('fnum','var')&&~isempty(fnum)
        if NCarriers>1
            disp_idC = ['idC',num2str(idC)];
        else
            disp_idC = [];
        end
        if exist('disp_title','var')&&~isempty(disp_title)
            NCarriersOfTitle=size(disp_title,1);
            if NCarriersOfTitle~=1
                disp_title_idC = cell2mat(disp_title(idC,:));
            else
                if iscell(disp_title)
                    disp_title_idC = cell2mat(disp_title);
                elseif ischar(disp_title)
                    disp_title_idC = (disp_title);
                end
            end
        else
            disp_title_idC=[];
        end
        if strcmp(flag_SUMCarrier,'SUM')
            fnum2 = [fnum(1),NCarriers,2,2*idC-1];
        else
            fnum2 = [fnum(1),NCarriers,1,idC];
        end
        PLOT_FFT_dB_g(waveform_NCO_idCcell{idC,:}, fs_idC, Nsamps, [disp_idC], 1, 'full', 'pwr', fnum2);
        title([disp_title_idC,' fNCO:',num2str(fNCO_idC/1e6),'MHz, fs:',num2str(fs_idC/1e6),'MHz'])
    end
    
    %% Combinate Multi-Carriers
    if strcmp(flag_SUMCarrier,'SUM')
        
        size_wf_idC=size(waveform_NCO_idCcell{idC,:});
        if idC==1
            waveform_NCO_Comb=waveform_NCO_idCcell{1,:};
            size_wfComb_idC=size(waveform_NCO_Comb);
        elseif size_wfComb_idC==size_wf_idC
            % combination
            waveform_NCO_Comb=waveform_NCO_Comb+waveform_NCO_idCcell{idC,:};
        elseif size_wf_idC(:,1)==size_wfComb_idC(:,1) % length of waveform the same
            if size_wf_idC(:,2)<size_wfComb_idC(:,2) && mod(size_wfComb_idC(:,2)/size_wf_idC(:,2),1)==0
                waveform_NCO_Comb=waveform_NCO_Comb+repmat(waveform_NCO_idCcell{idC,:},1,size_wfComb_idC(:,2)/size_wf_idC(:,2)); % repeat the matrix to multi-branches for combination
            elseif size_wf_idC(:,1)==size_wfComb_idC(:,1) && mod(size_wf_idC(:,2)/size_wfComb_idC(:,2),1)==0
                waveform_NCO_Comb=repmat(waveform_NCO_Comb,1,size_wf_idC(:,2)/size_wfComb_idC(:,2))+waveform_NCO_idCcell{idC,:};
            end
        else
            error('!')
        end
    end
    
end

%% 2020-04-28, Nsamps=307200, to reduce EVM calculation time
offset_mea = 0;

if Nsamps>307200 && (offset_mea~=0)
    NsampsEVM = 307200;
%     NsampsEVM = Nsamps;
else
    NsampsEVM = Nsamps;
end

disp_legend = [];
if ~isempty(bwInbandOut_idC)
    for idC=1:NCarriers
        mea = waveform_NCO_Comb;
        waveform_NCO_idC = waveform_NCO_idCcell{idC,:};
        ref_idC = waveform_NCO_idC;
        
        % evm calculation
        flag_EVMSync = 'off';
        flag_EVMSync= 1;
        flag_EVMFineSync = 2;
        [evmNCOComb(idC,:), indDelay] = dsp_evm_timexcorr_inband_g(ref_idC(1:NsampsEVM,:), circshift(mea(1:NsampsEVM,:),offset_mea ,DIMFFT), fs_idC, bwInbandOut_idC(idC,:),flag_EVMSync,flag_EVMFineSync);
%         [evmNCOComb(idC,:), indDelay] = dsp_evm_timexcorr_inband_g(ref_idC(:,:), circshift(mea(:,:),offset_mea ,DIMFFT), fs_idC, bwInbandOut_idC(idC,:),flag_EVMSync,flag_EVMFineSync);

        %         evmNCOComb_idC(idC,:)=evmNCOComb;
        if isempty(flag_EVMFineSync)
            disp_legend_EVM = ['idC',num2str(idC),', ','EVM:',num2str(round(evmNCOComb(idC,:),2))];
        elseif flag_EVMFineSync==2
            disp_legend_EVM = ['idC',num2str(idC),', ','evm:',num2str(round(evmNCOComb(idC,:),2))];
        end
        disp_legend = [disp_legend,disp_legend_EVM,newline];
    end
    % table
    idC = (1:NCarriers).';
    tableIOputNCO_idC(idC,:)=table(idC,fNCO_MHz,evmNCOComb);
    
else
    waveform_NCO_Combcell=[];
    evmNCOComb=[];
    
    % table
    idC = (1:NCarriers).';
    tableIOputNCO_idC(idC,:)=table(idC,fNCO_MHz);
    
end

if strcmp(flag_SUMCarrier,'SUM')
    % export waveform and plot
    waveform_NCO_Combcell = {waveform_NCO_Comb};
    fnum2 = [fnum(1),NCarriers,2,[2:2:2*NCarriers]];
    PLOT_FFT_dB_g(cell2mat(waveform_NCO_Combcell), fs_idC, Nsamps, [disp_legend], 1, 'full', 'pwr', fnum2);
    %     title([disp_title_idC,' NCO and Carrier Combination'])
    title([num2str(NCarriers), ' Carriers Combination'])
end

% end
tableIOputNCO_idC
if exist('fnum','var')&&~isempty(fnum)
    %% 2020-10-31, fnum_save_dir: save picture to folder
    if fnum_save_dir
        fnum_save_file = [fnum_save_dir,'\',num2str(fnum),'.fig']
        saveas(gcf,[fnum_save_file])
    end
end

% export
if ~flag_Format_Cell
    tmp_waveform_NCO_idCcell = waveform_NCO_idCcell;
    tmp_waveform_NCO_Combcell=waveform_NCO_Combcell;
    waveform_NCO_idCcell = zeros(Nsamps,NCarriers);
    waveform_NCO_Combcell = zeros(Nsamps,1);
    for idC=1:NCarriers
        waveform_NCO_idCcell(:,idC) = cell2mat(tmp_waveform_NCO_idCcell(idC));
    end
%     waveform_NCO_Combcell(:,1) = cell2mat(tmp_waveform_NCO_Combcell);
end
end
