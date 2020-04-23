function EMG_rm_main_group(FileBase,denoise_groups,varargin)
% EMG_rm_main_group(FileBases,denoise_groups,[savedir,cleanEMGch,...
%                        silence_periods,sp_loadingfuns,...
%                        rm_linenoise,line_thrd, ...
%                        numOfIC,hp_freq,nchunks, cmp_method,down_sample, ...
%                        nFFT, Winlength, ...
%                        save_together])
%
% The main function to perform denoising. 
% This function is working under the session directory
%
% Inputs:
%   FileBase: session name.
%   denoise_groups: the shank groups for denoising in cell. e.g.,
%                  {[shk_1,shank2],shk_n} 
%   Optional:
%       keep_former: keep the former denoising results
%       savedir: path to save the cleaned files.
%       cleanEMGch: the channels to detect EMG noise. defualt: 5 sampled in
%                   the first shank.
%       silence_periods: only compute and remove the EMG from non_silence 
%                       periods, give the periods you don't want to remove 
%                       anything or their filenames. default: false 
%       sp_loadingfuns: functions to load the discarded periods. 
%       rm_linenoise: if remove line noise. default: true
%       line_thrd: power ratio between line band and other band. 
%                  defualt: 1.8, I'm being conservative here.
%       Parameters used by the subfunctions: (computational details)
%           numOfIC: number of ICs.
%           hp_freq: high_pass_freq, defualt: 100 Hz
%           nchunks: number of chunks, defult: 6 or every chunk less than 15 mins,
%           cmp_method: methods to compute the EMG components. Whiten('w') or
%                       Highpass&Whiten('hw', a bit more stable). defualt: 'hw'
%           down_sample: down sample the data to compute the EMG
%                       components. default: 3.
%           nFFT: nfft to compute spectrum in EMG_rm_viewspec.m
%           WinLength: length of moving window in EMG_rm_viewspec.m
%           save_together: if save all the chunks together. defult: true
%           use_wb: whether to use the wide band signal when the algorithm
%               can't find a proper flat component at high frequency band
%               defualt: false NB: please check the report fig to see if you
%               really need this! 
%           PeriodLengthLimits:the maximum allowed length of your period.
%               default: 1e6
% 
% Related functions: 
% EMG_Cluster.m, EMG_rm_long.m, EMG_rm_pip.m, EMG_rm_viewspec.m, 
% EMG_rm_report.m, EMG_rm_viewnoise.m
%
% NB: please make sure you remove all the files generated previously when you
% try to redo the denoising.
%
% Error contact: chen at biologie.uni-muenchen.de
% 
[keep_former,savedir,cleanEMGch,silence_periods,sp_loadingfuns,rm_linenoise,line_thrd, numOfIC, hp_freq,nchunks, cmp_method,down_sample,nFFT, Winlength, save_together,use_wb,PeriodLengthLimits] = DefaultArgs(varargin, {false,pwd,[],false,[],true,1.8,[], 100,6, 'hw',3, 2^9,[],true,false,1e6});

if exist([FileBase,'.lfpinterp'],'file')
    LFPfile = [FileBase,'.lfpinterp'];
    try
        armodel = load([FileBase, '.WhitenSignal.ARmodel.lfpinterp.mat']);
    catch
        armodel = [];
    end
elseif exist([FileBase,'.lfp'],'file')
    LFPfile = [FileBase,'.lfp'];
    try
        armodel = load([FileBase, '.WhitenSignal.ARmodel.lfp.mat']);
    catch
        armodel = [];
    end
else
    fprintf('\n\nNo LFP file is found in this directory, exiting...\n\n')
    return
end

par=loadxml(FileBase);
nG = length(denoise_groups);
HPs = cell(nG,1);
Shk01 = cell(nG,1);
denoise_shank=[];
for  k = 1:nG
    nshk = length(denoise_groups{k});
    denoise_shank=[denoise_shank,denoise_groups{k}(:)'];
    Shk01{k} = zeros(nshk,2);
    last_channel = 0;
    for n  =1:nshk
        used_channels = ~par.AnatGrps(denoise_groups{k}(n)).Skip(:);
        tmp_channels = par.AnatGrps(denoise_groups{k}(n)).Channels(used_channels);
        HPs{k} = [HPs{k};tmp_channels(:)];
        Shk01{k}(n,:) = last_channel+[1 sum(used_channels)];
        last_channel=Shk01{k}(n,2);
    end
end

clear G_par
G_par.denoise_groups = denoise_groups;
G_par.HPs = HPs;
G_par.Shk01=Shk01;

if isempty(HPs{1})
    error('Please check your recording parameters.\n The Groups should not be empty!')
end

% use a long shank to detect the noisy periods.
[~,n]= max(par.GroupChannels);
HP = par.AnatGrps(n).Channels(~par.AnatGrps(n).Skip(:))+1;% assuming the chanels coding starts from 0. 
Fs = par.lfpSampleRate;
if isempty(cleanEMGch)
    cleanEMGch = fix(linspace(3,length(HP)-3,5));% channels to detect high coherence periods.
end
if isempty(numOfIC)
    numOfIC = max(length(HP)*.8,min(16, length(HP)));
end
if isempty(Winlength)
    Winlength = Fs;
end
% [Evts,~] = RecEvents([FileBase,'.cat.evt'],Fs);
% Evts = max(Evts,1);
% Automatically compute the number of samples. 
mlfp = memmapfile(LFPfile,'Format','int16');
Evts = [1 length(mlfp.Data)/par.nChannels];
clear mlfp

%% CREATING .LFPD FILE
myData = repmat(int16(0),1,par.nChannels*Evts(end));
FileName = [FileBase,'.lfpd'];
if ~exist(FileName,'file')
    fileID = fopen(FileName,'w');
    fwrite(fileID, myData,'int16');
    fclose(fileID);
    clear myData
end

% myData = int16(zeros(1,Evts(end)));
% EMGFileName = [FileBase,'.emg'];
% if ~exist(EMGFileName,'file')
%     fileID = fopen(EMGFileName,'w');
%     fwrite(fileID, myData,'int16');
%     fclose(fileID);
%     clear myData
% end

save_range = cell(2,1);
save_range{2} = [par.nChannels, Evts(end)];

%% DETECTING THE HIGH EMG PRIODS
if exist([FileBase, '.EMG_Cluster.mat'], 'file')
    load([FileBase, '.EMG_Cluster.mat'],'EMG_thrd', 'sug_period')
    fprintf('\n\n--------------------------------------------------\n                   Please!!! \n  Remove all the previously generated files \n     when you try to redo the denoising!\n--------------------------------------------------\n\n')
    if silence_periods
        load([FileBase, '.EMG_Cluster.mat'],'included_periods')
    else
        included_periods=[];
    end
else
    swin = 500;
    lwinms = 20;
    
    lfp = loadbinary(LFPfile,cleanEMGch,par.nChannels,2,[],[],Evts)';%
    lfp=bsxfun(@minus,lfp,mean(lfp,1));
    
    fprintf('\nDetecting EMG periods...\n')
    if silence_periods
        if ischar(silence_periods) % load from files
            included_periods = LoadAwake(silence_periods,Evts(end),sp_loadingfuns);
        elseif length(silence_periods)<2 && silence_periods % if true load from files
            included_periods = LoadAwake(sprintf('%s.sts.%s', FileBase, 'SWS'),Evts(end),sp_loadingfuns);
        elseif length(silence_periods)<Evts(end) % given periods
            included_periods = true(Evts(end),1);
            for k = 1:size(silence_periods,1)
                included_periods(silence_periods(k,1):silence_periods(k,2)) = false;
            end
        else % binary sample mask
            included_periods = ~silence_periods;
        end
        silence_periods = true;% a flag
        [EMG_thrd, ~, sug_period,included_periods] = EMG_Cluster_s(lfp,included_periods,hp_freq, Fs, lwinms, nchunks,...
            swin, true, FileBase, PeriodLengthLimits);
    else
        [EMG_thrd, ~, sug_period] = EMG_Cluster(lfp,hp_freq, Fs, lwinms, nchunks,...
            swin, true, FileBase);
        silence_periods = false;
        included_periods = [];
    end
    clear lfp
end

%% REMOVING THE ARTIFACTS
fprintf('\nEMG artifacts removing...\n')

nPeriod = size(sug_period,1);
HPses = [];
Ws = cell(nPeriod,nG);
As = cell(nPeriod,nG);
EMG_au = cell(nPeriod,nG);
AW  = cell(nPeriod,nG);
for n = 1:nG
    HP = HPs{n} +1;% - par.AnatGrps(1).Channels(1);
    HPses = [HPses;HP(:)];
    
    myData = int16(zeros(1,Evts(end)));
    EMGFileName = sprintf('%s%s.G%d.emg',savedir,FileName,n);
    if ~exist(EMGFileName,'file')
        fileID = fopen(EMGFileName,'w');
        fwrite(fileID, myData,'int16');
        fclose(fileID);
        clear myData
    end
    
    for k = 1:nPeriod
        fprintf('\rshank%d, period%d in %d...', n, k, nPeriod)
        tmp_Period = sug_period(k,:);
        save_range{1} = [tmp_Period;HP(1) HP(end)];
        if ~isempty(included_periods)
            tmp_included_periods=included_periods(tmp_Period(1):tmp_Period(2));
        else
            tmp_included_periods=[];
        end
        [~, Ws{k,n}, As{k,n}, EMG_au{k,n}, AW{k,n}, armodel,scaling_factor(k,n)] = EMG_rm_long(LoadBinary(LFPfile,HP,par.nChannels,2,[],[],tmp_Period)',   ...
            silence_periods,tmp_included_periods,...
            Fs, rm_linenoise,line_thrd,hp_freq, ...
            EMG_thrd(tmp_Period(1):tmp_Period(2)), true, ...
            armodel, cmp_method, down_sample,numOfIC,...
            true, save_range, FileBase, savedir, save_together,use_wb,denoise_groups{n},HP,Shk01{n});
    end
end

if save_together
    shank_names = sprintf('%d-',denoise_shank);
    save(sprintf('%s/%s.EMG_rm.G%s.mat',savedir,FileBase,shank_names(1:(end-1))), 'Ws','As','AW','EMG_au','armodel','sug_period','par','denoise_shank','LFPfile','scaling_factor')
end

%% COMPLETE OTHER CHANNELS
other_channels = setdiff(1:par.nChannels,HPses);
nothrch = length(other_channels);
m = memmapfile(FileName,'Format',{'int16',save_range{2},'x'},'Writable',true);
mlfp = memmapfile(LFPfile,'Format',{'int16',save_range{2},'x'});
for k = 1:nothrch
    m.Data.x(other_channels(k),:) = mlfp.Data.x(other_channels(k),:);
end
clear m mflp
fprintf('\nDone\n')
%% CHECK RESULTS:
PYR_Channel = 37;
EMG_rm_report();% ([],PYR_Channel);
EMG_rm_viewnoise();% (PYR_Channel,[])
if Evts(end)>1.6e7
    nFFT = 2^(32-ceil(log2(Evts(end))));
    fprintf('\nNB: we are using nFFT: %d, recompute EMG_rm_viewspec() if you want.',nFFT)
end
EMG_rm_viewspec([],[],[],[],nFFT,Winlength);
% EOF