function EMG_rm_main(FileBase,varargin)
% EMG_rm_main(FileBases,[savedir,denoise_shank,cleanEMGch,...
%                        rm_linenoise,line_thrd,...
%                        numOfIC,hp_freq,nchunks, cmp_method,down_sample, ...
%                        save_together])
%
% The main function to perform denoising. 
% This function is working under the session directory
%
% Inputs:
%   FileBase: session name.
%   Optional:
%       savedir: path to save the cleaned files.
%       denoise_shank: the shanks for denoising. defualt: shank number 1
%       cleanEMGch: the channels to detect EMG noise. defualt: 5 sampled in
%                   the first shank.
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
%           save_together: if save all the chunks together. defult: true
% 
% Related functions: 
% EMG_Cluster.m, EMG_rm_long.m, EMG_rm_pip.m
%
% Error contact: chen at biologie.uni-muenchen.de
% 
% Last Modified: 11.12.2019.

[savedir,denoise_shank,cleanEMGch,rm_linenoise,line_thrd, numOfIC, hp_freq,nchunks, cmp_method,down_sample, save_together] = DefaultArgs(varargin, {pwd,1,[],true,1.8,50, 100,6, 'hw',3, true});

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

par=LoadXml(FileBase);
HP = par.AnatGrps(denoise_shank(1)).Channels +1- par.AnatGrps(1).Channels(1);
Fs = par.lfpSampleRate;
if isempty(cleanEMGch)
    cleanEMGch = fix(linspace(3,length(HP)-3,5));% channels to detect high coherence periods.
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

save_range = cell(2,1);
save_range{2} = [par.nChannels, Evts(end)];

%% DETECTING THE HIGH EMG PRIODS
if exist([FileBase, '.EMG_Cluster.mat'], 'file')
    load([FileBase, '.EMG_Cluster.mat'],'EMG_thrd', 'sug_period')
else
    swin = 500;
    lwinms = 20;
    
    lfp = LoadBinary(LFPfile,cleanEMGch,par.nChannels,2,[],[],Evts)';%
    lfp=bsxfun(@minus,lfp,mean(lfp,1));
    
    fprintf('\nDetecting EMG periods...\n')
    
    [EMG_thrd, ~, sug_period] = EMG_Cluster(lfp,hp_freq, Fs, lwinms, nchunks,...
        swin, true, FileBase);
    clear lfp
end

%% REMOVING THE ARTIFACTS
fprintf('\nEMG artifacts removing...\n')

nPeriod = size(sug_period,1);
nshank = length(denoise_shank);
HPs = [];
Ws = cell(nPeriod,nshank);
As = cell(nPeriod,nshank);
EMG_au = cell(nPeriod,nshank);
AW  = cell(nPeriod,nshank);
for n = 1:nshank
    HP = par.AnatGrps(denoise_shank(n)).Channels +1- par.AnatGrps(1).Channels(1);
    HPs = [HPs;HP(:)];
    for k = 1:nPeriod
        fprintf('\rshank%d, period%d in %d...', n, k, nPeriod)
        tmp_Period = sug_period(k,:);
        save_range{1} = [tmp_Period;HP(1) HP(end)];
        [~, Ws{k,n}, As{k,n}, EMG_au{k,n}, AW{k,n}, armodel] = EMG_rm_long(LoadBinary(LFPfile,HP,par.nChannels,2,[],[],tmp_Period)',   ...
            Fs, rm_linenoise,line_thrd,hp_freq, ...
            EMG_thrd(tmp_Period(1):tmp_Period(2)), true, ...
            armodel, cmp_method, down_sample,numOfIC,...
            true, save_range, FileBase, savedir, save_together);
    end
end

if save_together
    shank_names = sprintf('%d-',denoise_shank);
    save(sprintf('%s/%s.EMG_rm.sh%s.mat',savedir,FileBase,shank_names(1:(end-1))), 'Ws','As','AW','EMG_au','armodel','sug_period','par','denoise_shank','LFPfile')
end

%% COMPLETE OTHER CHANNELS
other_channels = setdiff(1:par.nChannels,HPs);
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
EMG_rm_viewspec([],[],[],[],2^9,par.lfpSampleRate);

% EOF