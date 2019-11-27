function EMG_rm_main(FileBase,varargin)
% EMG_rm_main(FileBases,[savedir,denoise_shank])
%
% The main function to perform denoising. 
% This function is working under the session directory
%
% Inputs:
%   FileBase: session name.
%   Optional:
%       savedir: path to save the cleaned files.
%       denoise_shank: the shanks for denoising. defualt: shank number 1
% 
% Related functions: 
% EMG_Cluster.m, EMG_rm_long.m, EMG_rm_pip.m
%
% Error contact: chen at biologie.uni-muenchen.de
% 
% Last Modified: 27.11.2019.

[savedir,denoise_shank, hp_freq,ntrunks, cmp_method,down_sample] = DefaultArgs(varargin, {pwd,1, 100,6, 'hw',3});

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
cleanEMGch = fix(linspace(3,length(HP)-3,5));% channels to detect high coherence periods.
[Evts,~] = RecEvents([FileBase,'.cat.evt'],Fs);
Evts = max(Evts,1);

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
swin = 500;
lwinms = 20;

lfp = LoadBinary(LFPfile,cleanEMGch,par.nChannels,2,[],[],Evts')';% 
lfp=bsxfun(@minus,lfp,mean(lfp,1));

fprintf('\nDetecting EMG periods...\n')

[EMG_thrd, ~, sug_period] = EMG_Cluster(lfp,hp_freq, Fs, lwinms, ntrunks,...
swin, true, FileBase);
clear lfp

%% REMOVING THE ARTIFACTS
fprintf('\nEMG artifacts removing...\n')

nPeriod = size(sug_period,1);
HPs = [];
for n = 1:length(denoise_shank)
    HP = par.AnatGrps(denoise_shank(n)).Channels +1- par.AnatGrps(1).Channels(1);
    HPs = [HPs;HP(:)];
    for k = 1:nPeriod
        tmp_Period = sug_period(k,:);
        save_range{1} = [tmp_Period;HP(1) HP(end)];
        lfp = LoadBinary(LFPfile,HP,par.nChannels,2,[],[],tmp_Period)';%
        EMG_rm_long(lfp,   Fs, hp_freq, ...
            EMG_thrd(tmp_Period(1):tmp_Period(2)), true, ...
            armodel, cmp_method, down_sample,...
            true, save_range, FileBase, savedir);
        clear lfp
    end
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
