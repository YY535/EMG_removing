function EMG_rm_view(varargin)
% EMG_rm_view([FileBase, Period, Channels, WinLen, useT, sfactor, figure_handel,t_pause])
% 
% Viewing the cleaned data. 
% To use this function, please make sure that your .lfp and .lfpd files are
% within your searching directories.  
% 
% Inputs: 
%   FileBase: session name. 
%   Optional: 
%       Period: the periods you want to check, default: .sts.RUN.
%       Channels: starting from 1, default: data in the first shank/AnatGrps.
%       WinLen: shown Window length in samples, default: lfp sampling rate, aka. 1s. 
%       useT: if use time as unit of x axis, default: true.
%       sfactor: rescale factor, if not given, will be chosen by the program. 
%       figure_handel: the figure
%       t_pause: if play the data automatically. if t_pause>0, then play with
%               pause(t_pause), otherwise it's pause. 
%       ifEMGthr: if plot EMG_thrd
% 
% Related functions: 
% EMG_rm_main.m, EMG_rm_pip.m
%  
% Error contact: chen at biologie.uni-muenchen.de
% 
% Last Modified: 27.11.2019.


[FileBase,Period,Channels,WinLen,useT,sfactor,figure_handel,t_pause,ifEMGthr] = DefaultArgs(varargin, {[],[],[],[],true,-1,0,0,false});

if isempty(FileBase)
    cur_dir = pwd;
    FileBase = cur_dir((find(cur_dir == '/',1,'last')+1):end);
end

if exist([FileBase,'.lfpinterp'],'file')
    LFPfile = [FileBase,'.lfpinterp'];
elseif exist([FileBase,'.lfp'],'file')
    LFPfile = [FileBase,'.lfp'];
else
    fprintf('\n\nNo LFP file is found in this directory, exiting...\n\n')
    return
end
if exist([FileBase,'.lfpd'],'file')
    LFPdfile = [FileBase,'.lfpd'];
else
    fprintf('\n\nNo cleaned file is found on this directory, Please clean the data with EMG_rm_main first.\n\n')
    return
end
if figure_handel<1
    figure_handel = figure;
end
a = dir('*.EMG_Cluster.mat');
try
    load(a(1).name,'EMG_thrd')
catch
    warning('No EMG_Cluster file in current folder!!')
end
%% THE DATA TO SHOW
par=LoadXml(FileBase);
if isempty(Channels)
    Channels = par.AnatGrps(1).Channels +1- par.AnatGrps(1).Channels(1);
end
Fs = par.lfpSampleRate;
if isempty(WinLen)
    WinLen = Fs;
end
mlfp = memmapfile(LFPfile,'Format','int16');
Data_Range = [par.nChannels, length(mlfp.Data)/par.nChannels];

if isempty(Period)
    if exist([FileBase,'.sts.RUN'],'file')
        Period = load([FileBase,'.sts.RUN']);
    else
        a = dir('*.sts.*');
        if length(a)>0
            Period = load(a(1).name);
            fprintf('\nViewing the %s period...\n', a(1).name((find(a(1).name=='.',1,'last')+1):end))
        else
            Period = [1 12500];
        end
    end
end
% Smaller Chunks. 
nPeriod = size(Period,1);
tmp_Period = cell(nPeriod,1);
for k = 1:nPeriod
    tmp = Period(k,1):WinLen:Period(k,2);
    tmp_Period{k} = [tmp(1:(end-1))' tmp(2:end)'];
end
Period = cell2mat(tmp_Period);
nPeriod = size(Period,1);

% LOAD & PLOT DATA
mlfp = memmapfile(LFPfile,'Format',{'int16',Data_Range,'x'});
m = memmapfile(LFPdfile,'Format',{'int16',Data_Range,'x'});
if sfactor<0
    v = double(mlfp.Data.x(Channels(fix(end/2)),:));
    sfactor = 2^(ceil(log2(std(v)))); % heuristic
end
% report rescale factor, if the user is not happy with it, they can change
% that. 
fprintf('\nThe rescale factor is %d.\n', sfactor)

% SHOW THE PLOTS

opf1 = @(x)(bsxfun(@plus,x,-Channels(:)'));

if useT
    Xtitle = 'time (s)';
else
    Xtitle = 'Samples';
end
if ~t_pause
    fprintf('\nPlease press anything to continue.\n', sfactor)
end
for k = 1:nPeriod
    figure(figure_handel);clf
    
    % the original signal
    s1=subplot(2,1,1);
    
    tmp = Period(k,1):Period(k,2);
    if useT
        tmp_t = tmp/Fs;% unit: S
    else
        tmp_t = tmp;
    end
    tmp_x = double(mlfp.Data.x(Channels,tmp))';
    plot(tmp_t,opf1(tmp_x/sfactor),'k')
    if ifEMGthr
        hold on
        plot(tmp_t, EMG_thrd(tmp)*10,'g+')
    end
    ylabel('channels')
    title('Before')
    axis tight
    
    % the cleaned signal
    s2=subplot(2,1,2);
    
    tmp_x = double(m.Data.x(Channels,tmp))';
    plot(tmp_t,opf1(tmp_x/sfactor),'r');
    axis tight
    title('After')
    xlabel(Xtitle)
    
    linkaxes([s1 s2],'xy')
    if ~t_pause
        pause
    else
        pause(t_pause)
    end
end