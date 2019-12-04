function EMG_rm_viewspec(Channel,Periods)
% EMG_rm_viewspec(Channel,Periods)
% 
% view the effect of EMG removal in specragram.
% Inputs: 
%   Channel: the channels to compare. 
%   Periods: the periods to check.
% 
% Error contact: chen at biologie.uni-muenchen.de
% 
% Last Modified: 01.12.2019.

%%
a = dir('*.EMG_rm.sh*.mat');
FileName = a(1).name;
load(FileName)
FileBase = FileName(1:(find(FileName=='.',1,'first')-1));
a = dir('*.EMG_Cluster.mat');
load(a(1).name)
if isempty(Periods)
    a = dir('*.sts.*');
    nP = length(a);
    Periods = cell(nP,1);
    for k = 1:nP
        Periods{k} = load(a(k).name);
        P_title{k} = a(k).name((find(a(k).name=='.',1,'last')+1):end);
    end
elseif ~iscell(Periods)
    tmp = cell(1);
    tmp{1} = Periods;
    Periods = tmp;
    nP = length(Periods);
    P_title{1} = 'Period';
else
    nP = length(Periods);
end

%% DATA PROCESSING
mlfp = memmapfile(LFPfile,'Format','int16');
data_range = [par.nChannels, length(mlfp.Data)/par.nChannels];
mlfp = memmapfile(LFPfile,'Format',{'int16',data_range,'x'});
lfp = double(mlfp.Data.x(Channel,:))';
lfp = bsxfun(@minus,lfp,mean(lfp));%zscore(lfp);
mlfp = memmapfile([FileBase, '.lfpd'],'Format',{'int16',data_range,'x'});
dlfp = double(mlfp.Data.x(Channel,:))';
dlfp = bsxfun(@minus,dlfp,mean(dlfp));%zscore(dlfp);
EMG = lfp-dlfp;
if isfield(armodel,'ARmodel')
    lfp=WhitenSignal(lfp,[],[],armodel.ARmodel);
    dlfp=WhitenSignal(dlfp,[],[],armodel.ARmodel);
    EMG=WhitenSignal(EMG,[],[],armodel.ARmodel);
else
    lfp=WhitenSignal(lfp,[],[],armodel);
    dlfp=WhitenSignal(dlfp,[],[],armodel);
    EMG=WhitenSignal(EMG,[],[],armodel);
end
[y2, f, t, phi, ~,~] = mtcsdlong_E([EMG,lfp,dlfp],[],par.lfpSampleRate);
%% AVAERAGE CROSS-SPECTRAL DENSITY CROSS STATES
nt = length(t);
tt = 1:500;
nclm = 4;
nrows = ceil(nP/nclm)+3;
figure
for k = 1:nP
    subplot(nrows,nclm,k)
    tmp_t = false(nt,1);
    tmp_Period = Periods{k}/par.lfpSampleRate;
    for n = 1:size(tmp_Period,1)
        tmp_t(t<tmp_Period(n,2) & t>=tmp_Period(n,1)) = true;
    end
    plot(f, sq(mean(y2(tmp_t,:,[1 5 9 4 7]),1)))
    title(P_title{k})
end
legend('EMG','lfp','dlfp','EMG-lfp','EMG-dlfp')
%% THE SPECTROGRAM
subplot(nrows,1,nrows-2)
c = pcolor(t(tt),f(f<500),sq(y2(tt,f<500,1,1))');
c.EdgeColor='none';
ylabel('EMG')
subplot(nrows,1,nrows-1)
c = pcolor(t(tt),f,sq(y2(tt,:,2,2))');
c.EdgeColor='none';
ylabel('lfp')
subplot(nrows,1,nrows)
c = pcolor(t(tt),f,sq(y2(tt,:,3,3))');
c.EdgeColor='none';
ylabel('dlfp')
xlabel('time(s)')


