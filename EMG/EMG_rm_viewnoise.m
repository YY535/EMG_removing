function EMG_rm_viewnoise(varargin)
% EMG_rm_viewnoise(Periods,savedir,isnorm,nFFT,savespec)
% 
% view the effect of EMG removal in specragram.
% Inputs: 
%   Periods: the periods to check. or automatically select one.
%   useT: if use time as unit of x axis, default: true.
%   sfactor: rescale factor, if not given, will be chosen by the program. 
%   savedir: the path to save the figures.
% 
% Error contact: chen at biologie.uni-muenchen.de
% 
% Last Modified: 05.12.2019.

%%
[FileName,Period,useT,sfactor,savedir] = DefaultArgs(varargin, {[],[],true,-1,pwd});
if isempty(FileName)
    a = dir('*.EMG_rm.sh*.mat');
    FileName = a(1).name;
end
load(FileName);
FileBase = FileName(1:(find(FileName=='.',1,'first')-1));
a = dir([FileBase,'.EMG_Cluster.mat']);
load(a(1).name)


mlfp = memmapfile(LFPfile,'Format','int16');
data_range = [par.nChannels, length(mlfp.Data)/par.nChannels];
clear mlfp
if isempty(Period)
    EMG_Prd = StartEnd1d(EMG_thrd);
    use_prd = find(diff(EMG_Prd,1,2)>200 & EMG_Prd(:,1)>5000,1,'first');
    if isempty(use_prd)
        t_interval = 100;
        nnn = 0;
        while (nnn<50) && isempty(use_prd)
            use_prd = find(diff(EMG_Prd,1,2)>(t_interval*.9^nnn) & EMG_Prd(:,1)>5000,1,'first');
            nnn=nnn+1;
        end
    end
    Period = fix(mean(EMG_Prd(use_prd,:))) + [-par.lfpSampleRate  par.lfpSampleRate] ;
end
nsh = length(denoise_shank);

mlfp = memmapfile(LFPfile,'Format',{'int16',data_range,'x'});
m = memmapfile([FileBase, '.lfpd'],'Format',{'int16',data_range,'x'});
if useT
    Xtitle = 'time (s)';
else
    Xtitle = 'Samples';
end


for  ksh = 1:nsh
    nsh = figure;
    tmp_shank = denoise_shank(ksh);
    try
        emg = memmapfile([FileBase, '.sh',num2str(tmp_shank),'.emg'],'Format',{'int16',[1 data_range(2)],'x'});
    catch
        emg = memmapfile([FileBase, '.emg'],'Format',{'int16',[1 data_range(2)],'x'});
    end
    
    HP = par.AnatGrps(tmp_shank).Channels+1;
    
    opf1 = @(x)(bsxfun(@plus,x,-HP(:)'));
    if sfactor<0
        v = double(mlfp.Data.x(HP(fix(end/2)),:));
        sfactor = 2^(ceil(log2(std(v)))); % heuristic
    end
    % the original signal
    s1=subplot(2,1,1);
    
    tmp = Period(1):Period(2);
    if useT
        tmp_t = tmp/par.lfpSampleRate;% unit: S
    else
        tmp_t = tmp;
    end
    tmp_x = double(mlfp.Data.x(HP,tmp))';
    plot(tmp_t,opf1(tmp_x/sfactor),'k')
    hold on
    plot(tmp_t,EMG_thrd(tmp)*10-min(HP)+1,'g+')
    
    ylabel('channels')
    title('Before')
    axis tight
    
    % the cleaned signal
    s2=subplot(2,1,2);
    tmp_x = double(m.Data.x(HP,tmp))';
    EMGc = 5+zscore(double(emg.Data.x(tmp)))*10;
    p1 = plot(tmp_t,opf1(tmp_x/sfactor),'r');
    hold on
    p2 = plot(tmp_t,EMGc-min(HP)+1,'b');
    axis tight
    legend([p1(1), p2],{'dlfp', 'EMG'})
    title('After')
    xlabel(Xtitle)
    linkaxes([s1 s2],'xy')
    savefig(nsh,sprintf('%s/%s.sh.%d.EMG_viewnoise.fig',savedir,FileBase,tmp_shank))
end
