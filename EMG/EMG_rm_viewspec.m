function EMG_rm_viewspec(varargin)
% EMG_rm_viewspec(Channel,Periods,savedir,isnorm,nFFT,WinLength,savespec,plot_coh)
% 
% view the effect of EMG removal in specragram.
% Inputs: 
%   Channel: the channels to compare. If not given, will select the one
%            with the largest ripple band power.
%   Periods: the periods to check.
%   savedir: the path to save the figures.
%   isnorm: whether to plot the normalized spectrum for lfp and dlfp. 
%           defualt: true
%   nFFT: nfft to compute spectrum
%   WinLength: length of moving window
%   savespec: if one wants to save the computed spectrum
%
% Tips: play with nFFT and WinLength for the long files. 
% 
% Error contact: chen at biologie.uni-muenchen.de
% 
% Last Modified: 01.12.2019.

%%
[Channel,Periods,savedir,isnorm,nFFT,WinLength,savespec,plot_coh] = DefaultArgs(varargin, {[],[],pwd,false,[],[],true,true});

a = dir('*.EMG_rm.sh*.mat');
FileName = a(1).name;
load(FileName)
FileBase = FileName(1:(find(FileName=='.',1,'first')-1));
a = dir('*.EMG_Cluster.mat');
load(a(1).name)
if isnorm 
    NORM = 'normalized';
else
    NORM = [];
end
if isempty(Periods)
    a = dir('*.sts.*');
    nP = length(a);
    usePeriod = true(nP,1);
    Periods = cell(nP,1);
    for k = 1:nP
        Periods{k} = load(a(k).name);
        if size(Periods{k},2)~= 2
            Periods{k} = [];
            usePeriod(k) = false;
        end
        P_title{k} = a(k).name((find(a(k).name=='.',1,'last')+1):end);
    end
    Periods = Periods(usePeriod);
    P_title = P_title(usePeriod);
    nP = length(Periods);
elseif ~iscell(Periods)
    tmp = cell(1);
    tmp{1} = Periods;
    Periods = tmp;
    nP = length(Periods);
    P_title{1} = 'Period';
else
    nP = length(Periods);
end
if plot_coh
    CNAME = 'COH';
else
    CNAME = 'CSD';
end
mlfp = memmapfile(LFPfile,'Format','int16');
data_range = [par.nChannels, length(mlfp.Data)/par.nChannels];
clear mlfp
%% DATA PROCESSING
for ksh = 1:length(denoise_shank)
    tmp_shank = denoise_shank(ksh);
    HP = par.AnatGrps(tmp_shank).Channels-par.AnatGrps(1).Channels(1)+1;
    nch = length(HP);
    if exist(sprintf('%s/%s.sh.%d.EMG_Spec.mat',savedir,FileBase,tmp_shank),'file')
        load(sprintf('%s/%s.sh.%d.EMG_Spec.mat',savedir,FileBase,tmp_shank), 'y2','f','t')
    else
        y2 = cell(nch,1);
        emg = memmapfile([FileBase, '.emg'],'Format',{'int16',[1 data_range(2)],'x'});
        EMG = double(emg.Data.x);
        EMG = EMG(:) - mean(EMG);
        EMG=WhitenSignal(EMG,[],[],AW{1}.armodel);
        for n = 1:nch
            fprintf('\rShank:%d, Channel: %d, ...', ksh,n)
            mlfp = memmapfile(LFPfile,'Format',{'int16',data_range,'x'});
            lfp = double(mlfp.Data.x(HP(n),:))';
            lfp = bsxfun(@minus,lfp,mean(lfp));%zscore(lfp);
            mlfp = memmapfile([FileBase, '.lfpd'],'Format',{'int16',data_range,'x'});
            dlfp = double(mlfp.Data.x(HP(n),:))';
            dlfp = bsxfun(@minus,dlfp,mean(dlfp));%zscore(dlfp);
            
            if isfield(armodel,'ARmodel')
                lfp=WhitenSignal(lfp,[],[],armodel.ARmodel);
                dlfp=WhitenSignal(dlfp,[],[],armodel.ARmodel);
            else
                lfp=WhitenSignal(lfp,[],[],armodel);
                dlfp=WhitenSignal(dlfp,[],[],armodel);
            end
            [tmp_y2, f, t, ~, ~,~] = mtcsdlong_E([EMG,lfp,dlfp],nFFT,par.lfpSampleRate, WinLength);
            y2{n}=tmp_y2(:,f<400,[1 5 9 7]);
            clear tmp_y2
        end
        f = f(f<400);
        clear mlfp lfp dlfp EMG
        if savespec
            save(sprintf('%s/%s.sh.%d.EMG_Spec.mat',savedir,FileBase,tmp_shank),'y2','f','t','-v7.3')
        end
    end
    nt = length(t);
    if isempty(Channel)
        tmp = zeros(nch,1);
        for k  =1:nch
            tmp(k) = mean(mean(y2{k}(:,f>100 & f<300,2)));
        end
        [~,ch] = max(tmp);
        fprintf('\nSelect channel %d.\n',HP(ch))
    elseif length(Channel)>1
        ch = Channel(ksh);
    else
        ch = Channel;
    end
    %% AVAERAGE CROSS-SPECTRAL DENSITY CROSS STATES
    nclm = 4;
    nrows = nP+3;% ceil(nP/nclm)
    nf = length(f);
    EMGs = zeros(nf,nP);
    EMGch = zeros(nf,nP);
    kf = figure;
    if isnorm
        opf_p = @(x)bsxfun(@rdivide,x,sqrt(sum(x.^2)));
    else
        opf_p = @(x)(x);
    end
    ax1 = ax_subplots(nP, 3, [.05 .35 .65 .6],[.1 .1]);
    % For depth profile:
    for k = 1:nP
        tmp_t = false(nt,1);
        tmp_Period = Periods{k}/par.lfpSampleRate;
        for n = 1:size(tmp_Period,1)
            tmp_t(t<tmp_Period(n,2) & t>=tmp_Period(n,1)) = true;
        end
        % RAW DATA
        axes(ax1(k,1))
        % subplot(nrows,nclm,(k-1)*nclm +1)
        tmp = zeros(nf,nch);
        for n = 1:nch
            tmp(:,n) = mean(abs(y2{n}(tmp_t,:,2)),1);
        end
        c = pcolor(f, HP, opf_p(tmp)');
        c.EdgeColor='none';
        ylabel(P_title{k})
        colorbar
%         if k<nP
%             set(gca,'XTick',[])
%         end
        
        % CLEAN DATA
        axes(ax1(k,2))
%         subplot(nrows,nclm,(k-1)*nclm +2)
        tmp = zeros(nf,nch);
        for n = 1:nch
            tmp(:,n) = mean(abs(y2{n}(tmp_t,:,3)),1);
        end
        c = pcolor(f, HP, opf_p(tmp)');
        c.EdgeColor='none';   
        colorbar
        if k<nP
            set(gca,'XTick',[])
        end
        
        % EMG-CLEAN CSD
        axes(ax1(k,3))
%         subplot(nrows,nclm,(k-1)*nclm +3)
        tmp = zeros(nf,nch);
        for n = 1:nch
            tmp(:,n) = mean(abs(y2{n}(tmp_t,:,4)),1);
        end
        c = pcolor(f, HP, tmp');
        c.EdgeColor='none';
        colorbar
        if k<nP
            set(gca,'XTick',[])
        end
        
        EMGs(:,k) = mean(abs(y2{ch}(tmp_t,:,1)),1);
        tmp_ = mean(abs(y2{ch}(tmp_t,:,3)),1);
        if plot_coh
            EMGch(:,k) = mean(abs(y2{ch}(tmp_t,:,4)),1)./sqrt(EMGs(:,k)'.*tmp_);
        else
            EMGch(:,k) = mean(abs(y2{ch}(tmp_t,:,4)),1);
        end
        
    end
    erase_ticks(ax1);
    axes(ax1(nP,1))
%     subplot(nrows,nclm,1)
    title(['lfp ', NORM])
    axes(ax1(nP,2))
%     subplot(nrows,nclm,2)
    title(['dlfp ', NORM])
    axes(ax1(nP,3))
%     subplot(nrows,nclm,3)
    title('EMG-dlfp')
    linkaxes(ax1,'xy')
    
    ax2 = ax_subplots(1,2,[.75 .35 .25 .6]);
    axes(ax2(1))
    %     l(1) = subplot(nrows,nclm,nclm);
    plot(f,EMGch)
    lg = legend(P_title);
    if plot_coh
        title(['EMG-pyr ',CNAME])
    else
        title(['EMG-pyr ',CNAME])
    end
    
    axes(ax2(2))
%     l(2) = subplot(nrows,nclm,2*nclm);
    plot(f,EMGs)
    lg.Location = 'best';
    xlabel('frequency (Hz)')
    title('EMG Power')
    %     xlabel('frequency (Hz)')
    if ~plot_coh
        erase_ticks(ax2);
        linkaxes(ax2,'xy')
    end
    axis tight

% For single channel
%     for k = 1:nP
%         subplot(nrows,nclm,k)
%         tmp_t = false(nt,1);
%         tmp_Period = Periods{k}/par.lfpSampleRate;
%         for n = 1:size(tmp_Period,1)
%             tmp_t(t<tmp_Period(n,2) & t>=tmp_Period(n,1)) = true;
%         end
%         plot(f, sq(mean(y2(tmp_t,:,[1 5 9 4 7]),1)))
%         title(P_title{k})
%     end
%     legend('EMG','lfp','dlfp','EMG-lfp','EMG-dlfp')
    %% THE SPECTROGRAM
    nt = length(t);
    tt = 1:500;
    ax3 = ax_subplots(3,1,[.02 .03 .98 .3]);
    axes(ax3(1))
%     subplot(nrows,1,nrows-2)
    c = pcolor(t(tt),f(f<500),sq(y2{ch}(tt,f<500,1,1))');
    c.EdgeColor='none';
    ylabel('EMG')
    colorbar
    xlabel('time(s)')
    
    axes(ax3(2))
%     subplot(nrows,1,nrows-1)
    c = pcolor(t(tt),f,sq(y2{ch}(tt,:,2))');
    c.EdgeColor='none';
    ylabel('lfp')
    colorbar
    
    axes(ax3(3))
%     subplot(nrows,1,nrows)
    c = pcolor(t(tt),f,sq(y2{ch}(tt,:,3))');
    c.EdgeColor='none';
    ylabel('dlfp')
    colorbar
    erase_ticks(ax3);
    savefig(kf,sprintf('%s/%s.sh.%d.EMG_viewspec%s.%s.fig',savedir,FileBase,tmp_shank,NORM,CNAME))
end
