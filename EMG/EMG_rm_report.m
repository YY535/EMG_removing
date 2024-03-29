function [AW, Power,axs] = EMG_rm_report(varargin)
% [AW, Power] = EMG_rm_report(FileName,Channels,savedir)
% report the dynamics of EMG components.
% Inputs: 
%   FileName: EMG removing files you want to check. 
%             defualt: the .EMG_rm.sh* 
%   Channels: the channels to compare. 
%   savedir: the path you want to save the figures. 
% 
% Error contact: chen at biologie.uni-muenchen.de
% 
% Last Modified: 01.12.2019.
%%
[FileName,Channels,savedir,plot_coh] = DefaultArgs(varargin, {[],[],pwd,true});
if isempty(FileName)
    a = dir('*.EMG_rm.sh*.mat');
    nfile = length(a);
    FileName = cell(nfile,1);
    for k = 1:nfile
        FileName{k} = a(k).name;
    end
elseif ~iscell(FileName)
    tmp = cell(1);
    tmp{1} = FileName;
    FileName = tmp;
    nfile = length(FileName);
else
    nfile = length(FileName);
end
if plot_coh
    TITLE = 'COH';
else
    TITLE = 'CSD';
end
FileBase = FileName{1}(1:(find((FileName{1}(1:(end-1))=='.') & (FileName{1}(2:end)=='E'),1,'first')-1));
opf_nA = @(x)bsxfun(@rdivide,x,sqrt(sum(x.^2)));
opf_A = @(x)(bsxfun(@rdivide,x,SREaffineV(1:size(x,1),x)));
opf_C = @(x)abs(sum(opf_A(x)));
a = dir([FileBase, '.EMG_Cluster.mat']);
load(a(1).name)
nclm = 7;

for kk = 1:nfile
    %%
    tmp_File = FileName{kk};
    load(tmp_File)
    [nchunk,nshank] = size(AW);
    Power = cell(nshank,1);
    nshk = [];
    if isempty(Channels)
        ch = zeros(nshank,1);
        for n = 1:nshank
            ch(n) = par.AnatGrps(denoise_shank(n)).Channels(fix(end/3*2)) - par.AnatGrps(1).Channels(1) +1;
            chin(n) = par.AnatGrps(denoise_shank(n)).Channels(fix(end/3*2)) - par.AnatGrps(denoise_shank(n)).Channels(1) +1;
        end
    elseif length(Channels)<nshank
        ch = zeros(nshank,1);
        for n = 1:nshank
            ch(n) = repmat(Channels,nshank,1);% par.AnatGrps(denoise_shank(n)).Channels(Channels);
        end
        chin = ch;
    else
        ch = Channels;
        for n  =1:nshank
            if Channels(n)>par.AnatGrps(denoise_shank(n)).Channels(1)
                ch(n)=Channels(n)-par.AnatGrps(1).Channels(1)+1;
                chin(n) = Channels(n) - par.AnatGrps(denoise_shank(n)).Channels(1) +1;
            end
        end
    end
    %%
    mlfp = memmapfile(LFPfile,'Format','int16');
    data_range = [par.nChannels, length(mlfp.Data)/par.nChannels];
    mlfp = memmapfile(LFPfile,'Format',{'int16',data_range,'x'});
    lfp = double(mlfp.Data.x(ch,:))';
    lfp = bsxfun(@minus,lfp,mean(lfp));%zscore(lfp);
    mlfp = memmapfile([FileBase, '.lfpd'],'Format',{'int16',data_range,'x'});
    dlfp = double(mlfp.Data.x(ch,:))';
    dlfp = bsxfun(@minus,dlfp,mean(dlfp));%zscore(dlfp);
    if isfield(armodel,'ARmodel')
        lfp=whitensignal(lfp,[],[],armodel.ARmodel);
        dlfp=whitensignal(dlfp,[],[],armodel.ARmodel);
    else
        lfp=whitensignal(lfp,[],[],armodel);
        dlfp=whitensignal(dlfp,[],[],armodel);
    end
    
    %%
    % spatial criteria
    set_Position = [];% clear the ploting framwork for a new plot.
    FIRST_REMOVED_CHUNK = 0;
    for n=1:nshank
        tmp_sh = denoise_shank(n);
        nshk = figure;% (1);clf
        
        for k = 1:nchunk
            if isfield(AW{k},'A')
                EMG_REMOVED = true;
                if ~FIRST_REMOVED_CHUNK
                    FIRST_REMOVED_CHUNK = k;
                end
            else
                EMG_REMOVED = false;
            end
            % Ploting the Components profile 
            axs(k,1) = subplot(nchunk,nclm,1+nclm*(k-1));
            if k ==1 || isempty(find(set_Position(:,1), 1))
                set_Position(:,1) = axs(1,1).Position([1 3:4]);
            end
            if EMG_REMOVED
                plot(par.AnatGrps(tmp_sh).Channels,opf_nA(AW{k}.A),'Color',[.4 .4 .4])
                hold on
                plot(par.AnatGrps(tmp_sh).Channels,opf_nA(As{k}),'r', 'LineWidth',2)
                axis tight
                axs(k,1).Position([1 3:4]) = set_Position(:,1);
            end
           
            % Ploting the Components Flatness 
            axs(k,2) = subplot(nchunk,nclm,2+nclm*(k-1));
            if k ==1 || isempty(find(set_Position(:,2), 1))
                set_Position(:,2) = axs(1,2).Position([1 3:4]);
            end
            if EMG_REMOVED
                tmp_c = opf_C(AW{k}.A);
                [~, EMG_comp] = max(tmp_c);
                AW{k}.flatness = tmp_c;
                boxplot(axs(k,2),tmp_c);
                hold on
                plot(tmp_c(EMG_comp),'ro')
                axs(k,2).Position([1 3:4]) = set_Position(:,2);
                axs(k,2).Position(2) = axs(k,1).Position(2);
            end
            % Compare the High frequency power propotion in selected periods
            axs(k,3) = subplot(nchunk,nclm,3+nclm*(k-1));
            if k ==1 || isempty(find(set_Position(:,3), 1))
                set_Position(:,3) = axs(1,3).Position([1 3:4]);
            end
            if EMG_REMOVED
            hx = log(abs(hilbert(ButFilter(EMG_au{k},4,EMG_par.hp_freq/(EMG_par.lfpSamplingRate/2),'high'))));
            [h1, b1] = hist(hx(EMG_thrd(sug_period(k,1):sug_period(k,2))),100);
            [h2, b2] = hist(hx(~EMG_thrd(sug_period(k,1):sug_period(k,2))),100);
            plot(b1, h1/sum(h1))
            hold on
            plot(b2, h2/sum(h2))
            axis tight
            axs(k,3).Position([1 3:4]) = set_Position(:,3);
            end
            
            % Spectrum of the EMG component, original and denoised data
            axs(k,4) = subplot(nchunk,nclm,4+nclm*(k-1));
            if k ==1 || isempty(find(set_Position(:,4), 1))
                set_Position(:,4) = axs(1,4).Position([1 3:4]);
            end
            if EMG_REMOVED
                if isfield(armodel,'ARmodel')
                    [yo, fo] = mtcsdfast([whitensignal(EMG_au{k}*As{k}(chin(n)),[],[],armodel.ARmodel), lfp(sug_period(k,1):sug_period(k,2),n), dlfp(sug_period(k,1):sug_period(k,2),n)],...
                        [],par.lfpSampleRate);
                else
                    [yo, fo] = mtcsdfast([whitensignal(EMG_au{k}*As{k}(chin(n)),[],[],armodel), lfp(sug_period(k,1):sug_period(k,2),n), dlfp(sug_period(k,1):sug_period(k,2),n)],...
                        [],par.lfpSampleRate);
                end
                plot(fo,abs(sq([yo(:,1,1),yo(:,2,2),yo(:,3,3)])))
                axis tight
                axs(k,4).Position([1 3:4]) = set_Position(:,4);
            else
                [yo, fo] = mtcsdfast(lfp(sug_period(k,1):sug_period(k,2),n),...
                    [],par.lfpSampleRate);
                plot(fo,abs(sq(yo)))
                axis tight
                axs(k,4).Position([1 3:4]) = set_Position(:,4);
            end
            
            % Coherence of the EMG component vs original/denoised data
            axs(k,5) = subplot(nchunk,nclm,5+nclm*(k-1));
            if k ==1 || isempty(find(set_Position(:,5), 1))
                set_Position(:,5) = axs(1,5).Position([1 3:4]);
            end
            if EMG_REMOVED
                if plot_coh
                    plot(fo,abs(sq([yo(:,1,2)./sqrt(yo(:,1,1).*yo(:,2,2)),yo(:,1,3)./sqrt(yo(:,1,1).*yo(:,3,3)),yo(:,2,3)./sqrt(yo(:,2,2).*yo(:,3,3))])))
                else
                    plot(fo,abs(sq([yo(:,1,2),yo(:,1,3),yo(:,2,3)])))
                end
                axis tight
                Power{n}(:,:,k) = sq(yo(:,[1 5 9 4 7]));
                axs(k,5).Position([1 3:4]) = set_Position(:,5);
            end
            
            % LINE NOISE
            axs(k,6) = subplot(nchunk,nclm,6+nclm*(k-1));
            if k ==1 || isempty(find(set_Position(:,6), 1))
                set_Position(:,6) = axs(1,6).Position([1 3:4]);
            end
            if isfield(AW{k},'A_line')
                plot(par.AnatGrps(tmp_sh).Channels,opf_nA(AW{k}.A_line),'Color',[.4 .4 .4])
                hold on
                if isfield(AW{k},'A_rm_line')
                    if ~isempty(AW{k}.A_rm_line)
                        plot(par.AnatGrps(tmp_sh).Channels,opf_nA(AW{k}.A_rm_line),'r', 'LineWidth',2)
                    end
                end
                axis tight
            end
            axs(k,6).Position([1 3:4]) = set_Position(:,6);
            
            axs(k,7) = subplot(nchunk,nclm,7+nclm*(k-1));
            if k ==1 || isempty(find(set_Position(:,7), 1))
                set_Position(:,7) = axs(1,7).Position([1 3:4]);
            end
            if isfield(AW{k},'power_ratio')
                if ~isempty(AW{k}.power_ratio)
                    boxplot(axs(k,7),AW{k}.power_ratio);
                    hold on
                    plot(axs(k,7),[.5 1.5], AW{k}.thrd*[1 1],'b')
                    axis tight
                end
            end
            axs(k,7).Position([1 3:4]) = set_Position(:,7);
            if EMG_REMOVED
                axs(k,7).Position(2) = axs(k,1).Position(2);
            else
                axs(k,7).Position(2) = axs(k,6).Position(2);
            end
            
        end
        % subplot(nchunk,nclm,1)
        subplot(axs(1,1))
        title('norm. loading')
        % subplot(nchunk,nclm,(nchunk-1)*nclm +1)
        subplot(axs(nchunk,1))
        xlabel('Channel')
        % subplot(nchunk,nclm,2)
        subplot(axs(1,2))
        title('inv.SRE')
        % subplot(nchunk,nclm,3)
        subplot(axs(1,3))
        title('log(high frq power)')
        % subplot(nchunk,nclm,4)
        subplot(axs(1,4))
        title('Spectrum')
        
        % subplot(nchunk,nclm,5)
        subplot(axs(1,5))
        title(TITLE)
        
        % legends
        if FIRST_REMOVED_CHUNK
        subplot(axs(FIRST_REMOVED_CHUNK,3))
        ll = legend('higher EMG','otherwise');
        ll.Location = 'best';
        subplot(axs(FIRST_REMOVED_CHUNK,4))
        ll = legend('EMG','raw','clean');
        ll.Location = 'best';
        subplot(axs(FIRST_REMOVED_CHUNK,5))
        ll = legend('EMG-raw','EMG-clean','raw-clean');
        ll.Location = 'best';
        end
        subplot(axs(1,6))
        title('LineNoiseComp')
        subplot(axs(1,7))
        title('LN power ratio')
        for kn = 1:nclm
            tmp = nan(4,nchunk);
            for k = 1:nchunk
                tmp(:,k) = [axs(k,kn).XLim(1), axs(k,kn).XLim(2), ...
                       axs(k,kn).YLim(1), axs(k,kn).YLim(2)];
            end
            tmp(:,tmp(1,:)==0&tmp(2,:)==1) = [];
            tmp = [min(tmp(1,:)), max(tmp(2,:)),min(tmp(3,:)),max(tmp(4,:))];
            linkaxes(axs(:,kn),'xy');
            erase_ticks(axs(:,kn),'last');
            axes(axs(1,kn));
            axis(tmp);
        end
        subplot(axs(end))
        xlabel(sprintf('%s sh%d',FileBase,tmp_sh))
        savefig(nshk,sprintf('%s/%s.sh.%d.EMG_rm_report.%s.fig',savedir,FileBase,tmp_sh,TITLE))
    end
    save(sprintf('%s/%s.EMG_rm_report.power.mat',savedir,FileBase),'AW','Power','EMG_par')
end
