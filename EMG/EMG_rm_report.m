function EMG_rm_report(varargin)
% EMG_rm_report(FileName,Channels,savedir)
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
[FileName,Channels,savedir] = DefaultArgs(varargin, {[],[],pwd});
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
FileBase = FileName{1}(1:(find(FileName{1}=='.',1,'first')-1));
opf_nA = @(x)bsxfun(@rdivide,x,sqrt(sum(x.^2)));
opf_A = @(x)(bsxfun(@rdivide,x,SREaffineV(1:size(x,1),x)));
opf_C = @(x)abs(sum(opf_A(x)));
a = dir('*.EMG_Cluster.mat');
load(a(1).name)
nclm = 5;

for kk = 1:nfile
    %%
    tmp_File = FileName{kk};
    load(tmp_File)
    [nchunk,nshank] = size(AW);
    nshk = [];
    if isempty(Channels)
        ch = zeros(nshank,1);
        for n = 1:nshank
            ch(n) = par.AnatGrps(denoise_shank(n)).Channels(fix(end/3*2));
        end
    elseif length(Channels)<nshank
        ch = zeros(nshank,1);
        for n = 1:nshank
            ch(n) = par.AnatGrps(denoise_shank(n)).Channels(Channels);
        end
    else
        ch = Channels;
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
        lfp=WhitenSignal(lfp,[],[],armodel.ARmodel);
        dlfp=WhitenSignal(dlfp,[],[],armodel.ARmodel);
    else
        lfp=WhitenSignal(lfp,[],[],armodel);
        dlfp=WhitenSignal(dlfp,[],[],armodel);
    end
    
    %%
    % spatial criteriaa
    for n=1:nshank
        tmp_sh = denoise_shank(n);
        nshk = figure;% (1);clf
        for k = 1:nchunk
            subplot(nchunk,nclm,1+nclm*(k-1))
            plot(par.AnatGrps(tmp_sh).Channels,opf_nA(AW{k}.A),'Color',[.4 .4 .4])
            hold on
            plot(par.AnatGrps(tmp_sh).Channels,opf_nA(As{k}),'r', 'LineWidth',2)
            axis tight
            
            subplot(nchunk,nclm,2+nclm*(k-1))
            tmp_c = opf_C(AW{k}.A);
            [~, EMG_comp] = max(tmp_c);
            boxplot(tmp_c)
            hold on
            plot(tmp_c(EMG_comp),'ro')
            
            subplot(nchunk,nclm,3+nclm*(k-1))
            hx = log(abs(hilbert(ButFilter(EMG_au{k},4,EMG_par.hp_freq/(EMG_par.lfpSamplingRate/2),'high'))));
            [h1, b1] = hist(hx(EMG_thrd(sug_period(k,1):sug_period(k,2))),100);
            [h2, b2] = hist(hx(~EMG_thrd(sug_period(k,1):sug_period(k,2))),100);
            plot(b1, h1/sum(h1))
            hold on
            plot(b2, h2/sum(h2))
            axis tight
            
            subplot(nchunk,nclm,4+nclm*(k-1))
            if isfield(armodel,'ARmodel')
                [yo, fo] = mtcsdfast([WhitenSignal(EMG_au{k}*As{k}(ch(n)),[],[],armodel.ARmodel), lfp(sug_period(k,1):sug_period(k,2),n), dlfp(sug_period(k,1):sug_period(k,2),n)],...
                    [],par.lfpSampleRate);
            else
                [yo, fo] = mtcsdfast([WhitenSignal(EMG_au{k}*As{k}(ch(n)),[],[],armodel), lfp(sug_period(k,1):sug_period(k,2),n), dlfp(sug_period(k,1):sug_period(k,2),n)],...
                    [],par.lfpSampleRate);
            end
            plot(fo,abs(sq([yo(:,1,1),yo(:,2,2),yo(:,3,3)])))
            axis tight
            
            subplot(nchunk,nclm,5+nclm*(k-1))
            plot(fo,abs(sq([yo(:,1,2),yo(:,1,3),yo(:,2,3)])))
            axis tight
            
        end
        subplot(nchunk,nclm,1)
        title('norm. loading')
%         subplot(nchunk,nclm,2)
%         title('inv.SRE')
        subplot(nchunk,nclm,3)
        title('log(high frq power)')
        legend('higher EMG','otherwise')
        subplot(nchunk,nclm,4)
        title('Spectrum')
        legend('EMG','raw','clean')
        subplot(nchunk,nclm,5)
        title('Cross Specral Density')
        legend('EMG-raw','EMG-clean','raw-clean')
        savefig(nshk,sprintf('%s/%s.sh.%d.EMG_rm_report.fig',savedir,FileBase,tmp_sh))
    end
end
