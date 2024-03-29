function [EMG_thrd, EMG_heavy, sug_period] = EMG_Cluster(x,varargin)
% [EMG_thrd, EMG_heavy] = EMG_Cluster(x,[hp_freq, lfpSamplingRate, lwinms, nchunks,
%                                        swin, isave, FileName])
% function to detect large EMG activity period. 
% 
% Inputs: 
%   x: signal in nt x nch
%   Optional inputs:
%       hp_freq: high pass filter to detect high museual tune, defualt: 100
%       lfpSamplingRate: the sampling rate, defualt: 1250 Hz
%       lwinms: window length to calculate high coherence periods. default:
%               20 samples,~16ms
%       nchunks: chunk the data into n truncks, defualt:6. this directly
%           related to the window to compute the EMG artifact.
%           Heuristically should be less than 15 mins, aka. 1.125e6
%           samples.  
%       swin: the minimum length of the scilence window. default: 500 ms
%       isave:  if save to a file. defualt: false (when the entry is empty)
%               if the entry is not empty, the file is saved to a file
%               given by "Filename".
%       FileBase: the data is going to be aved in FileBase.EMG_Cluster.mat.
%                 defualt: current_dictionary_name
%                 
% Outputs: 
%   EMG_thrd: vector indicates the high EMG periods, used in te EMG removing.
%   EMG_heavy: vector indicates the highly coherent periods. The condition is
%           looser than the EMG_thrd. This trace is used to detect the
%           recommended working periods. 
%   sug_period: recommened working periods. Determined by the nchunks, nx2.
% 
% This is the perperocessing function of the EMG_rm_long.m to extract EMG
% activity condensed periods. 
% This function uses the StartEnd1d.m function in ./util/. 
% 
% Related functions: 
% EMG_rm_long, StartEnd1d.
% This function is a part of the EMG_removing toolbox.
% 
% errors contact: chen at biologie.uni-muenchen.de
% 
% last modified: 01.12.2019

%% ARGUEMENT COMPELETION

cwd_name = pwd;
[hp_freq, lfpSamplingRate, lwinms, nchunks, swin, isave, FileBase,savedir] = DefaultArgs(varargin, {100, 1250, 20, 6, [],false, cwd_name((find(cwd_name=='/',1,'last')+1):end),[]});
[nt,nch] = size(x);
if nt<nch
    x = x';
    [nt,~] = size(x);
end
nchunks = fix(nt/min(15*60*lfpSamplingRate, fix(nt/nchunks)));
if isempty(swin)
    swin = fix(lfpSamplingRate/10);
end

%% DATA PREPROCESSING

nch = size(x,2);
tmp_x = butfilter(x,4,hp_freq/(lfpSamplingRate/2),'high');
clear x
tmp_x = zscore(tmp_x); % zscore to compute the coherence instead of covariance. 
% I guess it's not a big deal to use the covariance. 

%% DETECT HIGH EMG PERIODS
% pairwise coherence with sliding window (here I just using the mean). 
cov_x = cell(nch);
lwin = fix(lwinms*lfpSamplingRate/1000);
Win = ones(lwin, 1)/lwin;
for k = 1:nch
    for n = k:nch
        cov_x{k,n} = conv(tmp_x(:,k).*tmp_x(:,n), Win, 'same');
    end
end
Xtrain = cmp_Xtrain(cov_x);
% Xtrain=[ttv_v,cov_v];

% cluster with ICA. 
[A, W] = fastica(Xtrain','verbose','off');
tmp_A = abs(A(1,:))./max(abs(A(2,:)),0);
if tmp_A(2)>tmp_A(1)
    A = A(:,[2 1]);
    W = W([2 1],:);
end
EMG_heavy = abs(Xtrain*W(2,:)')> abs(Xtrain*W(1,:)');
EMG_thrd =(abs(Xtrain*W(2,:)') > prctile(abs(Xtrain(~EMG_heavy,:)*W(2,:)'),99));

%% RECOMMENDED WORKING PERIOD SECTION

% suggested segmant points: 
sug_period = ceil(linspace(1,nt,nchunks+1));% rough section 

% the midpoints of low EMG long periods (LEP) closest to the rough segment
% points.

se = StartEnd1d(~EMG_heavy);
se_cand = se;
se_durations = diff(se,1,2);

se(diff(se,1,2)<min(swin,max(50,prctile(se_durations,50))),:) = [];
se_cand(diff(se_cand,1,2)<min(max(swin-10,1),max(20,prctile(se_durations,10))),:) = [];% 20 ms
for k  =2:(nchunks-1)
    tmp1 = find(se(:,1)<=sug_period(k),1,'last');
    tmp2 = find(se(:,2)>=sug_period(k),1,'first');
    if tmp1==tmp2 
        % when the rough point is covered by LEP
        sug_period(k) = fix(sum(se(tmp1,:))/2);
    elseif abs(sug_period(k)- se(tmp1,2))<abs(sug_period(k)- se(tmp2,1)) 
        % when the rough point is not covered by LEP
        sug_period(k) = fix(sum(se(tmp1,:))/2);
    else
        if ~isempty(tmp2)
            sug_period(k) = fix(sum(se(tmp2,:))/2);
        else
            tmp2 = find(se_cand(:,2)>=sug_period(k),1,'first');
            sug_period(k) = fix(sum(se_cand(tmp2,:))/2);
        end
    end
end
sug_period = unique(sug_period);

% DISCARD SHORT PERIODS
short_period_threshold = 20;
while sum(diff(sug_period)<(short_period_threshold*lfpSamplingRate))
    sug_period(find(diff(sug_period)<(short_period_threshold*lfpSamplingRate))+1)=[];
end

sug_period = [[1; sug_period(2:(end-1))'+1] [sug_period(2:(end-1))';nt]];

%% SAVE RESULTS

if isave
    EMG_par.hp_freq = hp_freq;
    EMG_par.lfpSamplingRate=lfpSamplingRate;
    EMG_par.lwinms=lwinms;
    EMG_par.nchunks=nchunks;
    EMG_par.swin=swin;
    EMG_par.FileBase=FileBase;
    EMG_par.datetime = datetime;
    EMG_par.A = A;
    EMG_par.W = W;
    save(sprintf('%s/%s.EMG_Cluster.mat',savedir,FileBase), 'EMG_heavy','EMG_thrd','sug_period','EMG_par')
end
function Xtrain = cmp_Xtrain(cov_x)
c_idx = find(triu(ones(size(cov_x))));
[~,rids] = sort(rand(length(c_idx),1));
rc_idx = c_idx(rids);
cov_v = cov_x{c_idx(1)}(:);
ttv_v = zeros(length(cov_x{1}),1);
for k = 2:length(c_idx)
    cov_v = cov_v+cov_x{c_idx(k)}(:);
    % total variance to garentee line shape, against the dipole structure.
    % ttv_v = sum(abs(diff(tmp_x(:,rids),1,2)),2);
    ttv_v = ttv_v+ abs(cov_x{rc_idx(k)}(:) - cov_x{rc_idx(k-1)}(:));
end
Xtrain=[ttv_v,cov_v];