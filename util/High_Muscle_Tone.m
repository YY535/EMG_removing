function selectedprd = High_Muscle_Tone(x,varargin)
% function selectedprd = High_Muscle_Tone(x, ... 
%                 (prctil_thr, high_pass_freq, sampling_frequency, Window))
% 
% To detect the period of higher muscle tone, based on high frequency
% synchronization for pairwise randomly selected channels. 
% Currently, the deection is based on the coherence between channels in
% high frequencies. 
% Contact: chen at biologie.uni-muenchen.de

[prctil_thr, hp_freq, FS, Win] = DefaultArgs(varargin, {[], 100, 1250, []});
    
[nt, nch] = size(x);
if nt < nch 
    x = x';
    [nt, nch] = size(x);
    if isempty(prctil_thr)
        if nt <4e6
            prctil_thr = 80;
        else
            prctil_thr = (2e6/nt)*100;
        end
    end
end
% Window to compute synchronization
if isempty(Win)
    % 20 ms window. 
    lwin = fix(20*FS/1000);
    Win = ones(lwin, 1)/lwin;
end
% prepare the data
if nch >5 
    % only need to randomly select 5 channels
    tmp_x = x(:,randperm(nch,5));
else
    tmp_x = x;
end
nch = size(tmp_x,2);
% Highpass
tmp_x = ButFilter(tmp_x,4,hp_freq/(FS/2),'high');
% homebrew cov, including instanteneous auto coherence.
tmp_x = zscore(tmp_x); % bsxfun(@minus, x, mean(x,1));
cov_x = nan(nt, nch, nch);
cmpr_n = (nch+1)*nch/2;
for k = 1:nch
    for n = k:nch
        cov_x(:,k,n) = conv(tmp_x(:,k).*tmp_x(:,n), Win, 'same');
    end
end
opf_r = @(x)reshape(x,nt,[]);
opf_c = @(x)bsxfun(@ge, x, prctile(x,prctil_thr));
selectedprd = (opf_c(opf_r(cov_x))*ones(nch*nch,1))>fix(cmpr_n*.8);
return % EOF