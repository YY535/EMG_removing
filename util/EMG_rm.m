function [x, Ws, As, EMG_au] = EMG_rm(x, varargin)
% function [x, Ws, As, EMG_au] = EMG_rm(x,[SamplingRate, ifrearrange, ...
%                                       prctil_thr, high_pass_freq, Window, 
%                                       Number_of_Sample_test_simi, if_rm_mean])
% function to remove the EMG noise. 
% EMG detected by high frequency (>high_pass_freq) correlated activity
% among channels. 
% Inputs: 
%   x: data, nt x nch
%   optional: 
%       SamplingRate: in Hz, default: 1000 Hz
%       ifrearrange: if rearrange the channels, if it is true, rearrange
%                   according to the similarity of channels. 
%       prctil_thr: threshold to detect high correlation period
%       high_pass_freq: beginning of frequency to detect the muscle tone
%                       default: 100 Hz
%       Window: window length to detect high corr
%       Number_of_Sample_test_simi: number of samples used to compute the
%                                   similarity between channels
%       if_rm_mean: if return the decenterized signal
% Outputs:
%   x: EMG removed signal
%   Ws: unmixing vector of the EMG component
%   As: loading of the EMG component
%   EMG_au: normalized EMG component in arbitrary unit.
% 
% ------- 25.02.2019 ------- % 
% Contact: chen at biologie.uni-muenchen.de

[LFPfs, rearrange, prctil_thr, high_pass_freq, Window, nsimi, if_rm_mean] = DefaultArgs(varargin, {1000, true, 99, 100, [], 1e4, true});
[nt, nch] = size(x);
if nt < nch 
    istr = true;
    x = x';
    nch = nt;
else
    istr = false;
end

if rearrange
    idx = ChannelAlignment(x,nsimi);
    selectedprd = High_Muscle_Tone(x(:,idx(1:fix(nch/5):end)), prctil_thr, high_pass_freq, LFPfs, Window);
else
    selectedprd = High_Muscle_Tone(x,  prctil_thr, high_pass_freq, LFPfs, Window);
end
% better select ~5 channels.
mx = mean(x);
x = bsxfun(@minus, x, mx);
x = ButFilter(x,4,high_pass_freq/(LFPfs/2),'high');
[A, W] = fastica(x(selectedprd,:)');
opf_A = @(x)(bsxfun(@rdivide,x,std(x)));
[~, EMG_comp] = max(abs(sum(opf_A(A)))); 
As = A(:,EMG_comp);
Ws = W(EMG_comp,:);
EMG_au = x*Ws';
x = x - EMG_au*As';
if ~if_rm_mean
    x = bsxfun(@plus, x, mx);
end
if istr
    x = x';
end