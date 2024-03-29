function [A_line,W_line,A,W,power_ratio,line_thrd] = EMG_rm_linenoise(x,varargin)
% [A_line,W_line,A_rm_line,W_rm_line] = EMG_rm_linenoise(wx,[line_thrd,lfpSampleRate])
% 
% This function is intend to remove the line noise here. To reconstruct the
% signal without EMG noise, do: x = x-x*W_line'*A_line' (x in nt*nch)
% 
% Inputs: 
%   x: data, nt x nch.
%   Optional: 
%   line_thrd: power ratio between line band and other band. 
%                  defualt: 1.8, I'm being conservative here.
%   lfpSampleRate: default: 1250
%   show_line_noise: default: false
% 
% Outputs:
%   A_line: the mixing vector of line noise component
%   W_line: the unmixing vector of line noise component
%   A: the mixing matrix of the line noise component
%   W: the unmixing matrix of the line noise component
%   power_ratio
%   line_thrd
% 
% Related functions: 
% EMG_rm_long.m, EMG_rm_main.m.
% This function is a part of the EMG_removing toolbox. But could serve
% general purpose. 
%  
% Error contact: chen at biologie.uni-muenchen.de
% 
% Last Modified: 11.12.2019.

[line_thrd,lfpSampleRate,show_line_noise] = DefaultArgs(varargin, {1.8, 1250,false});

[nt,nch] = size(x);
if nt<nch
    warning('Wrong dimension of the data matrix. Shifted. ')
    x = x';
    [nt,nch] = size(x);
end
[A, W] = fastica(x', 'verbose','off');

x = x*W';
yo = [];
for k = 1:size(x,2)
    [yo(:,k), fo] = mtcsdfast(x(:,k),[],lfpSampleRate);
end
%% LINE NOISE
line_noise_band = false(size(fo));
for k = 1:floor(lfpSampleRate/50)
    line_noise_band(fo<(k*50 +5) & fo>(k*50 -5)) = true;
end

power_ratio = mean(yo(line_noise_band,:))./mean(yo(~line_noise_band,:));
ids = find(power_ratio>line_thrd);
A_line = A(:,ids);
W_line = W(ids,:);

%% SHOW RESULTS
if show_line_noise
    opf_AA = @(x)(bsxfun(@rdivide,x,sqrt(sum(x.^2))));

    figure;
    clf;
    sp(1)=subplot(1,2,1);
    scaler = length(fo)/5/mean(sum(yo,2));
    pp1 = plot(fo,bsxfun(@plus, yo*scaler, 1:size(x,2)));
    hold on
    tmp_sum = sum(A.^2);
    pp2= plot(tmp_sum/max(tmp_sum)*fo(end),1:size(x,2));
    pp3 = plot(power_ratio/max(power_ratio)*100,1:length(power_ratio));
    if ~isempty(ids)
        pp4 = plot(power_ratio(ids)/max(power_ratio)*100, ids,'x');
        legend([pp2,pp3,pp4], {'tot power', 'power ratio', 'line noise component'})
    else
        legend([pp2,pp3], {'tot power', 'power ratio'})
    end
    sp(2)=subplot(1,2,2);
    plot(1:nch,bsxfun(@plus, opf_AA(A), 1:size(x,2)))
    axis tight
    grid on
    linkaxes(sp,'y')
end