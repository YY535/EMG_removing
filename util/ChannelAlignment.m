function idx = ChannelAlignment(x,varargin)
% function idx = ChannelAlignment(x,[nsample])
% Input: 
%   x: data nt x nch
%   nsample: numbers of samples to compute the similarity between channels,
%           defualt is 5000. 
if nargin>1
    nsample = varargin{1};
else
    nsample = 5000;
end
[nt, nch] = size(x);
if nt < nch 
    x = x';
    [nt, nch] = size(x);
end
z = linkage(x(randperm(nt, nsample),:)');
k = figure;
[~,~,idx] = dendrogram(z);
close(k)
end



% Not Necessary, you might need to condider this if you are really going to
% the sparse matrix realm. 
% % cov_x = cov(x);
% % % [~,idx] = sort(diag(cov_x));
% % idx = symrcm(cov_x);
% % cov_x = cov_x(:,idx);
% % cov_x = cov_x(idx,:);
% % 
% % cor_x = cov_x;
% % for k = 1:size(x,2)
% %     cor_x(k,:) = cor_x(k,:)/sqrt(cor_x(k,k));
% %     cor_x(:,k) = cor_x(:,k)/sqrt(cor_x(k,k));
% % end
