function [Anew,Wnew,dx]=Denoise_Matrix_Fitting(EMG_au,x,varargin)
% fit the EMG model in new data x according to a given EMG_au activity, 
% s.t. the component generating the same EMG_au activity would be removed.
% 
% Inputs:
%   EMG_au: unit variance EMG trace
%   x   : the new data
%   tol : noise tolerance. default: 1e-10
% Output:
%   Anew: the new mixing matrix
%   Wnew: the new unmixing matrix
%   dx  : the denoised data
% 
% Error contact: chen at biologie.uni-muenchen.de
% 


[tol] = DefaultArgs(varargin, {1e-10});
% To Garentee the Data Property
EMG_au = EMG_au(:);

[nch,T]=size(x);
if nch>T
    x=x';
    [nch,T]=size(x);
end

use_Period=EMG_au~=0;
EMG_au(use_Period) = EMG_au(use_Period)/sqrt(var(EMG_au(use_Period)));
cX = cov(x(:,use_Period)');
[u,s,~]=svd(cX);

Wnew=(x(:,use_Period)*EMG_au(use_Period))'/(cX+tol*eye(nch))/T;
Anew=u*s*u'*Wnew';
dx=x-Anew*EMG_au;