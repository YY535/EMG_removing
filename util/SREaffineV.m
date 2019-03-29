function r = SREaffineV(x,y,varargin)
% r = SREaffineV(x,y,varargin) 

if nargin > 2
    iscenter = varargin{1};
else
    iscenter = true;
end
[nt, ny] = size(y);
if nt < ny
    y = y';
    [nt, ny] = size(y);
end
[nx, nch] = size(x);
if nx ~= nt
    if nt == nch
        x = x';
        nch = nx;
    else
        fprintf('size mismatch')
    end
end
r = zeros(nch,ny);
% opf_x = @(x)(x'*x/(length(x)-1));
for n = 1:ny 
    tmp_y = y(:,n);
for k = 1:nch
    if iscenter
        tmp_x = [ones(nt, 1), x(:,k)];
        tmp_ch = 2;
    else
        tmp_x = x(:,k);
        tmp_ch = 1;
    end
    r(k,n) = norm(tmp_y - tmp_x*((tmp_x'*tmp_x+1e-10*eye(tmp_ch))\(tmp_x'*tmp_y)));
    % norm is optimized 
end
end