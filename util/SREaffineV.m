function varargout = SREaffineV(x,y,varargin)
% r = SREaffineV(x,y,varargin) 
% the "non-flateness" of a vector.
% fit linear model to the data then compute the SRE.
% error contact chen at biologie.uni-muenchen.de
[Shk01, iscenter,keepy] = DefaultArgs(varargin, {-1, true,false});
[nt, ny] = size(y);
if ~keepy
    if nt < ny
        y = y';
        [nt, ny] = size(y);
    end
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
alpha = zeros(nch,ny);
if Shk01(1)>0
    for k = 1:size(Shk01,1)
        tmp = Shk01(k,1):Shk01(k,2);
        r = r+SREaffineV(x(tmp,:),y(tmp,:),[],[],true);
    end
else
    % opf_x = @(x)(x'*x/(length(x)-1));
    noise_vec = nan(nt,ny);
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
            tmp = ((tmp_x'*tmp_x+1e-10*eye(tmp_ch))\(tmp_x'*tmp_y));
            alpha(k,n) = tmp(end);
            noise_vec(:,ny) = tmp_y - tmp_x*tmp;
            r(k,n) = norm(tmp_y - tmp_x*tmp);
            
            % norm is optimized
        end
    end
end
varargout{1} = r;
if nargout>1
    varargout{2} = alpha;
end
if nargout>2
    varargout{3} = noise_vec;
end