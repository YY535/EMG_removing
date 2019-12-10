function ax = erase_ticks(ax,varargin)
% ax = erase_ticks(ax,direction)
[direction] = DefaultArgs(varargin, {'first'});
[nr, nclm] = size(ax);
if isempty(direction)
    direction = 'first';
end
if ~ischar(direction(1))
    switch direction
        case 1
            direction = 'first';
        case -1
            direction = 'last';
        otherwise
    end
end
switch lower(direction)
    case 'first'
        for k = 1:nr
            for n = 1:nclm
                if k >1
                    ax(k,n).XTick = [];
                end
                if n >1
                    ax(k,n).YTick = [];
                end
            end
        end
    case 'last'
        for k = 1:nr
            for n = 1:nclm
                if k <nr
                    ax(k,n).XTick = [];
                end
                if n >1
                    ax(k,n).YTick = [];
                end
            end
        end
    otherwise
end