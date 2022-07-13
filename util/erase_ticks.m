function ax = erase_ticks(ax,varargin)
% ax = erase_ticks(ax,direction)
% This function erases the ticks of subplots when they share the same xaxis
% or y axis. It keeps the yaxis labels for the left most column and xaxis
% for the first or the last row given by the user. 
% 
% Inputs: 
%   ax: the axes objects or subplots that you want to erase the ticks.
%       Notice that the axes object here is given as a matrix. 
%   direction: the 'first'(or use 1) or the 'last'(or use -1) row. 
%              defult: 'first'.
% 
% Related functions:
% ax_subplots.m, EMG_rm_viewnoise.m, EMG_rm_viewspec.m
% 
% This function is a part of the EMG_removing toolbox but serves general
% purpose.
% 
% errors contact: chen at biologie.uni-muenchen.de
% 
% last modified: 10.12.2019

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