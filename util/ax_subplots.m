function ax = ax_subplots(nr, nclm, varargin)
% ax = ax_subplots(nr, nclm, tot_Position,margins,tiling)
% This function generate subplot objects with given position. 
% Inputs: 
%   nr: number of rows
%   nclm: number of columns 
%   tot_Position: size and location of the frame which you want to put all
%                 the subplots.four-element vector of the form 
%                 [left bottom width height].
%   margins: the ratio of the margins to subplots, e.g., [.1 .1]
%   tiling: the tiling methods: 1. filling in
%                               2. left some margin to the bottom. 
% 
% Related functions:
% erase_ticks.m, EMG_rm_viewnoise.m, EMG_rm_viewspec.m
% 
% This function is a part of the EMG_removing toolbox but serves general
% purpose.
% 
% errors contact: chen at biologie.uni-muenchen.de
% 
% last modified: 10.12.2019


[tot_Position,margins,tiling] = DefaultArgs(varargin, {[.05 .05 .9 .9], [.1 .1],2});
if length(margins)<2
    margins = repmat(margins, 1,2);
end

[x_ind,x_len] = set_ind(tot_Position(3),margins(1),nclm,tiling,tot_Position(1));
[y_ind,y_len] = set_ind(tot_Position(4),margins(2),nr,tiling,tot_Position(2));

for k = 1:nr 
    for n = 1:nclm
        ax(k,n) = axes('Position',[x_ind(n),y_ind(k), x_len,y_len]);
    end
end

function [x_ind,x_len] = set_ind(tot_len,margin,n,tiling_method,baseline)
switch tiling_method
    case 1
        x_len = tot_len/(margin*(n-1)+n);
        x_ind = [0:(n-1)]*x_len*(1:margin)+baseline;
    case 2
        x_ind = linspace(0,tot_len,n+1);
        x_len = (x_ind(2)-x_ind(1))/(1+margin);
        x_ind = x_ind(1:n)+(x_len*margin/2)+baseline;
    otherwise
end
