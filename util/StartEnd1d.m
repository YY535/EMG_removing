function se = StartEnd1d(x)
% se = StartEnd1d(x)
% detecting the startign and the endpoitn of a binary trace. 
% Input: x: the binary trace, or nonentry is consider as true. 
% Output: se: the starting and the endding points in samples. 
%           in shape of nsample*2, [start end]
%           we extrapolate the whole trace by 0s outside. 
% error contact: chen at biologie.uni-muenchen.de

se = diff([0;abs(x(:))>0;0]);
se = [find(se>0) find(se<0)-1];
