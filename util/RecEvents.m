function [t, EvtName] = RecEvents(EventName,Fs)
% [t, EvtName] = RecEvents(EventName,Fs)
% function to load the events saved in the .evt file.
% 
% Inputs: 
%       EventName: .evt file.
%       Fs: scalar, data sampling frequency (e.g. 1250 Hz).
% Outputs:
%       t: nx1 vector, event time in sample.
%       EvtName: cell, discription of the events.
% 
% Error contact: chen at biologie.uni-muenchen.de

evtrec = importdata(EventName);
nevt = length(evtrec);
t = nan(nevt,1);
EvtName = cell(nevt,1);
for k = 1:nevt
    first_space = find(evtrec{k}==' ',1,'first');
    t(k) = fix(str2double(evtrec{k}(1:(first_space-1)))*Fs/1000);%
    EvtName{k} = evtrec{k}((first_space+1):end);
end