function PeriodsAWAKE = LoadAwake(FileName,nSamples,varargin)
% PeriodsAWAKE = LoadAwake(FileName,nSamples,[loadmethods])
% loadmethods is any customer defined loading functions, e.g., 'my_load'. 
%            Please do make sure my_load.m is in your searching path. 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ADDED by Evgeny
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load SWS brain state periods
% modify the loadrangefiles to adapt to your own methods loadmethods
if nargin<3
    loadmethods = 0;
else
    loadmethods = varargin{1};
end

PeriodsSWS = loadrangefiles(FileName, loadmethods);
PeriodsAWAKE = true(nSamples,1);
for k = 1:size(PeriodsSWS,1)
    PeriodsAWAKE(max(PeriodsSWS(k,1),1):PeriodsSWS(k,2)) = false;
end
% %Compute all AWAKE (non-SWS) periods
% clear PeriodsAWAKE per1 per2 per3
% if ~isempty(PeriodsSWS)
%     
%     %first awake period
%     if PeriodsSWS(1,1)>0
%         per1 = [1 PeriodsSWS(1,1)-1];
%     else
%         per1 = [];
%     end
%     %middle awake periods
%     for k=1:size(PeriodsSWS,1)-1
%         per2(k,:) = [PeriodsSWS(k,2)+1  PeriodsSWS(k+1,1)-1];
%     end
%     %last awake period
%     if PeriodsSWS(end,2) < nSamples
%         per3 = [PeriodsSWS(end,2)+1 nSamples];
%     else
%         per3 = [];
%     end
%     %merge thhem all
%     PeriodsAWAKE = [per1; per2; per3];
%     
% else %if ~isempty(PeriodsSWS)
%     %No SWS periods, the whole session is AWAKE
%     PeriodsAWAKE = [1 nSamples];
% end %if ~isempty(PeriodsSWS)

% BrainStateDuration_sec =  sum(diff(PeriodsAWAKE,1,2)) / par.lfpSampleRate;

%DEBUGGING PLOT
% figure;
% subplot(211); cla; hold on
% PlotIntervals(PeriodsSWS)
% PlotIntervals(PeriodsAWAKE, 'bars')

% %Durations of individual AWAKE periods
% diff(PeriodsAWAKE,1,2) / par.lfpSampleRate/60


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PeriodsSWS = loadrangefiles(FileName,loadmethods)
if loadmethods
    eval(sprintf('PeriodsSWS = %s(FileName);',loadmethods));
else
    PeriodsSWS = load(FileName);
end