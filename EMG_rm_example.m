%  Example of EMG noise detection and removing: 
%% Data Loading
load lfp1.mat %  this is the first continuse period of LFP recording, at

%% EMG Priod Selection
tmp_t = 9e6:9.5e6; % example period with high muscle tone
[x, Ws, As, EMG_au] = EMG_rm(lfp.data, lfp.sr);
%% Validation 
% Ach Signals and movement related stuffs
% plot of spectrum, Ach recording and ripple power: 
load ACh_NREM.mat
opf1 = @(x)(bsxfun(@rdivide,x, std(x))); % not perticularily important.
opf = @(x)(opf1(bsxfun(@minus, x, min(x))));
end_t = length(lfp.data)/1000*5;
figure;
plot([1:end_t]/5, bsxfun(@plus, [10 20], opf([ACh_NREM.signals.ACh(1:end_t), ACh_NREM.signals.SWRpower(1:end_t)])));
hold on
plot([1:10:length(lfp.data)]/1000, opf(EMG_au(1:10:end)));
plot(ACh_NREM.SWRs/ACh_NREM.lfpSampRate, zeros(size(ACh_NREM.SWRs))-20,'r+')
plot(ACh_NREM.peaksACh/Achfs, zeros(size(ACh_NREM.peaksACh))-25,'k+')
legend('ripple peak power', 'Ach','IC', 'Ach peak', 'ripple peak')
xlim([0 1.4e4])