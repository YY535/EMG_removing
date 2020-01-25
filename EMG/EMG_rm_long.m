function varargout = EMG_rm_long(x, varargin)
% function [x, Ws, As, EMG_au, AW, armodel] = EMG_rm_long(x,[SamplingRate, 
%                                       high_pass_freq, EMG_thrd, if_rm_mean, 
%                                       armodel, cmp_method, down_sample,
%                                       isave, save_range, FileName, savedir])
% function to remove the EMG noise in a long period. 
% EMG detected by high frequency (>high_pass_freq) correlated activity
% among channels. 
% Inputs: 
%   x: data, nt x nch.
%   Optional: 
%       SamplingRate: in Hz, default: 1000 Hz.
%       rm_linenoise: if remove line noise. default: true
%       line_thrd: power ratio between line band and other band. 
%                  defualt: 1.8, I'm being conservative here.
%       high_pass_freq: beginning of frequency to detect the muscle tone
%                       default: 100 Hz.
%       EMG_thrd: the thresholded data of the high EMG periods. Precomputed
%       by the function of EMG_Cluster.m.
%       if_rm_mean: if return the decenterized signal, default: true.
%       armodel: the armodel or the name of the armodel file. if not given,
%               we woud use whiten signal to compute a corresponding
%               ARmodel. 
%       cmp_method: methods to compute the EMG components. Whiten('w') or
%                   Highpass&Whiten('hw', a bit more stable). defualt: 'hw'
%       down_sample:down sample the data to compute the periods.
%                   default:3.
%       numOfIC: adjusted according to number of channels.
%       isave: if save to a file. defualt: false (when the entry is empty)
%              if the entry is not empty, the file is saved to a file
%              given by "Filename".
%       save_range: the Period in the original whole dataset, used to write
%                   the cleaned data into a larger binary file. 
%                   formate: cell(2,1): 
%                           save_range{1}: [Start_Sample  End_Sample;
%                                           Start_Channel End_Channel];
%                           save_range{2}: [total_channels, total_samples],
%                   if not given, save as defualt: [1 nt;1 nch]. 
%       FileName: defualt:  current_dictionary_name.EMG_rm.mat and
%                           current_dictionary_name.lfpc.
%       savedir:    the directory that one is gonna save the data. the
%                   default is the current directory.
%       save_together: if save ICs from all the chunks together. 
%                       default: true
%       
% Outputs:
%   x: EMG removed signal.
%   Ws: unmixing vector of the EMG component.
%   As: loading of the EMG component.
%   EMG_au: normalized EMG component in arbitrary unit.
%   AW: everything about the ICA 
%   armodel: the armodel used. (use the same armodel for all the chunks)
%
%   Note:   function without output arguments saves the EMG component and
%           the binary lfp data to files. For large data I use memmapfile.
% 
% Related functions: 
% EMG_Cluster.m, SREaffineV.m, EMG_rm_main.m.
% This function is a part of the EMG_removing toolbox.
%  
% Error contact: chen at biologie.uni-muenchen.de
% 
% Last Modified: 01.12.2019.


%% COMPLETE VARIABLES

cwd=pwd; %
[LFPfs, rm_linenoise, line_thrd, high_pass_freq, EMG_thrd, if_rm_mean,armodel,cmp_method,down_sample,numOfIC,isave,save_range,FileName,savedir,save_together] = DefaultArgs(varargin, {1000, true, 2, 100, [], true,[],'hw',3,0,false,[],[],[],true});
[nt, nch] = size(x);
if nt < nch 
    istr = true;
    x = x';
    [nt, nch] = size(x);
else
    istr = false;
end
if isempty(EMG_thrd)
    error('Please Check the high EMG detection file.')
end
if ~numOfIC
    numOfIC = max(fix(.7*nch), min(10,nch));
end
if ~isempty(savedir)
    savedir = [savedir,'/'];
end
selectedprd = EMG_thrd;

%% PREPARE DATA

mx = mean(x);
x = bsxfun(@minus, x, mx);

% Whitening data
if isempty(armodel)
    armodel = dir('*.WhitenSignal.ARmodel.lfp.mat');
    if isempty(armodel)
        [wx,armodel]=WhitenSignal(x,[],[],[],1);
        tmp_ar = armodel;
    else
        armodel = load(armodel(1).name);
        wx=WhitenSignal(x,[],[],armodel.ARmodel);
        tmp_ar = armodel.ARmode;
    end
elseif isstr(armodel)
    armodel = load(armodel);
    tmp_ar = armodel.ARmodel;
    wx=WhitenSignal(x,[],[],tmp_ar);
elseif isfield(armodel,'ARmodel')
    tmp_ar = armodel.ARmodel;
    wx=WhitenSignal(x,[],[],tmp_ar);
else 
    try
        tmp_ar = armodel;
        wx=WhitenSignal(x,[],[],armodel);
    catch 
        warning('Something wrong with the AR model you give, try whiten without given model.')
        [wx, armodel]=WhitenSignal(x,[],[],[],1);
        tmp_ar = armodel;
    end
end
AW.armodel = tmp_ar;
% opf_A = @(x)(bsxfun(@rdivide,x,std(x)));
chmap = 1:nch;
opf_a = @(x)bsxfun(@rdivide,x,sqrt(sum(x.^2)));
opf_Aa = @(x)(bsxfun(@rdivide,x,SREaffineV(chmap,x)));
opf_A = @(x)opf_Aa(opf_a(x));
% SREaffineV use affine to fit, and compute the variance accordinglty.
% Accounting for the linear leaking from other areas.  
%% REMOVE LINE-NOISE
if rm_linenoise
    [A_rm_line,W_rm_line,A_line,W_line,power_ratio,thrd] = EMG_rm_linenoise(wx,line_thrd,LFPfs);
    AW.A_rm_line = A_rm_line;
    AW.W_rm_line = W_rm_line;
    AW.A_line = A_line;
    AW.W_line = W_line;
    AW.power_ratio = power_ratio;
    AW.thrd = thrd;
    if ~isempty(A_rm_line)
        signals = x*W_rm_line';
        for n = 1:size(signals,2)
            tmp_id = find(sign(signals(1:(end-1),n)).*sign(signals(2:end,n))<=0, 1,'first');
            signals(1:tmp_id,n) = 0;
            tmp_id = find(sign(signals(1:(end-1),n)).*sign(signals(2:end,n))<=0, 1,'last');
            signals(tmp_id:end,n) = 0;
        end
        x = x - signals*A_rm_line';
        wx = wx - wx*W_rm_line'*A_rm_line'; % won't affect the results. 
        fprintf('\n Line noise component removed...\n')
    end
end
%% EMG COMPONENTS AND ACTIVITIES
AW.usewb = false;
% Components from the high frequency.
switch lower(cmp_method) % 
    case 'hw'
        hx = ButFilter(x,4,high_pass_freq/(LFPfs/2),'high');
        [Ah, Wh] = fastica(hx(selectedprd,:)', 'numOfIC', numOfIC,'verbose','off');
        % [~, EMG_comp] = max(abs(sum(opf_A(Ah))));
        [~,rod] = sort(abs(sum(opf_A(Ah))),'descend');
        % EMG_au(:,1) = (x*Wh(EMG_comp,:)');
        % Components from all over
        [Ax, Wx] = fastica(Wh(rod,:)*wx(1:down_sample:end,:)','verbose','off');
        A = Ah(:,rod)*Ax;
        W = Wx*Wh(rod,:);
        if (sum(abs(sum(opf_A(A)))>nch)>1) || sum(abs(sum(opf_A(A)))>(2*nch))<1
            fprintf('recompute...\n')
            wx=WhitenSignal(wx,[],[],tmp_ar);
            nn = 1;
            while  ((sum(abs(sum(opf_A(A)))>nch)>1) || sum(abs(sum(opf_A(A)))>(2*nch))<1)&&(nn<5)
                % too many flat components.
                if ((sum(abs(sum(opf_A(Ah)))>nch)>1) || sum(abs(sum(opf_A(Ah)))>(2*nch))<1)
                    [Ah, Wh] = fastica(hx(selectedprd,:)','verbose','off'); % 'numOfIC', numOfIC,
                end
                [Ax, Wx] = fastica(Wh*wx(nn:down_sample:end,:)','verbose','off');
                A = Ah*Ax;
                W = Wx*Wh;
                nn = nn+1;
                fprintf('\r%d in %d...\n',nn,5)
            end
            
            fprintf('\n\n')
        end
        if sum(abs(sum(opf_A(A)))>nch)<1
            [Awb, Wwb] = fastica(x(1:down_sample:end,:)', 'verbose','off');% se
            AW.Awb = Awb;
            AW.Wwb = Wwb;
            if sum(abs(sum(opf_A(Awb)))>nch)>1
                A = Awb;
                W = Wwb;
                AW.usewb = true;
                fprintf('\r Using the wide band.\n')
            end
        end
        AW.Ah = Ah;
        AW.Ax = Ax;
        AW.Wh = Wh;
        AW.Wx = Wx;
    case 'w' % Use Whiten alone
        [A, W] = fastica(wx(1:down_sample:end,:)', 'numOfIC', numOfIC,'verbose','off');% selectedprd
        if sum(abs(sum(opf_A(A)))>nch)<1
            [Awb, Wwb] = fastica(x(1:down_sample:end,:)', 'verbose','off');% se
            AW.Awb = Awb;
            AW.Wwb = Wwb;
            if sum(abs(sum(opf_A(Awb)))>nch)>1
                A = Awb;
                W = Wwb;
                AW.usewb = true;
                fprintf('\r Using the wide band.\n')
            end
        end
    otherwise
        fprintf('Please using hw or w. ')
end

[~, EMG_comp] = max(abs(sum(opf_A(A))));
As = A(:,EMG_comp);
Ws = W(EMG_comp,:);
EMG_au = (x*Ws');%.*selectedprd; 
% (:,n) is to see the behavior of the highpassed resulted EMG.
% not need in the final varsion. 

AW.A = A;
AW.W = W;

%% SMOOTHING THE ENDS OF EMG TRACES.

% first_cross
if abs(EMG_au(1))<(median(abs(EMG_au))*1e-3)
    first_cross = 1;
else
    first_cross = find(sign(EMG_au(1:(end-1)).*EMG_au(2:end))<0,1,'first');
end
if first_cross>1250
    first_cross = find(EMG_au<=(median(abs(EMG_au))*1e-2),1,'first');
end
EMG_au(1:first_cross)=0;
AW.zero_first = first_cross;% zeros line at first.

% last_cross
if abs(EMG_au(end))<(median(abs(EMG_au))*1e-3)
    last_cross = nt;
else
    last_cross = find(sign(EMG_au(1:(end-1)).*EMG_au(2:end))<0,1,'last');
end
if (nt-last_cross)>1250
    last_cross = find(EMG_au<=(median(abs(EMG_au))*1e-2),1,'last');
end
EMG_au((last_cross+1):end)=0;
AW.zero_last = nt - last_cross+1;% zeros line at last.

x = x - EMG_au*As';% (:,end)
if ~if_rm_mean
    x = bsxfun(@plus, x, mx);
end
if istr % transpose back. 
    x = x';
end

%% OUTPUT OR SAVE DATA

if nargout>1
    varargout{1} = x;
    varargout{2} = Ws;
    varargout{3} = As;
    varargout{4} = EMG_au;
    varargout{5} = AW;
    varargout{6} = tmp_ar;
else 
    isave = true;
end

% function without output arguments saves the EMG component and the binary
% lfp data to files. 
if isempty(save_range)
    save_range =    [1 nt;...
                    1 nch];
end

if isave
    if isempty(FileName)
        if ~isempty(savedir)
            cwd = savedir;
        else
            cwd = pwd;
        end
        FileName =  cwd((find(cwd=='/',1,'last')+1):end);
        LFPfile = sprintf('%s%s.lfpd',savedir,FileName);
        EMGfile = sprintf('%s%s.emg',savedir,FileName);
        if exist(LFPfile,'file')
            warning(sprintf('%s already exist! Please check!\n Now saving to the %s.new files.', FileName, FileName))
            FileName = [FileName,'.new'];
            LFPfile = sprintf('%s%s.lfpd',savedir,FileName);
        end
        if ~save_together
            save(sprintf('%s.EMG_rm.t%d-%d.ch%d-%d.mat',FileName,save_range(1,:),save_range(2,:)), 'AW','EMG_au','armodel')
        end
        fileID = fopen(LFPfile,'w');
        fwrite(fileID, int16(x'),'int16');
        fclose(fileID);
        
        fileID = fopen(EMGfile,'w');
        fwrite(fileID, int16(EMG_au'),'int16');
        fclose(fileID);
    else
        if ~save_together
            save(sprintf('%s.EMG_rm.t%d-%d.ch%d-%d.mat',FileName,save_range{1}(1,:),save_range{1}(2,:)), 'AW','EMG_au','armodel')
        end
        LFPfile = sprintf('%s%s.lfpd',savedir,FileName);
        EMGfile = sprintf('%s%s.emg',savedir,FileName);
        m = memmapfile(LFPfile,'Format',{'int16',save_range{2},'x'},'Writable',true);
        m.Data.x(save_range{1}(2,1):save_range{1}(2,2), save_range{1}(1,1):save_range{1}(1,2)) = int16(x');
        clear m
        m = memmapfile(EMGfile,'Format',{'int16',[1 save_range{2}(2)],'x'},'Writable',true);
        m.Data.x(save_range{1}(1,1):save_range{1}(1,2)) = int16(EMG_au');
        clear m
    end
end