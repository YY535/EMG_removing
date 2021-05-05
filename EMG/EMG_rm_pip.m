function EMG_rm_pip(FileBases,denoise_shank,varargin)
% EMG_rm_pip(FileBases,denoise_shank,[rm_linenoise,line_thrd])
% The pipeline to process all the sessions on the current dictionary. 
% 
% Inputs: 
%   FileBases: a cell of all the session names, defualt: sessions in
%              the current dictionary. 
%   denoise_shank: the shanks for denoising. [shanks] or a cell define
%                   which shank to denoise in each session, respectively.
%       silence_periods: only compute and remove the EMG from non_silence 
%                       periods, give the periods you don't want to remove 
%                       anything or their filenames. default: false 
%       sp_loadingfuns: functions to load the discarded periods. 
%   Optional:
%       rm_linenoise: if remove line noise. default: true
%       line_thrd: power ratio between line band and other band. 
%                  defualt: 1.8, I'm being conservative here.

% Notice: all the other recording except for the channels in the
% denoise_shank will be directly copied to the .lfpd file.  
% 
% Related functions: 
% EMG_rm_main.m, EMG_rm_view.m, EMG_rm_viewspec.m, EMG_rm_viewnoise.m, 
% EMG_rm_report.m...
% 
% Error contact: chen at biologie.uni-muenchen.de
% 
% Last Modified: 12.12.2019.

[silence_periods,sp_loadingfuns,rm_linenoise,line_thrd,denoise_frequency_lowerbound] = DefaultArgs(varargin, {0,[],true,1.8,0});

HeadDir = pwd;
if isempty(FileBases)||nargin<1
    aa = dir;
    FileBases = cell(1);
    n = 1;
    for k = 3:length(aa)
        if aa(k).isdir
            FileBases{n} = aa(k).name;
            n=n+1;
        end
    end
        
end
nSession = length(FileBases);

if isempty(denoise_shank)
    Error('Please give the linear shank number.')
elseif length(denoise_shank)<nSession
    tmp = cell(nSession,1);
    for k = 1:nSession
        tmp{k} = denoise_shank;
    end
    denoise_shank = tmp;
end
for k  =1:nSession
    tmp_session = FileBases{k};
    fprintf('\n\nStart denoising %s...\n\n',tmp_session)
    tmp_dir = sprintf('%s/%s', HeadDir,tmp_session);
    cd(tmp_dir)
    EMG_rm_main(tmp_session,0,tmp_dir,denoise_shank{k},[],silence_periods,sp_loadingfuns,[],rm_linenoise,line_thrd)
    fprintf('\n\nFinished denoising %s...\n\n',tmp_session)
    cd(HeadDir)
end