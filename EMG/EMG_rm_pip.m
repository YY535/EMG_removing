function EMG_rm_pip(FileBases,denoise_shank)
% EMG_rm_pip(FileBases)
% The pipeline to process all the sessions on the current dictionary. 
% 
% Inputs: 
%   FileBases: a cell of all the session names, defualt: sessions in
%              the current dictionary. 
%   denoise_shank: the shanks for denoising. [shanks] or a cell define
%                   which shank to denoise in each session, respectively. 
% Notice: all the other recording except for the channels in the
% denoise_shank will be directly copied to the .lfpd file.  
% 
% Related functions: 
% EMG_rm_main.m, EMG_rm_view.m
% 
% Error contact: chen at biologie.uni-muenchen.de
% 
% Last Modified: 27.11.2019.

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
    EMG_rm_main(tmp_session,tmp_dir,denoise_shank)
    fprintf('\n\nFinished denoising %s...\n\n',tmp_session)
    cd(HeadDir)
end