function EMG_au_rescale(varargin)
% to fix the rescaling problem. 
% please go to the directory or give the directory and run the function.
% inputs: 
%   savedir: string indicate the directory that one is gonna save the data. 
%            the default is the current directory. 
%   EMG_rm_file: e.g., 'FileBase.EMG_rm.sh*.mat'
% 

owd = pwd;
[savedir,EMG_rm_file] = DefaultArgs(varargin, {[],[]});

if ~isempty(savedir)
    cd(savedir)
end
cwd = pwd;
savedir = cwd;
FileBase =  cwd((find(cwd=='/',1,'last')+1):end);
if isempty(EMG_rm_file)
    EMG_rm_files = dir('*EMG_rm.sh*.mat');
else
    EMG_rm_files.name=EMG_rm_file;
end
nfiles = length(EMG_rm_files);

for kk = 1:nfiles
    EMG_rm_file = EMG_rm_files(kk).name;
    load(EMG_rm_file, 'Ws','As','AW','EMG_au','armodel','sug_period','par','denoise_shank','LFPfile')
    
    nPeriod = size(sug_period,1);
    nshank = length(denoise_shank);
    nt = sug_period(end);
    
    try
        load(EMG_rm_file,'scaling_factor')
    catch
        scaling_factor = zeros(nPeriod,nshank);
    end
    
    EMGFileNames = cell(nshank,1);
    
    for n = 1:nshank
        myData = zeros(1,nt);
        
        for k = 1:nPeriod
            if ~isempty(As{k,n})
                nch = length(As{k,n});
                scaling_factor(k,n) = sign(sum(As{k,n}))*sqrt(As{k,n}'*As{k,n}/nch);
            end
            myData(sug_period(k,1):sug_period(k,2)) = scaling_factor(k,n)*EMG_au{k,n};
        end
        
        myData = int16(myData);
        EMGFileName = sprintf('%s%s.sh%d.emg',savedir,FileBase,denoise_shank(n));
        if ~exist(EMGFileName,'file')
            fileID = fopen(EMGFileName,'w');
            fwrite(fileID, myData,'int16');
            fclose(fileID);
            clear myData
        end
        EMGFileNames{n} = EMGFileName;
    end
    
    save(EMG_rm_file, 'Ws','As','AW','EMG_au','armodel','sug_period','par','denoise_shank','LFPfile',...
        'scaling_factor','EMGFileNames','savedir')
    
    fprintf('\n %s fixed',EMG_rm_file)
end
cd(owd)