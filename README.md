# EMG_removing
MATLAB toolbox for removing high frequency EMG artifacts from the multichannel extracellular recording with ICA.

This is generally an overcomplete ICA problem since the number of the potential signal sources is much larger than the recording sites. Lower frequency physiological signal is going to affect the separation.

The current version use the spectrum whitening to enphasize the high frequency EMG tone (`EMG_rm_main.m`).
In an early version we use highpassed (>100 Hz) data to find the EMG component. The related functions (`ÈMG_rm.m`) are still left there.

## Examples:

```matlab
% running the pipeline
addpath(genpath('/path/to/EMG_removing'))
cd('/path/to/the/head/dictionary/of/your/sessions')
denoise_shanks = 1;% or [shank1, shank2....], each shank is running separately 
EMG_rm_pip([],denoise_shanks)
```

or go to your session and try:

```matlab
% running for one session
addpath(genpath('/path/to/EMG_removing'))
cd('/path/to/the/session')
FileBase = session_name;
denoise_shanks = 1;
denoise_frequency_lowerbound = 10;
rm_line_noise = true; 
silence_periods = false;
sp_loadingfuns = [];% use load
EMG_rm_main(FileBase,denoise_frequency_lowerbound,[],denoise_shanks,[],silence_periods,sp_loadingfuns,rm_line_noise)
# Groups:
ngroups=1;
Grps=cell(ngroups,1);
Grps{1}=[1,2];
keep_old_lfpd = true;
EMG_rm_main_group(FileBase,Grps,denoise_frequency_lowerbound,keep_old_lfpd,[],[],silence_periods,sp_loadingfuns,rm_line_noise)
```

- Notice the function `EMG_rm_main.m` or `EMG_rm_long.m` by defualt automatically remove the line noise component. If you don't want to do this, set this parameter to `false`. 
- Sometimes people prefer to remove noise in the awake periods and left sleeping periods unchanged. To do this you would need to give the periods you don't want to touch in `silence_periods`. In this case, if you have any function other than `load` to import the periods, specify it in `sp_loadingfuns` with the function name and remember to add that to your path. 
- The EMG components are predominantly fitted in the higher frequency data. Removing EMG in the wide-band might affect the slower dynamics. You can choose to keep the original data up until `denoise_frequency_lowerbound` (Hz). EMG activity beyond this frequency would be removed. 
- `EMG_rm_main_group.m` allows one to fit the components simultaneously for all the shanks in a group, groups is given in cell. The `keep_old_lfpd` allows one to keep the old denoing results and only work on the new groups. 

## Check the Results:
The cleaned signals will be saved in `.lfpd` files and the EMG activity in `.emg`. The EMG signals (`EMG_au`) and the EMG components `AW.As` is saved in `FileBase.EMG_rm.mat`. To check the cleaned signal, use:

```matlab
cd('/path/to/your/sessions')
EMG_rm_view()
EMG_rm_view([], [t_beginning, t_end]) % to check arbitrary period
PYR_Channel = 37;% choose the channel to visualize the effect. 
EMG_rm_report([],PYR_Channel) % or we'll use the channel with the largest ripple power.
EMG_rm_viewspec()
% play with the nfft and the window length to compute the spectrum when you have a long file. 
EMG_rm_viewnoise()
```
When computing the properties of slower frequency signals, e.g., delta wave phase and power, it is recommended to use the original signal `x_orig`. The original signal could also be reconstructed by the EMG activities. In case of useing channel `ch_delta` (`x` is the cleaned data):

```
ch_delta ;% the selected channel
x_d = x(:,ch_delta);% the cleaned signal
x_r = AW.As(ch_delta)*EMG_au(:);% the EMG signal
x_orig = x_d + x_r;
```

NB: please make sure you've removed all the previously generated denoing files when you try to redo the denoise. 

This package is entirely based on matlab codes. Errors please contact: chen at biologie.uni-muenchen.de
