# EMG_removing
Removing high frequency EMG artifacts from the multichannel extracellular recording with ICA.

This is generally an overcomplete ICA problem since the number of the potential signal sources is much larger than the recording sites. Lower frequency physiological signal is going to affect the separation.

The current version use the spectrum whitening to enphasize the high frequency EMG tone (`EMG_rm_main.m`).
In an early version we use highpassed (>100 Hz) data to find the EMG component. The related functions (`ÈMG_rm.m`) are still left there.

example:

```matlab
% running the pipeline
addpath(genpath('/path/to/EMG_removing'))
cd('/path/to/the/head/dictionary/of/your/sessions')
denoise_shank = 1;
EMG_rm_pip([],denoise_shank)
```

or go to your session and try:

```matlab
% running for one session
addpath(genpath('/path/to/EMG_removing'))
cd('/path/to/the/session')
FileBase = session_name;
denoise_shank = 1;
rm_line_noise = true; 
EMG_rm_main(FileBase,[],denoise_shank,[],rm_line_noise)
```

- Notice the function `EMG_rm_main.m` or `EMG_rm_long.m` by defualt automatically remove the line noise component. If you don't want to do this, set this parameter to `false`. 

The cleaned signals will be saved in `.lfpd` files. The EMG signals (`EMG_au`) and the EMG components `AW.As` is saved in `FileBase.EMG_rm.mat`. To check the cleaned signal, use:

```matlab
cd('/path/to/your/sessions')
EMG_rm_view()
PYR_Channel = 37;% choose the channel to visualize the effect. 
EMG_rm_report([],PYR_Channel) % or we'll use the channel with the largest ripple power.
EMG_rm_viewspec()
% play with the nfft and the window length to compute the spectrum when you have a long file. 
EMG_rm_viewnoise()
```

This package is entirely based on matlab codes. Errors please contact: chen at biologie.uni-muenchen.de
