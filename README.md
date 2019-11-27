# EMG_removing
Remove high frequency EMG artifacts from the extracellular recording with ICA.

This is generally an overcomplete ICA problem since the number of the potential signal sources is much larger than the recording sites. Lower frequency physiological signal is going to affect the separation.

The current version use the spectrum whitening to enphasize the high frequency EMG tunes.
In an early version we use highpassed (>100 Hz) data to find the EMG component. The related functions are still left there.

example:

```matlab
% running the pipeline
addpath(genpath('/path/to/EMG_removing'))
cd('/path/to/the/head/dictionary/of/your/sessions')
denoise_shank = 1;
EMG_rm_pip([],denoise_shank)
```

The cleaned signals will be saved in .lfpd files.

This package is entirely based on matlab codes. Errors please contact: chen at biologie.uni-muenchen.de
