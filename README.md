# Welcome (i.e. what does this do?)

This is a MATLAB toolbox for Gaussian mixture model based spike-sorting.  The model is fit using a split-and-merge expectation maximization procedure. A more detailed README is on the way soon!

Depencies:  [markolab toolbox](https://github.com/jmarkow/markolab), Image Processing Toolbox, Statistics Toolbox.

This has been tested on Max OS X 10.9 using MATLAB R2010A (32-bit).

# Quickstart

This toolbox assumes that your data is organized into a samples x trials double-precision Matrix.  For example, if your matrix is 10e3 rows x 100 columns, then the script assumes there are 10e3 time samples and 100 trials (1 trial is fine).      To sort your data, try 

```
>>cluster=spikoclust_sort(mydata,'fs',24e3,'align_method','min','auto_clust',1,'freq_range',[400 6e3]);
```

This will attempt to spike sort the data organized into the matrix mydata (sampling rate 24000 Hz).  The waveforms will be filtered between 400 and 6e3 Hz, upsampled using cubic splines, aligned to their minima (use 'max' for maxima or 'com' for center-of-mass), noise whitened and sorted automatically. To sort using a GUI, set the option auto_clust to 0.  Note that all parameters are passed as parameter/value pairs (see the help documentation for spikoclust_sort for a detailed description of all options).

The output is a structure that contains: the spike waveforms, the spike times (in samples), the inter-spike-intervals, instantaneous firing rate and some clustering statistics.   

  















