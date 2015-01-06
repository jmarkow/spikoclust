# Welcome (i.e. what does this do?)

This is a MATLAB toolbox for Gaussian mixture model based spike-sorting.  The model is fit using a split-and-merge expectation maximization procedure. 

===

# Table of Contents

1. [Requirements](#requirements)
2. [Installation](#installation)
3. [Quick Start](#quick-start)
4. [Script options](#script-options)

###Requirements

This has been tested on Max OS X 10.10 using MATLAB R2010A (32-bit) and MATLAB R2012b (64-bit).  The following toolboxes must be installed: [markolab toolbox](https://github.com/jmarkow/markolab), Image Processing Toolbox, Statistics Toolbox.

###Installation

Download a release or simply clone the repository.  Add all directories in the repository (including sub-directories) to your MATLAB PATH.  To keep things clean, exclude all directories that start with .git (these are used internally by git, MATLAB doesn't need to see them).

### Quick start

This toolbox assumes that your data is organized into a samples x trials double-precision Matrix.  For example, if your matrix is 10e3 rows x 100 columns, then the script assumes there are 10e3 time samples and 100 trials (1 trial is fine).      To sort your data, try 

```
>>cluster=spikoclust_sort(mydata,24e3,'align_method','min','auto_clust',1,'freq_range',[400 6e3]);
```

This will attempt to spike sort the data organized into the matrix mydata (sampling rate 24000 Hz).  The waveforms will be filtered between 400 and 6e3 Hz, upsampled using cubic splines, aligned to their minima (use 'max' for maxima or 'com' for center-of-mass), noise whitened and sorted automatically. To sort using a GUI, set the option auto_clust to 0.  Note that all parameters (aside from the first two inputs, data and the sampling rate) are passed as parameter/value pairs (see the help documentation for spikoclust_sort for a detailed description of all options).

The output is a structure that contains: the spike waveforms, the spike times (in samples), the inter-spike-intervals, instantaneous firing rate and some clustering statistics.

To verify that everything is working correctly, download "simulation_1.mat" from [here](http://www2.le.ac.uk/departments/engineering/research/bioengineering/neuroengineering-lab/simulations/simulation-1.mat). This is some simulated data from Rodrigo Quiroga's lab.  Now run the following commands,

```
>>[cluster,spikeless]=spikoclust_sort(data,24e3,'auto_clust',1,'freq_range',[400 6e3],'align_method','max');
>>spikoclust_autostats(cluster,stats);
```

You should see the following three windows. 

![Stats window A](/spikoclust_demo_1.png?raw=true "Stats window A") ![Stats window B](/spikoclust_demo_2.png?raw=true "Stats window B") ![Stats window C](/spikoclust_demo_3.png?raw=true "Stats window C")

###Script options

Detailed list to be added soon, see the MATLAB help documentation for spikoclust_sort.m 

  















