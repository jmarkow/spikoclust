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

This toolbox assumes that your data is organized into a samples x trials x channels double-precision Matrix.  For example, if your matrix is 10e3 rows x 100 columns x 1, then the script assumes there are 10e3 time samples and 100 trials (1 trial is fine). To sort your data, try 

```
>>cluster=spikoclust_sort(mydata,24e3,'align_feature','min','freq_range',[400 6e3]);
```

This will attempt to spike sort the data organized into the matrix mydata (sampling rate 24000 Hz).  The waveforms will be filtered between 400 and 6e3 Hz, upsampled using cubic splines, aligned to their minima (use 'max' for maxima or 'com' for center-of-mass), noise whitened and sorted automatically. To sort using a GUI, set the option auto_clust to 0.  Note that all parameters (aside from the first two inputs, data and the sampling rate) are passed as parameter/value pairs (see the help documentation for spikoclust_sort for a detailed description of all options).

The output is a structure that contains: the spike waveforms, the spike times (in samples), the inter-spike-intervals, instantaneous firing rate and some clustering statistics.

To verify that everything is working correctly, download "simulation_1.mat" from [here](http://www2.le.ac.uk/departments/engineering/research/bioengineering/neuroengineering-lab/simulations/simulation-1.mat). This is some simulated data from Rodrigo Quiroga's lab.  Now run the following commands,

```
>>[cluster,spikeless]=spikoclust_sort(data,24e3,'freq_range',[400 6e3],'align_feature','max');
>>spikoclust_autostats(cluster,spikeless);
```

You should see the following three windows. 

![Stats window A](/spikoclust_demo_1.png?raw=true "Stats window A") ![Stats window B](/spikoclust_demo_2.png?raw=true "Stats window B") ![Stats window C](/spikoclust_demo_3.png?raw=true "Stats window C")

If you want to sort tetrode data, pass a matrix with additional channels (e.g. 10e3 x 100 x 4 would be 10e3 samples by 100 trials by 4 channels).  The first channel is the "master" channel, i.e. spikes will be detected on this channel, and those timepoints will be used to select windows from the other 3 channels.  Windows from the other 3 channels will be appended to the first window before computing the principal components for spike sorting.

###Script options

Yes there are many options, but many of the defaults should be "sensible" for typical spike sorting applications.  The options, with a description and their default settings, are listed in the following table:

| Parameter | Description | Format | Options | Default |
|-----------|-------------|--------|---------|---------|
| `noise_removal` | Method of noise rejection (common average) | string | `car`,`none` | `none` |
| `car_trim` | If using common average, trim N/2 percent of data on each end before taking mean | float | N/A | `40` |
| `freq_range` | Frequency range for filtering | 1-2 element vector of floats | N/A | [400] |
| `filt_type` | Filter class | `high`,`low`,`bandpass` | string | `high` |
| `filt_order` | Filter order | integer | N/A | `3` |
| `filt_name` | Filter type | string | 'butter','ellip','kaiser' | `ellip` |
| `sigma_t` | Detection threshold (multiple of robust standard deviation) | float | N/A | `4` |
| `detect_method` | Detect negative-going spikes, positive-going spikes, or both | string | `n`,`p`,`b` |
| `spike_window` | Size of spike window (in s before and after spike) | 2 element vector of floats (s) | N/A | [.0005 .0005] |
| `interp_f` | Interpolation factor (interpolate spikes by a factor of N) | integer | N/A | `8` | 
| `align_feature` | Feature used for spike re-alignment | string | `min`,`max`,`com` | `min` |
| `jitter` | Limit on number of samples a spike may be moved for re-alignment | integer | N/A | `10` |
| `gui_clust` | Use the GUI assistant | logic | N/A | `1` |
| `clust_check` | Number of spikes to check for | vector of integers | N/A | `[2:8]` |
| `pcs` | Number of principal components to compute/use | integer | N/A | `2` |
| `garbage` | Use a garbage component in the mixture to exclude outliers? | logic | N/A | `1` |
| `smem` | Use split-and-merge EM rather than standard EM? | logic | N/A | `1` |
| `modelselection` | Technique to use for selecting the number of neurons | string | `icl`,`bic`,`aic` | `icl` |
| `maxnoisetraces` | Maximum number of traces to use to compute noise covariance | integer | N/A | `1e6` |
| `noisewhiten` | Enable noise whitening? | logic | N/A | `1` |

*Note that all parameters are passed as parameter/value pairs.


  















