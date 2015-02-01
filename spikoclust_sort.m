function [cluster spikeless]=spikoclust_sort(EPHYS_DATA,FS,varargin)
%CLI-driven Gaussian mixture model-based spike sorting
%
%	[cluster spikeless]=ephys_spikesort(EPHYS_DATA,varargin)
%
%	EPHYS_DATA
%	data for spike sorting (samples x trials x channels)
%
%	FS
%	sampling rate of data
%	
%	the following may be specified as parameter/value pairs:
%
%		interpolate_f
%		interpolation factor used for spike alignment (cubic splines, default: 8, i.e. upsample by factor of 8)
%
%		sort_f
%		downsample factor after alignment for clustering (default: interpolate_f, i.e. spikes upsampled then downsampled to original fs)
%
%		tetrode_channels
%		channels may act as a putative n-trode (can specify an arbitrary number of channels),
%		defaults to manual cluster cutting
%
%		car_exclude
%		carelectrodes to exclude from noise estimate
%		
%		noise_removal
%		noise rejection method ('car' for common average 'nn' for nearest neighbor, or 'none',
%		default: 'none')
%
%		freq_range
%		vector with two elements to specify the frequency range (one element specifies low pass, default: [500 8e3])
%
%		filt_type
%		filter type (default: 'bandpass', options 'low','high','bandpass')
%
%		filt_order
%		filter order (Butterworth, default: 6)
%
%		gui_clust
%		gui assisted clustering? (default: 0)	
%
%		sigma_t
%		multiple of variance estimate for automatic threshold setting (uses the Quiroga formulation, default: 4)
%
%		modelselection
%		method of choosing the number of components for GMM clustering (default: icl, options, 
%		'BIC' for Bayes Information, 'mml' for minimum message length, 'icl' for BIC-ICL)
%
%		spike_window
%		seconds before and after threshold crossing to store (in seconds, default: [.0005 .0005])
%
%		clust_check
%		vector of number of components to use in GMM clustering (default: [1:10])
%
%		pcs
%		number of principal components to use in GMM clustering (default: 2)
%
%		garbage
%		include a garbage-collecting uniform density in GMM clustering (0 or 1, default: 1)
%
%		smem
%		use split-and-merge algorithm for GMM clustering (0, 1, or 2 for FREE smem, default: 1)
%	
%		maxnoisetraces
%		maximum number of noise traces to use in spike whitening (default: 1e6)
%
%		align_feature
%		method for spike alignment ('min','max', or 'com', default: 'min');
%		
%		noise
%		denoising method ('car' for common average, 'none' for none, default: 'none');
%
%		jitter
%		maximum jitter for spike realignment in samples (default: 10)
%
%		noisewhiten
%		whiten noise using Cholesky decomposition (0 or 1, default: 1)
%
%		wavelet_denoise
%		denoise using wavelets (0 or 1, default: 0)
%
%		decomp_level
%		parameter for wavelet denoising (default: 7)
%	
%	the following outputs are returned by the script:
%
%	cluster
%   	structure that contains clustering results, each field is a cell array, where cell n is the results for cluster n
%
%	spikeless
%	spike data with the spikes removed
%
%
% see also spikoclust_autosort.m,spikoclust_guisort.m,spikoclust_spike_detect.m


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETER COLLECTION %%%%%%%%%%%%%%%%%

% TODO:  remove multi-channel support

if nargin<1
	error('ephysPipeline:suavis:notenoughparams','Need 1 arguments to continue, see documentation');
end

if isvector(EPHYS_DATA)
	EPHYS_DATA=EPHYS_DATA(:);
end

nparams=length(varargin);

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

noise_removal='none'; % none, nn for nearest neighbor, or car for common average
car_exclude=[]; % exclude any channels for common average?

% 300 Hz E high-pass, see Quiroga et al. 2013

freq_range=[400]; % bandpassing <10e3 distorted results, reasoning that >800 Hz is fine for spikes < 1ms long
filt_type='high'; % high,low or bandpass
filt_order=3; % filter order
filt_name='e'; % filter type, e for elliptic and b for butterworth
gui_clust=0; % use GUI?
regularize=.01;

tetrode_channels=[];
sigma_t=4; % multiple of noise estimate for spike threshold (generally 3-4, using Quiroga's method)
jitter=10; % max jitter in samples for spike re-alignment (4 is reasonable
align_feature='min'; % how to align spike waveforms can be min, max or com for center-of-mass
interpolate_f=8; % interpolate factor (fs*interpolate_f)
sort_f=[]; % if empty, downsamples back to original fs (advised to leave empty)
detect_method='n'; % p for pos n for neg b for both

car_trim=40; % common average using the trimmed mean (car_trim/2 is percent cut off from either edge)
decomp_level=7; % wavelet decomposition level (not used unless wavelet_denoise is set to 1)
wavelet_denoise=0; % wavelet denoise?

spike_window=[.0005 .0005]; % window to the left and right of the detected spike time to use for sorting
clust_check=1:8; % number of neurons to check for 
pcs=2; % number of principal components to use
garbage=1; % use a garbage collecting uniform density in the mixture?
smem=1; % split-and-merge or standard EM?

spikeworkers=1; % only used in other packages, can safely ignore
modelselection='icl'; % how to select the number of neurons, 'bic', 'aic', or 'icl' (bic or icl recommended)
maxnoisetraces=1e6; % not recommended to change, upper bound on number of noise traces used for noise whitening
noisewhiten=1; % enable noies whitening?

% remove eps generation, too slow here...

for i=1:2:nparams
	switch lower(varargin{i})
		case 'noise_removal'
			noise_removal=varargin{i+1};
		case 'sigma_t'
			sigma_t=varargin{i+1};
		case 'car_exclude'
			car_exclude=varargin{i+1};
		case 'filt_type'
			filt_type=varargin{i+1};
		case 'filt_name'
			filt_name=varargin{i+1};
		case 'decomp_level'
			decomp_level=varargin{i+1};
		case 'freq_range'
			freq_range=varargin{i+1};
		case 'spikesort'
			spikesort=varargin{i+1};
		case 'gui_clust'
			gui_clust=varargin{i+1};
		case 'align_feature'
			align_feature=varargin{i+1};
		case 'jitter'
			jitter=varargin{i+1};
		case 'interpolate_f'
			interpolate_f=varargin{i+1};
		case 'spike_window'
			spike_window=varargin{i+1};
		case 'sort_f'
			sort_f=varargin{i+1};
		case 'maxnoisetraces'
			maxnoisetaces=varargin{i+1};
		case 'filt_order'
			filt_order=varargin{i+1};
		case 'clust_check'
			clust_check=varargin{i+1};
		case 'car_trim'
			car_trim=varargin{i+1};
		case 'pcs'
			pcs=varargin{i+1};
		case 'garbage'
			garbage=varargin{i+1};
		case 'smem'
			smem=varargin{i+1};
		case 'spikeworkers'
			spikeworkers=varargin{i+1};
		case 'modelselection'
			modelselection=varargin{i+1};
		case 'wavelet_denoise'
			wavelet_denoise=varargin{i+1};
		case 'noisewhiten'
			noisewhiten=varargin{i+1};
		case 'decomp_level'
			decomp_level=varargin{i+1};
		case 'detect_method'
			detect_method=varargin{i+1};
	end
end

if isempty(sort_f)
	sort_f=interpolate_f;
	disp(['Setting sort downsample factor to spike upsample factor:  ' num2str(sort_f)]);
end

interpolate_fs=FS*interpolate_f;
sort_fs=interpolate_fs/sort_f;

[samples,ntrials,ncarelectrodes]=size(EPHYS_DATA);
if ncarelectrodes==1 & strcmp(noise_removal,'car')
	disp('Turning off CAR, number of electrodes is 1');
	noise_removal='none';
	car_exclude=[];
end

channels=1:ncarelectrodes;
proc_data=zeros(samples,ntrials,length(channels));

TIME=[1:samples]./FS; % time vector for plotting

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIGNAL CONDITIONING %%%%%%%%%%%%%%%%

proc_data=spikoclust_denoise_signal(EPHYS_DATA,channels,channels,'method',noise_removal,'car_exclude',car_exclude,'car_trim',car_trim);
proc_data=spikoclust_condition_signal(proc_data,'s','freq_range',...
	freq_range,'filt_type',filt_type,'filt_order',filt_order,'filt_name',filt_name,...
	'wavelet_denoise',wavelet_denoise,'decomp_level',decomp_level);

if length(channels)>1
	tetrode_channels=channels(2:end);
end

if ~isempty(tetrode_channels)
	tetrode_data=spikoclust_denoise_signal(EPHYS_DATA,channels,tetrode_channels,'method',noise_removal,'car_exclude',car_exclude,'car_trim',car_trim);
	tetrode_data=spikoclust_condition_signal(tetrode_data,'s','freq_range',freq_range,'filt_type',filt_type,'filt_order',...
		filt_order,'filt_name',filt_name,'wavelet_denoise',wavelet_denoise,'decomp_level',decomp_level);
else
	tetrode_data=[];
end

clear EPHYS_DATA;

[samples,ntrials,newchannels]=size(proc_data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SPIKE DETECTION %%%%%%%%%%%%%%%%%%%%

disp('Entering spike detection...');
disp(['Alignment method:  ' align_feature]);
disp(['Interpolate FS:  ' num2str(interpolate_fs)]);
% need a noise cutoff...let's start with 3*std or Quiroga's measure

disp(['Electrode ' num2str(channels)])

if ~isempty(tetrode_channels)
	disp(['Will use tetrode channels ' num2str(tetrode_channels) ' for sorting']);
end

% collect spikes

sort_data=cat(3,proc_data(:,:,1),tetrode_data);
nchannels=size(sort_data,3);

spikethreshold=sigma_t*median(abs(sort_data(:,:,1))/.6745);
spikes=spikoclust_spike_detect(sort_data,spikethreshold,'fs',FS,'visualize','n','align_feature',align_feature,...
		'jitter',jitter,'window',spike_window,'method',detect_method);
spikeless{1}=spikoclust_spike_remove(sort_data(:,:,1),spikes);

totalspikes=length(spikes.times);

for i=2:nchannels
	tmp_thresh=sigma_t*median(abs(sort_data(:,:,i))/.6745);
	tmp_spikes=spikoclust_spike_detect(sort_data(:,:,i),tmp_thresh,'fs',FS,'visualize','n','align_feature',align_feature,...
		'jitter',jitter,'window',spike_window,'method',detect_method);
	spikeless{i}=spikoclust_spike_remove(sort_data(:,:,i),spikes)
end

disp([ num2str(totalspikes) ' total spikes']);
clear sort_data;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp(['Channel ' num2str(channels)]);

cluster.parameters.fs=FS;
cluster.parameters.interpolate_fs=interpolate_fs;
cluster.parameters.sort_fs=sort_fs;
cluster.parameters.threshold=spikethreshold;
cluster.parameters.tetrode_channels=tetrode_channels;
cluster.parameters.spike_window=spike_window;
cluster.parameters.align_feature=align_feature;

if noisewhiten
	disp('Noise-whitening spikes...');
	spikes=spikoclust_noisewhiten(spikes,spikeless,'maxnoisetraces',maxnoisetraces,'regularize',regularize);
end

% store unwhitened times and use the unwhitened spikes for spike times

spikes.storewindows=spikes.windows;

% upsample and align, then downsample and whiten!!!

spikes=spikoclust_upsample_align(spikes,'interpolate_fs',interpolate_fs,'align_feature',align_feature);	
[nsamples,ntrials,nchannels]=size(spikes.windows);

spikes.windows=reshape(permute(spikes.windows,[1 3 2]),[],ntrials);
spikes.storewindows=reshape(permute(spikes.storewindows,[1 3 2]),[],ntrials);

if ~gui_clust
	[cluster.windows cluster.times cluster.trials cluster.isi cluster.stats... 
		cluster.outliers cluster.spikedata cluster.model]=...
		spikoclust_autosort(spikes,...
			'fs',FS,'interpolate_fs',interpolate_fs,'proc_fs',sort_fs,...
			'clust_check',clust_check,'pcs',pcs,...
			'workers',spikeworkers,'garbage',garbage,'smem',smem,'modelselection',...
			modelselection,'align_feature',align_feature);
else
	[cluster.windows cluster.times cluster.trials cluster.isi cluster.stats...
		cluster.outliers cluster.spikedata cluster.model]=...
		spikoclust_guisort(spikes,spikeless,cluster.parameters,...
			'fs',FS,'interpolate_fs',interpolate_fs,'proc_fs',sort_fs,...
			'clust_check',clust_check,'pcs',pcs,...
			'workers',spikeworkers,'garbage',garbage,'smem',smem,'modelselection',...
			modelselection,'align_feature',align_feature);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IFR, SMOOTH RATE %%%%%%%%%%%%%%%%%%%

uniq_clusters=1:length(cluster.windows);
nclust=length(uniq_clusters);

if ~isempty(cluster.windows)

	% cycle through each cluster id

	for j=1:nclust

		IFR{j}=zeros(ntrials,samples);

		for k=1:ntrials

			clusterspikes=cluster.times{j}(cluster.trials{j}==k);

			% IFR will use the current sampling rate

			IFR{j}(k,:)=spikoclust_ifr(round(clusterspikes),samples,FS);

		end
	end
end

cluster.IFR=IFR;
