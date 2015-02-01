function [WINDOWS TIMES TRIALS ISI STATS OUTLIERS SPIKEDATA MODEL]=spikoclust_autosort(SPIKES,varargin)
%automated spike clustering using a GMM with split-and-merge EM
%

% spikewindows', rows x samples, each row is a windowed spike waveform

if nargin<2
	error('ephysPipeline:suavis:notenoughparams','Need 2 arguments to continue, see documentation');
end

nparams=length(varargin);

SPIKEDATA=[];
CLUSTERPOINTS=[];
LABELS=[];
TRIALS=[];
ISI=[];
WINDOWS=[];
OUTLIERS=[];
MODEL=[];

fs=25e3;
%interpolated_fs=200e3;
interpolate_fs=200e3;
proc_fs=25e3;
maxnoisetraces=1e6;
clust_check=10;
pcs=4;
workers=1;
garbage=1;
smem=1;
modelselection='icl';
align_feature='min';
regularize=.01;
noisewhiten=1;

%outlier_cutoff=.05; % posterior probability cutoff for outliers (.6-.8 work well) [0-1, high=more aggresive]

nfeatures=10; % number of features to use, ranked by dimreduction technique

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'fs'
			fs=varargin{i+1};
		case 'interpolate_fs'
			interpolate_fs=varargin{i+1};
		case 'proc_fs'
			proc_fs=varargin{i+1};
		case 'spikecut'
			spikecut=varargin{i+1};
		case 'merge'
			merge=varargin{i+1};
		case 'maxnoisetraces'
			maxnoisetaces=varargin{i+1};
		case 'ranklimit'
			ranklimit=varargin{i+1};
		case 'clust_check'
			clust_check=varargin{i+1};
		case 'pcs'
			pcs=varargin{i+1};
		case 'garbage'
			garbage=varargin{i+1};
		case 'smem'
			smem=varargin{i+1};
		case 'workers'
			workers=varargin{i+1};
		case 'modelselection'
			modelselection=varargin{i+1};
		case 'align_feature'
			align_feature=varargin{i+1};
		case 'noisewhiten'
			noisewhiten=varargin{i+1};
	end
end

% string the channels together for clustering
% get the covariance matrices for whitening

[nsamples,ntrials,nchannels]=size(SPIKES.windows);

[idx spikedata MODEL]=spikoclust_gmmsort(SPIKES.windows,...
	'proc_fs',proc_fs,'fs',fs,'interpolate_fs',interpolate_fs,...
	'smem',smem,'garbage',garbage,'maxnoisetraces',maxnoisetraces,...
	'clust_check',clust_check,'pcs',pcs,'workers',workers,'modelselection',...
	modelselection);

features=size(spikedata,2); % what's the dimensionality of the data used for sorting?

nclust=size(MODEL.mu,1);
clusters=1:nclust;

% number of spikes per cluster is simply the number of labels

nspikes=[];

for i=1:nclust
	nspikes(i)=sum(idx==clusters(i));
end

[val loc]=sort(nspikes,'descend');

% make the number contiguous and sort by number of spikes, descending

LABELS=idx;
LABELS=zeros(size(idx));

for i=1:length(clusters)
	LABELS(idx==clusters(loc(i)))=i;	
end

MODEL.R(:,1:nclust)=MODEL.R(:,loc);
MODEL.mixing(1:nclust)=MODEL.mixing(loc);
MODEL.sigma=MODEL.sigma(:,:,loc);
MODEL.mu=MODEL.mu(loc,:);

clusters=unique(LABELS(LABELS>0));
OUTLIERS=SPIKES.storewindows(:,LABELS==0);

% now assess the cluster quality ,
% take each cluster and check the FP and FN rate

[WINDOWS TIMES TRIALS SPIKEDATA ISI STATS]=...
	spikoclust_cluster_quality(SPIKES.storewindows,SPIKES.times,spikedata,LABELS,SPIKES.trial,MODEL);
