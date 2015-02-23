function [LABELS MODEL CLUSTER_DATA]=spikoclust_autosort(SPIKES,varargin)
%automated spike clustering using a GMM with split-and-merge EM
%

% spikewindows', rows x samples, each row is a windowed spike waveform

if nargin<2
	error('ephysPipeline:suavis:notenoughparams','Need 2 arguments to continue, see documentation');
end

nparams=length(varargin);

CLUSTER_DATA=[];
LABELS=[];
MODEL=[];

clust_check=10;
pcs=4;
workers=1;
garbage=1;
smem=1;
modelselection='icl';

%outlier_cutoff=.05; % posterior probability cutoff for outliers (.6-.8 work well) [0-1, high=more aggresive]

nfeatures=10; % number of features to use, ranked by dimreduction technique

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
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
	end
end

[nsamples,ntrials,nchannels]=size(SPIKES.windows);

[idx CLUSTER_DATA MODEL]=spikoclust_gmmsort(SPIKES.windows,...
	'smem',smem,'garbage',garbage,'clust_check',clust_check,...
	'pcs',pcs,'workers',workers,'modelselection',modelselection,...
	'sigma_fix',sigma_fix);

features=size(CLUSTER_DATA,2); % what's the dimensionality of the data used for sorting?

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


