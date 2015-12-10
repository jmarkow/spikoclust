function [LABELS SPIKE_DATA MODEL]=spikoclust_gmmsort(SPIKE_DATA,varargin)
% calculate merge and accept thresholds based on estimate of the noise
% match the noise sampling rate to the spike sampling rate and

nparams=length(varargin);

maxnoisetraces=1e6; % maximum number of noise traces to use for Cholesky decomposition
clust_check=1:6; % number of clusters to start with
clustreplicates=1; % replicates for clustering procedure
garbage=1; % use uniform garbage collection density
workers=1; % number of workers when deployed
modelselection='icl'; % icl, bic, aic
smem=1; % smem 1 uses split and merge, 0 is standard em, 2 for free split and merge
outlierpoints=[];

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'clust_check'
			clust_check=varargin{i+1};
		case 'clustreplicates'
			clustreplicates=varargin{i+1};
		case 'pcs'
			pcs=varargin{i+1};
		case 'garbage'
			garbage=varargin{i+1};
		case 'workers'
			workers=varargin{i+1};
		case 'modelselection'
			modelselection=varargin{i+1};
		case 'smem'
			smem=varargin{i+1};
		case 'outlierpoints'
			outlierpoints=varargin{i+1};
	end
end

%%% TODO: use new fields

[nsamples,ntrials]=size(SPIKE_DATA);

if isempty(outlierpoints)
	outlierpoints=false(nsamples,1);
end

% do we need to resample the noise data to match the sorting fs?

if isdeployed
	matlabpool('open',workers);
end

% these comprise the outliers before the projection... set to >1 to include all (default for now)

disp(['Robust PCA outliers:  ' num2str(sum(outlierpoints))]);

% only cluster the non-outliers

SPIKE_DATA=SPIKE_DATA(~outlierpoints,:);

% might consider removing outliers before we project to the new space (Sahani,1999)

% commented out, but need to consider later...
% outlierpoints=find(newmodel.R(:,2)>.98);

% outliers are any points with garbage components R>.9

oldstate1=warning('off','stats:gmdistribution:FailedToConverge');
oldstate2=warning('off','stats:kmeans:FailedToConvergeRep');
oldstate3=warning('off','stats:kmeans:FailedToConverge');

% for now free_gmem is deprecated, forcing standard gmem until it is fixed

if 1<0
	disp('Free GMEM...');
	tmpclustobj={};
	startmu=[];
	startcov=[];
	mixing=[];

	startobj=struct('mu',startmu,'sigma',startcov,'mixing',mixing);

	loglikelihood=zeros(1,clustreplicates);
	idx=kmeans(SPIKE_DATA,clust_check(1),'replicates',5);

	%% set up initial model

	mu=[];
	for i=1:clust_check(1)
		startmu(i,:)=mean(SPIKE_DATA(idx==i,:))';
		startcov(:,:,i)=diag(var(SPIKE_DATA));
	end

	startobj.mu=startmu;
	startobj.sigma=startcov;

	for i=1:clust_check(1)
		startobj.mixing(i)=sum(idx==i)/length(idx);
	end

	for i=1:clustreplicates
		tmpclustobj{i}=spikoclust_free_gmem(SPIKE_DATA,startobj,clust_check(1),...
			'garbage',garbage,'merge',smem,'debug',0);
		loglikelihood(i)=tmpclustobj{i}.likelihood;
	end

	[~,loc]=max(loglikelihood);
	clustermodel=tmpclustobj{loc(1)};
else
	disp('Standard GMEM...');
	for i=1:1:length(clust_check)

		tmpclustobj={};
		startmu=[];
		startcov=[];
		mixing=[];

		startobj=struct('mu',startmu,'sigma',startcov,'mixing',mixing);

		loglikelihood=zeros(1,clustreplicates);
		idx=kmeans(SPIKE_DATA,clust_check(i),'replicates',5);

		%% set up initial model

		mu=[];
		for j=1:clust_check(i)
			selection=find(idx==j);

			% bugfix, if only one sample then mean collapses to single
			% sample (fixed 8/8/2014)

			if length(selection)==1
				startmu(j,:)=SPIKE_DATA(selection,:);
			else
				startmu(j,:)=mean(SPIKE_DATA(selection,:));
			end

			%startmu(j,:)=mean(SPIKE_DATA(idx==j,1:rankcut))';
			startcov(:,:,j)=diag(var(SPIKE_DATA));
		end

		startobj.mu=startmu;
		startobj.sigma=startcov;

		for j=1:clust_check(i)
			startobj.mixing(j)=sum(idx==j)/length(idx);
		end

		for j=1:clustreplicates
			tmpclustobj{j}=spikoclust_gmem(SPIKE_DATA,startobj,clust_check(i),...
				'garbage',garbage,'merge',smem,'debug',0);
			loglikelihood(j)=tmpclustobj{j}.likelihood;
		end

		% only keep the clustobj with the best likelihood

		[~,loc]=max(loglikelihood);
		clustobj{i}=tmpclustobj{loc(1)};
		BIC(i)=clustobj{i}.BIC;
		MML(i)=clustobj{i}.MML;
		ICL(i)=clustobj{i}.ICL;

	end

	if isdeployed
		matlabpool('close');
	end

	warning(oldstate1);
	warning(oldstate2);
	warning(oldstate3);

	switch lower(modelselection(1))
		case 'b'
			[~,loc]=min(BIC);
		case 'm'
			[~,loc]=max(MML);
		case 'i'
			[~,loc]=min(ICL);
		otherwise
	end

	clustermodel=clustobj{loc(1)};

end

MODEL=clustermodel;

% get the labels from the responsibilities (probability of each cluster given each datapoint)

idx=[];
for i=1:size(clustermodel.R,1)
	posteriors=clustermodel.R;
	[~,idx(i)]=max(posteriors(i,:));
end

if garbage
	garbageidx=find(clustermodel.garbage);
	idx(idx==garbageidx)=NaN;
end

LABELS=zeros(ntrials,1);

% what did we label through clustering

LABELS(~outlierpoints)=idx;

% pre-pca outliers

LABELS(outlierpoints)=NaN;

grps=unique(LABELS(LABELS>0));
nclust=length(grps);

% ensure the labeling is contiguous

idx=LABELS;
for i=1:nclust
	idx(LABELS==grps(i))=i;
end

LABELS=idx;
clear idx;

disp(['N Clust:  ' num2str(length(unique(LABELS(LABELS>0))))]);
