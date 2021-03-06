function SPIKES=spikoclust_noisewhiten(SPIKES,NOISEDATA,varargin)
%
%
%
%
%

nparams=length(varargin);

maxnoisetraces=1e6;

if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'maxnoisetraces'
			maxnoisetaces=varargin{i+1};
	end
end

if ~iscell(NOISEDATA)
	tmp{1}=NOISEDATA;
	clear NOISEDATA;
	NOISEDATA=tmp;
end

[nsamples,ntrials,nchannels]=size(SPIKES(1).windows);
ntraces=length(SPIKES);

for i=1:length(NOISEDATA)

	len=length(NOISEDATA{i});
	noisetrials=floor(length(NOISEDATA{i})/nsamples);

	if noisetrials>maxnoisetraces
		noisetrials=maxnoisetraces;
	end

	disp(['Noise trials ' num2str(noisetrials)]);

	residual=len-(noisetrials*nsamples);

	NOISEDATA{i}=NOISEDATA{i}(1:len-residual);
	noisematrix=reshape(NOISEDATA{i},nsamples,[])';
	noisecov=cov(noisematrix);
    	noisecov=spikoclust_gmem_covcheck(noisecov);

	r=chol(noisecov);
	inv_r=inv(r);

	SPIKES.windows(:,:,i)=[SPIKES.windows(:,:,i)'*inv_r]';
end
