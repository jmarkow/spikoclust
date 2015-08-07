function [SPIKES]=spikoclust_spike_detect_mu(DATA,THRESH,FS,varargin)
%spikoclust_spike_detect.m performs spike detection on a vector with a pre-determined
%threshold
%
%	[SPIKES]=spikoclust_spike_detect(DATA,THRESH,FS)
%
%
%	DATA
%	sample x trial matrix of voltage recrordings (threshold crossings detected on the first column, the rest are slaved)
%
%	THRESH
%	threshold for detecting spikes (can pass a vector where each element corresponds to threshold for trial)
%

SPIKES=[];

if nargin<1
	error('Need the input data to continue!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input argument collection

censor=.75e-3; % minimum time between spikes, i.e. censor period
	       % per Hill, Mehta and Kleinfeld (2011), 750 microseconds
method='b'; % how to grab spikes, [p]os, [n]eg, or [b]oth (abs value)
visualize='y';
maxspikes=1e5;

%%%%%%%%%%%%%%%%%%%%%%%%%%

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'censor'
			censor=varargin{i+1};
		case 'method'
			method=varargin{i+1};
		case 'visualize'
			visualize=varargin{i+1};
	end
end

[nsamples,ntrials,nchannels]=size(DATA);
if nsamples==1
	warning('Either data is in wrong format or only contains 1 sample!');
end

% specify the spike window in terms of samples

% collect the pos-going and neg-going spikes

SPIKES.times=zeros(1,maxspikes);
SPIKES.trial=zeros(1,maxspikes);
SPIKES.threshold=zeros(1,maxspikes);
spike_counter=1;

for i=1:ntrials

	if lower(method(1))=='b' || lower(method(1))=='a'	
		spike_times=find(abs(DATA(:,i,1))>THRESH(i));
	elseif lower(method(1))=='p'
		spike_times=find(DATA(:,i,1)>THRESH(i));
	else
		spike_times=find(DATA(:,i,1)<-THRESH(i));
	end

	nspikes=length(spike_times);

	% censor period

	counter=2;
	while counter<=length(spike_times)
		dtime=spike_times(counter)-spike_times(counter-1);
		if dtime<censor*FS
			spike_times(counter)=[];
		else
			counter=counter+1;
		end
	end

	nspikes=length(spike_times);

	if nspikes==0
		continue;
	end

	idxs=spike_counter:spike_counter+nspikes-1;


	SPIKES.times(idxs)=spike_times;
	SPIKES.trial(idxs)=ones(size(spike_times))*i;
	SPIKES.threshold(idxs)=ones(size(spike_times))*THRESH(i);

	spike_counter=idxs(end)+1;

end

% how much of array is left unused...

SPIKES.times(spike_counter:maxspikes)=[];
SPIKES.trial(spike_counter:maxspikes)=[];
SPIKES.threshold(spike_counter:maxspikes)=[];

SPIKES.fs=FS;
SPIKES.censor=censor;

% make sure we haven't made any alignments that violate the censor period
% visualize the voltage trace, threshold(s) and spikes

if lower(visualize(1))=='y' & ntrials==1

	nsamples=length(DATA);
	figure();
	plot([1:nsamples]./FS,DATA(:,1),'b');hold on
	ylabel({'Voltage (in V)';['Threshold (in V):  ' num2str(THRESH)]},'FontSize',13,'FontName','Helvetica');
	xlabel('T (in s)','FontSize',13,'FontName','Helvetica');
	plot([1:nsamples]./FS,ones(nsamples,1).*THRESH,'r');
	plot([1:nsamples]./FS,ones(nsamples,1).*-THRESH,'r');

	plot(SPIKES.times/FS,DATA(SPIKES.times,1),'b*','markersize',10);
	set(gca,'FontSize',11,'FontName','Helvetica')
	box off
	axis tight;

end


