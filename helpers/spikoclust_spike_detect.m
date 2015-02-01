function [SPIKES]=ephys_spike_detect(DATA,THRESH,varargin)
%ephys_spike_detect.m performs spike detection on a vector with a pre-determined
%threshold
%
% [SPIKES]=spike_detect(DATA,fs,traces,THRESH)
%
%
% DATA
% sample x trace matrix of voltage recrordings (threshold crossings detected on the first column, the rest are slaved)
%
% fs
% sampling rate of the recording
%
% THRESH
% threshold for detecting spikes
%
% the following can be specified as a parameter/value pair
%
% censor
% length between spikes in S, i.e. the censor period (default .001)
%
% window
% two element vector specifying the distance before and after the spike to store
% (default [.001 .001])
%
% method
% string specifying whether to use (p)ositive threshold crossings, (n), or (b) 
%
% visualize
% generate a figure to visualize spike detection results
%
%
% align_feature
% ('min' for absolute minimum,
% 'max' for absolute maximum and 'com' for center of mass about the minimum)
%
% 

SPIKES=[];

if nargin<1
	error('Need the input data to continue!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input argument collection

censor=.75e-3; % minimum time between spikes, i.e. censor period
	       % per Hill, Mehta and Kleinfeld (2011), 750 microseconds
window=[.0004 .0004]; % how large of a window to grab, seconds before and after spike
method='b'; % how to grab spikes, [p]os, [n]eg, or [b]oth (abs value)
visualize='y';
fs=25e3;
jitter=4; % how much jitter do we allow before tossing out a spike (in samples of original fs)?
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
		case 'window'
			window=varargin{i+1};
		case 'method'
			method=varargin{i+1};
		case 'visualize'
			visualize=varargin{i+1};
		case 'fs'
			fs=varargin{i+1};
		case 'jitter'
			jitter=varargin{i+1};	
	end
end

[nsamples,ntrials,nchannels]=size(DATA);
if nsamples==1
	warning('Either data is in wrong format or only contains 1 sample!');
end

% specify the spike window in terms of samples

% the frame to grab (in the original sampling rate) around the threshold crossing
% collect jitter samples around the frame for upsample and realign_featurement

SPIKES.frame=window;
SPIKES.jitter=jitter;

frame=round(window*fs);
frame=frame+jitter;
frame_length=length([-frame(1):frame(2)]);
timepoints=-frame(1):frame(2);

% collect the pos-going and neg-going spikes

SPIKES.times=zeros(1,maxspikes);
SPIKES.windows=zeros(frame_length,maxspikes,nchannels);
SPIKES.storewindows=[];
SPIKES.storetimes=[];
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
		if dtime<censor*fs
			spike_times(counter)=[];
		else
			counter=counter+1;
		end
	end

	for j=1:length(spike_times)

		if spike_times(j)-frame(1)>0 && spike_times(j)+frame(2)<length(DATA(:,1))

			% find the absolute minimum (or max) and use as the spike peak for align_featurement

			tmp_time=spike_times(j);
			tmp_window=DATA(tmp_time-frame(1):tmp_time+frame(2),i,:);
			SPIKES.times(spike_counter)=tmp_time;
			SPIKES.windows(:,spike_counter,:)=tmp_window;
			SPIKES.threshold(spike_counter)=THRESH(i);
			SPIKES.trial(spike_counter)=i;
			spike_counter=spike_counter+1;

		end
	end

end

% how much of array is left unused...

SPIKES.times(spike_counter:maxspikes)=[];
SPIKES.windows(:,spike_counter:maxspikes,:)=[];
SPIKES.trial(spike_counter:maxspikes)=[];
SPIKES.threshold(spike_counter:maxspikes)=[];
SPIKES.fs=fs;
SPIKES.censor=censor;
SPIKES.window_time=timepoints./fs;

% make sure we haven't made any alignments that violate the censor period

for i=1:ntrials
	counter=2;
	while counter<=length(SPIKES.times)
		dtime=SPIKES.times(counter)-SPIKES.times(counter-1);
		if dtime<censor*fs && SPIKES.trial(counter)==SPIKES.trial(counter-1)
			SPIKES.times(counter)=[];
			SPIKES.windows(:,counter,:)=[];
			SPIKES.trial(counter)=[];
			SPIKES.threshold(counter)=[];
		else
			counter=counter+1;
		end
	end
end


% visualize the voltage trace, threshold(s) and spikes

if lower(visualize(1))=='y'

	nsamples=length(DATA);
	figure();
	plot([1:nsamples]./fs,DATA(:,1),'b');hold on
	ylabel({'Voltage (in V)';['Threshold (in V):  ' num2str(THRESH)]},'FontSize',13,'FontName','Helvetica');
	xlabel('T (in s)','FontSize',13,'FontName','Helvetica');
	plot([1:nsamples]./fs,ones(nsamples,1).*THRESH,'r');
	plot([1:nsamples]./fs,ones(nsamples,1).*-THRESH,'r');

	plot(SPIKES.times/fs,DATA(SPIKES.times,1),'b*','markersize',10);
	set(gca,'FontSize',11,'FontName','Helvetica')
	box off
	axis tight;

end

