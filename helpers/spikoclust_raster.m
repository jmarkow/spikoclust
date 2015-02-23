function spikoclust_raster(TIMES,TRIALS,varargin)
%spikoclust_raster generates spike rasters given a vector of spiketimes and trial IDs
%
%	spikoclust_raster(TIMES,TRIALS,FS,varargin)
%
%	TIMES
%	vector spike times
%
%	TRIALS
%	vector of spike trial
%
%	FS
%	sampling rate
%

nparams=length(varargin);
if mod(nparams,2)>0
	error('ephysPipeline:argChk','Parameters must be specified as parameter/value pairs!');
end

%%%
spike_height=.5;
max_time=[];
fs=[];
color='k';

for i=1:2:nparams
	switch lower(varargin{i})
		case 'spike_height'
			spike_height=varargin{i+1};	
		case 'max_time'
			max_time=varargin{i+1};
		case 'fs'
			fs=varargin{i+1};
		case 'color'
			color=varargin{i+1};
		otherwise

	end
end

% grab the spikes from a particular cluster and triasl

if iscell(TIMES)
	new_times=[];
	TRIALS=[];

	for i=1:length(TIMES)
		new_times=[new_times TIMES{i}];
		TRIALS=[TRIALS ones(1,length(TIMES{i}))*i];
	end

	TIMES=new_times;

end


if ~isempty(fs)
	TIMES=TIMES/fs;
end

trials=min(TRIALS):max(TRIALS);

maxt=-inf;
for i=trials
	
	% grab the spike times for each trial

	newspikes{i}=TIMES(TRIALS==i);

	maxspike=max(newspikes{i});
	if maxspike>maxt
		maxt=maxspike;
	end
end

newtrialvec=[];
newspikevec=[];

% generate spike matrices

for i=1:length(newspikes)
	newspikevec=[newspikevec [ newspikes{i};newspikes{i}] ];
	newtrialvec=[newtrialvec [ ones(1,length(newspikes{i})).*i+spike_height; ones(1,length(newspikes{i})).*i-spike_height] ];
end

if isempty(max_time)
	max_time=maxt;
end

plot(newspikevec,newtrialvec,'-','color',color);
set(gca,'ydir','rev');
xlim([0 max_time]);


