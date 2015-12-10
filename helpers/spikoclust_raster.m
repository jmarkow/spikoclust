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
object_style='h';
spike_width=.5;

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
        case 'spike_width'
            spike_width=varargin{i+1};
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

TIMES=TIMES(:)';
TRIALS=TRIALS(:)';

% create object for each trial, typically more timebins than trials

uniq_trials=unique(TRIALS);
uniq_times=unique(TIMES);

xpoints=[];
ypoints=[];

if strcmp(lower(object_style),'h')
	n=length(uniq_trials);
else
	n=length(uniq_times);
end

for i=1:n

	if strcmp(lower(object_style),'h')
		curr_times=TIMES(TRIALS==uniq_trials(i));
		curr_trials=uniq_trials(i)*ones(size(curr_times));
	else
		curr_trials=TRIALS(TIMES==uniq_times(i));
		curr_times=uniq_times(i)*ones(size(curr_trials));
	end
	
	trial_xpoints=[curr_times;curr_times;NaN(size(curr_times))];
	trial_ypoints=[curr_trials-spike_height;curr_trials+spike_height;NaN(size(curr_times))];
	trial_xpoints=trial_xpoints(:)';
	trial_ypoints=trial_ypoints(:)';

	xpoints=[xpoints trial_xpoints];
	ypoints=[ypoints trial_ypoints];

end

%xpoints=[TIMES;TIMES;NaN(size(TIMES))];
%ypoints=[TRIALS-spike_height;TRIALS+spike_height;NaN(size(TRIALS))];

spike_width
plot(xpoints,ypoints,'-','color',color,'linewidth',spike_width);
%set(gca,'ydir','rev');
%xlim([0 max_time]);


