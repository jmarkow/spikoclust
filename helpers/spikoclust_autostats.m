function CLUSTER=spikoclust_autostats(CLUSTER,SPIKELESS,varargin)
%script for visualizing the results of clustering
%

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

savemode=0;
savedir=pwd;
savefilename_stats='stats';
ntrials=[];
channelboundary=[];
spikecut=1;
snr_cutoff=6.6;
lratio_cutoff=1;
isod_cutoff=1;
isi_cutoff=.01;
channels=[];
clust_plot=1;
dim_labels=[];

for i=1:2:nparams
	switch lower(varargin{i})
		case 'savedir'
			savedir=varargin{i+1};
		case 'savefilename_stats'
			savefilename_stats=varargin{i+1};
		case 'ntrials'
			ntrials=varargin{i+1};
		case 'channelboundary'
			channelboundary=varargin{i+1};
		case 'spikecut'
			spikecut=varargin{i+1};
		case 'snr_cutoff'
			snr_cutoff=varargin{i+1};
		case 'lratio_cutoff'
			lratio_cutoff=varargin{i+1};
		case 'isi_cutoff'
			isi_cutoff=varargin{i+1};
		case 'isod_cutoff'
			isod_cutoff=varargin{i+1};
		case 'channels'
			channels=varargin{i+1};
		case 'savemode'
			savemode=varargin{i+1};
		case 'clust_plot'
			clust_plot=varargin{i+1};
		case 'dim_labels'
			dim_labels=varargin{i+1};
	end
end

disp(['L-ratio cutoff:  ' num2str(lratio_cutoff)]);
disp(['IsoD cutoff:  ' num2str(isod_cutoff)]);
disp(['ISI cutoff:  ' num2str(isi_cutoff)]);
disp(['SNR cutoff:  ' num2str(snr_cutoff)]);
nplots=length(CLUSTER.windows);

if isempty(ntrials)

	trialmax=-inf;

	for i=1:nplots
		tmpmax=max(CLUSTER.trials{i});
		if tmpmax>trialmax
			trialmax=tmpmax;
		end
	end

	ntrials=trialmax;
end

% estimate SNR from 6*std of spikeless trace

noise_p2p=std(SPIKELESS{1});

noise_p2p(isnan(noise_p2p))=[];
mean_noise_p2p=mean(noise_p2p);

for j=1:nplots
	mean_wave=mean(CLUSTER.windows{j},2);
	max_wave=max(mean_wave);
	min_wave=min(mean_wave);
	CLUSTER.stats.snr(j)=[ abs(max_wave-min_wave)/mean_noise_p2p ];
end

if savemode
	stats_fig=figure('Visible','off');
else
	stats_fig=figure('Visible','on');
end

note=[];

% TODO: modify note to show stats for all clusters!

for j=1:nplots
	note{j}=['L-ratio ' sprintf('%.2f',CLUSTER.stats.lratio(j)) ...
		' IsoD ' sprintf('%.2f',CLUSTER.stats.isod(j)) ];
end

if ~isempty(CLUSTER.parameters.tetrode_channels)
	channelboundary=round(cumsum(ones(1,length(CLUSTER.parameters.tetrode_channels)).*...
		sum(CLUSTER.parameters.spike_window*CLUSTER.parameters.interpolate_fs)));
end

stats_fig=spikoclust_visual_waveforms(CLUSTER,...
	'noise_p2p',mean_noise_p2p,'fig_num',stats_fig,'note',note,...
	'channelboundary',channelboundary);

set(stats_fig,'Position',[0 0 250*nplots 600]);
set(stats_fig,'PaperPositionMode','auto');

% label candidate units if they meet our criteria
% isi intervals < absolute refractory period

isi_violations=sum((CLUSTER.isi{j}./CLUSTER.parameters.fs)<.001);
isi_violations=isi_violations/length(CLUSTER.isi{j});

% are there enough spikes?

if savemode
	markolab_multi_fig_save(stats_fig,savedir,...
		[ savefilename_stats 'spikestats' ],'eps,png','res',100,'renderer','painters');
	close([stats_fig]);
end

for j=1:nplots
	if length(CLUSTER.times{j})>=spikecut*ntrials & savemode

		if CLUSTER.stats.snr(j)>=snr_cutoff && ...
				(isnan(CLUSTER.stats.lratio(j)) || ...
				CLUSTER.stats.lratio(j)<=lratio_cutoff) &&  ...
				isi_violations<isi_cutoff 
			if isnan(CLUSTER.stats.isod(j)) || CLUSTER.stats.isod(j)>=isod_cutoff
				fid=fopen(fullfile(savedir,['candidate_unit_ch' num2str(channels) ... 
					'_cl' num2str(j)]),'w');
				fclose(fid);
			end
		end
	end
end

if nplots<=10

	if savemode
		stats_fig=figure('Visible','off','renderer','painters');
	else
		stats_fig=figure('Visible','on','renderer','painters');
	end

	stats_fig=spikoclust_visual_clustering(CLUSTER,'fig_num',stats_fig);

	set(stats_fig,'Position',[0 0 250+200*nplots 250+200*nplots]);
	set(stats_fig,'PaperPositionMode','auto');

	if savemode
		markolab_multi_fig_save(stats_fig,savedir,...
			[ savefilename_stats 'cluststats' ],'eps,png','res',100,'renderer','painters');
		close([stats_fig]);
	end


	if size(CLUSTER.model.mu,2)>1 & clust_plot

		if savemode
			stats_fig=figure('Visible','off');
		else
			stats_fig=figure('Visible','on');
		end


		stats_fig=spikoclust_gaussvis(CLUSTER.model,CLUSTER.spikedata,'fig_num',stats_fig,'dim_labels',dim_labels);
		
		if savemode
			markolab_multi_fig_save(stats_fig,savedir,...
				[ savefilename_stats 'clustplot' ],'eps,png','res',100,'renderer','painters');
			close([stats_fig]);
		end
	end

end

