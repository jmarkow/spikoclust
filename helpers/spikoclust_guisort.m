function [LABELS MODEL CLUSTER_DATA]=spikoclust_guisort(SPIKES,varargin)
%GUI for spike cluster cutting
%
%

disp('Close main sorting window to use current clusters...');

%TODO: make use of the interface more explicit (label top, indicate that user needs to close when finished)

% spikewindows', rows x samples, each row is a windowed spike waveform

% The functions are all nested inside the main function but not each other
% thus by default all variables declared in the main function ARE GLOBAL

nparams=length(varargin);

% all features excluding IFR and SPIKES.times

features_all={'max','min','ne','^2','neo','wid','pgrad','ngrad','PC1','PC2','PC3','PC4'}; 
features={'PCA','pose','nege','posgrad','neggrad','min','max','width','ISI'}; 

% possible features include, min, max, PCA, width, energy and wavelet coefficients
 
channel_labels=[];
colors={'b','r','g','c','m','y','r','g','b'};
outliercolor='k';

ndims=3;

LABELS=[];
TRIALS=[];
ISI=[];
WINDOWS=[];
OUTLIERS=[];
MODEL=[];
TIMES=[];
STATS=[];

legend_labels={};
CLUSTER_DATA=[];
SPIKEDATA=[];

pcs=4;
workers=1;
garbage=1;
smem=1;
sigma_fix=1e-5;

% the template cutoff could be defined by the 95th prctile of the abs(noise) magnitude

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs!');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'features'
			features=varargin{i+1};
		case 'pcs'
			pcs=varargin{i+1};
		case 'garbage'
			garbage=varargin{i+1};
		case 'smem'
			smem=varargin{i+1};
		case 'sigma_fix'
			sigma_fix=varargin{i+1};
	end
end

spike_data=[];
property_names={};

[nsamples ntrials nchannels]=size(SPIKES.windows);

% cheap to compute standard features


if nchannels==1

	geom_features=spikoclust_shape_features(SPIKES.windows);
	
	if any(strcmp('max',lower(features)))
		spike_data=[spike_data geom_features(:,1)];
		property_names{end+1}='max';
	end

	if any(strcmp('min',lower(features)))
		spike_data=[spike_data geom_features(:,2)];
		property_names{end+1}='min';
	end

	if any(strcmp('pose',lower(features)))
		spike_data=[spike_data geom_features(:,3)];
		property_names{end+1}='pose';
	end

	if any(strcmp('nege',lower(features)))
		spike_data=[spike_data geom_features(:,4)];
		property_names{end+1}='nege';
	end

	if any(strcmp('tote',lower(features)))
		spike_data=[spike_data geom_features(:,5)];
		property_names{end+1}='tote';
	end

	if any(strcmp('neo',lower(features)))
		spike_data=[spike_data geom_features(:,6)];
		property_names{end+1}='neo';
	end

	if any(strcmp('width',lower(features)))
		spike_data=[spike_data geom_features(:,7)];
		property_names{end+1}='width';
	end

	if any(strcmp('posgrad',lower(features)))
		spike_data=[spike_data geom_features(:,8)];
		property_names{end+1}='posgrad';
	end

	if any(strcmp('neggrad',lower(features)))
		spike_data=[spike_data geom_features(:,9)];
		property_names{end+1}='neggrad';
	end
end

outlierpoints=[];
if any(strcmp('pca',lower(features)))
	newmodel=spikoclust_gmem(SPIKES.windows',[],1,'garbage',1,'merge',0,'debug',0,'sigma_fix',sigma_fix);
	[v,d]=eigs(newmodel.sigma(:,:,1));
	newscore=-SPIKES.windows'*v;
	spike_data=[spike_data newscore(:,1:pcs)];

	% these comprise the outliers before the projection... set to >1 to include all (default for now)

	outlierpoints=newmodel.R(:,2)>=2;

	for i=1:pcs
		property_names{end+1}=['PC ' num2str(i)];
	end
end

spike_data(isnan(spike_data))=0;
nfeatures=size(spike_data,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GUI setup

stats_window_a=figure('Visible','on','Position',[450,500,300,300],'Name','Stats Window A','NumberTitle','off');
stats_window_b=figure('Visible','on','Position',[450,500,300,300],'Name','Stats Window B','NumberTitle','off');

main_window=figure('Visible','on','Position',[360,500,700,600],'Name','Data Plotter','NumberTitle','off');

plot_axis=axes('Units','pixels','Position',[50,50,425,425]);

pop_up_x= uicontrol('Style','popupmenu',...
	'String',property_names,...
	'Position',[400,90,75,25],'call',@change_plot,...
	'Value',min(1,nfeatures));
pop_up_x_text= uicontrol('Style','text',...
	'String','X',...
	'Position',[405,130,50,45]);

pop_up_y= uicontrol('Style','popupmenu',...
	'String',property_names,...
	'Position',[490,90,75,25],'call',@change_plot,...
	'Value',min(2,nfeatures));
pop_up_y_text= uicontrol('Style','text',...
	'String','Y',...
	'Position',[495,130,50,45]);

pop_up_z= uicontrol('Style','popupmenu',...
	'String',property_names,...
	'Position',[580,90,75,25],'call',@change_plot,...
	'Value',min(3,nfeatures));
pop_up_z_text= uicontrol('Style','text',...
	'String','Z',...
	'Position',[585,130,50,45]);

pop_up_clusters= uicontrol('Style','popupmenu',...
	'String',{'1','2','3','4','5','6','7','8','9'},...
	'Position',[470,210,75,25]);
pop_up_clusters_text= uicontrol('Style','text',...
	'String','Number of Clusters',...
	'Position',[495,250,100,45]);

push_replot_save= uicontrol('Style','pushbutton',...
	'String','Show cluster stats',...
	'Position',[520,40,125,35],'call',@show_stats,'FontSize',11);
push_recluster= uicontrol('Style','pushbutton',...
	'String','Recluster',...
	'Position',[410,40,100,35],'value',0,...
	'Call',@change_cluster,'FontSize',11);

rows=ceil(length(property_names)/5);

i=1;
while i<=length(property_names)
	row=ceil(i/7);
	column=mod(i,7);
	if column==0, column=7; end
	cluster_data_check{i}=uicontrol('Style','checkbox',...
		'String',property_names{i},...
		'Value',i==1,'Position',[5+column*60,550-row*35,70,25]);
	set(cluster_data_check{i},'Units','Normalized')
	i=i+1;
end

cluster_data_text=uicontrol('Style','text',...
	'String','Cluster features',...
	'Position',[250 545 150 25]);
set(cluster_data_text,'backgroundcolor',get(main_window,'color'));
set(cluster_data_text,'units','normalized');

% now align everything and send the main_window handle to the output
% so we can use the gui with uiwait (requires the handle as a return value)

align([pop_up_clusters,pop_up_clusters_text],'Center','None');
align([pop_up_x,pop_up_x_text],'Center','None');
align([pop_up_y,pop_up_y_text],'Center','None');
align([pop_up_z,pop_up_z_text],'Center','None');

change_cluster();
change_plot();

% run change_plot, which updates the plot according to the defaults

set([main_window,plot_axis,pop_up_x,pop_up_x_text,pop_up_y,pop_up_y_text,pop_up_z,...
	pop_up_z_text,pop_up_clusters,pop_up_clusters_text,...
	push_replot_save,push_recluster],'Units','Normalized');
movegui(main_window,'center')

set(main_window,'Visible','On','DeleteFcn',@plot_close);
uiwait(main_window);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Callbacks

% this callback changes the plot and returns the sum of the distances
% from the centroid for each point in a cluster

% change the plot if we change any of our dimensions, DO NOT RECLUSTER!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot change

function change_plot(varargin)

% get the number of dimensions for the plot

% clear the axes and draw the points

cla;

viewdim(1)=get(pop_up_x,'value');
viewdim(2)=get(pop_up_y,'value');
viewdim(3)=get(pop_up_z,'value');

view_data=spike_data(:,viewdim);
clusters=unique(LABELS(LABELS>0));

if ndims==2
	for i=1:length(clusters)
		points=find(LABELS==clusters(i));
		h(:,i)=plot(view_data(points,1),view_data(points,2),...
			'o','markerfacecolor',colors{i},'markeredgecolor','none');hold on
	end

	points=find(LABELS==0);
	if ~isempty(points)
		h(:,length(clusters)+1)=plot(view_data(points,1),view_data(points,2),...
			'o','markerfacecolor',outliercolor,'markeredgecolor','none');hold on
	end
else
	for i=1:length(clusters)
		points=find(LABELS==clusters(i));
		h(:,i)=plot3(view_data(points,1),view_data(points,2),view_data(points,3),...
			'o','markerfacecolor',colors{i},'markeredgecolor','none');hold on

	end

	points=find(LABELS==0);
	if ~isempty(points)
		h(:,length(clusters)+1)=plot3(view_data(points,1),view_data(points,2),view_data(points,3),...
			'o','markerfacecolor',outliercolor,'markeredgecolor','none');hold on
	end


end

grid on
view(ndims)

xlabel(property_names{viewdim(1)});
ylabel(property_names{viewdim(2)});
zlabel(property_names{viewdim(3)});

L=legend(h,legend_labels,'Location','NorthEastOutside');legend boxoff
set(L,'FontSize',20,'FontName','Helvetica')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Recluster

function change_cluster(varargin)

% label everything
% use the dimensions ticked in the top box for clustering

dim=[];
for i=1:length(cluster_data_check)

	value=get(cluster_data_check{i},'Value');

	if value
		dim=[dim i];
	end

end

clusterchoices=get(pop_up_clusters,'string');
clusterselection=get(pop_up_clusters,'value');
clusterchoice=clusterchoices{clusterselection};

% perform the kmeans analysis and return the labels, centroid coordinates,
% sum of all points in each cluster from their respective centroid and
% the distance of all points from all centroids

% start with one cluster, go up to 10 and check the within distance for all clusters

CLUSTER_DATA=spike_data(:,dim);
[datapoints,features]=size(spike_data);

options=statset('Display','off');

clustnum=2:9;
if datapoints<=features
	disp('Too few spikes to fit');
	return;
end

% gaussian mixture seems to work better than fcm


startmu=[];
startcov=[];
mixing=[];

nclust=str2num(clusterchoice);

startobj=struct('mu',startmu,'sigma',startcov,'mixing',mixing);
idx=kmeans(CLUSTER_DATA,nclust,'replicates',5);

%% set up initial model

mu=[];
for j=1:nclust
	startmu(j,:)=mean(CLUSTER_DATA(idx==j,:))';
	startcov(:,:,j)=diag(var(CLUSTER_DATA));
end

startobj.mu=startmu;
startobj.sigma=startcov;

for j=1:nclust
	startobj.mixing(j)=sum(idx==j)/length(idx);
end

clustermodel=spikoclust_gmem(CLUSTER_DATA,startobj,nclust,...
		'garbage',garbage,'merge',smem,'debug',0,'sigma_fix',sigma_fix);

MODEL=clustermodel;
MODEL.pcs=v;
idx=[];

for i=1:size(clustermodel.R,1)
	posteriors=clustermodel.R;
	[~,idx(i)]=max(posteriors(i,:));
end

if garbage
	garbageidx=find(clustermodel.garbage);
	idx(idx==garbageidx)=NaN;
end

LABELS=zeros(size(CLUSTER_DATA,1),1);

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


OUTLIERS=SPIKES.storewindows(:,isnan(idx));

clusters=unique(idx(idx>0)); % how many clusters?
nclust=length(clusters);

% number of spikes per cluster is simply the number of labels

%%%%% rearrange code start

nspikes=[];

for i=1:nclust
	nspikes(i)=sum(idx==clusters(i));
end

[val loc]=sort(nspikes,'descend');

% make the number contiguous and sort by number of spikes, descending

LABELS=zeros(size(idx));

for i=1:nclust
	LABELS(idx==clusters(loc(i)))=i;	
end


clustermodel.R(:,1:nclust)=clustermodel.R(:,loc);
clustermodel.mixing(1:nclust)=clustermodel.mixing(loc);
clustermodel.sigma=clustermodel.sigma(:,:,loc);
clustermodel.mu=clustermodel.mu(loc,:);

%%%% rearrange code end (comment out to not sort clusters by nspikes)

% return labels, and windows and ISI sorted by cluster IDX

% clear the plot axis

% plot in either 2 or 3 dims
% turns out plot is MUCH faster than scatter, changed accordingly...

legend_labels={};
for i=1:nclust
	legend_labels{i}=['Cluster ' num2str(i)];
end

if garbage & any(isnan(idx))
	legend_labels{end+1}='Outliers';
end

% compute any other stats we want, ISI, etc...
[WINDOWS TIMES TRIALS SPIKEDATA ISI STATS]=...
	spikoclust_cluster_quality(SPIKES.storewindows,SPIKES.times,CLUSTER_DATA,LABELS,SPIKES.trial,MODEL);
change_plot();

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Show a window with spike stats

function show_stats(varargin)

% get the labels from the main_window

cluster=[];
cluster=struct('windows',{WINDOWS},'times',{TIMES},'trials',{TRIALS},'spikedata',{SPIKEDATA},'stats',{STATS},...
	'isi',{ISI},'model',MODEL);
cluster.parameters.interpolate_fs=SPIKES.fs;
cluster.parameters.fs=SPIKES.original_fs;

nclust=length(cluster.windows);

set(0,'CurrentFigure',stats_window_a);

pos=get(stats_window_a,'Position');
set(stats_window_a,'Position',[pos(1) pos(2) 250+200*nclust 250+200*nclust]);
spikoclust_visual_clustering(cluster,'fig_num',stats_window_a);

pos=get(stats_window_b,'Position');
set(0,'CurrentFigure',stats_window_b);
set(stats_window_b,'Position',[pos(1) pos(2) 250*nclust 600]);

spikoclust_visual_waveforms(cluster,'fig_num',stats_window_b);

end

function plot_close(varargin)

if ishandle(stats_window_a) && strcmp(get(stats_window_a,'type'),'figure')
	close(stats_window_a);
end

if ishandle(stats_window_b) && strcmp(get(stats_window_b,'type'),'figure')
	close(stats_window_b);
end

end

end
