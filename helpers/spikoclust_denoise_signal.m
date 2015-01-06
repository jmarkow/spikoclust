function [DATA,CAR]=spikoclust_denoise_signal(EPHYS_DATA,CHIN,CHOUT,varargin)
%Denoises signal using specified method
%
%
%
%

if nargin<3 | isempty(CHOUT), CHOUT=CHIN; end

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

%%%

car_exclude=[];
method='none';
car_trim=40;

for i=1:2:nparams
	switch lower(varargin{i})
		case 'car_exclude'
			car_exclude=varargin{i+1};
		case 'method'
			method=varargin{i+1};
		case 'car_trim'
			car_trim=varargin{i+1};
		otherwise

	end
end

CAR=[];

exclude_channels=[];

for i=1:length(car_exclude)
	exclude_channels(i)=find(CHIN==car_exclude(i));
end

car_electrodes=setdiff(1:length(CHIN),exclude_channels); % which electrodes are good for CAR?
ndims_ephys=ndims(EPHYS_DATA);

[samples,ntrials,nchannels]=size(EPHYS_DATA);

% map each channel appropriately

chmap=[];
for i=1:length(CHOUT)
	idx=find(CHOUT(i)==CHIN);
	if ~isempty(idx)
		chmap=[chmap idx];
	end
end

if ndims_ephys==3
	DATA=zeros(samples,ntrials,length(chmap));
elseif ndims_ephys==2
	DATA=zeros(samples,length(chmap));
elseif ndims_ephys==1
	DATA=zeros(samples,1);
else
	error('Data must contain 1-3 dimensions');
end

switch lower(method)

	case 'car'

		% trimmed mean to avoid subtracting in artifacts and spikes

		disp(['Using electrodes ' num2str(CHIN(car_electrodes)) ' for CAR']);
		disp(['Trimmed mean prctile ' num2str(car_trim)]);

		if ndims_ephys==3
			CAR=trimmean(EPHYS_DATA(:,:,car_electrodes),car_trim,'round',3);

			for i=1:length(chmap)
				DATA(:,:,i)=EPHYS_DATA(:,:,chmap(i))-CAR;
			end
		elseif ndims_ephys==2
			
			CAR=trimmean(EPHYS_DATA(:,car_electrodes),car_trim,'round',2);

			for i=1:length(chmap)
				DATA(:,i)=EPHYS_DATA(:,chmap(i))-CAR;
			end

		else
			error('Data must contain at least two dimensions!');
		end


    	otherwise


		if ndims_ephys==3

			for i=1:length(chmap)
				DATA(:,:,i)=EPHYS_DATA(:,:,chmap(i));
			end

		elseif ndims_ephys==2

			for i=1:length(chmap)
				DATA(:,i)=EPHYS_DATA(:,chmap(i));
			end
		else
			DATA=EPHYS_DATA;
		end


end



