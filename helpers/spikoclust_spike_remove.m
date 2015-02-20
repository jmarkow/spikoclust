function SPIKELESS=spikoclust_spike_remove(DATA,SPIKES)
%Removes spikes from data
%
%
%
%


nspikes=length(SPIKES.times);
edges=round(SPIKES.frame*SPIKES.fs);
coords=zeros(nspikes,3);

leftedge=1;
lasttrial=1;
to_del=[];

for i=1:nspikes

	% left edge of next spike

	rightedge=SPIKES.times(i)-edges(1);

	if lasttrial==SPIKES.trial(i)
		coords(i,:)=[ leftedge rightedge SPIKES.trial(i) ];
	else
		to_del=[to_del i];
	end

	leftedge=SPIKES.times(i)+edges(2);
	lasttrial=SPIKES.trial(i);

end

coords(to_del,:)=[];

skip=find(coords(:,1)<1);
coords(skip,:)=[];

ncoords=size(coords,1);
nsamples=sum(diff(coords(:,1:2),[],2))+ncoords;

SPIKELESS=zeros(nsamples,1);
counter=1;

for i=1:ncoords
	nsamples=(coords(i,2)-coords(i,1));

	if nsamples<0
		continue;
	end

	SPIKELESS(counter:counter+nsamples)=DATA(coords(i,1):coords(i,2),coords(i,3));
	counter=counter+nsamples+1;
end

