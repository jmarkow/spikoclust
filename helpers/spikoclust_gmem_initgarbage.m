function [P,mixing]=initgarbage(DATA,NCLUST,MODEL)

mixing=MODEL.mixing;
[datapoints,D]=size(DATA);

% set the uniform density

datarange=range(DATA);

% attempt to compute the volume of the space

p=prod(datarange);

if p==inf

	% issue warning here

	warning('gmem:volumetoolarge',...
		'Cannot compute volume, setting to : %e',p);
	p=1e30;
end

P=(1/p).*ones(datapoints,1);

% set the initial mixing proportion

mixing(end+1)=1/(NCLUST);

% renormalize the mixing proportions

mixing=mixing./sum(mixing);

end
