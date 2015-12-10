function NEWMODEL=randinit(DATA,NCLUST,regularize,display_mode)

if nargin<4
	display_mode=1;
end

[datapoints,D]=size(DATA);

% get the variance of all dimensions, set up
% a diagonal covariance

% randinit

if display_mode
	fprintf(1,'Initializing %g cluster(s) with random points\n',NCLUST);
end

initpoints=randsample(datapoints,NCLUST);

NEWMODEL.mu=DATA(initpoints,:);

% initialize all covariance matrices

datavar=var(DATA);
initsigma=diag(datavar);

for i=1:NCLUST
	NEWMODEL.sigma(:,:,i)=initsigma+eye(D).*regularize;
	NEWMODEL.mixing(i)=1/NCLUST;
end

end
