function NEWMODEL=randinit(DATA,NCLUST,regularize)

[datapoints,D]=size(DATA);

% get the variance of all dimensions, set up
% a diagonal covariance

% randinit

fprintf(1,'Initializing %g cluster(s) with random points\n',NCLUST);

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

