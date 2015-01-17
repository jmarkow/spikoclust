function CONTAMINATION=spikoclust_cluster_contamination(SPIKEDATA,MODEL)
%

nclust=size(MODEL.mu,1);
npoints=1e4; % point from each cluster
CONTAMINATION=ones(1,nclust).*NaN;

% or generate points from each cluster
%

total_points=nclust*npoints;
labels=zeros(total_points,1);
if iscell(SPIKEDATA)
	D=size(SPIKEDATA{1},2);
else
	D=size(SPIKEDATA,2);
end
points=zeros(total_points,D);

for i=1:nclust
	startpoint=((i-1)*npoints);
	labels(startpoint+1:startpoint+npoints)=i;

	% generate using uniform density spanning entire space
	points(startpoint+1:startpoint+npoints,:)=mvnrnd(MODEL.mu(i,:),MODEL.sigma(:,:,i),npoints);
end

for i=1:nclust
	prob(:,i)=mvnpdf(points,MODEL.mu(i,:),MODEL.sigma(:,:,i));
end

prob=prob.*repmat(MODEL.mixing(1:nclust),[total_points 1]);
[~,classification]=max(prob,[],2);

for i=1:size(MODEL.mu,1)

	% how many false positives and false negatives
	%
	fp=sum((classification==i)&(labels~=i));

	% false negatives
	%
	fn=sum((classification~=i)&(labels==i));
	CONTAMINATION(i)=(fp/npoints)+(fn/npoints);
end
