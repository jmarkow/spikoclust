function [PROJ,PCS,LAMBDA,MODEL]=cov_check(DATA,REPLICATES)
%
%
%
%
%

if nargin<2
	REPLICATES=10;
end

likelihood=zeros(1,REPLICATES);
for i=1:REPLICATES
	tmp_newmodel{i}=spikoclust_gmem(DATA,[],1,'garbage',1,'merge',0,'debug',0,'display_mode',0);
	likelihood(i)=tmp_newmodel{i}.likelihood;
end

% choose the model with the highest likelihood

[~,loc]=max(likelihood);
MODEL=tmp_newmodel{loc(1)};
[PCS,LAMBDA]=eigs(MODEL.sigma(:,:,1));
LAMBDA=diag(LAMBDA);
[~,idx]=sort(LAMBDA,'descend');
PCS=PCS(:,idx);
LAMBDA=LAMBDA(idx);

PROJ=-DATA*PCS;
