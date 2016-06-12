function [PROJ,PCS,LAMBDA,MODEL]=spikoclust_robpca(DATA,K,REPLICATES)
%
%
%
%
%


if nargin<3 | isempty(REPLICATES)
	REPLICATES=10;
end

if nargin<2 | isempty(K)
	K=6;
end

likelihood=zeros(1,REPLICATES);
for i=1:REPLICATES
	tmp_newmodel{i}=spikoclust_gmem(DATA,[],1,'garbage',1,'merge',0,'debug',0,'display_mode',0);
	likelihood(i)=tmp_newmodel{i}.likelihood;
end

% choose the model with the highest likelihood

[~,loc]=max(likelihood);
MODEL=tmp_newmodel{loc(1)};
[PCS,LAMBDA]=svds(MODEL.sigma(:,:,1),K,'largest');
LAMBDA=diag(LAMBDA);
[~,idx]=sort(LAMBDA,'descend');
PCS=PCS(:,idx);
LAMBDA=LAMBDA(idx);

PROJ=DATA*PCS;
