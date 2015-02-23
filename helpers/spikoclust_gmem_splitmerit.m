function [split,splitidx]=splitmerit(DATA,MODEL,unip)
%
%
%

[datapoints,D]=size(DATA);
clusterids=find(~MODEL.garbage);

R=spikoclust_gmem_estep(DATA,MODEL,unip);

% take the responsibilities, compare with the model PDFs

mu=MODEL.mu;
sigma=spikoclust_gmem_covcheck(MODEL.sigma);

for i=clusterids

	pxtheta=mvnpdf(DATA,mu(i,:),sigma(:,:,i)); % point x 1 vector
	
	% "empirical" density
	
	f=R(:,i)./sum(R(:,i));

	pxtheta=pxtheta+1e-5;
	idx=find(f>1e-5);

	% get KL divergence

	split(i)=sum(f(idx).*log(f(idx)./pxtheta(idx)));

end

[split splitidx]=sort(split,'descend');

end
