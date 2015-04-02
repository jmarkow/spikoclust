function SIGMA=cov_check(SIGMA)

[~,~,NCLUST]=size(SIGMA);

for i=1:NCLUST
	SIGMA(:,:,i)=spikoclust_gmem_covfix(SIGMA(:,:,i));
end
