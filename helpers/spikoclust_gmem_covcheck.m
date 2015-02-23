function SIGMA=cov_check(SIGMA)

[~,~,NCLUST]=size(SIGMA);

for i=1:NCLUST
	p=spikoclust_gmem_ispdm(SIGMA(:,:,i));

	if ~p
		disp('Fixing covariance...');
		SIGMA(:,:,i)=spikoclust_gmem_covfix(SIGMA(:,:,i));
	end
end

end
