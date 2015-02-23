function R=estep(DATA,MODEL,unip)
%
%
%
%
%
%
%
% performs the e step

[datapoints,D]=size(DATA);

garbage=0;
if ~isempty(unip)
	garbage=1;
end

mu=MODEL.mu;
sigma=spikoclust_gmem_covcheck(MODEL.sigma);
mixing=MODEL.mixing;

NCLUST=size(mu,1);

px=zeros(datapoints,1);
den=zeros(datapoints,1);

if ~garbage
	R=zeros(datapoints,NCLUST);
else
	R=zeros(datapoints,NCLUST+1);
end

for i=1:NCLUST

	pxtheta=mvnpdf(DATA,mu(i,:),sigma(:,:,i)); % point x 1 vector
	mixprob=mixing(i)*pxtheta;
	px=px+mixprob;
	den=den+mixprob;
	R(:,i)=mixprob;

end

if garbage

	mixprob=mixing(NCLUST+1)*unip;
	R(:,NCLUST+1)=mixprob;
	den=den+mixprob;
	px=px+mixprob;

end

for i=1:NCLUST
	R(:,i)=R(:,i)./(den+1e-300);
end

if garbage
	R(:,NCLUST+1)=R(:,NCLUST+1)./(den+1e-300);
end

end
