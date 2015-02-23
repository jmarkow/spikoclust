function [NEWMODEL]=partialem(DATA,MODEL,unip,maxiter,epsilon,lambda,idx)

%
%
%
%
%

garbage=0;
if ~isempty(unip)
	garbage=1;
end

% total responsibilities for normalization

mu=MODEL.mu;
sigma=MODEL.sigma;
mixing=MODEL.mixing;

R=spikoclust_gmem_estep(DATA,MODEL,unip);

[datapoints,D]=size(DATA);
NCLUST=size(mu,1);

if garbage
	normalizationR=sum(R(:,[idx NCLUST+1]),2);
else
	normalizationR=sum(R(:,[idx]),2);
end

prev_likelihood=1e-9;

for i=1:maxiter

	% compute likelihoods
	% responsibilities

	if garbage
		R=zeros(datapoints,NCLUST+1);
	else
		R=zeros(datapoints,NCLUST);
	end

	px=zeros(datapoints,1);
	den=zeros(datapoints,1);

	% e step, get the responsibilities

	sigma=spikoclust_gmem_covcheck(sigma);

	for j=idx

		pxtheta=mvnpdf(DATA,mu(j,:),sigma(:,:,j)); 
		mixprob=mixing(j)*pxtheta;
		px=px+mixprob;
		den=den+mixprob;
		R(:,j)=mixprob;

	end

	if garbage
		mixprob=mixing(NCLUST+1)*unip;
		R(:,NCLUST+1)=mixprob;
		den=den+mixprob;
		px=px+mixprob;
	end

	% likelihood

	den=normalizationR./(den+1e-300);

	for j=idx
		R(:,j)=R(:,j).*den;
	end

	if garbage
		R(:,NCLUST+1)=R(:,NCLUST+1).*den;
	end

	likelihood=sum(log(px+1e-300));
	deltalikelihood=(likelihood-prev_likelihood);

	% break if we've reached our stopping criterion

	if deltalikelihood>=0 && deltalikelihood<epsilon*abs(prev_likelihood)
		break;
	end

	prev_likelihood=likelihood;
	% update mu, sigma and mixing probabilities

	for j=idx

		% need the total r for normalization

		totalR=sum(R(:,j)); 

		% recompute mu, the inner product between all datapoints and their
		% responsibilities within the cluster, normalized by the totalR

		mu(j,:)=(DATA'*R(:,j))./totalR; 

		% get the deviation from the mean for each datapoint

		dx=(DATA-repmat(mu(j,:),[datapoints 1]))';

		% transpose so we have D x datapoints

		Rdx=repmat(R(:,j)',[D 1]).*dx;

		% now R for the cluster is repeated so we have D x datapoints

		% take the inner product between the mean deviation D x datapoints and datapoints x D 
		% for responsibilities

		% add the regularization constant and normalize

		sigma(:,:,j)=(Rdx*dx'+lambda*eye(D))/(totalR+lambda);
		mixing(j)=mean(R(:,j));
	end

	if garbage
		mixing(NCLUST+1)=mean(R(:,NCLUST+1));
	end



end

NEWMODEL.R=R;
NEWMODEL.sigma=sigma;
NEWMODEL.mu=mu;
NEWMODEL.mixing=mixing;
NEWMODEL.likelihood=likelihood;
