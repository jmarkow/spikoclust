
function [NEWMODEL]=fullem(DATA,MODEL,unip,maxiter,epsilon,lambda)

%
%
%
%
%

[datapoints,D]=size(DATA);

mu=MODEL.mu;
sigma=MODEL.sigma;
mixing=MODEL.mixing;
NCLUST=size(mu,1);

garbage=0;
if ~isempty(unip)
	garbage=1;
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

	for j=1:NCLUST
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
    
	for j=1:NCLUST
		R(:,j)=R(:,j)./(den+1e-300);
	end

	if garbage
		R(:,NCLUST+1)=R(:,NCLUST+1)./(den+1e-300);
	end

	mixing=mean(R);
	likelihood=sum(log(px+1e-300));

	deltalikelihood=(likelihood-prev_likelihood);

	% break if we've reached our stopping criterion

	if deltalikelihood>=0 && deltalikelihood<epsilon*abs(prev_likelihood)
		break;
	end

	prev_likelihood=likelihood;
	
	% update mu, sigma and mixing probabilities
    
	for j=1:NCLUST

		% need the total r for normalization

		totalR=sum(R(:,j)); 

		% recompute mu, the inner product between all datapoints and their
		% responsibilities within the cluster, normalized by the totalR

		mu(j,:)=(DATA'*R(:,j))./(totalR+1e-300); 

		% get the deviation from the mean for each datapoint

		dx=(DATA-repmat(mu(j,:),[datapoints 1]))';

		% transpose so we have D x datapoints

		Rdx=repmat(R(:,j)',[D 1]).*dx;

		% now R for the cluster is repeated so we have D x datapoints

		% take the inner product between the mean deviation D x datapoints and datapoints x D 
		% for responsibilities

		% add the regularization constant and normalize

		sigma(:,:,j)=(Rdx*dx'+lambda*eye(D))/(totalR+lambda);
	end

	% store the likelihood


end

NEWMODEL.R=R;
NEWMODEL.sigma=sigma;
NEWMODEL.mu=mu;
NEWMODEL.mixing=mixing;
NEWMODEL.likelihood=likelihood;
