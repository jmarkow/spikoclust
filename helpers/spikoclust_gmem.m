function [newmodel]=spikoclust_gmem(DATA,INIT,NCLUST,varargin)
%EM for Gaussian mixture with support for outlier distribution
%
% all credit to Daniel Wagenaar's code
% http://www.danielwagenaar.net/res/papers/00-Wage2.pdf
%
% algorithm derived from Maneesh Sahani's thesis

% data should be passed as a matrix observations x variables

if nargin<3 | isempty(NCLUST)
	NCLUST=2;
end

[datapoints,D]=size(DATA);
nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

%%%

maxiter=100;
regularize=1e-6;
epsilon=1e-10;
lambda=.05; % changed from .01 to .05 9/19/13
garbage=1;
merge=1;
splitepsi=1; % noise scale for splits
maxcand=5; % maximum number of smem candidates
smemiter=100; % maximum smem iterations
debug=0;
display_mode=1;

for i=1:2:nparams
	switch lower(varargin{i})
		case 'regularize'
			regularize=varargin{i+1};
		case 'maxiter'
			maxiter=varargin{i+1};
		case 'lambda'
			lambda=varargin{i+1};
		case 'epsilon'
			epsilon=varargin{i+1};
		case 'maxiter'
			maxiter=varargin{i+1};
		case 'garbage'
			garbage=varargin{i+1};
		case 'beta'
			beta=varargin{i+1};
		case 'merge'
			merge=varargin{i+1};
		case 'debug'
			debug=varargin{i+1};
		case 'display_mode'
			display_mode=varargin{i+1};
		otherwise
	end
end


if nargin<2 | isempty(INIT)
	INIT=spikoclust_gmem_randinit(DATA,NCLUST,regularize,display_mode);
end

mu=INIT.mu;
mixing=INIT.mixing;
sigma=INIT.sigma;

% make sure dimensions are correct

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% CHECK DIMENSIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(mu,2)~=D
	error('Mean dimensionality %g does not match data point dimensions %g',size(mu,2),D);
end

if length(mixing)~=NCLUST
	error('Wrong number of mixing proportions %g, NCLUST=%g',length(mixing),NCLUST);
end

if size(sigma,1)~=size(sigma,2)
	error('Covariance matrix must have the same dimensions...');
end

if size(sigma,3)~=NCLUST
	error('Must have covariance matrices for all clusters...');
end

if size(sigma,1)~=D
	error('Covariance matrix has wrong dimensionality %g, D=%g',size(sigma,1),D);
end

% number of parameters (full covariance)

% df per component is (1+d+d(d+1)/2)

nparams=(NCLUST*D*((D+1)/2));
nparams=nparams+NCLUST-1+NCLUST*D;

%nparams=NCLUST*(1+D+(D*(D+1))/2); % new parameter number

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% FULL EM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% return the parameters


if garbage
	[unip,mixing]=spikoclust_gmem_initgarbage(DATA,NCLUST,INIT);
	INIT.mixing=mixing;
else
	unip=[];
end

if debug
	fig=figure();
	gaussvis(INIT,DATA,'fig_num',fig);
	pause(.5);
end


[newmodel]=spikoclust_gmem_fullem(DATA,INIT,unip,maxiter,epsilon,lambda);

if debug
	spikoclust_gaussvis(newmodel,DATA,'fig_num',fig);
	drawnow;
	pause(.5);
end

newmodel.garbage=zeros(1,size(newmodel.R,2));
if garbage
	newmodel.garbage(end)=1;
end

% get the merge candidates
% perform SMEM

if merge & NCLUST>2

	% keep merging until BIC no longer improves

	breakflag=0;

	for ii=1:smemiter

		% get the merge candidates

		[merges,pairs]=spikoclust_gmem_mergemerit(newmodel);

		% get the split candidates

		[~,splits]=spikoclust_gmem_splitmerit(DATA,newmodel,unip);

		% go through each pair find a non-matching split
		% now for each pair we take the split candidates ~==pair


		triplet=[];
		for i=1:size(pairs,1)

			currpair=pairs(i,:);

			for j=1:length(splits)
				if ~any(splits(j)==currpair)
					triplet=[triplet;currpair splits(j)];
				end
			end
		end

		for i=1:maxcand

			if i>size(triplet,1)
				break;
			end

			currtrip=triplet(i,:);
			if display_mode
				fprintf(1,'Merging %g and %g, splitting %g\n',currtrip(1),currtrip(2),currtrip(3));
			end


			mergemodel1=spikoclust_gmem_mergeclust(newmodel,...
				currtrip(1),currtrip(2));
			mergemodel1=spikoclust_gmem_splitclust(mergemodel1,...
				currtrip(3),currtrip(2),splitepsi);

			% run partial em on our new merged cluster

			mergemodel1=spikoclust_gmem_partialem(DATA,mergemodel1,...
				unip,maxiter,epsilon,lambda,[currtrip]);
			mergemodel1=spikoclust_gmem_fullem(DATA,mergemodel1,...
				unip,maxiter,epsilon,lambda);

			if debug
				spikoclust_gaussvis(mergemodel1,DATA,'fig_num',fig);
				drawnow;
				pause(.5);
			end

			mergemodel1.garbage=newmodel.garbage;


			if mergemodel1.likelihood<newmodel.likelihood

				if display_mode
					fprintf(1,'No improvement in likelihood trying another candidate\n');
				end
				continue;

			else

				if display_mode
					fprintf(1,'Candidate improved, breaking out of inner SMEM loop\n');
				end
				newmodel=mergemodel1;

				% break out of the candidate loop if we've improved

				break;
			end

		end

		% if we've ended and there's no improvement, break out of the main loop

		if mergemodel1.likelihood<newmodel.likelihood

			if display_mode
				fprintf(1,'No improvements found in any candidates, breaking out of outer SMEM loop\n');
			end
			break;
		end

		% run partial em on the new cluster

	end

end

if debug
	spikoclust_gaussvis(newmodel,DATA,'fig_num',fig);
	title('Final');
	drawnow;
	pause(.5);
end


newmodel.BIC=-2*newmodel.likelihood+log(datapoints)*nparams;
% get the total entropy

newmodel.sigma=spikoclust_gmem_covcheck(newmodel.sigma);

for i=1:NCLUST
	pxtheta=mvnpdf(DATA,newmodel.mu(i,:),newmodel.sigma(:,:,i));
	tmp=newmodel.R(:,i).*pxtheta;
	tmp(tmp<0)=[];
	entropy(i)=sum(tmp.*log(tmp+1e-300));
end

newmodel.ICL=newmodel.BIC-2*sum(entropy);

% mml

nparams=1+D+.5*D*(D+1);

rightterm=-.5*nparams*sum(log((newmodel.mixing.*datapoints)/12))+...
	(NCLUST/2)*log(datapoints/12)+(NCLUST*(nparams+1))/2;
leftterm=newmodel.likelihood;
newmodel.MML=leftterm+rightterm;

%newmodel.MML2=.5*sum(log(newmodel.mixing))+((NCLUST*nparams+nparams)/2)*log(datapoints)-newmodel.likelihood;

if display_mode
	fprintf(1,'NComponents %g, Likelihood %5.4e, BIC %5.4e, MML %5.4e, ICL %5.4e\n',NCLUST,...
		newmodel.likelihood,newmodel.BIC,newmodel.MML,newmodel.ICL);
end
