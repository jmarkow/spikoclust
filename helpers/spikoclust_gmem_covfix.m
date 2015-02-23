function SIGMA=cov_fix(SIGMA)

D=size(SIGMA,1);

% much of this is cribbed from the gmmbayes package

maxloops=1e6; % to make sure we don't get stuck
nloops=0;
tolerance=eps*10; 
fix_per=.01; % percent correction along the diagonal per loop
cov_fix_mat = ones(D) + fix_per*eye(D);
% ensure all elements are finite

if ~all(isfinite(SIGMA(:)))
	error('Not all covariance entries are finite');
end

%%% check for positive-definiteness

% enforce symmetry

SIGMA=(SIGMA+SIGMA')/2;

while ~spikoclust_gmem_ispdm(SIGMA) & nloops<maxloops

	% continue to add elements to the diagonal 

	nloops=nloops+1;
	d=diag(SIGMA); 

	if any(d<=tolerance)
		m=max(abs(d))*fix_per;
		neg=min(d);
		if neg<0
			to_add=(m-neg)*eye(D);
		else
			if m<tolerance
				m=tolerance;
			end
			to_add=m*eye(D);
		end

		SIGMA=SIGMA+to_add;

	else
		SIGMA=SIGMA.*cov_fix_mat;
	end
end

if nloops>maxloops
	error('Could not fix covariance matrix.');
end

end

