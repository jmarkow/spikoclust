function [NEWMODEL]=mergeclust(MODEL,c1,c2)
%
%
%
%


% merge into c1


mu=MODEL.mu;
mixing=MODEL.mixing;
sigma=MODEL.sigma;

mu(c1,:)=(mixing(c1)*mu(c1,:)+mixing(c2)*mu(c2,:))./(mixing(c1)+mixing(c2));
mixing(c1)=mixing(c1)+mixing(c2);
sigma(:,:,c1)=(sigma(:,:,c1)+sigma(:,:,c2))./2;

% set the mixing proportion of the merged cluster to 0

mixing(c2)=0;

NEWMODEL.mu=mu;
NEWMODEL.sigma=sigma;
NEWMODEL.mixing=mixing;

end

