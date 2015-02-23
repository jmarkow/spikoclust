function [MODEL]=splitclust(MODEL,sc,sc1,splitepsi)

% take the largest eigenvalue of the covariance

[NCLUST,D]=size(MODEL.mu);

[U DD V]=svd(MODEL.sigma(:,:,sc));

% put the new cluster in the beginning

MODEL.mixing(sc)=MODEL.mixing(sc)/2;
MODEL.mixing(sc1)=MODEL.mixing(sc);

sd=sqrt(diag(D));

oldmu=MODEL.mu(sc,:);

MODEL.mu(sc,:)=oldmu+(splitepsi*U*(sd*randn(D,1)))';
MODEL.mu(sc1,:)=oldmu+(splitepsi*U*(sd*randn(D,1)))';

MODEL.sigma(:,:,sc)=DD(1,1)*eye(D);
MODEL.sigma(:,:,sc1)=MODEL.sigma(:,:,sc);

end
