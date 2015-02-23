function [merit,pairs]=mergemerit(MODEL)

%
%
%
%

clusterids=find(~MODEL.garbage);
pairs=nchoosek(clusterids,2);

for i=1:size(pairs,1)
	merit(i)=MODEL.R(:,pairs(i,1))'*MODEL.R(:,pairs(i,2));
	merit(i)=merit(i)./(norm(MODEL.R(:,pairs(i,1)))*norm(MODEL.R(:,pairs(i,2))));
end

[val,idx]=sort(merit,'descend');
merit=merit(idx);
pairs=pairs(idx,:);

