function CHECK=ispdm(MAT)

%%%%

[~,p]=chol(MAT);

if p==0
	CHECK=1;
else
	CHECK=0;
end


end
