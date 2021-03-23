function [RDec,tDec,nDec,alpha]=decomposeHomography(H)
RDec=[];
tDec=[];
nDec=[];

[U S V]=svd(H);
alpha=1/S(2,2);
H=H*alpha;

[U S V]=svd(H);


ro2=S(1,1);
ro3=S(3,3);

v1=V(:,2);
v2=V(:,1);
v3=V(:,3);


k=ro2-ro3;
p=ro2*ro3-1;


lambda1=S(1,1);
lambda2=S(2,2);
lambda3=S(3,3);

gamma=(-k+sqrt(k*k+4*(p+1)))/(2*k*(p+1));
theta=(-k-sqrt(k*k+4*(p+1)))/(2*k*(p+1));

%gamma=(-1+sqrt(1+4*(lambda1*lambda3)/((lambda1-lambda3)*(lambda1-lambda3))))/(2*lambda1*lambda3)
%theta=(-1-sqrt(1-4*(lambda1*lambda3)/((lambda1-lambda3)*(lambda1-lambda3))))/(2*lambda1*lambda3)

v2norm=gamma*gamma*(lambda1-lambda3)*(lambda1-lambda3)+2*gamma*(lambda1*lambda3-1)+1;
v3norm=theta*theta*(lambda1-lambda3)*(lambda1-lambda3)+2*theta*(lambda1*lambda3-1)+1;

v2=sqrt(v2norm)*v2;
v3=sqrt(v3norm)*v3;

t0=(v2+v3)/(gamma-theta);
n=(gamma*v3+theta*v2)/(gamma-theta);

R=H*inv(eye(3)+t0*n');

if (trace((R*R'-eye(3))'*(R*R'-eye(3)))<1e-5)
	RDec(1:3,size(RDec,2)+1:size(RDec,2)+3)=R;
	tDec(1:3,size(tDec,2)+1)=R*t0;
	nDec(1:3,size(nDec,2)+1)=n;
end

t0=(v2+v3)/(gamma-theta);
n=-(gamma*v3+theta*v2)/(gamma-theta);

R=H*inv(eye(3)+t0*n');

if (trace((R*R'-eye(3))'*(R*R'-eye(3)))<1e-5)
	RDec(1:3,size(RDec,2)+1:size(RDec,2)+3)=R;
	tDec(1:3,size(tDec,2)+1)=R*t0;
	nDec(1:3,size(nDec,2)+1)=n;
end


t0=-(v2+v3)/(gamma-theta);
n=-(gamma*v3+theta*v2)/(gamma-theta);


R=H*inv(eye(3)+t0*n');

if (trace((R*R'-eye(3))'*(R*R'-eye(3)))<1e-5)
	RDec(1:3,size(RDec,2)+1:size(RDec,2)+3)=R;
	tDec(1:3,size(tDec,2)+1)=R*t0;
	nDec(1:3,size(nDec,2)+1)=n;
end

t0=-(v2+v3)/(gamma-theta);
n=(gamma*v3+theta*v2)/(gamma-theta);


R=H*inv(eye(3)+t0*n');

if (trace((R*R'-eye(3))'*(R*R'-eye(3)))<1e-5)
	RDec(1:3,size(RDec,2)+1:size(RDec,2)+3)=R;
	tDec(1:3,size(tDec,2)+1)=R*t0;
	nDec(1:3,size(nDec,2)+1)=n;
end


t0=(v2-v3)/(gamma-theta);
n=(gamma*v3-theta*v2)/(gamma-theta);


R=H*inv(eye(3)+t0*n');
if (trace((R*R'-eye(3))'*(R*R'-eye(3)))<1e-5)
	RDec(1:3,size(RDec,2)+1:size(RDec,2)+3)=R;
	tDec(1:3,size(tDec,2)+1)=R*t0;
	nDec(1:3,size(nDec,2)+1)=n;
end

t0=(v2-v3)/(gamma-theta);
n=-(gamma*v3-theta*v2)/(gamma-theta);


R=H*inv(eye(3)+t0*n');
if (trace((R*R'-eye(3))'*(R*R'-eye(3)))<1e-5)
	RDec(1:3,size(RDec,2)+1:size(RDec,2)+3)=R;
	tDec(1:3,size(tDec,2)+1)=R*t0;
	nDec(1:3,size(nDec,2)+1)=n;
end


t0=-(v2-v3)/(gamma-theta);
n=-(gamma*v3-theta*v2)/(gamma-theta);

R=H*inv(eye(3)+t0*n');

if (trace((R*R'-eye(3))'*(R*R'-eye(3)))<1e-5)
	RDec(1:3,size(RDec,2)+1:size(RDec,2)+3)=R;
	tDec(1:3,size(tDec,2)+1)=R*t0;
	nDec(1:3,size(nDec,2)+1)=n;
end

t0=-(v2-v3)/(gamma-theta);
n=(gamma*v3-theta*v2)/(gamma-theta);

R=H*inv(eye(3)+t0*n');

if (trace((R*R'-eye(3))'*(R*R'-eye(3)))<1e-5)
	RDec(1:3,size(RDec,2)+1:size(RDec,2)+3)=R;
	tDec(1:3,size(tDec,2)+1)=R*t0;
	nDec(1:3,size(nDec,2)+1)=n;
end

end
