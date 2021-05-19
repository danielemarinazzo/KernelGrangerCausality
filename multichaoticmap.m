function xx=multichaoticmap(n,c,s,mu,eps)
% n = length of the series
% c is the connectivity matrix 
% Ex. c=zeros(3);c(1,3)=0.2;c(2,1)=0.2;
% s=noise variance (0.01)
% mu = 1.8 in our case
% eps=1 logistic map
nvar=length(c);
nt=1000; %points to remove to ensure independence from initial conditions
xx=zeros(nvar,n+nt);
cp=c+diag(1-sum(c));
test=1;
while test>0
    nn=s*randn(nvar,n+nt-1);
    xx(:,1)=rand(nvar,1);
    for i=2:n+nt
        xx(:,i)=cp'*(1-mu*abs(xx(:,i-1)).^(1+eps))+nn(:,i-1);
    end
    X=xx(:,nt:n+nt-1);
    x=xx(:,nt+1:n+nt)';
    test=isinf(max(max(abs(xx)))^2);
end
xx(:,1:nt)=[];