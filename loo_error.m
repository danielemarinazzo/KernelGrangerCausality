function [loo dim] =loo_err(X,x,type,par,f)
if nargin < 5
    f=1.e-6;
end
[nvar m n]=size(X);
for j=1:nvar
    x(:,j)=x(:,j)/sqrt(x(:,j)'*x(:,j));
end
Xr=reshape(X,nvar*m,n);
[V ifail]=filtro(Xr,type,par,f);
V(n+1,:)=[];
dim=size(V,2);
for i=1:n
    g(i)=sum(V(i,:).^2);
end
vx=V'*x;
xc =x-V*vx;
for j=1:m
    xc(:,j)=xc(:,j)./(1-g)';
end
loo=sum(sum(xc.^2));