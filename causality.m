function [cb rr pp thb ifail]=causality(X,x,type,par,Np,th,f)
%         nvar=number of time series,m= order of the model,n= number of samples
% Input:    X       : matrix (nvar x m x n) of input data;
%                       nvar = number of time series
%                       m    = order of the model
%                       n    = number of samples
%           x       : matrix (n x nvar) of target data;
%           par     : parameter of the kernel function;
%           th      :
%           f       :
% Output:   cb      :
%           ifail   :
if nargin<7
    f=1.e-6;
end;
if nargin<6
    th=0.05;
end;
[nvar m n]=size(X);
if nargin<5
    Np=ones(1,nvar);
end
Nt=[0 cumsum(Np)];
ngr=length(Np);
for j=1:size(x,2)
    x(:,j)=x(:,j)/sqrt(x(:,j)'*x(:,j));
end
Xr=reshape(X,nvar*m,n);
[VV ifail]=filtro(Xr,type,par,f);
if ~ifail
    VT=VV(1:n,:)*diag(VV(n+1,:).^0.5);
    kk=0;
    for i=1:ngr
        tic;
        XX=X;
        XX(Nt(i)+1:Nt(i+1),:,:)=[];
        Xr=reshape(XX,(nvar-Np(i))*m,n);
        ind=setdiff(1:ngr,i);
        for j=ind
            for k=Nt(j)+1:Nt(j+1);
                V=filtro(Xr,type,par,f);
                xv=V(1:n,:)*(V(1:n,:)'*x(:,k));
                [rrp ppp]=corrkk(VT,V(1:n,:),VV(1:n,:),x(:,k)-xv);
                nrr=length(rrp);
                rr(i,j,k-Nt(j),1:nrr)=rrp;
                pp(i,j,k-Nt(j),1:nrr)=ppp;
                kk=kk+nrr;
            end
        end
    end
    thb=th/kk;
    indpr=find(pp>thb);
    rn=rr;
    rn(indpr)=0;
    for i=1:ngr
        cb(:,i)=sum(sum(rn(:,i,:,:).^2,4),3)/Np(i);
    end
end