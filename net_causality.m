function [cb cr mmxx rxcum pxcum ifail rx px]=net_causality(X,x,type,par,thr,thb,th,f)
%         nvar=number of time series,m= order of the model,n= number of samples
% Input:    X       : matrix (nvar x m x n) of input data;
%                       nvar = number of time series
%                       m    = order of the model
%                       n    = number of samples
%           x       : matrix (n x nvar) of target data;
%           type    : type of kernel function to use: 'p' polynomial, 'g' gaussian;
%           par     : parameter of the kernel function;
%           thr     : 
%           thb     :
%           th      :
%           f       :
% Output:   cb      :
%           cr      :
%           rxcum   :
%           pxcum   :
%           ifail   :
%           rx      :
%           px      :
if nargin<8
    f=1.e-6;
end;
if nargin<7
    th=0.05;
end;
if nargin<6
    thb=1;
end;
if nargin<5
    thr=0:0.01:1;
end;
[nvar m n]=size(X);
cb=zeros(nvar);
c=zeros(nvar,nvar,length(thr));
for j=1:nvar
    x(:,j)=x(:,j)/sqrt(x(:,j)'*x(:,j));
end
Xr=reshape(X,nvar*m,n);
[VV ifail]=filtro(Xr,type,par,f);
if ~ifail
    VT=VV(1:n,:)*diag(VV(n+1,:).^0.5);
    kk=0;
    for i=1:nvar
        XX=X;
        XX(i,:,:)=[];
        Xr=reshape(XX,(nvar-1)*m,n);
        V=filtro(Xr,type,par,f);
        xv=V(1:n,:)*(V(1:n,:)'*x);
        [rx{i} px{i}]=corrkk(VT,V(1:n,:),VV(1:n,:),x-xv);
        rr=rx{i};
        pp=rx{i};
        rr(i,:)=[];
        pp(i,:)=[];
        [k1 k2]=size(rr);
        rxcum(kk+1:kk+k1*k2)=reshape(rr,k1*k2,1);
        pxcum(kk+1:kk+k1*k2)=reshape(pp,k1*k2,1);
        kk=kk+k1*k2;
        mmxx(i)=size(rx{i},2);
    end
    if thb==1
        thb=2*th/(sum(mmxx)*(nvar-1));
    end
    cr=net_recov(rx,thr);
    for i=1:nvar
        for l=1:nvar
            if l~=i
                % Bonferroni test
                indxb=find(px{i}(l,:)<thb);
                cb(i,l)=sum(rx{i}(l,indxb).^2);
           end
        end
    end
end