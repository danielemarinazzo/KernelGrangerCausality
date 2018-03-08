function [X x]=prepare(XX,n,m,tau,nt)
if nargin<5
    nt=1;
end
if nargin<4
    tau=1;
end
if nargin<3
    m=1;
end
delta=0;
nvar=size(XX,2);
nt0=(m-1)*tau;
for i=1:nvar;
    XX(:,i)=(XX(:,i)-mean(XX(:,i)))/std(XX(:,i));
end
X=zeros(nvar,m,n);
x=zeros(n,nvar);
i=0;
% disp([nt nt0]);
for i=1:n
    for k=1:nvar
        for j=1:m
            X(k,j,i)=XX(nt+nt0-(j-1)*tau+i-1-delta,k);
        end
        x(i,k)=XX(nt+nt0+tau+i-1-delta,k);
    end
end
% media nulla
Xm=reshape(repmat(mean(X,3),1,n),nvar,m,n);
X=X-Xm;
xm=repmat(mean(x,1),n,1);
x=x-xm;
x=x./repmat(sqrt(diag(x'*x)'),n,1);
