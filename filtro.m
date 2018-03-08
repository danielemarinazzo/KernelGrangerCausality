function [V ifail]=filtro(X,type,p,f)
if nargin<4
    f=1.e-6;
end
if type=='f'
    [V ifail]=fourier(X,p,f);
else
    [m n]=size(X);
    if type=='g'
        fast=false;
        rmax=n;
    else
        rmax=min(n,nchoosek(p+m,m));
        fast=nchoosek(p+m,m)<rmax+1;
    end
    %disp(sprintf('toolbox:filtro fast=%d m=%d p=%d dim=%d',fast,m,p,nchoosek(p+m,m)));
    if fast
        %polynomial kernel (p+m)!/(p!m!)-1< rmax'
        term=polypower('X',m,p);
        k=length(term);
        for i=1:k
            L(:,i)=eval(term{i});
            L(:,i)=(L(:,i)-mean(L(:,i)));
        end
        [VN D]=eig(L'*L);
        V=L*VN;
        ifail=false;
    else
        [L P ifail]=cholesky(X,type,p,rmax,f);
        [VN D]=eig(L'*L);
        V=L(P,:)*VN;
    end
    xnorm=repmat(sqrt(dot(V,V)),n,1);
    V=V./xnorm;
    [s ind]=sort(abs(diag(D)),'descend');
    r=sum(s>f*s(1));
    ind=ind(1:r);
    V=[V(:,ind);diag(D(ind,ind))'];
end