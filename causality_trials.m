function [cb ifail rr pp]=causality_trials(Y,type,par,m)
% Input:    Y       : matrix  (ntrials x n x nvar) of input data;
%                       ntrials = number of trials
%                       nvar = number of time series
%                       n    = number of samples
%           type    : kernel function 'p' polinomial 'g' gaussian ;
%           par     : parameter of the kernel function;
%           m       : order of the model   % default m=1;
% Output:   cb      : cb(i,j) = i->j
%           ifail   :  0 no error
%                      1 cholesky algoritm fail
%                      2 complex eigenvalues
%                      3 complex eigensystem
%                      4 error in polypower
%Last version (16/07/2012)
if nargin<4
    m=1;
end
[ntrials ntot nvar]=size(Y);
X=zeros(nvar,m,ntot);        % size of the final training matrix.
x=zeros(ntot,nvar);          % size of the final test matrix.    
for it=1:ntrials
    Yt=squeeze(Y(it,:,:));
    [Xt xt]=init_causality(Yt,m);
    X(:,:,(it-1)*(ntot-m)+1:it*(ntot-m))=Xt;
    x((it-1)*(ntot-m)+1:it*(ntot-m),:)=xt;
end
[nvar m n]=size(X);
f=1.e-6;
th=0.05;
Xr=reshape(X,nvar*m,n);
[VV D ifail]=filtro(Xr,type,par,f,true);
cb=zeros(nvar,nvar);
VT=VV*D.^0.5;
polycall=true;
kk=0;
for i=1:nvar
    XX=X;
    XX(i,:,:)=[];
    Xr=reshape(XX,(nvar-1)*m,n);
    V=filtro(Xr,type,par,f,polycall);
    polycall=false;
    [VN ifail]=vnorma(VT,V,VV);
    if ifail>0
        rr=0;
        pp=0;
        return
    end
    xv=x-V*V'*x;
    [rrt ppt]=corr(xv,VN);
    [nr nc]=size(rrt);
    rr(i,1:nr,1:nc)=rrt;
    pp(i,1:nr,1:nc)=ppt;
    rr(i,i,1:nc)=0;
    kk=kk+(nr-1)*nc;
end
rn=rr.^2;
thb=th/kk;
indpr=find(pp>thb);
rn(indpr)=0;
cb=sum(rn,3);
