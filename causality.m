function [cb, ifail, rr, pp]=causality(Y,type,par,m,corrtype)
% Input:    Y       : matrix  (n x nvar) of input data;
%                       nvar = number of time series
%                       n    = number of samples
%           type    : kernel function 'p' polinomial 'g' gaussian ;
%           par     : parameter of the kernel function;
%           m       : order of the model   % default m=1;
%           corrtype: type of multiple comparison correction ('Bonferroni'(default) or 'FDR');
% Output:   cb      : cb(i,j) = i->j
%           ifail   :  0 no error
%                      1 cholesky algoritm fail
%                      2 complex eigenvalues
%                      3 complex eigensystem
%                      4 error in polypower
%Last versione (11/01/2011)
if nargin<4
    m=1;
    corrtype='Bonferroni';
end
if nargin<5
    corrtype='Bonferroni';
end
[X, x]=init_causality(Y,m);
[nvar, m, n]=size(X);
f=1.e-6;
th=0.05;
Xr=reshape(X,nvar*m,n);
[VV, D, ifail]=filtro(Xr,type,par,f,true);
cb=zeros(nvar,nvar);
VT=VV*D.^0.5;
polycall=true;
kk=0;

%%% this loop just to initialize the rr matrix, but it fails sometimes, removed for now %%%
% XX=X;
% XX(1,:,:)=[];
% Xr=reshape(XX,(nvar-1)*m,n);
% V=filtro(Xr,type,par,f,polycall);
% polycall=false;
% [VN, ifail]=vnorma(VT,V,VV);
% if ifail>0
%     rr=0;
%     pp=0;
%     return
% end
% xv=x-V*V'*x;
% [rrt, ppt]=corr(xv,VN);
% [nr, nc]=size(rrt);
% rr=zeros(nvar,nr,nc);pp=rr;
%%%%%%
for i=1:nvar
    XX=X;
    XX(i,:,:)=[];
    Xr=reshape(XX,(nvar-1)*m,n);
    V=filtro(Xr,type,par,f,polycall);
    polycall=false;
    [VN, ifail]=vnorma(VT,V,VV);
    if ifail>0
        rr=0;
        pp=0;
        disp('vnorma failed')
        return
    end
    xv=x-V*V'*x;
    [rrt, ppt]=corr(xv,VN);
    [nr, nc]=size(rrt);
    rr(i,1:nr,1:nc)=rrt;
    pp(i,1:nr,1:nc)=ppt;
    rr(i,i,1:nc)=0;
    kk=kk+(nr-1)*nc;
end
rn=rr.^2;
if corrtype=="Bonferroni"
    thb=th/kk;
    indpr=pp>thb;
    rn(indpr)=0;
    cb=sum(rn,3);
elseif corrtype=="FDR"
    h=fdr_bh(pp,th); % the output is 0 (nonsignificant) or 1 (significant)
    rn=rn.*h;
    cb=sum(rn,3);
else
    error('corrtype must be "Bonferroni" or "FDR"')
end
