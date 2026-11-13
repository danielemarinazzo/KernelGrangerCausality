function [lam, pl, nk]=lambda_test_bartlett(xv,VN,c,U,N,nsamp)
% Input:    xv      : matrix (n x nr) of target residuals (given V);
%           VN      : matrix (n x nc) of new directions, orthonormal and orthogonal to V;
%           c       : dimension of the conditioning space, size(V,2);
%           U       : maximum lag of the autocorrelation   % default floor(N/5)
%           N       : length of each trial   % default N=n
%           nsamp   : number of Monte Carlo samples   % default 10000
% Output:   lam     : lam(i) = prod_k (1-rho_k^2) for target xv(:,i);
%           pl      : p-values P(Lambda<=lam) (1 if nk<=0);
%           nk      : matrix (nc x nr) of effective degrees of freedom
% Modified Lambda-test, null lam ~ prod_k Beta(nk/2,1/2):
% Cliff et al., Phys. Rev. Research 3, 013145 (2021), Eqs. 30, 44-50.
%Last version (2026)
[n, nr]=size(xv);
nc=size(VN,2);
if nargin<5 || isempty(N)
    N=n;
end
if nargin<4 || isempty(U)
    U=floor(N/5);
end
if nargin<6
    nsamp=10000;
end
if mod(n,N)~=0
    error('n must be a multiple of N')
end
U=min(U,N-1);
u=(1:U)';
w=(N-u)/N;
ay=acf_trials(VN,U,N);
lam=ones(1,nr);
pl=ones(1,nr);
nk=zeros(nc,nr);
nund=0;
for i=1:nr
    e=xv(:,i);
    ex=zeros(n,nc);
    for k=1:nc
        ex(:,k)=e;
        e=e-VN(:,k)*(VN(:,k)'*e);
    end
    lam(i)=sum(e.^2)/sum(xv(:,i).^2);
    if nc==0
        continue
    end
    ax=acf_trials(ex,U,N);
    s2=(1+2*sum(ax.*ay.*repmat(w,1,nc),1))/n;
    nk(:,i)=(1+1./s2)'-(c+(0:nc-1)')-2; % Check if the dof are 0 or negative
    if any(nk(:,i)<=0)
        nund=nund+1;
        continue
    end
    logL=sum(log(betarnd(repmat(nk(:,i)'/2,nsamp,1),0.5)),2);
    pl(i)=(sum(logL<=log(lam(i)))+1)/(nsamp+1);
end
if nund>0
    warning('lambda_test_bartlett: %d targets with nk<=0, p set to 1',nund)
end
