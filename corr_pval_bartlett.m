function [rr, pp, nu]=corr_pval_bartlett(xv,VN,c,U,N)
% Input:    xv      : matrix (n x nr) of target residuals;
%           VN      : matrix (n x nc) of new directions (orthogonal to V);
%           c       : dimension of the conditioning space, size(V,2);
%           U       : maximum lag of the autocorrelation   % default floor(N/5)
%           N       : length of each trial   % default N=n
% Output:   rr      : rr(i,k) correlation between xv(:,i) and VN(:,k);
%           pp      : two-tailed p-values (1 if nu<=0);
%           nu      : effective degrees of freedom
% Effective dof from Bartlett's formula, minus conditioning dimension:
% Cliff et al., Phys. Rev. Research 3, 013145 (2021), Eqs. 8, 9, 25.
%Last version (2026)
[n nr]=size(xv);
nc=size(VN,2);
if nargin<5
    N=n;
end
if nargin<4 || isempty(U)
    U=floor(N/5);
end
if mod(n,N)~=0
    error('n must be a multiple of N')
end
U=min(U,N-1);
rr=corr(xv,VN);
ax=acf_trials(xv,U,N);
ay=acf_trials(VN,U,N);
u=(1:U)';
w=(N-u)/N;
s2=(1+2*(ax.*w)'*ay)/n;
nu=1+1./s2-c-2;
pp=ones(nr,nc);
ind=nu>0;
t=rr(ind).*sqrt(nu(ind)./(1-rr(ind).^2));
pp(ind)=2*tcdf(-abs(t),nu(ind));
if any(~ind(:))
    warning('corr_pval_bartlett: %d tests with nu<=0, p set to 1',sum(~ind(:)))
end
