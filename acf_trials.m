function a=acf_trials(Z,U,N)
% Input:    Z       : matrix (n x nz), trials concatenated along the rows;
%           U       : maximum lag;
%           N       : length of each trial   % default N=n
% Output:   a       : a(k,j) autocorrelation at lag k of Z(:,j), averaged over trials
%Last version (2026)
[n nz]=size(Z);
if nargin<3
    N=n;
end
M=n/N;
a=zeros(U,nz);
for im=1:M
    Zc=Z((im-1)*N+1:im*N,:);
    Zc=Zc-repmat(mean(Zc,1),N,1);
    d=sum(Zc.^2,1);
    for k=1:U
        a(k,:)=a(k,:)+sum(Zc(1:N-k,:).*Zc(1+k:N,:),1)./d;
    end
end
a=a/M;
