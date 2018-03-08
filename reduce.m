function VN=reduce(VT,V,VV)
f=1.e-8;
n=size(V,1)-1;
fast=size(VV,2)<n-2;
if fast
    A=VV'*VT;
    B=(VV'*V)*(V'*VT);
    KKN=A'*A-B*A-A'*B'+B*B';
    [VVN D]=eig(KKN);
    d=diag(D);
    ind=find(d>f*max(d));
    VN=VV*VVN(:,ind);
else
    K=VT*VT';
    P=V*V';
    KT=(K-P*K-K*P+P*K*P);
    [VN D]=eig(KT);
    d=diag(D);
    ind=find(lam>f*max(lam));
    VN=VN(:,ind);
end