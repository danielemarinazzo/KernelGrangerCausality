function [cb cr rxcum pxcum ifail rxx pxx]=bivariate(X,x,type,par,thr,th,f);
%disp('toolbox:bivariate');
[nvar m n]=size(X);
rxcum=[];
pxcum=[];
thb=2*th/(nvar*(nvar-1));
for i=1:nvar
    for j=i+1:nvar
        Y=[X(i,:,:); X(j,:,:)];
        y=[x(:,i) x(:,j)];
        [cbb crr dimx rxc pxc ifail rx px]=net_causality(Y,y,type,par,thr,thb,th,f);
        np=length(px);
        for k=1:np
            [nr nk]=size(rx{k});
            pxx(i,j,k,1:nr,1:nk)=px{k};
            rxx(i,j,k,1:nr,1:nk)=rx{k};
        end
        cb(i,j)=cbb(1,2);
        cb(j,i)=cbb(2,1);
        cr(i,j,:)=crr(1,2,:);
        cr(j,i,:)=crr(2,1,:);
        rxcum=[rxcum rxc];
        pxcum=[pxcum pxc];
    end
end
