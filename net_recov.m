function cr=net_recov(rxx,thr)
nvar=length(rxx);
cr=zeros(nvar,nvar,length(thr));
for k=1:length(thr)
    for i=1:nvar
        rx=rxx{i};
        [rxo indxo]=sort(abs(rx),2);
        for l=1:nvar
            if l~=i
                indxb=find(abs(rx(l,:))>thr(k));
                cr(i,l,k)=sum(rx(l,indxb).^2);
            end
        end
    end
end
cr=squeeze(cr);