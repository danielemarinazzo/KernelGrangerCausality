function [nlink err eff imp]=efficienza(cr,ct);
nvar=size(cr,1);
nth=size(cr,3);
nlink=zeros(nth,1);
err=zeros(nth,1);
eff=zeros(nth,1);
imp=zeros(nth,1);
if nargin==1
    for i=1:nth
        nlink(i)=length(find(abs(cr(:,:,i))));
    end
else
    indt=find(abs(ct));
    nl=length(indt);
    indt0=find(abs(ct)==0);
    eff=zeros(nth,1);
    imp=zeros(nth,1);
    for i=1:nth
        indr=find(abs(cr(:,:,i)));
        indr0=find(abs(cr(:,:,i))==0);
        nl0=length(indt0);
        tp=length(intersect(indt,indr));
        tn=length(intersect(indt0,indr0));
        fp=length(indr)-tp;
        fn=nl-tp;
        err(i)=nvar^2-tp-tn;
        if nl*nl0>0
            eff(i)=100*tp/nl;
            imp(i)=100*(1-tn/nl0);
        end
        nlink(i)=length(indr);
    end
end