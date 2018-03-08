function [term r]=polypower(varname,m,p)
%disp(['toolbox:polypower:    p=' num2str(p)]);
sym 1;
mystr='syms';
command='1';
for i=1:m
    mystr=strcat(mystr,' a',num2str(i));
    command=(strcat(command,'+a',num2str(i)));
end
eval(mystr);
s=eval(command);
r=char(expand(s^p-1));
for i=m:-1:1
    ch=num2str(i);
    r=strrep(r,['a' ch],[varname '(' ch ',:)']);
end
r=strrep(r, '^','.^');
r=strrep(r, '*','.*');
ind=[0 strfind(r,'+') length(r)+1];
k=length(ind)-1;
term=cell(k,1);
for i=1:k
    rr=r(ind(i)+1:ind(i+1)-1);
    iend=length(rr);
    if isstrprop(rr(1),'digit')
        kk=strfind(rr,'*');
        rr=strcat('sqrt(',rr(1:kk(1)-1),')',rr(kk(1):iend));
    end
    term{i}=rr;
end
