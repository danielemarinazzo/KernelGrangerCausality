function [V ifail]=fourier(X,par,f)
[m n]=size(X);
nn=m*(m-1)/2;
switch par
    case 1
        for i=1:m
            L(:,i)=cos(X(i,:));
            L(:,i)=(L(:,i)-mean(L(:,i)))/sqrt((L(:,i)-mean(L(:,i)))'*(L(:,i)-mean(L(:,i))));
            L(:,i+m)=sin(X(i,:));
            L(:,i+m)=(L(:,i+m)-mean(L(:,i+m)))/sqrt((L(:,i+m)-mean(L(:,i+m)))'*(L(:,i+m)-mean(L(:,i+m))));
        end
    case 2
        k=0;
        for i=1:m
            for j=i+1:m
                k=k+1;
                L(:,k)=cos(X(i,:)-X(j,:));
                L(:,k)=(L(:,k)-mean(L(:,k)))/sqrt((L(:,k)-mean(L(:,k)))'*(L(:,k)-mean(L(:,k))));
                L(:,k+nn)=sin(X(i,:)-X(j,:));
                L(:,k+nn)=(L(:,k+nn)-mean(L(:,k+nn)))/sqrt((L(:,k+nn)-mean(L(:,k+nn)))'*(L(:,k+nn)-mean(L(:,k+nn))));
            end
        end
    case 3
        for i=1:m
            L(:,i)=cos(X(i,:));
            L(:,i)=(L(:,i)-mean(L(:,i)))/sqrt((L(:,i)-mean(L(:,i)))'*(L(:,i)-mean(L(:,i))));
            L(:,i+m)=sin(X(i,:));
            L(:,i+m)=(L(:,i+m)-mean(L(:,i+m)))/sqrt((L(:,i+m)-mean(L(:,i+m)))'*(L(:,i+m)-mean(L(:,i+m))));
        end
        k=2*m;
        for i=1:m
            for j=i+1:m
                k=k+1;
                L(:,k)=cos(X(i,:)-X(j,:));
                L(:,k)=(L(:,k)-mean(L(:,k)))/sqrt((L(:,k)-mean(L(:,k)))'*(L(:,k)-mean(L(:,k))));
                L(:,k+nn)=sin(X(i,:)-X(j,:));
                L(:,k+nn)=(L(:,k+nn)-mean(L(:,k+nn)))/sqrt((L(:,k+nn)-mean(L(:,k+nn)))'*(L(:,k+nn)-mean(L(:,k+nn))));
            end
        end
        case 4
            for i=1:m
            L(:,i)=cos(X(i,:));
            L(:,i)=(L(:,i)-mean(L(:,i)))/sqrt((L(:,i)-mean(L(:,i)))'*(L(:,i)-mean(L(:,i))));
            L(:,i+m)=sin(X(i,:));
            L(:,i+m)=(L(:,i+m)-mean(L(:,i+m)))/sqrt((L(:,i+m)-mean(L(:,i+m)))'*(L(:,i+m)-mean(L(:,i+m))));
            L(:,i+2*m)=cos(2*X(i,:));
            L(:,i+2*m)=(L(:,i+2*m)-mean(L(:,i+2*m)))/sqrt((L(:,i+2*m)-mean(L(:,i+2*m)))'*(L(:,i+2*m)-mean(L(:,i+2*m))));
            L(:,i+3*m)=sin(2*X(i,:));
            L(:,i+3*m)=(L(:,i+3*m)-mean(L(:,i+3*m)))/sqrt((L(:,i+3*m)-mean(L(:,i+3*m)))'*(L(:,i+3*m)-mean(L(:,i+3*m))));
            k=4*m;
            for j=i+1:m
                k=k+1;
                L(:,k)=cos(X(i,:)-X(j,:));
                L(:,k)=(L(:,k)-mean(L(:,k)))/sqrt((L(:,k)-mean(L(:,k)))'*(L(:,k)-mean(L(:,k))));
                L(:,k+nn)=sin(X(i,:)-X(j,:));
                L(:,k+nn)=(L(:,k+nn)-mean(L(:,k+nn)))/sqrt((L(:,k+nn)-mean(L(:,k+nn)))'*(L(:,k+nn)-mean(L(:,k+nn))));
                L(:,k+2*nn)=cos(X(i,:)+X(j,:));
                L(:,k+2*nn)=(L(:,k+2*nn)-mean(L(:,k+2*nn)))/sqrt((L(:,k+2*nn)-mean(L(:,k+2*nn)))'*(L(:,k+2*nn)-mean(L(:,k+2*nn))));
                L(:,k+3*nn)=sin(X(i,:)+X(j,:));
                L(:,k+3*nn)=(L(:,k+3*nn)-mean(L(:,k+3*nn)))/sqrt((L(:,k+3*nn)-mean(L(:,k+3*nn)))'*(L(:,k+3*nn)-mean(L(:,k+3*nn))));
            end
        end
end
[VN D]=eig(L'*L);
V=L*VN;
xnorm=repmat(sqrt(dot(V,V)),n,1);
%disp(sprintf('%e ',sqrt(dot(V,V))));
V=V./xnorm;
[s ind]=sort(abs(diag(D)),'descend');
r=sum(s>f*s(1));
ind=ind(1:r);
V=[V(:,ind);diag(D(ind,ind))'];
ifail=false;
