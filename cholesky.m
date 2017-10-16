function [A]=cholesky(A)

% n=length(A);
% temp1=0;
% temp2=0;
% for i=1:n
%     for k=1:i
%         temp1=temp1+A(i,k)^2;
%     end
%     A(i,i)=sqrt(A(i,i)-temp1);
%     for j=i+1:n
%         for k=1:i
%             temp2=temp2+A(i,k)*A(j,k);
%         end
%         A(j,i)=(A(i,j)-temp2)/A(i,i);
%     end
% end
% end

  [n nn]=size(A);
for k=1:n
    A(k,k)=sqrt(A(k,k));
    A(k+1:n,k)=A(k+1:n,k)/A(k,k);
    for j=k+1:n
        A(j:n,j)=A(j:n,j)-A(j,k)*A(j:n,k);
    end
end
end
