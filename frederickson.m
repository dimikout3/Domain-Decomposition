function[M]=frederickson(A)
n=length(A);
M=A;
e1=speye(n,n);
for j=1:n
    Q=A(:,M(:,j)~=0);
    temp=sum(abs(Q),2);
    Q=Q(temp~=0,:);
    M(M(:,j)~=0,j)=(pinv(full(Q))*e1(temp~=0,j))';
end
  