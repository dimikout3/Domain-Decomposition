function[M]=kolotolina(A,fill)
n=length(A);
M=A^fill;
e1=speye(n,n);
for j=1:n
    Q=A(M(:,j)~=0,M(:,j)~=0);

    M(M(:,j)~=0,j)=(pinv(full(Q))*e1(M(:,j)~=0,j))';
end
  