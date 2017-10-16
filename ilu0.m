function [ L,U ] = ilu0( A )
    n=length(A);
    L=speye(n,n);
    for i=2:n
        l=find(A(i,1:i-1)); 
        for k=1:length(l) 
            A(i,l(k))=A(i,l(k))/A(l(k),l(k));
            u=find(A(i,k+1:n));
            for j=1:length(u)
                A(i,u(j))=A(i,u(j))-A(i,l(k))*A(l(k),u(j));
            end
        end
        L(i,1:i-1)=A(i,1:i-1);
        U(i,i:n)=A(i,i:n);
    end
end

