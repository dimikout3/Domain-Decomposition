clear;
clc;
%c=1; % Source coefficient
m=50; % Number of Squares per dim
h=2/m; % Mesh size
ne=m^2+m*m/2; % Number of elements
non=(m/2)*(m/2+1)+(m/2+1)*(m+1); % Number of nodes
coo=zeros(non,2); % Initializing the matrix of coordinates
nn=zeros(ne,3); % Initializing the matrix of elements
t1=(-1:h:1)';
t2=(0:h:1)';
ax=1;
ay=1;
c=0;

bounds=[];
bounds=[bounds 1:m+1  (m+2):(m+1):(m+1)*(m/2+1)+1   ]; %katw arstera
bounds=[bounds ((m+1)*(m/2+1)+m/2+2):(m/2+1):(non-m/2)]; %panw aristera
bounds=[bounds 2*(m+1):(m+1):(m+1)*(m/2+1)]; %deksia katw
bounds=[bounds ((m+2)*(m/2+1)):(m/2+1):non]; %deksia panw
bounds=[bounds (non-1):(-1):(non-m/2+1) ((m/2)*(m+1)+2):(m/2)*(m+1)+m/2+1]; %panw

for i=1:m/2+1
    coo(((i-1)*(m+1)+1):i*(m+1),1)=t1;
    coo(((i-1)*(m+1)+1):i*(m+1),2)=(i-1)*h-1;
end

for i=1:m/2
    coo(((i-1)*(m/2+1)+1+(m/2+1)*(m+1)):(m/2+1)*(m+2)+(i-1)*(m/2+1),1)=t2;
    coo(((i-1)*(m/2+1)+1+(m/2+1)*(m+1)):(m/2+1)*(m+2)+(i-1)*(m/2+1),2)=i*h;
end



k=1;
for i=1:m/2
    for j=1:m
        nn(k,:)=[j j+1 j+m+1]+(i-1)*(m+1);
        k=k+1;
        nn(k,:)=[j+m+2 j+m+1 j+1]+(i-1)*(m+1);
        k=k+1;
    end
end

for i=1:m/2
    for j=1:m/2
        nn(k,:)=[j j+1 j+m/2+1]+(i-1)*(m/2+1)+(m/2)*(m+2);
        k=k+1;
        nn(k,:)=[j+m/2+2 j+m/2+1 j+1]+(i-1)*(m/2+1)+(m/2)*(m+2);
        k=k+1;
    end
end


K=spalloc(non,non,6*non);
b=zeros(non,1);


for i=1:ne
    x1=coo(nn(i,1),1);
    x2=coo(nn(i,2),1);
    x3=coo(nn(i,3),1);
    y1=coo(nn(i,1),2);
    y2=coo(nn(i,2),2);
    y3=coo(nn(i,3),2);
    
    A=((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))/2;
    M(1,1)=(y2-y3)^2/(4*A)+(x3-x2)^2/(4*A);
    M(1,2)=(y2-y3)*(y3-y1)/(4*A)+(x3-x2)*(x1-x3)/(4*A);
    M(2,1)=M(1,2);
    M(1,3)=(y2-y3)*(y1-y2)/(4*A)+(x3-x2)*(x2-x1)/(4*A);
    M(3,1)=M(1,3);
    M(2,3)=(y3-y1)*(y1-y2)/(4*A)+(x1-x3)*(x2-x1)/(4*A);
    M(3,2)=M(2,3);
    M(2,2)=(y3-y1)^2/(4*A)+(x3-x1)^2/(4*A);
    M(3,3)=(y1-y2)^2/(4*A)+(x2-x1)^2/(4*A);
    T(1,1)=c*A/6;
    T(1,2)=c*A/12;
    T(2,1)=c*A/12;
    T(1,3)=c*A/12;
    T(3,1)=c*A/12;
    T(2,3)=c*A/12;
    T(3,2)=c*A/12;
    T(2,2)=c*A/6;
    T(3,3)=c*A/6;
    F(1,1)=A/3;
    F(2,1)=A/3;
    F(3,1)=A/3;
    

    for j=1:3
        for k=1:3
            K(nn(i,j),nn(i,k))=K(nn(i,j),nn(i,k))+M(j,k)+T(j,k);
        end
        b(nn(i,j))=b(nn(i,j))+F(j);
    end
end
K(bounds,:)=[];
K(:,bounds)=[];
b(bounds)=[];



%sunartisi

% P=[1:(m-1)*(m/2-1) m*(m/2-1)+1:m*(m/2-1)+(m/2-1)^2 (m-1)*(m/2-1)+1:m*(m/2-1)];
% A=K(P,P);
A1=K(1:(m-1)*(m/2-1),1:(m-1)*(m/2-1));
A2=K(m*(m/2-1)+1:m*(m/2-1)+(m/2-1)^2,m*(m/2-1)+1:m*(m/2-1)+(m/2-1)^2);
D=K((m-1)*(m/2-1)+1:m*(m/2-1),(m-1)*(m/2-1)+1:m*(m/2-1));
B1=K(1:(m-1)*(m/2-1),(m-1)*(m/2-1)+1:m*(m/2-1));
B2=K(m*(m/2-1)+1:m*(m/2-1)+(m/2-1)^2,(m-1)*(m/2-1)+1:m*(m/2-1));
C1=K((m-1)*(m/2-1)+1:m*(m/2-1),1:(m-1)*(m/2-1));
C2=K((m-1)*(m/2-1)+1:m*(m/2-1),m*(m/2-1)+1:m*(m/2-1)+(m/2-1)^2);


S=D - (C1*(A1\B1)) - (C2*(A2\B2));

b1=b(1:(m-1)*(m/2-1));
b2=b(m*(m/2-1)+1:m*(m/2-1)+(m/2-1)^2);
b3=b((m-1)*(m/2-1)+1:m*(m/2-1)) - (C1*(A1\b1)) - (C2*(A2\b2));

% sol3=S\b3;
% sol2=A2\(b2-B2*sol3);
% sol1=A1\(b1-B1*sol3);
% 
% M1=frederickson(A1);
% M2=frederickson( A2);
% M3=frederickson( S);

% M1=ilu0(A1);
% M2=ilu0( A2);
% M3=ilu0( S);

M1=SPAI(A1,2);
M2=SPAI( A2,2);
M3=SPAI( S,2);
% 
% M1=kolotolina(A1,2);
% M2=kolotolina( A2,2);
% M3=kolotolina( S,2);

% M1=ichol(A1);
% M2=ichol( A2);
% M3=ichol( S);

% tic
% sol3=GMRES(S,b3,1e-10,500,20,M3);
% sol2=GMRES(A2,b2-B2*sol3,1e-10,500,20,M2);
% sol1=GMRES(A1,b1-B1*sol3,1e-10,500,20,M1);
% toc

tic
sol3=bicgstabPreco(S,b3,500,M3);
sol2=bicgstabPreco(A2,b2-B2*sol3,500,M2);
sol1=bicgstabPreco(A1,b1-B1*sol3,500,M1);
toc




x=[sol1; sol3; sol2];

%reshape

xa=x(1:(m-1)*(m/2-1));
ya=reshape(xa,m-1,m/2-1);

xb=x((m-1)*(m/2-1)+1:size(x,1));
yb=reshape(xb,m/2-1,[]);

yc=[zeros(m/2,m/2); yb];
yall=[ya yc];
y=reshape(yall,m-1,m-1);
mesh(y);





