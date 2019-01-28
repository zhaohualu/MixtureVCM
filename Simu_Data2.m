function[Y,X,U,Z,V,B]=Simu_Data2()
n=250;
maxLong=10;
P=2;
Q=1;
S=1;
X=zeros(P,maxLong,n);
Y=zeros(n,maxLong);
U=zeros(n,maxLong);
Z=zeros(Q,maxLong,n);
V=zeros(S,maxLong,n);
B=zeros(n,S);
sigma=[0.3,0.7,0.5];
Bsigma=[0.3,0.7,0.5];
for i=1:n
    tem=rand(1,1);
    miss=rand(maxLong,1);
    miss(1,1)=1;
    miss(2,1)=1;
    if tem<1/3
        G=normrnd(0,Bsigma(1,1),1,1);
    elseif tem>=1/3&&tem<2/3
        G=normrnd(0,Bsigma(1,2),1,1);
    else
        G=normrnd(0,Bsigma(1,3),1,1);
    end
    B(i,:)=G;
    for j=1:maxLong
        
        %random miss
        
        U(i,j)=4*rand(1,1);
        X(2,j,i)=rand(1,1)+0.25*U(i,j);
        X(1,j,i)=1;
        Z(1,j,i)=normrnd(0,1);
        V(1,j,i)=normrnd(0,4);
        if tem<1/3
            Y(i,j)=X(:,j,i)'*[2*sin(pi*(U(i,j)));-1*U(i,j)]+Z(:,j,i)'*-1+V(1,j,i)*G+normrnd(0,sigma(1));
        elseif tem>=1/3 && tem<=2/3
            Y(i,j)=X(:,j,i)'*[(-1*exp(U(i,j)));(1*exp(U(i,j)))]+Z(:,j,i)'*0+V(1,j,i)*G+normrnd(0,sigma(2));
        else
            Y(i,j)=X(:,j,i)'*[2*cos(pi*U(i,j));1*U(i,j)^2]+Z(:,j,i)'*1+V(1,j,i)*G+normrnd(0,sigma(3));
        end
        
        %random miss
        if miss(j,1)<0.5
            Y(i,j)=inf;
        end
    end
end 