function[Y,X,U,Z,V,B]=Simu_Data1()
n=150;
maxLong=2;
P=1;
Q=1;
S=1;
X=zeros(P,maxLong,n);
Y=zeros(n,maxLong);
U=zeros(n,maxLong);
Z=zeros(Q,maxLong,n);
V=zeros(S,maxLong,n);
B=zeros(n,S);
sigma=[0.6,0.6];
Bsigma=[0.6,0.6];
for i=1:n
    tem=rand(1,1);
    miss=rand(maxLong,1);
    miss(1,1)=1;
    miss(2,1)=1;
    if tem<1/2
        G=normrnd(0,Bsigma(1,1),1,1);
    else
        G=normrnd(0,Bsigma(1,2),1,1);
    end
    B(i,:)=G;
    for j=1:maxLong
        
        U(i,j)=4*(rand(1,1)-0.5);
        X(1,j,i)=rand(1,1);
        Z(1,j,i)=normrnd(0,1);
        Z(2,j,i)=normrnd(0,1);
        V(1,j,i)=normrnd(0,4);
        if tem<1/2
            Y(i,j)=X(:,j,i)*2*[-1+U(i,j)^2]+Z(:,j,i)'*[1;1]+V(1,j,i)*G+normrnd(0,sigma(1));
        else
            Y(i,j)=X(:,j,i)*2*[cos(0.5*pi*U(i,j))]+Z(:,j,i)'*[-1;-1]+V(1,j,i)*G+normrnd(0,sigma(2));
        end
        
        if miss(j,1)<0.5
        end
    end
end 