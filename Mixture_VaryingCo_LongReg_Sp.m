function[A,B,sigma,C,Loglikeli]=Mixture_VaryingCo_LongReg_Sp(Y,X,U,Z,Component,Knot,Int,lambda)
[n,maxLong]=size(Y);
[P,~,~]=size(X);
[Q,~,~]=size(Z);
C=Component;
KL=sum(Knot)+4*P;
if isempty(Int)
    Int=[min(min(U)),max(max(U))];
end
peneps=(n^-0.6)/log(n);
rep=3;
penalty=2;
ResLoglikeli=-inf;
Long_Dic=zeros(n,maxLong);
for i=1:n
    for j=1:maxLong
        if Y(i,j)~=inf
            Long_Dic(i,j)=1;
        else
            Y(i,j)=0;
            X(:,j,i)=zeros(P,1);
            Z(:,j,i)=zeros(Q,1);
            U(i,j)=0;
        end
    end
end
D=zeros(KL,maxLong,n);
for i=1:n
    for j=1:maxLong
        if Long_Dic(i,j)==1
            BU=zeros(P,KL);
            for p=1:P
                BU(p,(4*(p-1)+sum(Knot(1:(p-1),1))+1):(4*(p-1)+sum(Knot(1:(p-1),1))+Knot(p,1)+4))=Trun_PowB(Int,Knot(p,1),3,U(i,j));
            end
            D(:,j,i)=(X(:,j,i)'*BU)';
        end
    end
end

Design=zeros(KL+Q,maxLong,n);
for i=1:n
    for j=1:maxLong
        if Long_Dic(i,j)==1
            Design(:,j,i)=[D(:,j,i);Z(:,j,i)];
        end
    end
end

while (rep>0)
    rep=rep-1;
    C=Component;
    A=zeros(KL,C);
    B=zeros(Q,C);
    sigma=zeros(1,C);
    Pi=1/C*ones(1,C);
    R=zeros(C,n);

    Looptime=0;
    MaxLooptime=15;
    Epis=1;
    preLoglikeli=-inf;
    DecT=0.01;
    Df=P+Q+1;

    temY=zeros(n,1);
    temDesign=zeros(n,KL+Q);
    temX=zeros(n,P);
    temZ=zeros(n,Q);
    temU=zeros(n,1);
    for i=1:n
        Randj=randperm(maxLong);
        for j=1:maxLong
            if Long_Dic(i,Randj(j))==1
                temY(i,1)=Y(i,Randj(j));
                temX(i,:)=X(:,Randj(j),i)';
                temZ(i,:)=Z(:,Randj(j),i)';
                temU(i,1)=U(i,Randj(j));
                temDesign(i,:)=[D(:,Randj(j),i);Z(:,Randj(j),i)]';
                break;
            end
        end
    end

    list=kmeans([temY,temX,temU,temZ],C);
    for c=1:C
        index=find(list==c);
        Coef=(temDesign(index,:)'*temDesign(index,:)+0.2*diag(ones(KL+Q,1)))\(temDesign(index,:)'*temY(index,1));
        A(:,c)=Coef(1:KL,1);
        B(:,c)=Coef((KL+1):(KL+Q),1);
        sigma(1,c)=(temY(index,1)-temDesign(index,:)*[A(:,c);B(:,c)])'*(temY(index,1)-temDesign(index,:)*[A(:,c);B(:,c)])/length(index);
    end

    while(Looptime<MaxLooptime)
        Looptime=Looptime+1;
    
        for i=1:n
            for c=1:C
                R(c,i)=Pi(1,c);
                for j=1:maxLong
                    if Long_Dic(i,j)==1
                        R(c,i)=R(c,i)*1/sqrt(2*pi*sigma(1,c))*exp(-0.5/sigma(1,c)*(Y(i,j)-D(:,j,i)'*A(:,c)-Z(:,j,i)'*B(:,c))^2);
                    end
                end
            end
            tem=sum(R(:,i));
            for c=1:C
                R(c,i)=R(c,i)/tem;
            end
        end

        if penalty==2
            temPi=Pi;
            DOWN=0;
            for c=1:C
                DOWN=DOWN-n*lambda*Df*SCAD(temPi(1,c),lambda,1)*temPi(1,c)/(peneps+SCAD(temPi(1,c),lambda,0));
            end
        end
        for c=1:C
            if penalty==0
                Pi(1,c)=mean(R(c,:));
            elseif penalty==1
                Pi(1,c)=max(0,(mean(R(c,:))-lambda*Df)/(1-C*lambda*Df));
            elseif penalty==2
                Pi(1,c)=max(0,(sum(R(c,:))/(n+DOWN+n*lambda*Df*SCAD(temPi(1,c),lambda,1)/(peneps+SCAD(temPi(1,c),lambda,0)))));
            end
        end
    
        Delec=[];
        for c=1:C
            if Pi(1,c)<DecT && penalty~=0
                Delec=[Delec,c];
            end
        end
    
        tem=0;
        for c=1:C
            if sum(find(Delec==c))==0
                tem=tem+Pi(1,c);
            end
        end
        Pi(1,:)=Pi(1,:)/tem;
    
        for c=1:C
            if sum(find(Delec==c))==0
                UP=zeros(KL+Q,1);
                DOWN=zeros(KL+Q,KL+Q);
                for i=1:n
                    temD=[];
                    temZ=[];
                    temY=[];
                    for j=1:maxLong
                        if Long_Dic(i,j)==1
                            temD=[temD;D(:,j,i)'];
                            temZ=[temZ;Z(:,j,i)'];
                            temY=[temY;Y(i,j)];
                        end
                    end
                    UP=UP+[temD,temZ]'*(diag(R(c,i)*ones(sum(Long_Dic(i,:)),1)))*temY;
                    DOWN=DOWN+[temD,temZ]'*(diag(R(c,i)*ones(sum(Long_Dic(i,:)),1)))*[temD,temZ];
                end

                Coef=DOWN\UP;
                A(:,c)=Coef(1:KL,1);
                B(:,c)=Coef((KL+1):(KL+Q),1);
            end
        end
    
        for c=1:C
            if sum(find(Delec==c))==0
                UP=zeros(1,1);
                DOWN=zeros(1,1);
                for i=1:n
                    temD=[];
                    temZ=[];
                    temY=[];
                    for j=1:maxLong
                        if Long_Dic(i,j)==1
                            temD=[temD;D(:,j,i)'];
                            temZ=[temZ;Z(:,j,i)'];
                            temY=[temY;Y(i,j)];
                        end
                    end
                    UP=UP+sum((temY-[temD,temZ]*[A(:,c);B(:,c)]).^2)*R(c,i);
                    DOWN=DOWN+R(c,i)*sum(Long_Dic(i,:));
                end
                sigma(1,c)=UP/DOWN;
            end
        end
    
        A(:,Delec)=[];
        B(:,Delec)=[];
        R(Delec,:)=[];
        sigma(Delec)=[];
        Pi(Delec)=[];
        C=C-length(Delec);
        Delec=[];
    
        Loglikeli=0;
        for i=1:n
            tem=0;
            for c=1:C
                ttem=1;
                for j=1:maxLong
                    if Long_Dic(i,j)==1
                        ttem=ttem*1/sqrt(2*pi*sigma(1,c))*exp(-0.5/sigma(1,c)*(Y(i,j)-[D(:,j,i);Z(:,j,i)]'*[A(:,c);B(:,c)])^2);
                    end
                end
                tem=tem+Pi(1,c)*ttem;
            end
            Loglikeli=Loglikeli+log(tem);
        end
        Loglikeli;
        if Loglikeli>=preLoglikeli && Loglikeli-preLoglikeli<Epis
            break;
        else
            preLoglikeli=Loglikeli;
        end
    end
    if Loglikeli>ResLoglikeli
        ResLoglikeli=Loglikeli;
        ResA=A;
        ResB=B;
        Ressigma=sigma;
        ResC=C;
    end
end

    
        
    
    


