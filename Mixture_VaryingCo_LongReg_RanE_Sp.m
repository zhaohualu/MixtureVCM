function[ResA,ResB,Ressigma,ResBsigma,ResPi,ResLoglikeli,ResC]=Mixture_VaryingCo_LongReg_RanE_Sp(Y,X,U,Z,V,Component,Knot,Int,lambda,isplot)
[n,maxLong]=size(Y);
[P,~,~]=size(X);
[Q,~,~]=size(Z);
[S,~,~]=size(V);
KL=sum(Knot)+4*P;
if isempty(Int)
    Int=[min(min(U)),max(max(U))];
end
Long_Dic=zeros(maxLong,n);
for i=1:n
    for j=1:maxLong
        if Y(i,j)~=inf
            Long_Dic(i,j)=1;
        end
    end
end
lambda=lambda*(n^-0.25);

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

MaxLooptime=100;
MaxInLooptime=100;
Epis=0.08;
InEpis=0.08;
DecT=0.001;
Df=KL+Q+1;
peneps=(n^-0.6)/log(n);
rep=2;
ResLoglikeli=-inf;
penalty=2;

%--------------------------------------------------------------------------

while(rep>=0)
    Looptime=0;
    preLoglikeli=-inf;
    
    [A,B,sigma,C]=Mixture_VaryingCo_LongReg_Sp(Y,X,U,Z,Component,Knot,Int,lambda);
    Pi=1/C*ones(1,C);
    Bsigma=zeros(S,S,C);
    for c=1:C
        Bsigma(:,:,c)=0.5*eye(S);
    end
    sigma=0.5*ones(1,C);
    
    R=zeros(C,n);
    Ei=zeros(S,n,C);
    EResi=zeros(C,n);
    EResbi=zeros(S,S,n,C);
       
    while(Looptime<MaxLooptime)
        Looptime=Looptime+1;
    
        for i=1:n
            for c=1:C
                index=find(Long_Dic(i,:)==1);
                temY=Y(i,index)'-D(:,index,i)'*A(:,c)-Z(:,index,i)'*B(:,c);
                temW=(V(:,index,i)'*Bsigma(:,:,c)*V(:,index,i)+eye(sum(Long_Dic(i,:)))*sigma(1,c));
                R(c,i)=exp(-0.5*temY'*(temW\temY))*Pi(1,c)/sqrt(det(temW));
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
            if penalty==1
                if Pi(1,c)==0
                    Delec=[Delec,c];
                end
            elseif penalty==2
                if Pi(1,c)<DecT
                    Delec=[Delec,c];
                end
            end
        end
    
        tem=0;
        for c=1:C
            if sum(find(Delec==c))==0
                tem=tem+Pi(1,c);
            end
        end
        Pi(1,:)=Pi(1,:)/tem;
    
        preInLoglikeli=-inf;
        InLooptime=0;
        while(InLooptime<MaxInLooptime)
            InLooptime=InLooptime+1;
            for c=1:C
                if sum(find(Delec==c))==0
                    for i=1:n
                        index=find(Long_Dic(i,:)==1);
                        ni=sum(Long_Dic(i,:));
                        temW=(V(:,index,i)'*Bsigma(:,:,c)*V(:,index,i)+eye(ni)*sigma(1,c));
                        temInvW=temW^-1;
                        temY=(Y(i,index)'-D(:,index,i)'*A(:,c)-Z(:,index,i)'*B(:,c));
                        Ei(:,i,c)=Bsigma(:,:,c)*V(:,index,i)*temInvW*temY;
                        EResi(c,i)=trace(eye(ni)*sigma(1,c)-sigma(1,c)^2*temInvW)+sigma(1,c)^2*temY'*temInvW^2*temY;
                        EResbi(:,:,i,c)=Bsigma(:,:,c)-Bsigma(:,:,c)*V(:,index,i)*temInvW*V(:,index,i)'*Bsigma(:,:,c);
                    end
                end
            end
            
            for c=1:C
                if sum(find(Delec==c))==0
                    UP=zeros(KL+Q,1);
                    DOWM=zeros(KL+Q,KL+Q);
                    for i=1:n
                        index=find(Long_Dic(i,:)==1);
                        UP=UP+R(c,i)*[D(:,index,i);Z(:,index,i)]*(Y(i,index)'-V(:,index,i)'*Ei(:,i,c));
                        DOWM=DOWM+R(c,i)*[D(:,index,i);Z(:,index,i)]*[D(:,index,i);Z(:,index,i)]';
                    end
                    Coef=DOWM\UP;
                    A(:,c)=Coef(1:KL,1);
                    B(:,c)=Coef((KL+1):(Q+KL),1);
                end
            end

            for c=1:C
                if sum(find(Delec==c))==0
                    UP=0;
                    DOWM=0;
                    for i=1:n
                        UP=UP+R(c,i)*EResi(c,i);
                        DOWM=DOWM+R(c,i)*sum(Long_Dic(i,:));
                    end
                    sigma(1,c)=UP/DOWM;
                end
            end
        
            for c=1:C
                if sum(find(Delec==c))==0
                    UP=zeros(S,S);
                    for i=1:n
                        UP=UP+R(c,i)*(Ei(:,i,c)*Ei(:,i,c)'+EResbi(:,:,i,c));
                    end
                    Bsigma(:,:,c)=UP/sum(R(c,:));
                end
            end
            
            InLoglikeli=0;
            for c=1:C
                if sum(find(Delec==c))==0
                    for i=1:n
                        index=find(Long_Dic(i,:)==1);
                        ni=sum(Long_Dic(i,:));
                        temW=(V(:,index,i)'*Bsigma(:,:,c)*V(:,index,i)+eye(ni)*sigma(1,c));
                        temInvW=temW^-1;
                        temY=(Y(i,index)'-D(:,index,i)'*A(:,c)-Z(:,index,i)'*B(:,c));
                        
                        InLoglikeli=InLoglikeli+R(c,i)*(-0.5*temY'*temInvW*temY-0.5*log(det(temW))-ni/2*log(2*pi));
                    end
                end
            end
            
            if InLoglikeli>preInLoglikeli && InLoglikeli-preInLoglikeli<InEpis
                break;
            else
                preInLoglikeli=InLoglikeli;
            end
        end
                        
        A(:,Delec)=[];
        B(:,Delec)=[];
        sigma(Delec)=[];
        Pi(Delec)=[];
        Bsigma(:,:,Delec)=[];
        
        R(Delec,:)=[];
        Ei(:,:,Delec)=[];
        EResi(Delec,:)=[];
        EResbi(:,:,:,Delec)=[];
        
        C=C-length(Delec);
        Delec=[];
    
        Loglikeli=0;
        for i=1:n
            index=find(Long_Dic(i,:)==1);
            ni=sum(Long_Dic(i,:));
            tem=0;
            for c=1:C
                temW=(V(:,index,i)'*Bsigma(:,:,c)*V(:,index,i)+eye(ni)*sigma(1,c));
                temInvW=temW^-1;
                temY=(Y(i,index)'-D(:,index,i)'*A(:,c)-Z(:,index,i)'*B(:,c));
                tem=tem+Pi(1,c)*exp(-0.5*temY'*temInvW*temY)/(sqrt(det(temW))*(2*pi)^(0.5*ni));
            end
            Loglikeli=Loglikeli+log(tem);
        end
    
        if (Loglikeli>=preLoglikeli && Loglikeli-preLoglikeli<Epis)
            break;
        else
            preLoglikeli=Loglikeli;
        end
    end
    
    if Loglikeli~=-inf
        rep=rep-1;
    end
        
    if Loglikeli>ResLoglikeli
        ResA=A;
        ResB=B;
        ResC=C;
        Ressigma=sigma;
        ResBsigma=Bsigma;
        ResPi=Pi;
        ResLoglikeli=Loglikeli;
    end
end

if (isplot)
    test=min(min(U)):0.01:max(max(U));
    testy=zeros(1,length(test));
    for c=1:C
        for p=1:P
            for i=1:length(test)
                testy(1,i)=Trun_PowB(Int,Knot(p,1),3,test(1,i))*ResA((sum(Knot(1:(p-1),1))+(p-1)*4+1):(sum(Knot(1:p,1))+p*4),c);
            end
            plot(test,testy);       
            title(strcat('The',{32},num2str(p),'-th varying coefficient of the',{32},num2str(c),'-th component'))
            pause();
        end
    end
end
    
        
    
    


