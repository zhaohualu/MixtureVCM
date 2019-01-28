function[AICResC,BICResC,AICResKnot,BICResKnot]=Mix_VarCo_LLRR_PenBIC(Y,X,U,Z,V,Component,Knot,Int,lambda)
[S,~,~]=size(V);
[P,~,~]=size(X);
[Q,~,~]=size(Z);
[~,SK]=size(Knot);
[~,SL]=size(lambda);
N=sum(sum(Y~=inf));
AICResPenLoglikeli=-inf;
BICResPenLoglikeli=-inf;
for sk=1:SK
    temKnot=Knot(:,sk);
    KL=sum(temKnot)+4*P;
    for sl=1:SL
        temlambda=lambda(1,sl);
        [~,~,~,~,~,Loglikeli,temC]=Mixture_VaryingCo_LongReg_RanE_Sp(Y,X,U,Z,V,Component,temKnot,Int,temlambda,0);
        Degree=temC*(KL+Q+(S+1)*S/2+1+1)-1;
        temAICPenLoglikeli=Loglikeli-0.5*Degree*2;
        temBICPenLoglikeli=Loglikeli-0.5*Degree*log(N);
        if temAICPenLoglikeli>AICResPenLoglikeli
            AICReslambda=temlambda;
            AICResC=temC;
            AICResKnot=temKnot;
            AICResPenLoglikeli=temAICPenLoglikeli;
        end
        if temBICPenLoglikeli>BICResPenLoglikeli
            BICReslambda=temlambda
            BICResC=temC
            BICResKnot=temKnot
            BICResPenLoglikeli=temBICPenLoglikeli
        end
    end
end