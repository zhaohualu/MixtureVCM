function[Res]=Trun_PowB(Int,NumK,Degree,x)
Interval=Int(1,1):(Int(1,2)-Int(1,1))/(NumK+1):Int(1,2);
Res=zeros(1,Degree+NumK+1);
for i=1:(Degree+1)
    Res(1,i)=x^(i-1);
end
for i=1:NumK
    Res(1,i+Degree+1)=max(0,(x-Interval(1+i))^Degree);
end

    

