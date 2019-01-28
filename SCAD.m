function[res]=SCAD(x,lambda,deri)

if x<0
    disp('x should be non-nagetive');
end

a=3.7;

if deri==1
    if x<=lambda
        res=1;
    else
        res=max(0,a*lambda-x)/((a-1)*lambda);
    end
end

if deri==0
    if x<=lambda
        res=x;
    elseif x>lambda && x<(a*lambda)
        res=-0.5*lambda/(a-1)+a*lambda/((a-1)*lambda)*x-0.5/((a-1)*lambda)*x^2;
    elseif x>=a*lambda
        res=0.5*(a+1)*lambda;
    end
end