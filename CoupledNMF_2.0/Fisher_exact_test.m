function [P,F,E,ka]=Fisher_exact_test(r1,r2)
t=length(r1);
t1=sum(r1);
t2=sum(r2);

F(1,1)=sum(r1.*r2);
F(1,2)=sum(r1)-F(1,1);
F(2,1)=abs(sum((r1-1).*r2));
F(2,2)=abs(sum(r1-1))-F(2,1);

E(1,1)=(F(1,1)+F(2,1))*t1/t;
E(1,2)=(F(1,2)+F(2,2))*t1/t;
E(2,1)=(F(1,1)+F(2,1))*(t-t1)/t;
E(2,2)=(F(1,2)+F(2,2))*(t-t1)/t;
ka=sum(sum((F-E).*(F-E)./E));
P=1-chi2cdf(ka,1);
if F(1,1)<E(1,1)
    P=1;
end
%H=P<0.05;