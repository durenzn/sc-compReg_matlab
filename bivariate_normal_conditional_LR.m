function [L_M0,L_M1,pValue,stat]=bivariate_normal_conditional_LR(X1,X2) 
sigma_min_pseudocount=0;
XX=[X1;X2];
beta=pinv([ones(size(XX,1),1) XX(:,2)])*XX(:,1);
Est_X=beta(1)+beta(2)*XX(:,2);
Condi_var=var(XX(:,1)-Est_X);
L_M0=sum(log(normpdf(XX(:,1),Est_X,sqrt(Condi_var)+sigma_min_pseudocount)));
%%%
beta1=pinv([ones(size(X1,1),1) X1(:,2)])*X1(:,1);
Est_X1=beta1(1)+beta1(2)*X1(:,2);
Condi_var1=var(X1(:,1)-Est_X1);
beta2=pinv([ones(size(X2,1),1) X2(:,2)])*X2(:,1);
Est_X2=beta2(1)+beta2(2)*X2(:,2);
Condi_var2=var(X2(:,1)-Est_X2);
L_M1_1=sum(log(normpdf(X1(:,1),Est_X1,sqrt(Condi_var1)+sigma_min_pseudocount)));
L_M1_2=sum(log(normpdf(X2(:,1),Est_X2,sqrt(Condi_var2)+sigma_min_pseudocount)));
L_M1=L_M1_1+L_M1_2;
[h,pValue,stat,cValue] = lratiotest(L_M1,L_M0,3);