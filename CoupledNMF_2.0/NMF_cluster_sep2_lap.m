function [W1,H1,W2,H2,Score,detr]=NMF_cluster_sep2_lap(PeakO,X,D,K,maxiter,lambda1,lambda2,W10,H10,W20,H20)
tolx=1e-4;
tolfun=1e-6;
sqrteps = sqrt(eps);
D1=D>0;
beta=2;
[n,m]=size(PeakO);
[n1,m1]=size(X);
%s=eye(K)-(1/(K-1))+diag(1*(1/(K-1)*ones(1,K)));
%s=2*eye(K)-1;
s=eye(K);
dnorm0=norm(PeakO-W10*H10,'fro')^2+lambda1*norm(X-W20*H20,'fro')^2;
for iter=1:maxiter
    S1=0.5*lambda2*D'*W20*s;
    numer = W10'*PeakO;
    H1 = max(0,H10.* (numer./((W10'*W10)*H10+eps(numer))));
    numer = PeakO*H1'+S1;
    mu11=diag(sum(0.5*W10.*max(0,numer)));
    mu12=diag(sum(0.5*W10.*(W10*(H1*H1'))));
    W1=max(0,W10.*((numer+2*(W10.^(beta-1))*mu12)./(W10*(H1*H1')+2*(W10.^(beta-1))*mu11+eps(numer))));
    S2=0.5*(lambda2/(lambda1+eps))*(D*W1*s);
    numer = W20'*X;
    H2 = max(0,H20.* (numer./((W20'*W20)*H20+eps(numer)))); 
    numer = X*H2'+S2;
    mu21=diag(sum(0.5*W20.*max(numer,0)));
    mu22=diag(sum(0.5*W20.*(W20*(H2*H2'))));
    W2=max(0,W20.*((numer+2*(W20.^(beta-1))*mu22)./(W20*(H2*H2')+2*(W20.^(beta-1))*mu21+eps(numer))));
    %sqrt([sum(W1.^2) sum(W2.^2)])
   % dnorm=norm(PeakO-W1*H1,'fro')^2+lambda1*norm(X-W2*H2,'fro')^2-lambda2*trace(W2'*D*W1)+mu*(norm(W1,'fro')^2+norm(W2,'fro')^2);
    dnorm=norm(PeakO-W1*H1,'fro')^2+lambda1*norm(X-W2*H2,'fro')^2;
    dw1 = max(max(abs(W1-W10) / (sqrteps+max(max(abs(W10))))));
    dh1 = max(max(abs(H1-H10) / (sqrteps+max(max(abs(H10))))));
    dw2 = max(max(abs(W2-W20) / (sqrteps+max(max(abs(W20))))));
    dh2 = max(max(abs(H2-H20) / (sqrteps+max(max(abs(H20))))));
    delta = max([dw1,dh1,dw2,dh2]);
    if iter>1
        if delta <= tolx
            disp(['delta ',num2str(delta),' is small'])
            break;
        elseif dnorm0-dnorm <= tolfun*max(1,dnorm0)
            disp(['dnorm0-dnorm ',num2str(dnorm0-dnorm),' is small'])
            break;
        elseif iter==maxiter
            break
        end
    end
    W10 = W1;
    H10 = H1;
    W20 = W2;
    H20 = H2;   
        if mod(iter,50)==0
                [d S20]=max(H20);
                [d S10]=max(H10);
                for j=1:K
                    FC1(:,j)=sum(PeakO(:,S10==j)'>0)'./(sum(PeakO(:,S10~=j)'>0)'/sum(S10~=j)*sum(S10==j)+1);
                    FC2(:,j)=sum(X(:,S20==j)'>0)'./(sum(X(:,S20~=j)'>0)'/sum(S20~=j)*sum(S20==j)+1);
                end
        WP1=quantilenorm(FC1,'MEDIAN',true);
        WP2=quantilenorm(FC2,'MEDIAN',true);
	S=(WP2'*D*WP1)';
	[assignment,cost] = munkres(-S);
 	W20=W20(:,assignment);
 	H20=H20(assignment,:);
        assignment
        -cost
        [W1'*W1 W2'*W2]
        (dnorm0-dnorm)/dnorm0
        end
        dnorm0 = dnorm;
end
    W1 = W10;
    H1 = H10;
    W2 = W20;
    H2 = H20;   
%%%%%
Score=norm(PeakO-W1*H1,'fro')^2+lambda1*norm(X-W2*H2,'fro')^2-lambda2*trace(W2'*D*W1);
if nargin <13
detr=0;
else
[d S20]=max(H2);
[d S10]=max(H1);
for j=1:K
    FC1(:,j)=sum(PeakO(:,S10==j)'>0)'./(sum(PeakO(:,S10~=j)'>0)'/sum(S10~=j)*sum(S10==j)+1);
    FC2(:,j)=sum(X(:,S20==j)'>0)'./(sum(X(:,S20~=j)'>0)'/sum(S20~=j)*sum(S20==j)+1);
end
WP1=quantilenorm(FC1,'MEDIAN',true);
WP2=quantilenorm(FC2,'MEDIAN',true);
T=WP2'*D*WP1;
T1=sum(sum(T))*diag(1./sum(T,2))*T*diag(1./sum(T));
detr=trace(T);
%Score=trace(T1);
end
%%%
[d S20]=max(H2);
[d S10]=max(H1);
for j=1:K
    FC1(:,j)=sum(PeakO(:,S10==j)'>0)'./(sum(PeakO(:,S10~=j)'>0)'/sum(S10~=j)*sum(S10==j)+1);
    FC2(:,j)=sum(X(:,S20==j)'>0)'./(sum(X(:,S20~=j)'>0)'/sum(S20~=j)*sum(S20==j)+1);
end
WP1=quantilenorm(FC1,'MEDIAN',true);
WP2=quantilenorm(FC2,'MEDIAN',true);
S=(WP2'*D*WP1)';
[assignment,cost] = munkres(-S);
W2=W2(:,assignment);
H2=H2(assignment,:);
assignment
-cost
