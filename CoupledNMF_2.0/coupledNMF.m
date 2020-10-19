function [W1,H1,W2,H2,lambda1,lambda2]=coupledNMF(PeakO,X,D,K,arfa,betamax_scale,betamin)
if nargin < 5
	arfa=0.5;
end

if nargin < 6
	betamax_scale=4;
end

if nargin > 6
	beta=betamin*10.^[betamax_scale:-1:0];
else 
	beta=10.^[betamax_scale:-1:-3];
end
[w2 h2]=nnmf(X,K);
[w1 h1]=nnmf(PeakO,K);
W10=rand(size(PeakO,1),K);
H10=rand(K,size(PeakO,2));
%W20=rand(size(X,1),K);
W20=D*W10;
H20=rand(K,size(X,2));
H1=diag(sum(W10.^2))*H10;
H2=diag(sum(W20.^2))*H20;
W1=W10*diag(1./sqrt(sum(W10.^2)));
W2=W20*diag(1./sqrt(sum(W20.^2)));
for j=1:length(beta)
[lambda1 lambda2]=defaultpar_CoupledNMF(PeakO,w1,h1,X,w2,h2,D,arfa,beta(j))
[W1,H1,W2,H2,Score,detr]=NMF_cluster_sep2_lap(PeakO,X,D,K,300,lambda1,lambda2,W1,H1,W2,H2);
WW1{1,j}=sparse(W1);
WW2{1,j}=sparse(W2);
HH1{1,j}=sparse(H1);
HH2{1,j}=sparse(H2);
H1=diag(1./sqrt(sum(H1.^2,2)))*H1;
H2=diag(1./sqrt(sum(H2.^2,2)))*H2;
W1=W1*diag(1./sqrt(sum(W1.^2)));
W2=W2*diag(1./sqrt(sum(W2.^2)));
end
