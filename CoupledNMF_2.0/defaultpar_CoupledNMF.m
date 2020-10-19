function [lambda1 lambda2]=defaultpar_CoupledNMF(PeakO,W1,H1,X,W2,H2,D,beta,arfa)
if nargin < 8
arfa=0.001;
beta=0.5;
end
if nargin < 9
arfa=0.001;
end
r1=mean(mean(PeakO*H1'))/mean(mean(D'*W2));
lambda2=2*arfa*beta*r1;
r2=mean(mean(X*H2'))/mean(mean(D*W1));
lambda1=beta/(1-beta+eps)*r1/r2;