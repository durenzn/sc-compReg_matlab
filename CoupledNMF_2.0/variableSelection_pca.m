function [X1, Symbol1]=variableSelection_pca(X,Symbol,cut,numPC)
if nargin <4
numPC=10;
end
h0=sum(X);
h0=h0/sqrt(sum(h0.^2));
[w h]=nnmf(X,1,'h0',h0,'algorithm','mult');
X2=X-w*h;
[COEFF_Exp, SCORE_Exp, LATENT_Exp] = fastpca(X2',numPC);
s1=abs(COEFF_Exp(:,1:numPC));
s2=sum(s1,2);
[d f]=sort(s2,'descend');
X1=quantilenorm(X(f(1:cut),:),'MEDIAN',true);
Symbol1=Symbol(f(1:cut));