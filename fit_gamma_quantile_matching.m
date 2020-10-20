function x=fit_gamma_quantile_matching(X,eita)
eita=eita';
cut=prctile(X,eita*100);
fun=@(arfa) gamcdf(cut,arfa(1),arfa(2))-eita;
x0=gamfit(min(X,cut(end)));
x = fsolve(fun,x0');