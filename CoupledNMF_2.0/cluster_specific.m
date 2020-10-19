function [p2 FC2 Score2 Specific]=cluster_specific(X,S2,Symbol,method,datatype,Outdir)
%%%%%%%%Cluster specific gene
K=length(unique(S2));
if method ==1
    for i=1:size(X,1)
    for j=1:K
    [p2(i,j),F,E,ka]=Fisher_exact_test(full(X(i,:)>0),S2==j);
    FC2(i,j)=F(1,1)/(1+E(1,1));
    end
    end
else
    for j=1:K
        FC2(:,j)=2.^(mean(X(:,S2==j)')'-mean(X(:,S2~=j)')');
    end
    for j=1:K
    [h p2(j,:)]=ttest2(X(:,S2==j)',X(:,S2~=j)','tail','right');
    end
    p2=p2';
end
p2(isnan(p2))=1;
p2(p2<10^(-10))=10^(-10);
Score2=-log10(p2).*FC2.*(FC2>1);
for j=1:K
    [d idx(:,j)]=sort(Score2(:,j),'descend');
    Specific{1,j}=Symbol(idx(1:min(1000,sum(Score2(:,j)>2)),j));
end
if nargin > 4
for j=1:K
    filename=[Outdir,'/cluster_',int2str(j),'_specific_',datatype,'.txt'];
    fid=fopen(filename,'wt');
    for i=1:size(Specific{1,j},1)
        fprintf(fid, '%s\t',Specific{1,j}{i,1});
	fprintf(fid, '%g\t',p2(idx(i,j),j));
	fprintf(fid, '%f\t',FC2(idx(i,j),j));
        fprintf(fid, '%f\n',Score2(idx(i,j),j));
    end
     fclose(fid);
end
end