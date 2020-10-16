function [diffNet,Hub_TF]=compReg(TF_binding,Match,E1,E1_idx,E2,E2_idx,O1_mean,O2_mean,Symbol,TFName,Element_name,peak_gene_prior)
%%TF_binding
a=full(maxk(TF_binding,5000,2));
TF_binding(TF_binding-a(:,end)<0)=0;
%%%%beta
fileID = fopen(peak_gene_prior);
C = textscan(fileID,'%s %s %f32 %f32');
fclose(fileID);
[d f]=ismember(C{1,1},Element_name);
[d1 f1]=ismember(C{1,2},Symbol);
[f2 ia ic]=unique([f(d.*d1==1) f1(d.*d1==1)],'rows');
c3=accumarray(ic,C{1,3}(d.*d1==1),[],@min);
c4=accumarray(ic,C{1,4}(d.*d1==1),[],@min);
c4(c4<0.2)=0;
d0=500000;
c=double(exp(-1*c3/d0).*c4);
beta=sparse(f2(:,2),f2(:,1),c,length(Symbol),length(Element_name));
%%pairwise
for ii=1:size(Match,1)
    i1=Match(ii,1);
    i2=Match(ii,2);
    BO1=(TF_binding.*O1_mean(:,i1)')*beta';
    BO2=(TF_binding.*O2_mean(:,i2)')*beta';
    TG1=E1(:,E1_idx==i1);
    TG2=E2(:,E2_idx==i2);
    TG2=E2(:,E2_idx==i2)*mean(mean(TG1))/mean(mean(TG2));
    [d f]=ismember(TFName,Symbol);
    TF=[TG1(f,:) TG2(f,:)];
    n1=size(TG1,2);
    n2=size(TG2,2)
    [h p]=ttest2(TG1',TG2');
    [h crit_p adj_p]=fdr_bh(p);
    diffGene=find(adj_p<0.1);
    [r1 p1]=corr(TF(:,1:n1)',TG1');
    [r2 p2]=corr(TF(:,1+n1:n1+n2)',TG2');
    p_combine=min(p1,p2);
    clear net_idx
    [net_idx(:,1) net_idx(:,2)]=find(p_combine<0.05);
    %%%
    LR_summary=[];
    %j=1;% gene
    for j=1:length(diffGene)
    clear stat
    OTF1=BO1(:,diffGene(j)).*TF(:,1:n1);
    OTF2=BO2(:,diffGene(j)).*TF(:,1+n1:n1+n2);
    id1=net_idx(net_idx(:,2)==diffGene(j),1);
    %id=intersect(id1,find(BO1(:,diffGene(j))+BO2(:,diffGene(j))>0));
    %id=find(BO1(:,diffGene(j))+BO2(:,diffGene(j))>0);
    id=find(sum(abs(OTF1)')+sum(abs(OTF2)')>0);
    id=intersect(id,id1);
    for i=1:length(id)  
        X1=full([OTF1(id(i),:)' TG1(diffGene(j),:)']);
        X2=full([OTF2(id(i),:)' TG2(diffGene(j),:)']);
        [L_M0,L_M1,pValue(i,1),stat(i,1)]=bivariate_normal_conditional_LR(X1,X2);
        LR_summary=[LR_summary;[id(i) diffGene(j) stat(i,1) pValue(i,1)]];
    end
    j
    end
    LR_summary(isnan(LR_summary(:,3)),:)=[];
    LR_summary(isinf(LR_summary(:,3)),:)=[];
    LR_summary(LR_summary(:,3)<10^(-16),3)=10^(-16);
    phat=fit_gamma_quantile_matching(LR_summary(:,3),[0.01:0.01:0.2]);
    %load('phat_gamma.mat')
    p=1-gamcdf(LR_summary(:,3),phat(1,1),phat(2,1));
    [h crit_p adj_p]=fdr_bh(p);
    LR_summary=[LR_summary p adj_p];
    %idx=find(adj_p<=0.1);
    idx=find(adj_p<0.1);
    [d f]=sort(LR_summary(idx,3),'descend');
    idx1=idx(f);
    diffNet{1,ii}=[TFName(LR_summary(idx1,1)) Symbol(LR_summary(idx1,2)) num2cell(LR_summary(idx1,3:6))];
    u=unique(diffNet{1,ii}(:,1));
    [d f]=ismember(diffNet{1,ii}(:,1),u);
    clear u_count
    for i=1:length(u)
        u_count(i,1)=sum(f==i);
    end
    [d1 f1]=sort(u_count,'descend');
    Hub_TF{1,ii}=[u(f1) num2cell(d1)];
end