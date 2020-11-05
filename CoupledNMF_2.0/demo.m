%%%%load data
load('../example/sample1.mat')
load('../example/coupling_matrix_sample1.mat')
Symbol=Symbol1;
X=E1;
a=sum(O1'>0);
PeakName=PeakName1(a>100);
PeakO=O1(a>100,:);
%%%parametersd
K=7;
arfa=0.5
betamax_scale=4
atac_binarize=0
Outdir='./.';
%%%%%run coupldNMF
PeakO=full(PeakO);
if atac_binarize >0
PeakO = 1*(PeakO>0);
end
[W1,H1,W2,H2,lambda1,lambda2]=coupledNMF(PeakO,X,D,K,arfa,betamax_scale);
[d S1]=max(H1);
[d S2]=max(H2);
save([Outdir,'/Paramaters_',num2str(lambda1),'_',num2str(lambda2),'.mat'],'W1','W2','H1','H2','lambda1','lambda2','S1','S2')
dlmwrite([Outdir,'/scATAC_label.txt'],S1,'\t')
dlmwrite([Outdir,'/scRNA_label.txt'],S2,'\t')
%%%cluster specific 
H1_norm=H1./sum(H1);
H2_norm=H2./sum(H2);
H1_max=max(H1_norm);
H2_max=max(H2_norm);
S1_assign=S1;
S2_assign=S2;
for i=1:K
	S1_assign(S1==i)=i*(H1_max(S1==i)>median(H1_max(S1==i)));
	S2_assign(S2==i)=i*(H2_max(S2==i)>median(H2_max(S2==i)));
end
dlmwrite([Outdir,'/scATAC_label_assign.txt'],S1_assign,'\t')
dlmwrite([Outdir,'/scRNA_label_assign.txt'],S2_assign,'\t')
[p1 FC1 Score1 Specific_peak]=cluster_specific(PeakO(:,S1_assign>0),S1(S1_assign>0),PeakName,2,'scATAC',Outdir);
[p2 FC2 Score2 Specific_gene]=cluster_specific(X(:,S2_assign>0),S2(S2_assign>0),Symbol,2,'scRNA',Outdir);
%%%%PCA based embedding
rmPC=1;
mdl = fitlsa(PeakO',20);
tSNE1=tsne(mdl.DocumentScores,'Distance','spearman','Standardize',true,'LearnRate',250,'Perplexity',50);
[COEFF_Exp, SCORE_Exp, LATENT_Exp] = fastpca(X',20+rmPC);
dep2=sum(X>0);
r2=corr(SCORE_Exp,dep2');
[d f]=sort(abs(r2),'descend');
tSNE2=tsne(SCORE_Exp(:,setdiff([1:20+rmPC],f(1:rmPC))),'Distance','spearman','Standardize',true,'LearnRate',250,'Perplexity',50);
figure
gscatter(tSNE1(:,1),tSNE1(:,2),S1,[],[],10)
set(gca,'FontSize',12);
set(gcf, 'Position', [0, 0, 600 500])
saveas(gcf,[Outdir,'/tSNE_opn'],'pdf')
figure
gscatter(tSNE2(:,1),tSNE2(:,2),S2,[],[],10)
set(gca,'FontSize',12);
set(gcf, 'Position', [0, 0, 600 500])
saveas(gcf,[Outdir,'/tSNE_exp'],'pdf')
%%%%CoupledNMF based embedding
tSNE1_H=tsne(H1(:,S1_assign>0)','Distance','cosine','Standardize',true,'LearnRate',250,'Perplexity',50);
tSNE2_H=tsne(H2(:,S2_assign>0)','Distance','cosine','Standardize',true,'LearnRate',250,'Perplexity',50);
figure
gscatter(tSNE1_H(:,1),tSNE1_H(:,2),S1(S1_assign>0),[],[],10)
set(gca,'FontSize',12);
set(gcf, 'Position', [0, 0, 600 500])
saveas(gcf,[Outdir,'/tSNE_opn_H'],'pdf')
figure
gscatter(tSNE2_H(:,1),tSNE2_H(:,2),S2(S2_assign>0),[],[],10)
set(gca,'FontSize',12);
set(gcf, 'Position', [0, 0, 600 500])
saveas(gcf,[Outdir,'/tSNE_exp_H'],'pdf')
%%%%%joint embedding
HH=[H1(:,S1_assign>0) H2(:,S2_assign>0)];
HH=HH./sqrt(sum(HH.^2));
tSNE_joint=tsne(HH','Distance','cosine','Standardize',true,'LearnRate',250,'Perplexity',50);
figure
gscatter(tSNE_joint(:,1),tSNE_joint(:,2),[S1(S1_assign>0) S2(S2_assign>0)],[],[],10)
set(gca,'FontSize',12);
set(gcf, 'Position', [0, 0, 600 500])
saveas(gcf,[Outdir,'/tSNE_joint_byCluster'],'pdf')
c={'scATAC';'scRNA'};
figure
gscatter(tSNE_joint(:,1),tSNE_joint(:,2),c([ones(sum(S1_assign>0),1);2*ones(sum(S2_assign>0),1)]),[],[],10)
set(gca,'FontSize',12);
set(gcf, 'Position', [0, 0, 600 500])
saveas(gcf,[Outdir,'/tSNE_joint_byDataType'],'pdf')
dlmwrite([Outdir,'/tSNE_atac.txt'],tSNE1,'\t')
dlmwrite([Outdir,'/tSNE_rna.txt'],tSNE2,'\t')
dlmwrite([Outdir,'/tSNE_atac_H.txt'],tSNE1_H,'\t')
dlmwrite([Outdir,'/tSNE_rna_H.txt'],tSNE2_H,'\t')
dlmwrite([Outdir,'/tSNE_joint.txt'],tSNE_joint,'\t')

