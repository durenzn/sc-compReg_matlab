
%%%%%%%%%%%%%%%%%%
if exist('arfa','var')==0
arfa=0.5
end
if exist('betamax_scale','var')==0
betamax_scale=4
end
if exist('atac_binarize','var')==0
atac_binarize=0
end

[X, Symbol, PeakO, PeakName, D] = preprocess_coupledNMF(bin_by_cell_count_file,bin_name_file,atac_barcode,rnaseq_genebycell_folder,species,Outdir);
PeakO=sparse(PeakO);
save([Outdir,'/CouplingData.mat'],'X','Symbol','PeakO','PeakName','D')
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
addpath ./umap_matlab/umap
addpath ./umap_matlab/util
javaaddpath('./umap_matlab/umap/umap.jar');
[COEFF_Opn, SCORE_Opn, LATENT_Opn] = fastpca(PeakO',20+rmPC);
dep1=sum(PeakO>0);
r1=corr(SCORE_Opn,dep1');
[d f]=sort(abs(r1),'descend');
Umap1 = run_umap(SCORE_Opn(:,setdiff([1:20+rmPC],f(1:rmPC))),'verbose','none');
[COEFF_Exp, SCORE_Exp, LATENT_Exp] = fastpca(X',20+rmPC);
dep2=sum(X>0);
r2=corr(SCORE_Exp,dep2');
[d f]=sort(abs(r2),'descend');
Umap2 = run_umap(SCORE_Exp(:,setdiff([1:20+rmPC],f(1:rmPC))),'verbose','none');
figure
gscatter(Umap1(:,1),Umap1(:,2),S1,[],[],10)
set(gca,'FontSize',12);
set(gcf, 'Position', [0, 0, 600 500])
saveas(gcf,[Outdir,'/Umap_opn'],'pdf')
figure
gscatter(Umap2(:,1),Umap2(:,2),S2,[],[],10)
set(gca,'FontSize',12);
set(gcf, 'Position', [0, 0, 600 500])
saveas(gcf,[Outdir,'/Umap_exp'],'pdf')
%%%%CoupledNMF based embedding
Umap1_H = run_umap(H1(:,S1_assign>0)','verbose','none');
Umap2_H = run_umap(H2(:,S2_assign>0)','verbose','none');
figure
gscatter(Umap1_H(:,1),Umap1_H(:,2),S1(S1_assign>0),[],[],10)
set(gca,'FontSize',12);
set(gcf, 'Position', [0, 0, 600 500])
saveas(gcf,[Outdir,'/Umap_opn_H'],'pdf')
figure
gscatter(Umap2_H(:,1),Umap2_H(:,2),S2(S2_assign>0),[],[],10)
set(gca,'FontSize',12);
set(gcf, 'Position', [0, 0, 600 500])
saveas(gcf,[Outdir,'/Umap_exp_H'],'pdf')
%%%%%joint embedding
HH=[H1(:,S1_assign>0) H2(:,S2_assign>0)];
HH=HH./sqrt(sum(HH.^2));
Umap_joint = run_umap(HH','verbose','none');
figure
gscatter(Umap_joint(:,1),Umap_joint(:,2),[S1(S1_assign>0) S2(S2_assign>0)],[],[],10)
set(gca,'FontSize',12);
set(gcf, 'Position', [0, 0, 600 500])
saveas(gcf,[Outdir,'/Umap_joint_byCluster'],'pdf')
c={'scATAC';'scRNA'};
figure
gscatter(Umap_joint(:,1),Umap_joint(:,2),c([ones(sum(S1_assign>0),1);2*ones(sum(S2_assign>0),1)]),[],[],10)
set(gca,'FontSize',12);
set(gcf, 'Position', [0, 0, 600 500])
saveas(gcf,[Outdir,'/Umap_joint_byDataType'],'pdf')
dlmwrite([Outdir,'/Umap_atac.txt'],Umap1,'\t')
dlmwrite([Outdir,'/Umap_rna.txt'],Umap2,'\t')
dlmwrite([Outdir,'/Umap_atac_H.txt'],Umap1_H,'\t')
dlmwrite([Outdir,'/Umap_rna_H.txt'],Umap2_H,'\t')
dlmwrite([Outdir,'/Umap_joint.txt'],Umap_joint,'\t')

