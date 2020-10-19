%%step 1, CoupledNMF
% step 2, subpopulation linking
load('./example/sample1.mat')
load('./example/sample2.mat')
fileID = fopen('./example/PeakName_intersect.txt');
C = textscan(fileID,'%s %s');
fclose(fileID);
PeakName_intersect=[C{1,1} C{1,2}];
[E1_mean,E2_mean,Symbol,O1_mean,O2_mean,Element_name]=cluster_profile(O1,E1,O1_idx,E1_idx,O2,E2,O2_idx,E2_idx,Symbol1,Symbol2,PeakName1,PeakName2,PeakName_intersect);
%load('Opn_bulk.mat')
Match=subpopulation_link(E1_mean,E2_mean,O1_mean,O2_mean);
%step 3, compReg
%load('MotifMatch_human_rmdup.mat')
%TFName=intersect(Symbol,unique(Match2(:,2)));
%TF_binding=mfbs(TFName,Element_name,motifName,motifWeight,Match2);
load('./example/TFbinding.mat')
[diffNet,Hub_TF]=compReg(TF_binding,Match,E1,E1_idx,E2,E2_idx,O1_mean,O2_mean,Symbol,TFName,Element_name,'peak_gene_prior_intersect.bed');
