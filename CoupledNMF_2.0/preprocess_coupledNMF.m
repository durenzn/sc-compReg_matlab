function [X, Symbol, PeakO, PeakName, D] = preprocess_coupledNMF(bin_by_cell_count_file,bin_name_file,atac_barcode,rnaseq_genebycell_folder,species,Outdir)
%bin_by_cell_count_file bin cell count
%bin_name_file, chr str end chr_str_end
%rnaseq_genebycell_folder
%species, human or mouse
%%%%%%%%%%atac-seq process
fileID = fopen(bin_name_file);
C = textscan(fileID,'%s %s %s %s');
fclose(fileID);
PeakName=C{1,4};
opn=[];
opn_count=[];
barcode=[];
for i=1:size(bin_by_cell_count_file,1)
	a=dlmread(bin_by_cell_count_file{i,1},'\t');
	opn=[opn sparse(a(:,1),a(:,2),log10(1+a(:,3)),length(PeakName),max(a(:,2)))];
	opn_count=[opn_count sparse(a(:,1),a(:,2),a(:,3),length(PeakName),max(a(:,2)))];
	barcode=[barcode;importdata(atac_barcode{i,1})];
end
Bin_cell_count=sum(opn'>0)';
Bin_cell_count_z=zscore(log10(1+Bin_cell_count));
%[d f]=sort(Bin_cell_count,'descend');
opn=opn(abs(Bin_cell_count_z)<2,:);
PeakName=PeakName(abs(Bin_cell_count_z)<2,:);
%dep=sum(opn);
PeakO=opn;
%PeakO=PeakO./dep * median(dep);
PeakO=full(PeakO);
if size(PeakO,1)>50000
[PeakO1, PeakName1]=variableSelection_pca(PeakO,PeakName,50000);
else
PeakO1=PeakO;
PeakName1=PeakName;
end
%%%%%%%%%%%%RNA-seq process
A=[];
coding=importdata(['protein_coding_gene_',species]);
for i=1:size(rnaseq_genebycell_folder,1)
	[E Feature]=getmatrix([rnaseq_genebycell_folder{i,1},'/matrix.mtx'],[rnaseq_genebycell_folder{i,1},'/features.tsv'],3);
	[d f]=ismember(coding,Feature(:,2));
        E1=zeros(length(coding),size(E,2));
        E1(d,:)=E(f(d==1),:);
	A=[A E1];
end
Symbol=coding;
[d f]=ismember(Symbol,coding);
Symbol=Symbol(d==1,:);
A=full(A(d==1,:));
Bin_cell_count=sum(A'>0)';
Bin_cell_count_z=zscore(log10(1+Bin_cell_count));
[d f]=sort(Bin_cell_count,'descend');
%A=A(Bin_cell_count<d(100),:);
%Symbol=Symbol(Bin_cell_count<d(100),:);
%A=A(abs(Bin_cell_count_z)<2,:);
%Symbol=Symbol(abs(Bin_cell_count_z)<2,:);
%A=log2(1+quantilenorm(A,'MEDIAN',true));
A=log2(1+A);
[X, Symbol]=variableSelection_pca(A,Symbol,5000);
%X=A;
%%%%%%%%%%%%%coupling matrix
fileID = fopen([Outdir,'/peak_gene_100k_corr.bed']);
C = textscan(fileID,'%s %s %f32 %f32');
fclose(fileID);
[d f]=ismember(C{1,1},PeakName1);
[d1 f1]=ismember(C{1,2},Symbol);
C{1,4}(C{1,4}<0)=0;
d0=30000;
c=double(exp(-1*C{1,3}(d.*d1==1)/d0).*C{1,4}(d.*d1==1));
%c(c<0.02)=0;
D=sparse(f1(d.*d1==1),f(d.*d1==1),c,length(Symbol),length(PeakName1));
PeakO=PeakO1;PeakName=PeakName1;