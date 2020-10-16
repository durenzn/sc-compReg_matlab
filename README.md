# sc-compReg_matlab
 The comparison of gene regulatory networks between disease versus healthy or between two different treatments is an important scientific problem. We propose sc-compReg as a method for the single cell comparative regulatory analysis based on single cell gene expression (scRNA-seq) and single cell chromatin accessibility data (scATAC-seq) from two different conditions. Our software sc-compReg can be used as a stand-alone package that provides joint clustering and embedding of the cells from both scRNA-seq and scATAC-seq, and construction of differential regulatory networks across two conditions.

## install



## Run sc-compReg:

Run sc-compReg by following two steps:

### Step 1: get cell label (here we recomend CoupledNMF method)

script

### Step 2: comparative Regulatory analysis
    1) pre-processing: 
         edit genome version in line 4 of processing_data_sc_compReg.sh, we support hg19, hg39, mm9 and mm10
         provide Peak name file of sample 1 and sample 2. example: PeakName1.txt and PeakName2.txt
         bash processing_data_sc_compReg.sh
         This step may take long time becasue it scan motif across entire genome (to generate the TFbinding matrix used in next step).
    2) run sc-compReg:
         See example in demo.m file. Please dowmload the data from the following link: www.***
         note 1: the input is the  sample1.mat file and sample2.mat file. Plesae see example.
         note 2: Generatation of TFbinding.mat file is time consuming, so we have provided this data in this example. To generate your TFbinding on your owen data, please use 3 line codes start with %:
            load('./prior/MotifMatch_human_rmdup.mat')
            TFName=intersect(Symbol,unique(Match2(:,2)));
            TF_binding=mfbs(TFName,Element_name,motifName,motifWeight,Match2);

