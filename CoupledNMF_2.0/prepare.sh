genome=mm10
sort -k1,1 -k2,2n Peak.bed >a
mv a Peak.bed
bedtools intersect -a ../prior/Promoter_100k_${genome}.bed  -b Peak.bed -wa -wb -sorted|awk 'BEGIN{OFS="\t"}{print $5,$6,$7,$4,$3-100000-$6}'|sed 's/-//g'|sort -k1,1 -k2,2n >peak_gene_100k.bed
bedtools intersect -a peak_gene_100k.bed -b ../prior/RE_gene_corr_${genome}.bed -wa -wb -sorted|awk 'BEGIN{OFS="\t"}{if ($4==$9) print $1,$2,$3,$4,$5,$10}'|sed 's/\t/\_/1'|sed 's/\t/\_/1'>peak_gene_100k_corr.bed
