##how to filter for sites with > 90% of samples having DP > 10 
##https://github.com/mhguo1/TRAPD
gnomAD_coverage=""   ## path to the gnomAD covarege file

zcat $gnomAD_coverage | tail -n+2 | awk '$7>0.9 {print $1"\t"($2-1)"\t"$2}' | bedtools merge -i stdin > gnomad.dp10.bed

