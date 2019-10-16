control_cov=""
case_cov=""

bedtools intersect -a $control_cov -b $case_cov | sort -V -k1,1n -k2,2n | bedtools merge -i stdin > combined.dp10.bed


