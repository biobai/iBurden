


## extract vcf by using the combined bed file
#java -jar $gatk3 -T SelectVariants -R $hg19 -V /scratch/cqs/baiy7/tools/gnomad/gnomad.exomes.r2.1.1.sites.rminfo.vcf.gz -o /scratch/cqs/baiy7/Loyd_proj/analysis/FM_exome_result/Burden_TEST/exome_control/gnomad.exomes.r2.1.1.sites.rminfo.combined.vcf.gz -L /scratch/cqs/baiy7/tools/iBurden/src/coverage/result/combine_case_gnomAD_exome.dp10.bed
java -jar $gatk3 -T SelectVariants -R $hg19 -V gnomad.genomes.r2.1.1.sites.rminfo.vcf.gz -o /scratch/cqs/baiy7/Loyd_proj/analysis/FM_exome_result/Burden_TEST/genome_control/gnomad.genomes.r2.1.1.sites.rminfo.combined.vcf.gz -L /scratch/cqs/baiy7/tools/iBurden/src/coverage/result/combine_case_gnomAD_genome.dp10.bed