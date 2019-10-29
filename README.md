# iBurden
easy burden test 


# 1. download gnomAD data and redo the vep annotation based on your vep version.
# exomes v2.1.1 and coverage file
wget -c https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz
wget -c https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz.tbi
wget -c https://storage.googleapis.com/gnomad-public/release/2.1/coverage/exomes/gnomad.exomes.coverage.summary.tsv.bgz

# genomes v2.1.1 and coverage file
wget -c https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.vcf.bgz
wget -c https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.vcf.bgz.tbi
wget -c https://storage.googleapis.com/gnomad-public/release/2.1/coverage/genomes/gnomad.genomes.coverage.summary.tsv.bgz


# Here use exomes as example
# get the coverage for at least dp>10 as suggested at TRAPD (https://github.com/mhguo1/TRAPD)
zcat gnomad.exomes.coverage.summary.tsv.bgz | tail -n+2 | awk '$7>0.9 {print $1"\t"($2-1)"\t"$2}' | bedtools merge -i stdin > gnomad.dp10.bed

bedtools intersect -a gnomad.dp10.bed -b cases.dp10.bed | sort -V -k1,1n -k2,2n | bedtools merge -i stdin > combined.dp10.bed

# redo VEP annotation
`/scratch/cqs/softwares/ensembl-vep/vep -i gnomad.exomes.r2.1.1.sites.vcf.bgz --cache --dir_cache /scratch/cqs/references/vep_data -o gnomAD.exome.vep.vcf --fork 10 --assembly GRCh37 --offline  --variant_class –-vcf --canonical --symbol`
`bgzip -@ 10 gnomAD.exome.vep.vcf`
`tabix -p vcf gnomAD.exome.vep.vcf.gz`



# 2. filter by the coverage file at both contola and case side.
`java -jar $gatk3 -T SelectVariants -R $hg19 -V gnomad.genomes.r2.1.1.sites.reDoVEP.vcf.gz -o ./gnomad.genomes.r2.1.1.sites.combined.reDoVEP.vcf.gz -L /scratch/cqs/baiy7/Loyd_proj/analysis/FM_exome_result/Burden_TEST/coverage_FM_exome/combine_gnomAD_genome.dp10.bed`
`java -jar $gatk3 -T SelectVariants -R $hg19 -V FM_exome.maf_filtered.vep.vcf.gz -o FM_exome.maf_filtered.vep.combined.vcf.gz -L /scratch/cqs/baiy7/Loyd_proj/analysis/FM_exome_result/Burden_TEST/coverage_FM_exome/combine_gnomAD_genome.dp10.bed`

`python /scratch/cqs/baiy7/tools/TRAPD/code/make_snp_file.py --vcffile gnomad.genomes.r2.1.1.sites.combined.reDoVEP.vcf.gz --outfile control.snpfile.txt --genecolname SYMBOL --includeinfo "controls_AF_nfe[<=]0.001" --includeinfo "controls_AC_nfe[>]0" --includevep "Consequence[=]missense_variant" --vep --snponly`

`java -jar $gatk3 -T SelectVariants -R $hg19 -V FM_exome.maf_filtered.vep.combined.vcf.gz -o FM_exome.maf_filtered.vep.combined.novel.vcf.gz --discordance gnomad.genomes.r2.1.1.sites.combined.reDoVEP.vcf.gz -selectType SNP`

/scratch/cqs/baiy7/tools/TRAPD/code/make_snp_file.py --vcffile FM_exome.maf_filtered.vep.combined.novel.vcf.gz --outfile FM.novel.snpfile.txt --genecolname SYMBOL --includevep "Consequence[=]missense_variant" --vep --snponly
python /scratch/cqs/baiy7/tools/TRAPD/code/merge_snp_file.py -s control.snpfile.txt,FM.novel.snpfile.txt -o case.snpfile.txt

# counts
python /scratch/cqs/baiy7/tools/TRAPD/code/count_controls_md.py --vcffile gnomad.genomes.r2.1.1.sites.combined.reDoVEP.vcf.gz --snpfile control.snpfile.txt --outfile control_counts.txt --database gnomad
python /scratch/cqs/baiy7/tools/TRAPD/code/count_cases.py --vcffile FM_exome.maf_filtered.vep.combined.vcf.gz -s case.snpfile.txt -o  case_counts.txt

Rscript /scratch/cqs/baiy7/tools/TRAPD/code/burden.R --casefile case_counts.txt --casesize 72 --controlfile control_counts.txt --controlsize 2762 --outfile burden.out.txt


# 3. filter by AF or population

# 4. filter by CCDC or other Score