# iBurden
easy burden test 


## 1. download gnomAD data and redo the vep annotation based on your vep version.
### exomes v2.1.1 and coverage file

`wget -c https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz`

`wget -c https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz.tbi`

`wget -c https://storage.googleapis.com/gnomad-public/release/2.1/coverage/exomes/gnomad.exomes.coverage.summary.tsv.bgz`

### genomes v2.1.1 and coverage file

`wget -c https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.vcf.bgz`

`wget -c https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.vcf.bgz.tbi`

`wget -c https://storage.googleapis.com/gnomad-public/release/2.1/coverage/genomes/gnomad.genomes.coverage.summary.tsv.bgz`

## 2. Here use exomes as example. Get coverage for base and control.

* get the coverage for sites with at least dp>10 in 90% samples as suggested at TRAPD (https://github.com/mhguo1/TRAPD).

`zcat gnomad.exomes.coverage.summary.tsv.bgz | tail -n+2 | awk '$7>0.9 {print $1"\t"($2-1)"\t"$2}' | bedtools merge -i stdin > gnomad.dp10.bed`

* That's how to calculate case coverage in the same with with gnomAD. `cases.dp10.bed` is the sites with DP >= 10.

`bash src/coverage/case_cov.sh -p paht/to/bam.files -b exome_bed.file -d 10 -o cases.dp10.bed

* combine the coverage bed file using bedtools. Here you proboly want to intersect case.coverage with gnomAD.exomes.dp10.bed or gnomAD.genomes.dp10.bed

`bedtools intersect -a gnomad.dp10.bed -b cases.dp10.bed | sort -V -k1,1n -k2,2n | bedtools merge -i stdin > combined.dp10.bed`



## 3. redo VEP annotation

Because the default VEP annotation of gnomAD is based on different VEP vesion from the VEP that we have, it's better to re-do VEP annotation for gnomAD data to keep all genes annotation up with case.

`/scratch/cqs/softwares/ensembl-vep/vep -i gnomad.exomes.r2.1.1.sites.vcf.bgz --cache --dir_cache /scratch/cqs/references/vep_data -o gnomAD.exome.vep.vcf --fork 10 --assembly GRCh37 --offline --variant_class –-vcf --canonical --symbol --max_af`


`bgzip -@ 10 gnomAD.exome.vep.vcf`

`tabix -p vcf gnomAD.exome.vep.vcf.gz`



## 4. filter by the coverage file at both contola and case side.

`java -jar $gatk3 -T SelectVariants -R $hg19 -V gnomad.genomes.r2.1.1.sites.reDoVEP.vcf.gz -o ./gnomad.genomes.r2.1.1.sites.combined.reDoVEP.vcf.gz -L /scratch/cqs/baiy7/Loyd_proj/analysis/FM_exome_result/Burden_TEST/coverage_FM_exome/combine_gnomAD_genome.dp10.bed`

`java -jar $gatk3 -T SelectVariants -R $hg19 -V FM_exome.maf_filtered.vep.vcf.gz -o FM_exome.maf_filtered.vep.combined.vcf.gz -L /scratch/cqs/baiy7/Loyd_proj/analysis/FM_exome_result/Burden_TEST/coverage_FM_exome/combine_gnomAD_genome.dp10.bed`



`python src/TRAPD_md/code/make_snp_file.py --vcffile gnomad.genomes.r2.1.1.sites.combined.reDoVEP.vcf.gz --outfile control.snpfile.txt --genecolname SYMBOL --includeinfo "controls_AF_nfe[<=]0.001" --includeinfo "controls_AC_nfe[>]0" --includevep "Consequence[=]missense_variant" --vep --snponly`

`java -jar $gatk3 -T SelectVariants -R $hg19 -V FM_exome.maf_filtered.vep.combined.vcf.gz -o FM_exome.maf_filtered.vep.combined.novel.vcf.gz --discordance gnomad.genomes.r2.1.1.sites.combined.reDoVEP.vcf.gz -selectType SNP`

`src/TRAPD_md/code/make_snp_file.py --vcffile FM_exome.maf_filtered.vep.combined.novel.vcf.gz --outfile FM.novel.snpfile.txt --genecolname SYMBOL --includevep "Consequence[=]missense_variant" --vep --snponly`

`python src/TRAPD_md/code/merge_snp_file.py -s control.snpfile.txt,FM.novel.snpfile.txt -o case.snpfile.txt`

## 5. counts

get information from control/nfe population by using modified `count_controls_md.py`.

`python src/TRAPD_md/code/count_controls_md.py --vcffile gnomad.genomes.r2.1.1.sites.combined.reDoVEP.vcf.gz --snpfile control.snpfile.txt --outfile control_counts.txt --database gnomad`

`python src/TRAPD_md/code/count_cases.py --vcffile FM_exome.maf_filtered.vep.combined.vcf.gz -s case.snpfile.txt -o  case_counts.txt`

`Rscript src/TRAPD_md/code/burden.R --casefile case_counts.txt --casesize 72 --controlfile control_counts.txt --controlsize 2762 --outfile burden.out.txt`

`Rscript src/TRAPD_md/code/QQ.R --pvalfile burden.out.txt --plotfile burden.out.png`

# 3. filter by AF or population

# 4. filter by CCDC or other Score

# FAQ
* how to extract some variation given information from snptable file?
`grep "MADCAM1" case.snpfile.txt|cut -f2|tr "," "\n"|cut -d ":" -f 1,2|awk -F ':' '{ print $1":"$2"-"$2}'|xargs -L 1 -I {} tabix project.vep.vcf.gz {} >>test.vcf`