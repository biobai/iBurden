#! /usr/bin/bash

set -e
set -o pipefail
version=0.1
LC_ALL=C
##  usage
function usage() {
    echo "
Version: $version
usage:   vt_rminfo.sh [options]

options: 
        -h       show this message
        -b STR   bed file for exome data [null]
        -d INT   threshold for depth [10]
        -o STR   output file
        -p STR   path to bam files
        -t INT   threads [8]
"
}


if test -z "$1"
then
    usage
exit 1
fi

while getopts ":hb:t:o:p:d:" OPTION
do
    case "${OPTION}" in
        h)
            usage
            exit 1
            ;;
        b)
            exome_target_file="$OPTARG"
            ;;
        t)
            threads="$OPTARG"
            ;;
        o)
            output="$OPTARG"
            ;;
        p)
            bam_path="$OPTARG"
            ;;
        d)
            depth="$OPTARG"
            ;;
    esac
done

if [[ -z $threads ]]; then
    threads="8"
fi
if [[ -z $exome_target_file ]]; then
    exome_target_file=""
fi
if [[ -z $output ]]; then
    output="case.dp10.bed"
fi
if [[ -z $output ]]; then
    output="case.dp10.bed"
fi
if [[ -z $depth ]]; then
    depth=10
fi

if [ -x "samtools --version" ]
then
  echo 'Error: samtools is not installed.' >&2
  exit 1
fi

bcftools view -f PASS -Ov gnomad.genomes.r2.1.1.sites.vcf.bgz | \
parallel -j 2 -N1 \
"vt rminfo -t \
ab_hist_alt_bin_freq,age_hist_het_bin_freq,age_hist_het_n_larger,age_hist_het_n_smaller,age_hist_hom_bin_freq,age_hist_hom_n_larger,age_hist_hom_n_smaller,allele_type,BaseQRankSum,ClippingRankSum,\
controls_AC,controls_AC_afr,controls_AC_afr_female,controls_AC_afr_male,controls_AC_amr,controls_AC_amr_female,controls_AC_amr_male,controls_AC_asj,controls_AC_asj_female,controls_AC_asj_male,controls_AC_eas,controls_AC_eas_female,controls_AC_eas_jpn,controls_AC_eas_kor,controls_AC_eas_male,controls_AC_eas_oea,controls_AC_female,controls_AC_fin,controls_AC_fin_female,controls_AC_fin_male,controls_AC_male,controls_AC_nfe,controls_AC_nfe_bgr,controls_AC_nfe_est,controls_AC_nfe_female,controls_AC_nfe_male,controls_AC_nfe_nwe,controls_AC_nfe_onf,controls_AC_nfe_seu,controls_AC_nfe_swe,controls_AC_oth,controls_AC_oth_female,controls_AC_oth_male,controls_AC_raw,controls_AC_sas,controls_AC_sas_female,controls_AC_sas_male,controls_AF,controls_AF_afr,controls_AF_afr_female,controls_AF_afr_male,controls_AF_amr,controls_AF_amr_female,controls_AF_amr_male,controls_AF_asj,controls_AF_asj_female,controls_AF_asj_male,controls_AF_eas,controls_AF_eas_female,controls_AF_eas_jpn,controls_AF_eas_kor,controls_AF_eas_male,controls_AF_eas_oea,controls_AF_female,controls_AF_fin,controls_AF_fin_female,controls_AF_fin_male,controls_AF_male,controls_AF_nfe,controls_AF_nfe_bgr,controls_AF_nfe_est,controls_AF_nfe_female,controls_AF_nfe_male,controls_AF_nfe_nwe,controls_AF_nfe_onf,controls_AF_nfe_seu,controls_AF_nfe_swe,controls_AF_oth,controls_AF_oth_female,controls_AF_oth_male,controls_AF_raw,controls_AF_sas,controls_AF_sas_female,controls_AF_sas_male,controls_AN,controls_AN_afr,controls_AN_afr_female,controls_AN_afr_male,controls_AN_amr,controls_AN_amr_female,controls_AN_amr_male,controls_AN_asj,controls_AN_asj_female,controls_AN_asj_male,controls_AN_eas,controls_AN_eas_female,controls_AN_eas_jpn,controls_AN_eas_kor,controls_AN_eas_male,controls_AN_eas_oea,controls_AN_female,controls_AN_fin,controls_AN_fin_female,controls_AN_fin_male,controls_AN_male,controls_AN_nfe,controls_AN_nfe_bgr,controls_AN_nfe_est,controls_AN_nfe_female,controls_AN_nfe_male,controls_AN_nfe_nwe,controls_AN_nfe_onf,controls_AN_nfe_seu,controls_AN_nfe_swe,controls_AN_oth,controls_AN_oth_female,controls_AN_oth_male,controls_AN_raw,controls_AN_sas,controls_AN_sas_female,controls_AN_sas_male,controls_faf95,controls_faf95_afr,controls_faf95_amr,controls_faf95_eas,controls_faf95_nfe,controls_faf95_sas,controls_faf99,controls_faf99_afr,controls_faf99_amr,controls_faf99_eas,controls_faf99_nfe,controls_faf99_sas,controls_nhomalt,controls_nhomalt_afr,controls_nhomalt_afr_female,controls_nhomalt_afr_male,controls_nhomalt_amr,controls_nhomalt_amr_female,controls_nhomalt_amr_male,controls_nhomalt_asj,controls_nhomalt_asj_female,controls_nhomalt_asj_male,controls_nhomalt_eas,controls_nhomalt_eas_female,controls_nhomalt_eas_jpn,controls_nhomalt_eas_kor,controls_nhomalt_eas_male,controls_nhomalt_eas_oea,controls_nhomalt_female,controls_nhomalt_fin,controls_nhomalt_fin_female,controls_nhomalt_fin_male,controls_nhomalt_male,controls_nhomalt_nfe,controls_nhomalt_nfe_bgr,controls_nhomalt_nfe_est,controls_nhomalt_nfe_female,controls_nhomalt_nfe_male,controls_nhomalt_nfe_nwe,controls_nhomalt_nfe_onf,controls_nhomalt_nfe_seu,controls_nhomalt_nfe_swe,controls_nhomalt_oth,controls_nhomalt_oth_female,controls_nhomalt_oth_male,controls_nhomalt_raw,controls_nhomalt_sas,controls_nhomalt_sas_female,controls_nhomalt_sas_male,\
decoy,DP,dp_hist_all_bin_freq,dp_hist_all_n_larger,dp_hist_alt_bin_freq,dp_hist_alt_n_larger,\
faf95_afr,faf95_amr,faf95_eas,faf95_nfe,faf95_sas,faf99,faf99_afr,faf99_amr,faf99_eas,faf99_nfe,faf99_sas,\
FS,gq_hist_all_bin_freq,gq_hist_alt_bin_freq,has_star,InbreedingCoeff,lcr,MQ,MQRankSum,n_alt_alleles,\
non_cancer_AC,non_cancer_AC_afr,non_cancer_AC_afr_female,non_cancer_AC_afr_male,non_cancer_AC_amr,non_cancer_AC_amr_female,non_cancer_AC_amr_male,non_cancer_AC_asj,non_cancer_AC_asj_female,non_cancer_AC_asj_male,non_cancer_AC_eas,non_cancer_AC_eas_female,non_cancer_AC_eas_jpn,non_cancer_AC_eas_kor,non_cancer_AC_eas_male,non_cancer_AC_eas_oea,non_cancer_AC_female,non_cancer_AC_fin,non_cancer_AC_fin_female,non_cancer_AC_fin_male,non_cancer_AC_male,non_cancer_AC_nfe,non_cancer_AC_nfe_bgr,non_cancer_AC_nfe_est,non_cancer_AC_nfe_female,non_cancer_AC_nfe_male,non_cancer_AC_nfe_nwe,non_cancer_AC_nfe_onf,non_cancer_AC_nfe_seu,non_cancer_AC_nfe_swe,non_cancer_AC_oth,non_cancer_AC_oth_female,non_cancer_AC_oth_male,non_cancer_AC_raw,non_cancer_AC_sas,non_cancer_AC_sas_female,non_cancer_AC_sas_male,non_cancer_AF,non_cancer_AF_afr,non_cancer_AF_afr_female,non_cancer_AF_afr_male,non_cancer_AF_amr,non_cancer_AF_amr_female,non_cancer_AF_amr_male,non_cancer_AF_asj,non_cancer_AF_asj_female,non_cancer_AF_asj_male,non_cancer_AF_eas,non_cancer_AF_eas_female,non_cancer_AF_eas_jpn,non_cancer_AF_eas_kor,non_cancer_AF_eas_male,non_cancer_AF_eas_oea,non_cancer_AF_female,non_cancer_AF_fin,non_cancer_AF_fin_female,non_cancer_AF_fin_male,non_cancer_AF_male,non_cancer_AF_nfe,non_cancer_AF_nfe_bgr,non_cancer_AF_nfe_est,non_cancer_AF_nfe_female,non_cancer_AF_nfe_male,non_cancer_AF_nfe_nwe,non_cancer_AF_nfe_onf,non_cancer_AF_nfe_seu,non_cancer_AF_nfe_swe,non_cancer_AF_oth,non_cancer_AF_oth_female,non_cancer_AF_oth_male,non_cancer_AF_raw,non_cancer_AF_sas,non_cancer_AF_sas_female,non_cancer_AF_sas_male,non_cancer_AN,non_cancer_AN_afr,non_cancer_AN_afr_female,non_cancer_AN_afr_male,non_cancer_AN_amr,non_cancer_AN_amr_female,non_cancer_AN_amr_male,non_cancer_AN_asj,non_cancer_AN_asj_female,non_cancer_AN_asj_male,non_cancer_AN_eas,non_cancer_AN_eas_female,non_cancer_AN_eas_jpn,non_cancer_AN_eas_kor,non_cancer_AN_eas_male,non_cancer_AN_eas_oea,non_cancer_AN_female,non_cancer_AN_fin,non_cancer_AN_fin_female,non_cancer_AN_fin_male,non_cancer_AN_male,non_cancer_AN_nfe,non_cancer_AN_nfe_bgr,non_cancer_AN_nfe_est,non_cancer_AN_nfe_female,non_cancer_AN_nfe_male,non_cancer_AN_nfe_nwe,non_cancer_AN_nfe_onf,non_cancer_AN_nfe_seu,non_cancer_AN_nfe_swe,non_cancer_AN_oth,non_cancer_AN_oth_female,non_cancer_AN_oth_male,non_cancer_AN_raw,non_cancer_AN_sas,non_cancer_AN_sas_female,non_cancer_AN_sas_male,non_cancer_faf95,non_cancer_faf95_afr,non_cancer_faf95_amr,non_cancer_faf95_eas,non_cancer_faf95_nfe,non_cancer_faf95_sas,non_cancer_faf99,non_cancer_faf99_afr,non_cancer_faf99_amr,non_cancer_faf99_eas,non_cancer_faf99_nfe,non_cancer_faf99_sas,non_cancer_nhomalt,non_cancer_nhomalt_afr,non_cancer_nhomalt_afr_female,non_cancer_nhomalt_afr_male,non_cancer_nhomalt_amr,non_cancer_nhomalt_amr_female,non_cancer_nhomalt_amr_male,non_cancer_nhomalt_asj,non_cancer_nhomalt_asj_female,non_cancer_nhomalt_asj_male,non_cancer_nhomalt_eas,non_cancer_nhomalt_eas_female,non_cancer_nhomalt_eas_jpn,non_cancer_nhomalt_eas_kor,non_cancer_nhomalt_eas_male,non_cancer_nhomalt_eas_oea,non_cancer_nhomalt_female,non_cancer_nhomalt_fin,non_cancer_nhomalt_fin_female,non_cancer_nhomalt_fin_male,non_cancer_nhomalt_male,non_cancer_nhomalt_nfe,non_cancer_nhomalt_nfe_bgr,non_cancer_nhomalt_nfe_est,non_cancer_nhomalt_nfe_female,non_cancer_nhomalt_nfe_male,non_cancer_nhomalt_nfe_nwe,non_cancer_nhomalt_nfe_onf,non_cancer_nhomalt_nfe_seu,non_cancer_nhomalt_nfe_swe,non_cancer_nhomalt_oth,non_cancer_nhomalt_oth_female,non_cancer_nhomalt_oth_male,non_cancer_nhomalt_raw,non_cancer_nhomalt_sas,non_cancer_nhomalt_sas_female,non_cancer_nhomalt_sas_male,\
non_neuro_AC,non_neuro_AC_afr,non_neuro_AC_afr_female,non_neuro_AC_afr_male,non_neuro_AC_amr,non_neuro_AC_amr_female,non_neuro_AC_amr_male,non_neuro_AC_asj,non_neuro_AC_asj_female,non_neuro_AC_asj_male,non_neuro_AC_eas,non_neuro_AC_eas_female,non_neuro_AC_eas_jpn,non_neuro_AC_eas_kor,non_neuro_AC_eas_male,non_neuro_AC_eas_oea,non_neuro_AC_female,non_neuro_AC_fin,non_neuro_AC_fin_female,non_neuro_AC_fin_male,non_neuro_AC_male,non_neuro_AC_nfe,non_neuro_AC_nfe_bgr,non_neuro_AC_nfe_est,non_neuro_AC_nfe_female,non_neuro_AC_nfe_male,non_neuro_AC_nfe_nwe,non_neuro_AC_nfe_onf,non_neuro_AC_nfe_seu,non_neuro_AC_nfe_swe,non_neuro_AC_oth,non_neuro_AC_oth_female,non_neuro_AC_oth_male,non_neuro_AC_popmax,non_neuro_AC_raw,non_neuro_AC_sas,non_neuro_AC_sas_female,non_neuro_AC_sas_male,non_neuro_AF,non_neuro_AF_afr,non_neuro_AF_afr_female,non_neuro_AF_afr_male,non_neuro_AF_amr,non_neuro_AF_amr_female,non_neuro_AF_amr_male,non_neuro_AF_asj,non_neuro_AF_asj_female,non_neuro_AF_asj_male,non_neuro_AF_eas,non_neuro_AF_eas_female,non_neuro_AF_eas_jpn,non_neuro_AF_eas_kor,non_neuro_AF_eas_male,non_neuro_AF_eas_oea,non_neuro_AF_female,non_neuro_AF_fin,non_neuro_AF_fin_female,non_neuro_AF_fin_male,non_neuro_AF_male,non_neuro_AF_nfe,non_neuro_AF_nfe_bgr,non_neuro_AF_nfe_est,non_neuro_AF_nfe_female,non_neuro_AF_nfe_male,non_neuro_AF_nfe_nwe,non_neuro_AF_nfe_onf,non_neuro_AF_nfe_seu,non_neuro_AF_nfe_swe,non_neuro_AF_oth,non_neuro_AF_oth_female,non_neuro_AF_oth_male,non_neuro_AF_raw,non_neuro_AF_sas,non_neuro_AF_sas_female,non_neuro_AF_sas_male,non_neuro_AN,non_neuro_AN_afr,non_neuro_AN_afr_female,non_neuro_AN_afr_male,non_neuro_AN_amr,non_neuro_AN_amr_female,non_neuro_AN_amr_male,non_neuro_AN_asj,non_neuro_AN_asj_female,non_neuro_AN_asj_male,non_neuro_AN_eas,non_neuro_AN_eas_female,non_neuro_AN_eas_jpn,non_neuro_AN_eas_kor,non_neuro_AN_eas_male,non_neuro_AN_eas_oea,non_neuro_AN_female,non_neuro_AN_fin,non_neuro_AN_fin_female,non_neuro_AN_fin_male,non_neuro_AN_male,non_neuro_AN_nfe,non_neuro_AN_nfe_bgr,non_neuro_AN_nfe_est,non_neuro_AN_nfe_female,non_neuro_AN_nfe_male,non_neuro_AN_nfe_nwe,non_neuro_AN_nfe_onf,non_neuro_AN_nfe_seu,non_neuro_AN_nfe_swe,non_neuro_AN_oth,non_neuro_AN_oth_female,non_neuro_AN_oth_male,non_neuro_AN_raw,non_neuro_AN_sas,non_neuro_AN_sas_female,non_neuro_AN_sas_male,non_neuro_faf95,non_neuro_faf95_afr,non_neuro_faf95_amr,non_neuro_faf95_eas,non_neuro_faf95_nfe,non_neuro_faf95_sas,non_neuro_faf99,non_neuro_faf99_afr,non_neuro_faf99_amr,non_neuro_faf99_eas,non_neuro_faf99_nfe,non_neuro_faf99_sas,non_neuro_nhomalt,non_neuro_nhomalt_afr,non_neuro_nhomalt_afr_female,non_neuro_nhomalt_afr_male,non_neuro_nhomalt_amr,non_neuro_nhomalt_amr_female,non_neuro_nhomalt_amr_male,non_neuro_nhomalt_asj,non_neuro_nhomalt_asj_female,non_neuro_nhomalt_asj_male,non_neuro_nhomalt_eas,non_neuro_nhomalt_eas_female,non_neuro_nhomalt_eas_jpn,non_neuro_nhomalt_eas_kor,non_neuro_nhomalt_eas_male,non_neuro_nhomalt_eas_oea,non_neuro_nhomalt_female,non_neuro_nhomalt_fin,non_neuro_nhomalt_fin_female,non_neuro_nhomalt_fin_male,non_neuro_nhomalt_male,non_neuro_nhomalt_nfe,non_neuro_nhomalt_nfe_bgr,non_neuro_nhomalt_nfe_est,non_neuro_nhomalt_nfe_female,non_neuro_nhomalt_nfe_male,non_neuro_nhomalt_nfe_nwe,non_neuro_nhomalt_nfe_onf,non_neuro_nhomalt_nfe_seu,non_neuro_nhomalt_nfe_swe,non_neuro_nhomalt_oth,non_neuro_nhomalt_oth_female,non_neuro_nhomalt_oth_male,non_neuro_nhomalt_raw,non_neuro_nhomalt_sas,non_neuro_nhomalt_sas_female,non_neuro_nhomalt_sas_male,non_neuro_popmax,\
non_topmed_AC,non_topmed_AC_afr,non_topmed_AC_afr_female,non_topmed_AC_afr_male,non_topmed_AC_amr,non_topmed_AC_amr_female,non_topmed_AC_amr_male,non_topmed_AC_asj,non_topmed_AC_asj_female,non_topmed_AC_asj_male,non_topmed_AC_eas,non_topmed_AC_eas_female,non_topmed_AC_eas_jpn,non_topmed_AC_eas_kor,non_topmed_AC_eas_male,non_topmed_AC_eas_oea,non_topmed_AC_female,non_topmed_AC_fin,non_topmed_AC_fin_female,non_topmed_AC_fin_male,non_topmed_AC_male,non_topmed_AC_nfe,non_topmed_AC_nfe_bgr,non_topmed_AC_nfe_est,non_topmed_AC_nfe_female,non_topmed_AC_nfe_male,non_topmed_AC_nfe_nwe,non_topmed_AC_nfe_onf,non_topmed_AC_nfe_seu,non_topmed_AC_nfe_swe,non_topmed_AC_oth,non_topmed_AC_oth_female,non_topmed_AC_oth_male,non_topmed_AC_popmax,non_topmed_AC_raw,non_topmed_AC_sas,non_topmed_AC_sas_female,non_topmed_AC_sas_male,non_topmed_AF,non_topmed_AF_afr,non_topmed_AF_afr_female,non_topmed_AF_afr_male,non_topmed_AF_amr,non_topmed_AF_amr_female,non_topmed_AF_amr_male,non_topmed_AF_asj,non_topmed_AF_asj_female,non_topmed_AF_asj_male,non_topmed_AF_eas,non_topmed_AF_eas_female,non_topmed_AF_eas_jpn,non_topmed_AF_eas_kor,non_topmed_AF_eas_male,non_topmed_AF_eas_oea,non_topmed_AF_female,non_topmed_AF_fin,non_topmed_AF_fin_female,non_topmed_AF_fin_male,non_topmed_AF_male,non_topmed_AF_nfe,non_topmed_AF_nfe_bgr,non_topmed_AF_nfe_est,non_topmed_AF_nfe_female,non_topmed_AF_nfe_male,non_topmed_AF_nfe_nwe,non_topmed_AF_nfe_onf,non_topmed_AF_nfe_seu,non_topmed_AF_nfe_swe,non_topmed_AF_oth,non_topmed_AF_oth_female,non_topmed_AF_oth_male,non_topmed_AF_popmax,non_topmed_AF_raw,non_topmed_AF_sas,non_topmed_AF_sas_female,non_topmed_AF_sas_male,non_topmed_AN,non_topmed_AN_afr,non_topmed_AN_afr_female,non_topmed_AN_afr_male,non_topmed_AN_amr,non_topmed_AN_amr_female,non_topmed_AN_amr_male,non_topmed_AN_asj,non_topmed_AN_asj_female,non_topmed_AN_asj_male,non_topmed_AN_eas,non_topmed_AN_eas_female,non_topmed_AN_eas_jpn,non_topmed_AN_eas_kor,non_topmed_AN_eas_male,non_topmed_AN_eas_oea,non_topmed_AN_female,non_topmed_AN_fin,non_topmed_AN_fin_female,non_topmed_AN_fin_male,non_topmed_AN_male,non_topmed_AN_nfe,non_topmed_AN_nfe_bgr,non_topmed_AN_nfe_est,non_topmed_AN_nfe_female,non_topmed_AN_nfe_male,non_topmed_AN_nfe_nwe,non_topmed_AN_nfe_onf,non_topmed_AN_nfe_seu,non_topmed_AN_nfe_swe,non_topmed_AN_oth,non_topmed_AN_oth_female,non_topmed_AN_oth_male,non_topmed_AN_popmax,non_topmed_AN_raw,non_topmed_AN_sas,non_topmed_AN_sas_female,non_topmed_AN_sas_male,non_topmed_faf95,non_topmed_faf95_afr,non_topmed_faf95_amr,non_topmed_faf95_eas,non_topmed_faf95_nfe,non_topmed_faf95_sas,non_topmed_faf99,non_topmed_faf99_afr,non_topmed_faf99_amr,non_topmed_faf99_eas,non_topmed_faf99_nfe,non_topmed_faf99_sas,non_topmed_nhomalt,non_topmed_nhomalt_afr,non_topmed_nhomalt_afr_female,non_topmed_nhomalt_afr_male,non_topmed_nhomalt_amr,non_topmed_nhomalt_amr_female,non_topmed_nhomalt_amr_male,non_topmed_nhomalt_asj,non_topmed_nhomalt_asj_female,non_topmed_nhomalt_asj_male,non_topmed_nhomalt_eas,non_topmed_nhomalt_eas_female,non_topmed_nhomalt_eas_jpn,non_topmed_nhomalt_eas_kor,non_topmed_nhomalt_eas_male,non_topmed_nhomalt_eas_oea,non_topmed_nhomalt_female,non_topmed_nhomalt_fin,non_topmed_nhomalt_fin_female,non_topmed_nhomalt_fin_male,non_topmed_nhomalt_male,non_topmed_nhomalt_nfe,non_topmed_nhomalt_nfe_bgr,non_topmed_nhomalt_nfe_est,non_topmed_nhomalt_nfe_female,non_topmed_nhomalt_nfe_male,non_topmed_nhomalt_nfe_nwe,non_topmed_nhomalt_nfe_onf,non_topmed_nhomalt_nfe_seu,non_topmed_nhomalt_nfe_swe,non_topmed_nhomalt_oth,non_topmed_nhomalt_oth_female,non_topmed_nhomalt_oth_male,non_topmed_nhomalt_popmax,non_topmed_nhomalt_raw,non_topmed_nhomalt_sas,non_topmed_nhomalt_sas_female,non_topmed_nhomalt_sas_male,non_topmed_popmax,\
nonpar,pab_max,QD,ReadPosRankSum,rf_label,rf_negative_label,rf_positive_label,rf_tp_probability,rf_train,segdup,SOR,transmitted_singleton,variant_type,VQSLOD,VQSR_culprit,VQSR_NEGATIVE_TRAIN_SITE,VQSR_POSITIVE_TRAIN_SITE,was_mixed /data1/baiy7/TMP/gatk_{}.vcf.gz |\
bgzip -c >/data1/baiy7/TMP/gatk_{}.rm.vcf.gz " ::: `seq 22 -1 1` X Y