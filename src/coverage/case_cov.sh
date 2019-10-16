#!/usr/bin/bash

##  Coverage was calculated separately for exomes and genomes on a ~10% subset of the samples using the samtools depth tool. The base quality threshold was set to 10 for the -q option and the mapping quality threshold set to 20 for the -Q option. It is calculated per base of the respective calling intervals, includes sites with zero depth (-a flag), and is capped at 100x for a given sample and base pair. Mean coverage is then plotted on the browser. The numbers in columns 1, 5, 10, etc of our downloadable coverage files refer to the fraction of samples with a depth of coverage of at least 1 read, 5 reads, 10 reads, etc. for the given chromosome and position.
set -e
set -o pipefail
version=0.1
LC_ALL=C
##  usage
function usage() {
    echo "
Version: $version
usage:   case_cov.sh [options]

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

## xgen target file
# exome_target_file="/scratch/cqs/references/exomeseq/IDT/xgen-exome-research-panel-targetsae255a1532796e2eaa53ff00001c1b3c.slop50.nochr.bed" 
## the path of the 'bwa_refine' result
# bam_path="/scratch/cqs/baiy7/Loyd_proj/analysis/FM_exome_result/bwa_refine/result" 
##                          
# threads=8                             



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


sample_num=`ls ${bam_path}/*bam|wc -l`

##  samtools depth calulation
parallel -j $threads -N1 samtools depth -a -b $exome_target_file -d 100 -q 10 -Q 20 -r 1:1-100000 {} '>' {/}.cov ::: `ls ${bam_path}/*bam`

##  merge coverage into one file
parallel -j $threads -N1 "cut -f3 {} >{}_tmp" ::: `ls -v *cov |tail -n+2`
ls -v *_tmp|xargs paste >right_column
paste `ls -v *cov |head -n 1 ` right_column >all.coverage
rm *_tmp *.cov right_column

cat all.coverage |awk -v var="$depth" -v num="sample_num" '{count=0} {for(i=3; i<=$num; i++) if($i>var) count++} {if(count/num>0.9) print $1"\t"($2-1)"\t"$2}'|bedtools merge -i stdin > $output

##  save the total coverage file
#pigz -p $threads all.coverage


