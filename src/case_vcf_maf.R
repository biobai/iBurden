#! /usr/bin/env Rscript
'extract snp counts in gene level for case datasets
Usage:
    case_vcf_maf.R [--vcfSource] [--snvType] <input> <output>
    
Options:
    -h --help  Show this screen.
    -v --version  Show version.
    --vcfSource=<DB>   contral vcf file type, could be ExAC, gnomAD or else. [default: gnomAD]
    --snvType=<missense> missense or other "missense_variant|start_lost|stop_lost|protein_altering_variant|stop_gained" [default: missense]
Arguments:
    input  contral vcf file, could be download from ExAC, genomAD or else
    output  output filename
' -> doc


suppressMessages(library(maftools))
suppressMessages(library(data.table))
suppressMessages(library(docopt))
arguments <- docopt(doc, version = 'control_table v0.1\n\n')
print(str(arguments))

mafFile=arguments$input
maf=read.maf(maf =mafFile)

#Double Check needed for this
#which(maf@data$Tumor_Seq_Allele1==maf@data$Tumor_Seq_Allele2)
#integer(0)
#which(maf@data$Tumor_Seq_Allele1!=maf@data$Reference_Allele)
#integer(0)

#table(maf@data$Variant_Classification)
#table(maf@data$Variant_Classification)
#table(maf@data$Variant_Classification,maf@data$Variant_Classification)


makeMafToGeneBurdenCount<-function(maf,SELECTED_VAR_TYPE="Missense_Mutation",varTypeCol="Variant_Classification") {	
	totalSampleNum=as.integer(maf@summary[3,"summary"])
	selectedVarInd=grep(SELECTED_VAR_TYPE,maf@data[[varTypeCol]])
	if (length(selectedVarInd)>0) {
		mafSubForCount=maf@data[selectedVarInd,]
	} else { #Not found SELECTED_VAR_TYPE, maybe in slient data
		selectedVarInd=grep(SELECTED_VAR_TYPE,maf@maf.silent[[varTypeCol]])
		mafSubForCount=maf@maf.silent[selectedVarInd,]
	}
	if (length(selectedVarInd)==0) { #Still can't find SELECTED_VAR_TYPE, Stop
		stop(paste0("Can't find ",SELECTED_VAR_TYPE," in MAF"))
	}
	
	
	mafSubForCount$HetOrHomoVar="None"
	mafSubForCount$HetOrHomoVar[which(mafSubForCount$Tumor_Seq_Allele2!=mafSubForCount$Reference_Allele)]="COUNT_HET"
	mafSubForCount$HetOrHomoVar[which(mafSubForCount$Tumor_Seq_Allele2!=mafSubForCount$Reference_Allele & mafSubForCount$Tumor_Seq_Allele1!=mafSubForCount$Reference_Allele)]="COUNT_HOM"
	
	#remove more than one variant in same gene on same sample. More robust
	mafSubForCount=unique(mafSubForCount[,c("Hugo_Symbol","Tumor_Sample_Barcode","HetOrHomoVar")])
	
	geneBurdenCount=mafSubForCount[, .N, , c("Hugo_Symbol","HetOrHomoVar")]
	geneBurdenCountTable=dcast(geneBurdenCount,Hugo_Symbol~HetOrHomoVar,value.var="N")
	if (! ("COUNT_HOM" %in% colnames(geneBurdenCountTable))) {
		geneBurdenCountTable$COUNT_HOM=0
	}
	geneBurdenCountTable[is.na(geneBurdenCountTable)]=0
	
	geneBurdenCountTable$TOTAL_AC=geneBurdenCountTable$COUNT_HET+2*geneBurdenCountTable$COUNT_HOM
	geneBurdenCountTable$TOTAL_AN=totalSampleNum*2
	
	colnames(geneBurdenCountTable)[1]="Gene"
	return(geneBurdenCountTable)
}

#SELECTED_VAR_TYPE="missense_variant"	
#geneBurdenCountTable=makeMafToGeneBurdenCount(maf,SELECTED_VAR_TYPE=SELECTED_VAR_TYPE)
#write.csv(geneBurdenCountTable,paste0(basename(mafFile),".",gsub("\\|",".",SELECTED_VAR_TYPE),".csv"),row.names=FALSE)
#SELECTED_VAR_TYPE="synonymous_variant"	
#geneBurdenCountTable=makeMafToGeneBurdenCount(maf,SELECTED_VAR_TYPE=SELECTED_VAR_TYPE)
#write.csv(geneBurdenCountTable,paste0(basename(mafFile),".",gsub("\\|",".",SELECTED_VAR_TYPE),".csv"),row.names=FALSE)

SELECTED_VAR_TYPE=arguments$snvType	
geneBurdenCountTable=makeMafToGeneBurdenCount(maf,SELECTED_VAR_TYPE=SELECTED_VAR_TYPE,varTypeCol="Variant_Classification")
write.csv(geneBurdenCountTable,paste0(basename(mafFile),".",gsub("\\|",".",SELECTED_VAR_TYPE),".csv"),row.names=FALSE)
