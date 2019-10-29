#! /usr/bin/env Rscript
'extract snp counts in gene level for control datasets
Usage:
    control_vcf_vep.R [--vcfSource] [--maxAF] [--snvType] <input> <output>
    
Options:
    -h --help  Show this screen.
    -v --version  Show version.
    --vcfSource=<DB>   contral vcf file type, could be ExAC, gnomAD or else. [default: gnomAD]
    --maxAF=<AF>  AF maximum for filtering. [default: 0.01]
    --snvType=<missense> missense or other "missense_variant|start_lost|stop_lost|protein_altering_variant|stop_gained" [default: missense]
Arguments:
    input  contral vcf file, could be download from ExAC, genomAD or else
    output  output filename
' -> doc


suppressMessages(library(VariantAnnotation))
suppressMessages(library(data.table))
suppressMessages(library(docopt))
arguments <- docopt(doc, version = 'control_table v0.1\n\n')
print(str(arguments))

controlVcfFile=arguments$input
#SELECTED_VAR_TYPE="missense_variant|start_lost|stop_lost|protein_altering_variant|stop_gained"
#SELECTED_VAR_TYPE="synonymous_variant"

SELECTED_VAR_TYPE=arguments$snvType
maxAF=arguments$maxAF


if ( arguments$vcfSource == "ExAC"){
	controlVcfFileType="ExAC"
} else if (arguments$vcfSource == "gnomAD") {
	controlVcfFileType="gnomAD"
} else {
	if (grepl("ExAC",basename(controlVcfFile),ignore.case=T)) {
	controlVcfFileType="ExAC"
} else if (grepl("gnomAD",basename(controlVcfFile),ignore.case=T)) {
	controlVcfFileType="gnomAD"
} else {
	stop(paste0(basename(controlVcfFile)," doesn't match ExAC or gnomAD"))
}
}

AF_max <- function(x) {
	AF <- unlist(info(x)$AF)
	as.vector(AF<=maxAF)
}


filters <- FilterRules(list(AF=AF_max))

updateGeneSnpCountTable<-function(geneSnpCountTable=NULL,selectedSnp,selectedSnpGene,selectedSnpNhomalt,selectedSnpAC,selectedSnpAN) {
	selectedSnpGeneUnlist=unlist(selectedSnpGene) #one variant maybe on more genes, need to be careful here 
	
	selectedSnpAC=unlist(selectedSnpAC)
	selectedSnpAN=unlist(selectedSnpAN)
	selectedSnpNhomalt=unlist(selectedSnpNhomalt)
	
	if (is.null(geneSnpCountTable)) {
		geneSnpCountTable=data.table(Gene=unique(selectedSnpGeneUnlist),Snp="",matrix(0,nrow=length(unique(selectedSnpGeneUnlist)),ncol=4))
		colnames(geneSnpCountTable)[-c(1,2)]=c("COUNT_HET","COUNT_HOM","TOTAL_AC","TOTAL_AN")
		setkey(geneSnpCountTable,Gene)
	} 
	
	selectedSnpGeneNew=setdiff(unique(selectedSnpGeneUnlist),geneSnpCountTable[,Gene])
	if (length(selectedSnpGeneNew)>0) {
		temp=data.table(Gene=unique(selectedSnpGeneNew),Snp="",matrix(0,nrow=length(unique(selectedSnpGeneNew)),ncol=4))
		colnames(temp)=colnames(geneSnpCountTable)
		geneSnpCountTable<-rbind(geneSnpCountTable,temp)
		setkey(geneSnpCountTable,Gene)
	}
	
	
	for (selectedSnpGeneOne in unique(selectedSnpGeneUnlist)) {
#		i=which(selectedSnpGene==selectedSnpGeneOne)
		i=which(sapply(selectedSnpGene,function(x) any(x==selectedSnpGeneOne)))
		
		selectedSnpOne=selectedSnp[i]
		selectedSnpOne=unique(c(selectedSnpOne,strsplit(geneSnpCountTable[selectedSnpGeneOne,Snp],";")[[1]]))
		geneSnpCountTable[selectedSnpGeneOne,"Snp"]=paste(selectedSnpOne,collapse=";")
		
		selectedSnpANOne=max(selectedSnpAN[i]) #using max AN, not sum
		geneSnpCountTable[selectedSnpGeneOne,"TOTAL_AN"]=max(geneSnpCountTable[selectedSnpGeneOne,TOTAL_AN],selectedSnpANOne)
		
		selectedSnpACOne=sum(selectedSnpAC[i])
		selectedSnpNhomaltOne=sum(selectedSnpNhomalt[i])
		selectedSnpNhetaltOne=selectedSnpACOne-2*selectedSnpNhomaltOne
		geneSnpCountTable[selectedSnpGeneOne,c("COUNT_HET","COUNT_HOM","TOTAL_AC")]=geneSnpCountTable[selectedSnpGeneOne,c("COUNT_HET","COUNT_HOM","TOTAL_AC")]+
				c(selectedSnpNhetaltOne,selectedSnpNhomaltOne,selectedSnpACOne)
	}
	
	return(geneSnpCountTable)
}

if (controlVcfFileType=="ExAC") {
	vepCol="CSQ"
	acCol="AC_Adj"
	anCol="AN_Adj"
	nhomaltCol="AC_Hom"
	ScanVcfParamInfo=c(acCol,anCol,"QD",nhomaltCol,"DP_HIST",vepCol) #Don't need AF here as we will use adjusted AF which is AC_Adj/AN_Adj
	ScanVcfParamFixed=c("ALT","FILTER") #ALT is needed to expand the vcf to one alt allele per row
	filters <- FilterRules(list(AF=AF_max_ExAC,QD=QD_min,DP=DP_hist_min_ExAC,PASS=isPassFilter))
	
} else if (controlVcfFileType=="gnomAD") {
	vepCol="vep"
	acCol="AC"
	anCol="AN"
	nhomaltCol="nhomalt"
	ScanVcfParamInfo=c(acCol,anCol,"AF","QD",nhomaltCol,"dp_hist_all_bin_freq",vepCol)
	ScanVcfParamFixed=c("FILTER")
	filters <- FilterRules(list(AF=AF_max))
	
}
readVcfParam <- ScanVcfParam(info=ScanVcfParamInfo, geno=NA,fixed=ScanVcfParamFixed)



geneSnpCountTable=NULL
tab <- TabixFile(controlVcfFile, yieldSize=50000)
open(tab)
vcfProgress=0
while (nrow(vcfYield <- readVcf(tab, "hg19", param=readVcfParam))) {

	vcfProgress=vcfProgress+length(vcfYield)
	print(paste0("Reading vcf: ", vcfProgress))
	
	if (controlVcfFileType=="ExAC") {
		vcfYield <- VariantAnnotation::expand(x = vcfYield, row.names = TRUE)
	}
	vcfChunk <- subsetByFilter(vcfYield, filters)
	
	vcfChunkVepSelected=lapply(info(vcfChunk)[,vepCol],function(x) grep(SELECTED_VAR_TYPE,x,value =TRUE))
	vcfChunkVepSelectedInd=which(lapply(vcfChunkVepSelected,length)>0)
	
	if (length(vcfChunkVepSelectedInd)>0) {
		selectedSnpAC=info(vcfChunk)[,acCol][vcfChunkVepSelectedInd]
		selectedSnpAN=info(vcfChunk)[,anCol][vcfChunkVepSelectedInd]
		selectedSnpNhomalt=info(vcfChunk)[,nhomaltCol][vcfChunkVepSelectedInd]
		
		selectedSnp=names(vcfChunk)[vcfChunkVepSelectedInd]
		selectedSnpGene=lapply(vcfChunkVepSelected[vcfChunkVepSelectedInd],function(x) unique(sapply(strsplit(x,"\\|"),function(y) y[4])))
		
		geneSnpCountTable=updateGeneSnpCountTable(geneSnpCountTable,selectedSnp,selectedSnpGene,selectedSnpNhomalt,selectedSnpAC,selectedSnpAN)
	}
}
close(tab)
#if( arguments$output ){
#	out<-arguments$output
#}else{
#	out<-paste0(controlVcfFile,".",gsub("\\|",".",SELECTED_VAR_TYPE),".csv")
#}
write.csv(geneSnpCountTable,paste0(arguments$output,".",gsub("\\|",".",SELECTED_VAR_TYPE),".csv"),row.names=FALSE)

