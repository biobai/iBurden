#! /usr/bin/env Rscript
'do burden test using control and case table
Usage:
    doBurdenTest.R <control> <case>

Arguments:
    control	contral table
    case	case table
' -> doc


suppressMessages(library(data.table))
suppressMessages(library(docopt))
arguments <- docopt(doc, version = 'doBurdenTest v0.1\n\n')
print(str(arguments))

controlBurdenFile=arguments$control
caseBurdenFile=arguments$case

controlBurden=fread(controlBurdenFile)
setkey(controlBurden,"Gene")
colnames(controlBurden)[-1]=paste0("Control_",colnames(controlBurden)[-1])
caseBurden=fread(caseBurdenFile)
colnames(caseBurden)[1]="Gene"
setkey(caseBurden,"Gene")

allBruden=merge(caseBurden,controlBurden,by="Gene")


burdenTest<-function(allBruden,testCol=c("Gene","TOTAL_AC","TOTAL_AN","Control_TOTAL_AC","Control_TOTAL_AN")) {
	testP=apply(allBruden[,..testCol],1,function(x) fisher.test(matrix(as.integer(x[-1]),ncol=2,nrow=2))$p.value)
	relativeRatio=(allBruden[[testCol[2]]]/allBruden[[testCol[3]]])/(allBruden[[testCol[4]]]/allBruden[[testCol[5]]])
	resultOut=data.table(allBruden[,..testCol],relativeRatio,testP)
	return(resultOut[order(resultOut$testP),])
}

burdenTestResult=burdenTest(allBruden)

write.csv(burdenTestResult,paste0(basename(caseBurdenFile),".GeneBurdenTest.csv"),row.names=FALSE)

