args <- commandArgs(TRUE) 
dataSource <- args[1]
RSDthresh <- as.numeric(args[2])
metaFile <- args[3]
qcLabel <- args[4]
qcTrail <- as.numeric(args[5])
outputName <- args[6]

DF <- read.csv(dataSource,header=T,row.names=1)

suppressMessages(library(XML))

meta = xmlInternalTreeParse(metaFile)

Classes = xpathSApply(meta,"//sampleID",xmlValue)
SampleNames = xpathSApply(meta,"//ID",xmlValue)

if(length(Classes)!=length(SampleNames))
	stop("Classes and Sample Names must match")

if(length(SampleNames)!=nrow(DF))
	stop("Matrix rows and MetaData do not match")

M = match(rownames(DF),SampleNames)
Classes = Classes[M]

if(qcTrail>0){
X2 <- DF[-(1:qcTrail),]
}else{
X2 <- DF
}

X2 <- X2[Classes==qcLabel,]

RSD <- function(x) sd(x,na.rm=T)/mean(x,na.rm=T)

vecRSD <- apply(X2,2,RSD)
X <- DF[,-(which(vecRSD>RSDthresh))]

write.csv(X,outputName)
