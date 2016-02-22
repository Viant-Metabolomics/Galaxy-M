###Peak Filtering FOR GALAXY 0.0.1

args <- commandArgs(TRUE)
dataSource <- args[1]
noiseFl <- as.logical(args[2])
metaFile = args[3]
Feat <- as.numeric(args[4])
Samp <- as.numeric(args[5])
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

PeakFilt <- function(X,noiseFill,cls,FeatureFilt,SampleFilt)
{
	if(noiseFill){
		fac <- factor(cls)
		noise <- min(X,na.rm=T)/2

		for(i in 1:nlevels(fac))
		{
			X2 <- X[cls==i,]	
			respo <- apply(X2,2,function(r) all(is.na(r)))
			peks <- which(respo)
			X[(cls==i),peks] <- noise
		}
	}

	#####FEATURE FILTER

	FUN <- function(irr) return(length(which(is.na(irr)))/length(irr))
	vef <- apply(X,2,FUN)
	vef <- which(vef>=FeatureFilt)
	if(length(vef)>=1)
		X <- X[,-vef]

	#####SAMPLE FILTER

	vef <- apply(X,1,FUN)
	vef <- which(vef>=SampleFilt)
	if(length(vef)>=1)
		X <- X[-vef,]

	return(X)
}

newDtst <- PeakFilt(DF,noiseFl,Classes,Feat,Samp)
write.csv(newDtst,outputName)
