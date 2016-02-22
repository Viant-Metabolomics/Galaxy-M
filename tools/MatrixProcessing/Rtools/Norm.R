###NORMALIZATION FOR GALAXY 0.1

args <- commandArgs(TRUE)
dataSource <- args[1]
md <- args[2]
metaFile = args[3]
qcLabel <- args[4]
outputName <- args[5]

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


Norm <- function(dtst,mode="SUM",QC_Label=1)
{	
	if(mode=="SUM"){
		dtst <- sweep(dtst,1,rowSums(dtst,na.rm=T)/100,FUN="/")
		return (dtst)
	}

	else if(mode=="PQN"){	
		ref <- dtst[Classes==QC_Label,]

		if(nrow(ref)<2){
			stop("No QC found with selected label")
		}

		ref <- apply(ref,2,mean,na.rm=T)

		coef <- vector()

		for (i in 1:dim(dtst)[1])
		{
			tempMat <- rbind(ref,dtst[i,])

 			vecelim <- which(apply(tempMat,2,function(x) any(is.na(x))))

			if(length(vecelim)!=0){			
			tempMat <- tempMat[,-vecelim]
			}

			coef[i] <- median(as.numeric(tempMat[2,]/tempMat[1,],na.rm=T))
		}

		dtst <- apply(dtst,2,function(x) x/coef)
		return (dtst)
	}
	
	else{
		print("Error. Invalid normalisation mode")
	}
}

X <- Norm(DF,md,qcLabel)

write.csv(X,outputName)
