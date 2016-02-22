###MISSING VALUE IMPUTATION


args <- commandArgs(TRUE)
dataSource <- args[1]
MVImode <- args[2]
outputName <- args[3]

X <- read.csv(dataSource,header=T,row.names=1)

if(MVImode=="SV"){
	X[is.na(X)] <- min(X,na.rm=T)/2 ##SMV

}else if(MVImode=="MN"){
	meanrep <- function(mat) apply(mat,2,mean,na.rm=T) ###MEAN REP
	meanVec <- meanrep(X)
	for(i in 1:(ncol(X)))
		X[,i][is.na(X[,i])] <- meanVec[i]	

}else if(MVImode=="MD"){
	medianrep <- function(mat) apply(mat,2,median,na.rm=T)
	medianVec <- medianrep(X)
	for(i in 1:(ncol(X)))
		X[,i][is.na(X[,i])] <- medianVec[i]

}else if(MVImode=="KNN"){
	suppressMessages(library(impute))
	obj <- suppressWarnings(impute.knn(as.matrix(t(X))))
	X <- as.data.frame(t(obj$data))
}else if(MVImode=="BPCA"){
	suppressMessages(library(pcaMethods))
	pcaOb <- pcaMethods::pca(X,method="bpca",scale="none")
	X <- pcaOb@completeObs
	X[X<0] <- min(X [X>0])	##GUARD AGAINST NEGATIVE VALUES

}else {
	stop("Error occurred. MVI was not selected")
}

write.csv(X,outputName)


