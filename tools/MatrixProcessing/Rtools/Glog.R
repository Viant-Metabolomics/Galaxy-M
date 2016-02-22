args <- commandArgs(TRUE)
dataSource <- args[1]
metaFile <- args[2]
QC_Label <- args[3]
outputName <- args[4]

X <- read.csv(dataSource,header=T,row.names=1)

suppressMessages(library(XML))

meta = xmlInternalTreeParse(metaFile)

Classes = xpathSApply(meta,"//sampleID",xmlValue)
SampleNames = xpathSApply(meta,"//ID",xmlValue)

if(length(Classes)!=length(SampleNames))
	stop("Classes and Sample Names must match")

if(length(SampleNames)!=nrow(X))
	stop("Matrix rows and MetaData do not match")

M = match(rownames(X),SampleNames)
Classes = Classes[M]

glog <- function(y,alpha,lambda)
	z <- log((y-alpha)+sqrt((y-alpha)^2+lambda))

jglog <- function(y,y0,lambda)
{
	z <- glog(y,y0,lambda)
	D <- log(sqrt((y-y0)^2+lambda))
	
	gmn <- exp(apply(D,2,mean,na.rm=T))
	zj <- z*gmn
	return(zj)
}
SSE <- function(lambda,alpha,y)
{
	N <- dim(y)[2]
	len <- dim(y)[1]
	
	z <- jglog(y,alpha,lambda)
	s <- 0
	mean_spec <- apply(z,1,mean,na.rm=T)
	
	s <- sum((z-mean_spec)^2,na.rm=T)

	#cat(lambda,"\t",s,"\n")

	return(s)
}


x <- X[Classes==QC_Label,]
x <- t(x)

y0 <- 0
N <- ncol(x)
L <- max(dim(x))

scal_fact <- 1
pow_fact <- 1
offset <- min(x,na.rm=T)
x <- x-offset

step_threshold <- 1e-16

small <- min(x,na.rm=T) #which is 0

if(min(apply(t(x),2,var,na.rm=T))==0)
{
	varbs <- apply(t(x),2,var,na.rm=T)
	newminVar <- sort(unique(varbs))[2]
}else{
	newminVar <- min(apply(t(x),2,var,na.rm=T))
}

low_lim <- -small^2
upper_lim <- max(pmax(apply(t(x),2,var,na.rm=T),max(apply(t(x),2,var,na.rm=T))/newminVar))

lambda <- optimize(SSE,interval=c(low_lim,upper_lim),y0,x,tol=step_threshold)

lambda <- as.numeric(lambda[[1]])

lambda_std <- 5.0278*10^(-09)

error_flag=F

if(abs(upper_lim-lambda)<=1e-5)
{
	cat("Error!Lambda tending to infinity!Using standard\n")
	error_flag=T
} else if(abs(low_lim-lambda)<=1e-5)
{
	cat("Error!Lambda tending to -infinity!Using standard\n")
	error_flag=T
}

x <- X
x <- t(x)	
N <- dim(x)[2]

if(error_flag)
{
	lambda <- lambda_std
	scal_fact <- apply(x,2,sum,na.rm=T)
	scal_fact <- mean(scal_fact)
	scal_fact <- 1/scal_fact
}

x <- x*scal_fact
x <- x^pow_fact
x <- x-min(x,na.rm=T)

Z <- glog(x,y0,lambda)

X <- as.data.frame(t(Z))

write.csv(X,outputName)
