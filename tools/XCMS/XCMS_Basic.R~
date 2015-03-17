suppressMessages(library(Hmisc)) 
suppressMessages(library(xcms)) 
suppressMessages(library(CAMERA)) 
args <- commandArgs(TRUE) 

cdfDir = args[1] 
stepvar=as.numeric(args[2]) 
snthreshvar=as.numeric(args[3]) 
mzdiffvar=as.numeric(args[4]) 
bwvar=as.numeric(args[5]) 
mzwidvar=as.numeric(args[6]) 
cor_eic_var=as.numeric(args[7]) 
outfile = args[8]
outfile_ch = args[9]
outfile_rh = args[10]
outfile_pt = args[11]


setwd(cdfDir) 
xset <- xcmsSet(step=stepvar,snthresh=snthreshvar,mzdiff = mzdiffvar) 
grp <- group(xset,bw=bwvar,mzwid = mzwidvar) 
pc.tmp = "" 
an <- annotate(grp, cor_eic_th=cor_eic_var) 
peaklist<-getPeaklist(an) 
write.table(peaklist, outfile_pt)

mz_vals<-peaklist$mz
write.table(mz_vals, outfile_ch, row.names=FALSE, col.names=FALSE)

peaklist<-subset(peaklist, select=-c(mz,mzmin,mzmax,rt,rtmin,rtmax,npeaks,isotopes, adduct, pcgroup))
peaklist<-peaklist[c(-1)]
peaklist[is.na(peaklist)]<-0

sample_names<-colnames(peaklist)
write.table(sample_names, outfile_rh, row.names=FALSE, col.names=FALSE)

write.table(t(peaklist),outfile, row.names=FALSE, col.names=FALSE, sep=",")  
