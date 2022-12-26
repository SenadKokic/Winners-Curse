# Master run for all genes
rm(list=ls(all=TRUE))
# This script will anlyze GAW17 data.
#  gene under consideration
# 
name = 'BCHE' # Can change to other CAUSAL GENES
DD = read.csv('/Users/davids/Documents/Work/Thesis/GAW17/GAW17 CD/gene_info',header=T,sep=',')
ii = which(DD[,1]==name)
start = DD[ii,4]
last  = DD[ii,5]
chr = DD[ii,2]
#
#   Read data from the ch3
name_to_read =  paste('/Users/davids/Documents/Work/Thesis/GAW17/GAW17 CD/c',chr,'_snps.unr',sep='')
A = read.csv(name_to_read,header=T,sep=',')
# helper functions
# helper functions

source('/Users/davids/Documents/Work/Thesis/Winners Curse/R code/gawhelp.r')
source('/Users/davids/Documents/Work/Thesis/Winners Curse/R code/gawpermM.r')
source('/Users/davids/Documents/Work/Thesis/Winners Curse/R code/methodsX.r')



ID = A[,1]
subID = getPop(ID) 

B = read.csv('/Users/davids/Documents/Work/Thesis/GAW17/GAW17 CD/snp_info',header=T,sep=',')
subSNP = subSet(B,names(A),start,last)
X = convertX(A[subID,subSNP])


# Get MAF

MAF = getmaf(X)
SD = MAF<0.05 & MAF>0
Xsd  =(X[,SD]) 
#Xsd  =matrix(X[,SD]) 
dim(Xsd)
xx = names(A)
xx[subSNP[MAF>0]]
P = NULL

for (j in 1:200){
# get Y

file = paste('/Users/davids/Documents/Work/Thesis/GAW17/GAW17 CD/unr_phen.',j,sep='')
Y = getY(ID[subID],file)
p = pvalCalc(Y,Xsd,10^4)
P = cbind(P,p)
cat(p,'\n')
}

powerLW = sum(P[1,]<0.01)/200
powerLS = sum(P[2,]<0.01)/200
powerC = sum(P[3,]<0.01)/200
powerH = sum(P[4,]<0.01)/200
powerM = sum(P[5,]<0.01)/200
powerF = sum(P[6,]<0.01)/200

powerLW
powerLS
powerC
powerH
powerM
powerF



 