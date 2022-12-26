# Master run for all genes
rm(list=ls(all=TRUE))
# This script will anlyze GAW17 data.
#  gene under consideration
# 
name = 'BCHE' # Can change to other CAUSAL GENES
name = 'PDGFD'
name = 'PLAT'
name = 'GCKR'
name = 'LPL'

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
xx[subSNP[MAF<0.05]]
P = NULL
j=1
for (j in 1:200){
# get Y

file = paste('/Users/davids/Documents/Work/Thesis/GAW17/GAW17 CD/unr_phen.',j,sep='')
Y = getY(ID[subID],file)
p = pvalCalc(Y,Xsd,10^4)
P = cbind(P,p)
cat(p,'\n')
print(j)
}

powerLW = sum(P[1,]<0.01)/200
powerLS = sum(P[2,]<0.01)/200
powerC = sum(P[3,]<0.01)/200
powerH = sum(P[4,]<0.01)/200
powerM = sum(P[5,]<0.01)/200
powerF = sum(P[6,]<0.01)/200

powerLW = sum(P[1,]<0.05)/200  #Linear weighted     ###table
powerLS = sum(P[2,]<0.05)/200  #Linear unweighted
powerC = sum(P[3,]<0.05)/200   #Quadratic unweighted 
powerH = sum(P[4,]<0.05)/200   #Quadratic           ###table
powerM = sum(P[5,]<0.05)/200   #Min-p
powerF = sum(P[6,]<0.05)/200   #Fisher

powerLW
powerLS
powerC
powerH
powerM
powerF





### Demonstrate the winner's curse

for (j in 1:200){
file = paste('/Users/davids/Documents/Work/Thesis/GAW17/GAW17 CD/unr_phen.',j,sep='')
Y = getY(ID[subID],file)
Y
X=Xsd

L = length(X[,1])
J = length(X[1,])
maf = colSums(X)/L
SW[j] = as.numeric(testQTL(Y,X,1/sqrt(maf*(1-maf))))
SL[j] = as.numeric(testQTL(Y,X,rep(1,J)))
SC[j] = as.numeric(testQTC(Y,X))
SQ[j] = as.numeric(testQTH(Y,X))
print(j)
}
SW
SL
SC
SQ

hist(SW)
hist(SL)
hist(SC)
hist(SQ)

mean(SL)
mean(SQ)




 