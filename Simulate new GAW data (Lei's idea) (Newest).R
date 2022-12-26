# load lattice library for data vis
library(lattice)

#save.image("/Users/davids/Documents/Work/R stuff/Workspaces/WinnersGAW17.RData")
#load("/Users/davids/Documents/Work/R stuff/Workspaces/WinnersGAW17.RData")

# Master run for all genes
#rm(list=ls(all=TRUE))
# This script will anlyze GAW17 data.
#  gene under consideration
#

##### Simulate the new data based on GAW17 parameters
# select SIRT1 gene
name='SIRT1'

print(name)
# read gene information from file
DD = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/gene_info',header=T,sep=',')
# find which indexes match the desired gene name
#   (DD[,1] gives gene names from the csv)
ii = which(DD[,1]==name)
# not quite sure what these mean... seem to be positions?
start = DD[ii,4]
last  = DD[ii,5]
chr = DD[ii,2]

#
#   Read data from the ch3
name_to_read =  paste('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/c',chr,'_snps.unr',sep='')
A = read.csv(name_to_read,header=T,sep=',')
# helper functions
# helper functions

source('/Users/senadkokic/Desktop/NSERC/R code/gawhelp.r')
source('/Users/senadkokic/Desktop/NSERC/R code/gawpermM.r')
source('/Users/senadkokic/Desktop/NSERC/R code/methodsX-winners.r')

# Headers of snp files
ID = A[,1]
# call function that gets the Han Chinese data from the file
subID = getPop(ID) 

B = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/snp_info',header=T,sep=',')
# function that gets subset of snps, based on their positions
subSNP = subSet(B,names(A),start,last)
# converts the call matrix to a numerical matrix
X = convertX(A[subID,subSNP])
# Get MAF

MAF = getmaf(X)
MAF
# generates a matrix of boolean values according to conditions
SD = MAF<0.05 & MAF>0
# gives a lot of 0s and 1s... based on the entries of the SD matrix
Xsd  =(X[,SD])
apply(Xsd,2,sum) 
apply(Xsd,2,sum)/321 
if(name=='GCKR'){
Xsd  =matrix(X[,SD]) 
}
# Indexing where the MAFs of the SNPs match conditions, and printing such 
xx = names(A)
xx[subSNP[MAF>0]]
xx[subSNP[MAF<0.05]]
xx[subSNP[MAF<0.05 & MAF>0]]

P = NULL

# Applies the maximum function to the SD boolean matrix
XCMC=apply(Xsd,1,FUN=max)
# I think this is just where everything is defined?
XCMC[XCMC>0]=XCMC[XCMC>0]^0
# not quite sure why this line produces logical(0)
qj[go]==mean(XCMC)
# set parameters
Betaj=NULL
pvalj=NULL
JJ=2000

qj=NULL
Betaj=NULL
tBetaj=NULL
pvalj=NULL
Bj=NULL
Tj=NULL
Pj=NULL

SW= NULL
SL= NULL
SWwin= NULL
SLwin= NULL
SC= NULL
SQ= NULL
NBetas=11
j=5
IBetas=matrix(rep(NA,JJ*NBetas),ncol= NBetas)
TIBetas=matrix(rep(NA,JJ*NBetas),ncol= NBetas)
PIBetas=matrix(rep(NA,JJ*NBetas),ncol= NBetas)
Phenos=matrix(rep(NA,JJ*321),ncol= JJ)

# This runs about 2000 times by the way
for(j in 1:JJ){
	print(j)
	vec=c(0,0,0.83224,0.97060,0,0,0,0,0,0.93459,0.53073)
	Design=apply(sweep(Xsd,MARGIN=2,vec,`*`),1,sum)
	# add normal simulations to the Design matrix
	Y=Design+rnorm(length(Design),0,1)
	# column entry changed to Y
	Phenos[,j]=Y
	# Linear model 
	summary(lm(Y~Xsd))$coeff[,1]
#	IBetas[j,]=summary(lm(Y~Xsd))$coeff[,1]
#	TIBetas[j,]=summary(lm(Y~Xsd))$coeff[,3]
#	PIBetas[j,]=summary(lm(Y~Xsd))$coeff[,4]
	# Assigns the XCMC of the estimate
	Betaj[j]=summary(lm(Y~XCMC))$coefficients[2,1]
	# XCMC of t value
	tBetaj[j]=summary(lm(Y~XCMC))$coefficients[2,3]
	# XCMC of Pr(>|t|)
	pvalj[j]=summary(lm(Y~XCMC))$coefficients[2,4]
	
	for(i in 1:11){
	  # estimate
		IBetas[j,i]=summary(lm(Y~Xsd[,i]))$coeff[2,1]
		# t value
		TIBetas[j,i]=summary(lm(Y~Xsd[,i]))$coeff[2,3]
		# Pr(>|t|)
		PIBetas[j,i]=summary(lm(Y~Xsd[,i]))$coeff[2,4]
	}
	
	X=Xsd
L = length(X[,1])
J = length(X[1,])
# taking the maf as the # of occurrences divided by the sample size of 321
maf = colSums(X)/L
# run tests to calculate estimates of values to place in matrix
SWwin[j] = as.numeric(testQTLwin(Y,X,1/sqrt(maf*(1-maf))))
SW[j] =    as.numeric(testQTL(Y,X,1/sqrt(maf*(1-maf))))
SL[j] = as.numeric(testQTL(Y,X,rep(1,J)))
SLwin[j] = as.numeric(testQTLwin(Y,X,rep(1,J))) 
SC[j] = as.numeric(testQTC(Y,X))
SQ[j] = as.numeric(testQTH(Y,X))

# p value calculations
	p = pvalCalc(Y,Xsd,10^3)
	P = cbind(P,p)
	cat(p,'\n')
	
}
dim(SIRT1dat)
# combine vectors/matrices by the columns
Bj=cbind(Bj,Betaj)
Tj=cbind(Tj,tBetaj)
Pj=cbind(Pj,pvalj)

# transposes the P matrix 
t(P)
# frame combines the various estimates
frame=cbind(t(P),SW,SL,SWwin,SLwin,SC,SQ,Bj,Tj,Pj,IBetas)
colnames(frame)[c(1,2,3,4,5,6,13:26)]=c('pLW','pLU','pQU','pQH','pMin','pFisher','CMCBeta','TCMC','pBeta','B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11')
#colnames(frame)[c(1,2,3,4,5,6,13:26)]=c('pLW','pLU','pQU','pQH','pMin','pFisher','CMCBeta','TCMC','pBeta','B0','B1','B2','B3','B4','B5','B6','B7','B8','B9','B10')
head(frame)
dim(frame)
write.csv(frame,'/Users/senadkokic/Desktop/NSERC/GAW17/newsimSIRT1')
newSIRT1=read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/newsimSIRT1')
head(newSIRT1)
#colnames(TIBetas)=c('B0','B1','B2','B3','B4','B5','B6','B7','B8','B9','B10')
colnames(IBetas)=c('B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11')
head(IBetas)
colnames(TIBetas)=c('B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11')
head(TIBetas)
colnames(PIBetas)=c('B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11')
head(PIBetas)
PIBetas=data.frame(PIBetas)
TIBetas=data.frame(TIBetas)
# # calculates correlation between TCMC and stats
with(newSIRT1,cor(TCMC,SW))
with(newSIRT1,cor(TCMC,SL))
with(newSIRT1,cor(TCMC,SC))
with(newSIRT1,cor(TCMC,SQ))

 IBetas1=IBetas
 TIBetas1=TIBetas
 PIBetas1=PIBetas
 newSIRT1.1=newSIRT1


# IBetas2=IBetas
# TIBetas2=TIBetas
# PIBetas2=PIBetas
# newSIRT1.2=newSIRT1


##### vec=c(0,0,0.83224,0.97060,0,0,0,0,0,-0.93459,0.53073)
# IBetas3=IBetas
# TIBetas3=TIBetas
# PIBetas3=PIBetas
# newSIRT1.3=newSIRT1


########### Now the analysis etc.

SIRT1dat=newSIRT1.1
#SIRT1dat=newSIRT1.2
#SIRT1dat=newSIRT1.3



sum(SIRT1dat$pQH<0.05)/JJ
sum(SIRT1dat$pQU<0.05)/JJ
sum(SIRT1dat$pLW<0.05)/JJ
sum(SIRT1dat$pBeta<0.05)/JJ
sum(PIBetas$B1<0.05)/JJ
sum(PIBetas$B3<0.05)/JJ
sum(PIBetas$B4<0.05)/JJ
sum(PIBetas$B10<0.05)/JJ
sum(PIBetas$B11<0.05)/JJ



mean(SIRT1dat$CMCBeta) # est of True AGE over the RV (causal and non causal)
mean(SIRT1dat$B1)
mean(SIRT1dat$B3)
mean(SIRT1dat$B4)
mean(SIRT1dat$B10)
mean(SIRT1dat$B11)

sd(SIRT1dat$CMCBeta) # est of True AGE over the RV (causal and non causal)
sd(SIRT1dat$B1)
sd(SIRT1dat$B3)
sd(SIRT1dat$B4)
sd(SIRT1dat$B10)
sd(SIRT1dat$B11)





###################################################################################################
# Obtain the bias for each of the signif replicates
###### Bootstrapping Algorithm (1)
n=length(Y)
J=10000
Bk=matrix(rep(0,n*J),ncol=J)
dim(Bk)
for(k in 1:J){
sam=1:321
Bk[,k]=sample(sam,321,replace=T)
}
#head(Bk)
#############################




############## Bootstrap replication under Wald Beta test
###### SIRT1

which(SIRT1dat$pBeta<0.05&SIRT1dat$pQH<0.05)
which(SIRT1dat$pQH<0.05)
which(SIRT1dat$pBeta<0.05)
SIRT1dat=newSIRT1.1
#SIRT1dat=newSIRT1.2
#SIRT1dat=newSIRT1.3

SIRT1dat$muBkmed=NA
SIRT1dat$muPkmed=NA
SIRT1dat$muBkmea=NA
SIRT1dat$muPkmea=NA


print(name)
DD = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/gene_info',header=T,sep=',')
# see which entry in gene column matches the desired gene (SIRT1)
ii = which(DD[,1]==name)
# get gene start and gene end, as well as chromosome
start = DD[ii,4]
last  = DD[ii,5]
chr = DD[ii,2]
# takes the file name based on the chromosome 
name_to_read =  paste('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/c',chr,'_snps.unr',sep='')
# read the file containing info about the SNPs associated with the chromosome
A = read.csv(name_to_read,header=T,sep=',')
ID = A[,1]
subID = getPop(ID) 
subID
B = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/snp_info',header=T,sep=',')
subSNP = subSet(B,names(A),start,last)
X = convertX(A[subID,subSNP])


MAF = getmaf(X)
SD = MAF<0.05 & MAF>0
Xsd  =(X[,SD])
if(name=='GCKR'){
Xsd  =matrix(X[,SD]) 
}
apply(Xsd,2,FUN=sum) 

#MAF = getmaf(X[BK,])
#SD = MAF<0.05 & MAF>0
#Xsd  =(X[,SD])
#apply(Xsd[BK,],2,FUN=sum) 

XCMC=apply(Xsd,1,FUN=max)

dim(Xsd)
xx = names(A)
#xx[subSNP[MAF>0]]
#xx[subSNP[MAF<0.05]]

P = NULL

c=4
for(c in which(SIRT1dat$pBeta<0.05)) {

print(paste(c,name))

BetaBk=NA
pBk=NA
BetaCk=NA
pCk=NA
pH=NULL
pt=NULL
sigK=0



for(k in 1:J){
	print(paste(k,"-",sigK))
	if(sigK==200) break
	BK=Bk[,k]
	CK=sam[-BK]	
	print(apply(Xsd[BK,],2,FUN=sum))
#	file = paste('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/unr_phen.',c,sep='')
#	Y = getY(ID[subID],file)
	Y = Phenos[,c]

#	SQ = as.numeric(testQTHboot(Y[BK],Xsd[BK,],Xsd))
#	perm=10^3
#	storeQH = rep(0,perm)
#	for (i in 1:perm){
#		Ynew = sample(Y[BK])
#		storeQH[i] = as.numeric(testQTHboot(Ynew,Xsd[BK,],Xsd))
#	}
#	pH[k] = sum(abs(storeQH)>=abs(SQ))/perm
#	pH[k]
	
	
    

	pt[k]=summary(lm(Y[BK]~XCMC[BK]))$coefficients[2,4]
	pt[k]
	if(pt[k]<0.05){
		
		sigK=sigK+1
		
		pBk[k]=mean(XCMC[BK])
		if(pBk[k]>0){
		BetaBk[k]=summary(lm(Y[BK]~XCMC[BK]))$coefficients[2,1]
		}
		pCk[k]=mean(XCMC[CK])
		if(pCk[k]>0){
		BetaCk[k]=summary(lm(Y[CK]~XCMC[CK]))$coefficients[2,1]
		}
	}
	
}

#cbind(pH,pt,BetaBk,pBk,BetaCk,pCk)

muBkmed=median(BetaBk-BetaCk,na.rm = T)
muPkmed=median(pBk-pCk,na.rm = T)
SIRT1dat$muBkmed[c]=muBkmed
SIRT1dat$muPkmed[c]=muPkmed

muBkmea=mean(BetaBk-BetaCk,na.rm = T)
muPkmea=mean(pBk-pCk,na.rm = T)
SIRT1dat$muBkmea[c]=muBkmea
SIRT1dat$muPkmea[c]=muPkmea

}


# write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot',name,sep=''))
# bootSIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot',name,sep=''))

# write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot1',name,sep=''))
# boot1SIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot1',name,sep=''))

# write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot2',name,sep=''))
# boot2SIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot2',name,sep=''))

head(boot1SIRT1,100)


###Apply Ghosh et. al. to the significant replicates
SIRT1dat$BT1=NA
SIRT1dat$BT2=NA
SIRT1dat$BT3=NA

source("/Users/senadkokic/Desktop/NSERC/R code/GW-functions.R")# YOU READ IN THE FUNCTIONS FOR OUR METHOD


cee=abs(qnorm(.5*0.05)) # Bonferroni threshold for achieving study-wide significance = 0.1
                          # we use "cee" so R does not get confused with the function 'c'
cc=1
for(cc in which(SIRT1dat$pBeta<0.05)) {
print(cc)

betahat=SIRT1dat[cc,]$CMCBeta # Reported OR = 1.54 
z=sign(betahat)*abs(qnorm(0.5*SIRT1dat[cc,]$pBeta)) # Reported p-value = 5.7e-4, which we convert to a z-value


###################################################
#              THE PROPOSED APPROACH              #
###################################################

se=betahat/z # standard error of betahat

mutilde1=optimize(f=conditional.like,c(-20,20),maximum=T,z=z,cee=cee)$maximum # the conditional mle
betatilde1=mutilde1*se # the conditional mle on beta scale
SIRT1dat[cc,]$BT1=betatilde1# the condition mle OR estimate
SIRT1dat[cc,]$BT1 # print out OR_TILDE_1

a=sum(apply(as.matrix(seq(-50,50,by=0.01)),1,conditional.like.z,cee=cee,z=z))*0.01 # numeric integration
b=sum(apply(as.matrix(seq(-50,50,by=0.01)),1,conditional.like,cee=cee,z=z))*0.01
betatilde2=(a/b)*se
SIRT1dat[cc,]$BT2=betatilde2
SIRT1dat[cc,]$BT2 # print out OR_TILDE_2

betatilde3=(betatilde1+betatilde2)/2
SIRT1dat[cc,]$BT3=betatilde3
SIRT1dat[cc,]$BT3 # OR_TILDE_3

cimu=getCI(z=z,cee=cee,alpha=0.05)
cibeta=c(cimu$lowermu*(betahat/z),cimu$uppermu*(betahat/z))
cibeta # BIAS-CORRECTED 95% CI FOR OR
}

write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot',name,sep=''))
bootSIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot',name,sep=''))

write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot1',name,sep=''))
boot1SIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot1',name,sep=''))
head(boot1SIRT1,100)

write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot2',name,sep=''))
boot2SIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot2',name,sep=''))


##### Histogram comparisons:

SIRT1dat=bootSIRT1



IBetas=IBetas1
TIBetas=TIBetas1
PIBetas=PIBetas1


head(SIRT1dat)
tail(SIRT1dat)

hist(SIRT1dat$CMCBeta,xlab="CMC Beta",main=NULL,freq=T,breaks=10)
hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=2,breaks=10)
hist(with(SIRT1dat,CMCBeta-muBkmed),xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=3,breaks=10)
abline(v=mean(SIRT1dat$CMCBeta),lwd=4,col=1)
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta),lwd=4,col=2)
abline(v=mean(with(SIRT1dat,CMCBeta-muBkmed),na.rm=T),lwd=4,col=3)
title('SIRT1')
title('SIRT1 - Bootstrap correction')

length(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta);length(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta)/JJ
mean(SIRT1dat$CMCBeta);sd(SIRT1dat$CMCBeta)
mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta);sd(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta)
mean(with(SIRT1dat,CMCBeta-muBkmed),na.rm=T);sd(with(SIRT1dat,CMCBeta-muBkmed),na.rm=T)

hist(SIRT1dat$CMCBeta,xlab="CMC Beta",main=NULL,freq=T,breaks=10)
hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=2,breaks=10)
hist(with(SIRT1dat,BT1),xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=3,breaks=10)
abline(v=mean(SIRT1dat$CMCBeta),lwd=4,col=1)
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta),lwd=4,col=2)
abline(v=mean(with(SIRT1dat,BT1),na.rm=T),lwd=4,col=3)
title('SIRT1')
title('SIRT1 - Likelihood Adjusted')
title('SIRT1 - Likelihood Correction 1')



length(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta);length(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta)/JJ
mean(SIRT1dat$CMCBeta);sd(SIRT1dat$CMCBeta)
mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta);sd(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta)
mean(with(SIRT1dat,BT1),na.rm=T);sd(with(SIRT1dat, BT1),na.rm=T)

hist(SIRT1dat$CMCBeta,xlab="CMC Beta",main=NULL,freq=T,breaks=10)
hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=2,breaks=10)
hist(with(SIRT1dat,BT2),xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=3,breaks=10)
abline(v=mean(SIRT1dat$CMCBeta),lwd=4,col=1)
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta),lwd=4,col=2)
abline(v=mean(with(SIRT1dat,BT2),na.rm=T),lwd=4,col=3)
title('SIRT1 - Likelihood Correction 2')

length(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta)
mean(SIRT1dat$CMCBeta);sd(SIRT1dat$CMCBeta)
mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta);sd(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta)
mean(with(SIRT1dat,BT2),na.rm=T);sd(with(SIRT1dat,BT2),na.rm=T)

hist(SIRT1dat$CMCBeta,xlab="CMC Beta",main=NULL,freq=T,breaks=10)
hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=2,breaks=5)
hist(with(SIRT1dat,BT3),xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=3,breaks=10)
abline(v=mean(SIRT1dat$CMCBeta),lwd=4,col=1)
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta),lwd=4,col=2)
abline(v=mean(with(SIRT1dat,BT3),na.rm=T),lwd=4,col=3)
title('SIRT1 - Likelihood Correction 3')

length(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta)
mean(SIRT1dat$CMCBeta);sd(SIRT1dat$CMCBeta)
mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta);sd(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta)
mean(with(SIRT1dat,BT3),na.rm=T);sd(with(SIRT1dat,BT3),na.rm=T)


#Individual Betas (3,4,10,11 causal)

hist(SIRT1dat$B3,xlab="CMC Beta",main=NULL,freq=T,breaks=10)

p1=hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$B3,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=2,breaks=10)

hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$B3,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$B4,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=3,breaks=10)
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B3),lwd=4,col=2)
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B4),lwd=4,col=3)

p1=hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$B3,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
p2=hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$B4,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
p3=hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$B10,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
p4=hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$B11,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=5)

plot( p1, col=rgb(0,0,1,1/4), xlim=c(-2,5),ylim=c(0,200))  
plot( p2, col=rgb(1,0,0,1/4), add=T)        # The variability in the corrected estimates is dependant of the MAF -- smaller MAF = larger varaince in adjusted estimate
plot( p3, col=rgb(1,1,0,1/4), add=T)  
plot( p4, col=rgb(1,1,1,1/4), add=T)

vec1=c(0,0,0.83224,0.97060,0,0,0,0,0,0.93459,0.53073)

p1=hist(SIRT1dat$B1,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
p2=hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$B1,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
plot(p1)  
plot( p2, col=2, add=T) 
abline(v=mean(SIRT1dat$B1),lwd=4)
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B1),lwd=4,col=2)

par(mfrow=c(1,1))

p1=hist(SIRT1dat$B3,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
p2=hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$B3,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10,add=T)
plot(p1)  
plot( p2, col=2, add=T) 
abline(v=mean(SIRT1dat$B3),lwd=4)
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B3),lwd=4,col=2)

p1=hist(SIRT1dat$B4,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
p2=hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$B4,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
plot(p1)  
plot( p2, col=2, add=T) 
abline(v=mean(SIRT1dat$B4),lwd=4)
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B4),lwd=4,col=2)

p1=hist(SIRT1dat$B10,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
p2=hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$B10,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
plot(p1)  
plot( p2, col=2, add=T) 
abline(v=mean(SIRT1dat$B10),lwd=4)
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B10),lwd=4,col=2)

p1=hist(SIRT1dat$B11,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
p2=hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$B11,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
plot(p1)  
plot( p2, col=2, add=T) 
abline(v=mean(SIRT1dat$B11),lwd=4)
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B11),lwd=4,col=2)


abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B3),lwd=4,col=rgb(0,0,1,1/4))
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B4),lwd=4,col=rgb(1,0,0,1/4))
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B10),lwd=4,col=rgb(1,0,0,1/4))
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B11),lwd=4,col=rgb(1,0,1,1/4))
vec=vec1
##Truth
vec[c(3,4,10,11)]
##Truth over replicates (JJ=2000)
mean(SIRT1dat$B3);mean(SIRT1dat$B4);mean(SIRT1dat$B10);mean(SIRT1dat$B11);
##Bias-Winner's Curse
mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B3)
mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B4)
mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B10)
mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B11)
##Absolute Bias
(mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B3)-vec[3])
(mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B4)-vec[4])
(mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B10)-vec[10])
(mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B11)-vec[11])
##Relative Bias
(mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B3)-vec[3])/vec[3]
(mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B4)-vec[4])/vec[4]
(mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B10)-vec[10])/vec[10]
(mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B11)-vec[11])/vec[11]

vec=c(0,0,0.83224,0.97060,0,0,0,0,0,0.93459,0.53073)
apply(Xsd,2,FUN=sum)[c(3,4,10,11)]/321


hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=2,breaks=5)
hist(with(SIRT1dat,CMCBeta-muBkmed),xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=3,breaks=10)
abline(v=mean(SIRT1dat$CMCBeta),lwd=4,col=1)
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta),lwd=4,col=2)
abline(v=mean(with(SIRT1dat,CMCBeta-muBkmed),na.rm=T),lwd=4,col=3)
title('SIRT1')



##############################################################################################

SIRT1dat=boot1SIRT1

IBetas=IBetas2
TIBetas=TIBetas2
PIBetas=PIBetas2

par(mfrow = c(1, 1))
hist(SIRT1dat$CMCBeta,xlab="CMC Beta",main=NULL,freq=T,breaks=10)
hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=2,breaks=10)
hist(with(SIRT1dat,CMCBeta-muBkmed),xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=3,breaks=10)
abline(v=mean(SIRT1dat$CMCBeta),lwd=4,col=1)
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta),lwd=4,col=2)
abline(v=mean(with(SIRT1dat,CMCBeta-muBkmed),na.rm=T),lwd=4,col=3)
title('SIRT1')

length(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta)
mean(SIRT1dat$CMCBeta);sd(SIRT1dat$CMCBeta)
mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta);sd(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta)
mean(with(SIRT1dat,CMCBeta-muBkmed),na.rm=T);sd(with(SIRT1dat,CMCBeta-muBkmed),na.rm=T)

hist(SIRT1dat$CMCBeta,xlab="CMC Beta",main=NULL,freq=T,breaks=10)
hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=2,breaks=10)
hist(with(SIRT1dat,BT1),xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=3,breaks=10)
abline(v=mean(SIRT1dat$CMCBeta),lwd=4,col=1)
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta),lwd=4,col=2)
abline(v=mean(with(SIRT1dat,BT1),na.rm=T),lwd=4,col=3)
title('SIRT1')
title('SIRT1 - Likelihood Adjusted')
title('SIRT1 - Likelihood Correction 1')


length(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta)
mean(SIRT1dat$CMCBeta);sd(SIRT1dat$CMCBeta)
mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta);sd(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta)
mean(with(SIRT1dat,BT1),na.rm=T);sd(with(SIRT1dat, BT1),na.rm=T)





##############################################################################################



SIRT1dat=boot2SIRT1



IBetas=IBetas3
TIBetas=TIBetas3
PIBetas=PIBetas3


head(SIRT1dat,20)

hist(SIRT1dat$CMCBeta,xlab="CMC Beta",main=NULL,freq=T,breaks=10)
hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=2,breaks=10)
hist(with(SIRT1dat,CMCBeta-muBkmed),xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=3,breaks=5)
abline(v=mean(SIRT1dat$CMCBeta),lwd=4,col=1)
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta),lwd=4,col=2)
abline(v=mean(with(SIRT1dat,CMCBeta-muBkmed),na.rm=T),lwd=4,col=3)
title('SIRT1')
title('SIRT1 - Bootstrap correction')

length(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta);length(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta)/JJ
mean(SIRT1dat$CMCBeta);sd(SIRT1dat$CMCBeta)
mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta);sd(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta)
mean(with(SIRT1dat,CMCBeta-muBkmed),na.rm=T);sd(with(SIRT1dat,CMCBeta-muBkmed),na.rm=T)

hist(SIRT1dat$CMCBeta,xlab="CMC Beta",main=NULL,freq=T,breaks=10)
hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=2,breaks=10)
hist(with(SIRT1dat,BT1),xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=3,breaks=5)
abline(v=mean(SIRT1dat$CMCBeta),lwd=4,col=1)
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta),lwd=4,col=2)
abline(v=mean(with(SIRT1dat,BT1),na.rm=T),lwd=4,col=3)
title('SIRT1')
title('SIRT1 - Likelihood Adjusted')
title('SIRT1 - Likelihood Correction 1')



length(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta)
mean(SIRT1dat$CMCBeta);sd(SIRT1dat$CMCBeta)
mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta);sd(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta)
mean(with(SIRT1dat,BT1),na.rm=T);sd(with(SIRT1dat, BT1),na.rm=T)

hist(SIRT1dat$CMCBeta,xlab="CMC Beta",main=NULL,freq=T,breaks=10)
hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=2,breaks=10)
hist(with(SIRT1dat,BT2),xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=3,breaks=10)
abline(v=mean(SIRT1dat$CMCBeta),lwd=4,col=1)
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta),lwd=4,col=2)
abline(v=mean(with(SIRT1dat,BT2),na.rm=T),lwd=4,col=3)
title('SIRT1 - Likelihood Correction 2')

length(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta)
mean(SIRT1dat$CMCBeta);sd(SIRT1dat$CMCBeta)
mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta);sd(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta)
mean(with(SIRT1dat,BT2),na.rm=T);sd(with(SIRT1dat,BT2),na.rm=T)

hist(SIRT1dat$CMCBeta,xlab="CMC Beta",main=NULL,freq=T,breaks=10)
hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=2,breaks=10)
hist(with(SIRT1dat,BT3),xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=3,breaks=10)
abline(v=mean(SIRT1dat$CMCBeta),lwd=4,col=1)
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta),lwd=4,col=2)
abline(v=mean(with(SIRT1dat,BT3),na.rm=T),lwd=4,col=3)
title('SIRT1 - Likelihood Correction 3')

length(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta)
mean(SIRT1dat$CMCBeta);sd(SIRT1dat$CMCBeta)
mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta);sd(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta)
mean(with(SIRT1dat,BT3),na.rm=T);sd(with(SIRT1dat,BT3),na.rm=T)


#Individual Betas (3,4,10,11 causal)

hist(SIRT1dat$B3,xlab="CMC Beta",main=NULL,freq=T,breaks=10)

p1=hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$B3,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=2,breaks=10)

hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$B3,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$B4,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=3,breaks=10)
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B3),lwd=4,col=2)
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B4),lwd=4,col=3)

p1=hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$B3,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
p2=hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$B4,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
p3=hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$B10,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
p4=hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$B11,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)

plot( p1, col=rgb(0,0,1,1/4), xlim=c(-3,5),ylim=c(0,80))  
plot( p2, col=rgb(1,0,0,1/4), add=T)        # The variability in the corrected estimates is dependant of the MAF -- smaller MAF = larger varaince in adjusted estimate
plot( p3, col=rgb(1,1,0,1/4), add=T)  
plot( p4, col=rgb(1,1,1,1/4), add=T)

  

vec1=c(0,0,0.83224,0.97060,0,0,0,0,0,-0.93459,0.53073)

p1=hist(SIRT1dat$B1,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
p2=hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$B1,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
plot(p1)  
plot( p2, col=2, add=T) 
abline(v=mean(SIRT1dat$B1),lwd=4)
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B1),lwd=4,col=2)

par(mfrow=c(1,1))

p1=hist(SIRT1dat$B3,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
p2=hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$B3,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10,add=T)
plot(p1)  
plot( p2, col=2, add=T) 
abline(v=mean(SIRT1dat$B3),lwd=4)
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B3),lwd=4,col=2)

p1=hist(SIRT1dat$B4,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
p2=hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$B4,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=5)
plot(p1)  
plot( p2, col=2, add=T) 
abline(v=mean(SIRT1dat$B4),lwd=4)
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B4),lwd=4,col=2)

p1=hist(SIRT1dat$B10,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
p2=hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$B10,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=5)
plot(p1)  
plot( p2, col=2, add=T) 
abline(v=mean(SIRT1dat$B10),lwd=4)
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B10),lwd=4,col=2)

p1=hist(SIRT1dat$B11,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
p2=hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$B11,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
plot(p1)  
plot( p2, col=2, add=T) 
abline(v=mean(SIRT1dat$B11),lwd=4)
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B11),lwd=4,col=2)


abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B3),lwd=4,col=rgb(0,0,1,1/4))
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B4),lwd=4,col=rgb(1,0,0,1/4))
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B10),lwd=4,col=rgb(1,0,0,1/4))
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B11),lwd=4,col=rgb(1,0,1,1/4))
vec=vec1
##Truth
vec[c(3,4,10,11)]
##Truth over replicates (JJ=2000)
mean(SIRT1dat$B3);mean(SIRT1dat$B4);mean(SIRT1dat$B10);mean(SIRT1dat$B11);
##Bias-Winner's Curse
mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B3)
mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B4)
mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B10)
mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B11)
##Absolute Bias
(mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B3)-vec[3])
(mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B4)-vec[4])
(mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B10)-vec[10])
(mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B11)-vec[11])
##Relative Bias
(mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B3)-vec[3])/vec[3]
(mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B4)-vec[4])/vec[4]
(mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B10)-vec[10])/vec[10]
(mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B11)-vec[11])/vec[11]

apply(Xsd,2,FUN=sum)[c(3,4,10,11)]/321


hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=2,breaks=5)
hist(with(SIRT1dat,CMCBeta-muBkmed),xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=3,breaks=10)
abline(v=mean(SIRT1dat$CMCBeta),lwd=4,col=1)
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta),lwd=4,col=2)
abline(v=mean(with(SIRT1dat,CMCBeta-muBkmed),na.rm=T),lwd=4,col=3)
title('SIRT1')



##############################################################################################






##### HistogrpQHam comparisons:

SIRT1dat=bootSIRT1
vec1=c(0,0,0.83224,0.97060,0,0,0,0,0,0.93459,0.53073)
IBetas=IBetas1
TIBetas=TIBetas1
PIBetas=PIBetas1

length(SIRT1dat[SIRT1dat$pQH<0.05,]$CMCBeta);length(SIRT1dat[SIRT1dat$pQH<0.05,]$CMCBeta)/JJ
length(SIRT1dat[SIRT1dat$pQU<0.05,]$CMCBeta);length(SIRT1dat[SIRT1dat$pQU<0.05,]$CMCBeta)/JJ
length(SIRT1dat[SIRT1dat$pLU<0.05,]$CMCBeta);length(SIRT1dat[SIRT1dat$pLU<0.05,]$CMCBeta)/JJ


SIRT1dat=boot2SIRT1
vec1=c(0,0,0.83224,0.97060,0,0,0,0,0,-0.93459,0.53073)
IBetas=IBetas3
TIBetas=TIBetas3
PIBetas=PIBetas3


head(SIRT1dat)
tail(SIRT1dat)

hist(SIRT1dat$CMCBeta,xlab="CMC Beta",main=NULL,freq=T,breaks=10)
hist(SIRT1dat[SIRT1dat$pQH<0.05,]$CMCBeta,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=2,breaks=10)
hist(with(SIRT1dat,CMCBeta-muBkmed),xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=3,breaks=10)
abline(v=mean(SIRT1dat$CMCBeta),lwd=4,col=1)
abline(v=mean(SIRT1dat[SIRT1dat$pQH<0.05,]$CMCBeta),lwd=4,col=2)
abline(v=mean(with(SIRT1dat,CMCBeta-muBkmed),na.rm=T),lwd=4,col=3)
title('SIRT1')
title('SIRT1 - Bootstrap correction')

length(SIRT1dat[SIRT1dat$pQH<0.05,]$CMCBeta);length(SIRT1dat[SIRT1dat$pQH<0.05,]$CMCBeta)/JJ
mean(SIRT1dat$CMCBeta);sd(SIRT1dat$CMCBeta)
mean(SIRT1dat[SIRT1dat$pQH<0.05,]$CMCBeta);sd(SIRT1dat[SIRT1dat$pQH<0.05,]$CMCBeta)
mean(with(SIRT1dat,CMCBeta-muBkmed),na.rm=T);sd(with(SIRT1dat,CMCBeta-muBkmed),na.rm=T)


hist(SIRT1dat$CMCBeta,xlab="CMC Beta",main=NULL,freq=T,breaks=10)
hist(SIRT1dat[SIRT1dat$pQH<0.05,]$CMCBeta,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=2,breaks=10)
hist(with(SIRT1dat,BT1),xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=3,breaks=10)
abline(v=mean(SIRT1dat$CMCBeta),lwd=4,col=1)
abline(v=mean(SIRT1dat[SIRT1dat$pQH<0.05,]$CMCBeta),lwd=4,col=2)
abline(v=mean(with(SIRT1dat,BT1),na.rm=T),lwd=4,col=3)
title('SIRT1')
title('SIRT1 - Likelihood Adjusted')
title('SIRT1 - Likelihood Correction 1')



length(SIRT1dat[SIRT1dat$pQH<0.05,]$CMCBeta);length(SIRT1dat[SIRT1dat$pQH<0.05,]$CMCBeta)/JJ
mean(SIRT1dat$CMCBeta);sd(SIRT1dat$CMCBeta)
mean(SIRT1dat[SIRT1dat$pQH<0.05,]$CMCBeta);sd(SIRT1dat[SIRT1dat$pQH<0.05,]$CMCBeta)
mean(with(SIRT1dat,BT1),na.rm=T);sd(with(SIRT1dat, BT1),na.rm=T)

hist(SIRT1dat$CMCBeta,xlab="CMC Beta",main=NULL,freq=T,breaks=10)
hist(SIRT1dat[SIRT1dat$pQH<0.05,]$CMCBeta,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=2,breaks=10)
hist(with(SIRT1dat,BT2),xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=3,breaks=10)
abline(v=mean(SIRT1dat$CMCBeta),lwd=4,col=1)
abline(v=mean(SIRT1dat[SIRT1dat$pQH<0.05,]$CMCBeta),lwd=4,col=2)
abline(v=mean(with(SIRT1dat,BT2),na.rm=T),lwd=4,col=3)
title('SIRT1 - Likelihood Correction 2')

length(SIRT1dat[SIRT1dat$pQH<0.05,]$CMCBeta)
mean(SIRT1dat$CMCBeta);sd(SIRT1dat$CMCBeta)
mean(SIRT1dat[SIRT1dat$pQH<0.05,]$CMCBeta);sd(SIRT1dat[SIRT1dat$pQH<0.05,]$CMCBeta)
mean(with(SIRT1dat,BT2),na.rm=T);sd(with(SIRT1dat,BT2),na.rm=T)

hist(SIRT1dat$CMCBeta,xlab="CMC Beta",main=NULL,freq=T,breaks=10)
hist(SIRT1dat[SIRT1dat$pQH<0.05,]$CMCBeta,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=2,breaks=5)
hist(with(SIRT1dat,BT3),xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=3,breaks=10)
abline(v=mean(SIRT1dat$CMCBeta),lwd=4,col=1)
abline(v=mean(SIRT1dat[SIRT1dat$pQH<0.05,]$CMCBeta),lwd=4,col=2)
abline(v=mean(with(SIRT1dat,BT3),na.rm=T),lwd=4,col=3)
title('SIRT1 - Likelihood Correction 3')

length(SIRT1dat[SIRT1dat$pQH<0.05,]$CMCBeta)
mean(SIRT1dat$CMCBeta);sd(SIRT1dat$CMCBeta)
mean(SIRT1dat[SIRT1dat$pQH<0.05,]$CMCBeta);sd(SIRT1dat[SIRT1dat$pQH<0.05,]$CMCBeta)
mean(with(SIRT1dat,BT3),na.rm=T);sd(with(SIRT1dat,BT3),na.rm=T)


#Individual Betas (3,4,10,11 causal)

hist(SIRT1dat$B3,xlab="CMC Beta",main=NULL,freq=T,breaks=10)

p1=hist(SIRT1dat[SIRT1dat$pQH<0.05,]$B3,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=2,breaks=10)

hist(SIRT1dat[SIRT1dat$pQH<0.05,]$B3,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
hist(SIRT1dat[SIRT1dat$pQH<0.05,]$B4,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=3,breaks=10)
abline(v=mean(SIRT1dat[SIRT1dat$pQH<0.05,]$B3),lwd=4,col=2)
abline(v=mean(SIRT1dat[SIRT1dat$pQH<0.05,]$B4),lwd=4,col=3)

p1=hist(SIRT1dat[SIRT1dat$pQH<0.05,]$B3,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
p2=hist(SIRT1dat[SIRT1dat$pQH<0.05,]$B4,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
p3=hist(SIRT1dat[SIRT1dat$pQH<0.05,]$B10,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
p4=hist(SIRT1dat[SIRT1dat$pQH<0.05,]$B11,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=5)

plot( p1, col=rgb(0,0,1,1/4), xlim=c(-2,5),ylim=c(0,200))  
plot( p2, col=rgb(1,0,0,1/4), add=T)        # The variability in the corrected estimates is dependant of the MAF -- smaller MAF = larger varaince in adjusted estimate
plot( p3, col=rgb(1,1,0,1/4), add=T)  
plot( p4, col=rgb(1,1,1,1/4), add=T)


p1=hist(SIRT1dat$B1,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
p2=hist(SIRT1dat[SIRT1dat$pQH<0.05,]$B1,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
plot(p1)  
plot( p2, col=2, add=T) 
abline(v=mean(SIRT1dat$B1),lwd=4)
abline(v=mean(SIRT1dat[SIRT1dat$pQH<0.05,]$B1),lwd=4,col=2)

par(mfrow=c(1,1))

p1=hist(SIRT1dat$B3,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
p2=hist(SIRT1dat[SIRT1dat$pQH<0.05,]$B3,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10,add=T)
plot(p1)  
plot( p2, col=2, add=T) 
abline(v=mean(SIRT1dat$B3),lwd=4)
abline(v=mean(SIRT1dat[SIRT1dat$pQH<0.05,]$B3),lwd=4,col=2)

p1=hist(SIRT1dat$B4,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
p2=hist(SIRT1dat[SIRT1dat$pQH<0.05,]$B4,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
plot(p1)  
plot( p2, col=2, add=T) 
abline(v=mean(SIRT1dat$B4),lwd=4)
abline(v=mean(SIRT1dat[SIRT1dat$pQH<0.05,]$B4),lwd=4,col=2)

p1=hist(SIRT1dat$B10,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
p2=hist(SIRT1dat[SIRT1dat$pQH<0.05,]$B10,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
plot(p1)  
plot( p2, col=2, add=T) 
abline(v=mean(SIRT1dat$B10),lwd=4)
abline(v=mean(SIRT1dat[SIRT1dat$pQH<0.05,]$B10),lwd=4,col=2)

p1=hist(SIRT1dat$B11,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
p2=hist(SIRT1dat[SIRT1dat$pQH<0.05,]$B11,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
plot(p1)  
plot( p2, col=2, add=T) 
abline(v=mean(SIRT1dat$B11),lwd=4)
abline(v=mean(SIRT1dat[SIRT1dat$pQH<0.05,]$B11),lwd=4,col=2)


abline(v=mean(SIRT1dat[SIRT1dat$pQH<0.05,]$B3),lwd=4,col=rgb(0,0,1,1/4))
abline(v=mean(SIRT1dat[SIRT1dat$pQH<0.05,]$B4),lwd=4,col=rgb(1,0,0,1/4))
abline(v=mean(SIRT1dat[SIRT1dat$pQH<0.05,]$B10),lwd=4,col=rgb(1,0,0,1/4))
abline(v=mean(SIRT1dat[SIRT1dat$pQH<0.05,]$B11),lwd=4,col=rgb(1,0,1,1/4))

pQU
##Truth
vec=vec1
vec[c(3,4,10,11)]
mean(vec[c(3,4,10,11)])
##Truth over replicates (JJ=2000)
mean(SIRT1dat$B3);mean(SIRT1dat$B4);mean(SIRT1dat$B10);mean(SIRT1dat$B11);
##Bias-Winner's Curse
mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B3)
mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B4)
mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B10)
mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B11)
##Absolute Bias
(mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B3)-vec[3])
(mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B4)-vec[4])
(mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B10)-vec[10])
(mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B11)-vec[11])
##Relative Bias
(mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B3)-vec[3])/vec[3]
(mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B4)-vec[4])/vec[4]
(mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B10)-vec[10])/vec[10]
(mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B11)-vec[11])/vec[11]

vec=c(0,0,0.83224,0.97060,0,0,0,0,0,0.93459,0.53073)
apply(Xsd,2,FUN=sum)[c(3,4,10,11)]/321

pQU








##### HistogrpLUam comparisons:
head(SIRT1dat)
SIRT1dat=bootSIRT1
vec1=c(0,0,0.83224,0.97060,0,0,0,0,0,0.93459,0.53073)
IBetas=IBetas1
TIBetas=TIBetas1
PIBetas=PIBetas1


SIRT1dat=boot2SIRT1
vec1=c(0,0,0.83224,0.97060,0,0,0,0,0,-0.93459,0.53073)
IBetas=IBetas3
TIBetas=TIBetas3
PIBetas=PIBetas3


head(SIRT1dat)
tail(SIRT1dat)

hist(SIRT1dat$CMCBeta,xlab="CMC Beta",main=NULL,freq=T,breaks=10)
hist(SIRT1dat[SIRT1dat$pLU<0.05,]$CMCBeta,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=2,breaks=5)
hist(with(SIRT1dat,CMCBeta-muBkmed),xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=3,breaks=10)
abline(v=mean(SIRT1dat$CMCBeta),lwd=4,col=1)
abline(v=mean(SIRT1dat[SIRT1dat$pLU<0.05,]$CMCBeta),lwd=4,col=2)
abline(v=mean(with(SIRT1dat,CMCBeta-muBkmed),na.rm=T),lwd=4,col=3)
title('SIRT1')
title('SIRT1 - Bootstrap correction')

length(SIRT1dat[SIRT1dat$pLU<0.05,]$CMCBeta);length(SIRT1dat[SIRT1dat$pLU<0.05,]$CMCBeta)/JJ
mean(SIRT1dat$CMCBeta);sd(SIRT1dat$CMCBeta)
mean(SIRT1dat[SIRT1dat$pLU<0.05,]$CMCBeta);sd(SIRT1dat[SIRT1dat$pLU<0.05,]$CMCBeta)
mean(with(SIRT1dat,CMCBeta-muBkmed),na.rm=T);sd(with(SIRT1dat,CMCBeta-muBkmed),na.rm=T)

hist(SIRT1dat$CMCBeta,xlab="CMC Beta",main=NULL,freq=T,breaks=10)
hist(SIRT1dat[SIRT1dat$pLU<0.05,]$CMCBeta,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=2,breaks=10)
hist(with(SIRT1dat,BT1),xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=3,breaks=10)
abline(v=mean(SIRT1dat$CMCBeta),lwd=4,col=1)
abline(v=mean(SIRT1dat[SIRT1dat$pLU<0.05,]$CMCBeta),lwd=4,col=2)
abline(v=mean(with(SIRT1dat,BT1),na.rm=T),lwd=4,col=3)
title('SIRT1')
title('SIRT1 - Likelihood Adjusted')
title('SIRT1 - Likelihood Correction 1')



length(SIRT1dat[SIRT1dat$pLU<0.05,]$CMCBeta);length(SIRT1dat[SIRT1dat$pLU<0.05,]$CMCBeta)/JJ
mean(SIRT1dat$CMCBeta);sd(SIRT1dat$CMCBeta)
mean(SIRT1dat[SIRT1dat$pLU<0.05,]$CMCBeta);sd(SIRT1dat[SIRT1dat$pLU<0.05,]$CMCBeta)
mean(with(SIRT1dat,BT1),na.rm=T);sd(with(SIRT1dat, BT1),na.rm=T)

hist(SIRT1dat$CMCBeta,xlab="CMC Beta",main=NULL,freq=T,breaks=10)
hist(SIRT1dat[SIRT1dat$pLU<0.05,]$CMCBeta,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=2,breaks=10)
hist(with(SIRT1dat,BT2),xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=3,breaks=10)
abline(v=mean(SIRT1dat$CMCBeta),lwd=4,col=1)
abline(v=mean(SIRT1dat[SIRT1dat$pLU<0.05,]$CMCBeta),lwd=4,col=2)
abline(v=mean(with(SIRT1dat,BT2),na.rm=T),lwd=4,col=3)
title('SIRT1 - Likelihood Correction 2')

length(SIRT1dat[SIRT1dat$pLU<0.05,]$CMCBeta)
mean(SIRT1dat$CMCBeta);sd(SIRT1dat$CMCBeta)
mean(SIRT1dat[SIRT1dat$pLU<0.05,]$CMCBeta);sd(SIRT1dat[SIRT1dat$pLU<0.05,]$CMCBeta)
mean(with(SIRT1dat,BT2),na.rm=T);sd(with(SIRT1dat,BT2),na.rm=T)

hist(SIRT1dat$CMCBeta,xlab="CMC Beta",main=NULL,freq=T,breaks=10)
hist(SIRT1dat[SIRT1dat$pLU<0.05,]$CMCBeta,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=2,breaks=5)
hist(with(SIRT1dat,BT3),xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=3,breaks=10)
abline(v=mean(SIRT1dat$CMCBeta),lwd=4,col=1)
abline(v=mean(SIRT1dat[SIRT1dat$pLU<0.05,]$CMCBeta),lwd=4,col=2)
abline(v=mean(with(SIRT1dat,BT3),na.rm=T),lwd=4,col=3)
title('SIRT1 - Likelihood Correction 3')

length(SIRT1dat[SIRT1dat$pLU<0.05,]$CMCBeta)
mean(SIRT1dat$CMCBeta);sd(SIRT1dat$CMCBeta)
mean(SIRT1dat[SIRT1dat$pLU<0.05,]$CMCBeta);sd(SIRT1dat[SIRT1dat$pLU<0.05,]$CMCBeta)
mean(with(SIRT1dat,BT3),na.rm=T);sd(with(SIRT1dat,BT3),na.rm=T)


#Individual Betas (3,4,10,11 causal)

hist(SIRT1dat$B3,xlab="CMC Beta",main=NULL,freq=T,breaks=10)

p1=hist(SIRT1dat[SIRT1dat$pLU<0.05,]$B3,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=2,breaks=10)

hist(SIRT1dat[SIRT1dat$pLU<0.05,]$B3,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
hist(SIRT1dat[SIRT1dat$pLU<0.05,]$B4,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=3,breaks=10)
abline(v=mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B3),lwd=4,col=2)
abline(v=mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B4),lwd=4,col=3)

p1=hist(SIRT1dat[SIRT1dat$pLU<0.05,]$B3,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
p2=hist(SIRT1dat[SIRT1dat$pLU<0.05,]$B4,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
p3=hist(SIRT1dat[SIRT1dat$pLU<0.05,]$B10,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
p4=hist(SIRT1dat[SIRT1dat$pLU<0.05,]$B11,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=5)

plot( p1, col=rgb(0,0,1,1/4), xlim=c(-2,5),ylim=c(0,200))  
plot( p2, col=rgb(1,0,0,1/4), add=T)        # The variability in the corrected estimates is dependant of the MAF -- smaller MAF = larger varaince in adjusted estimate
plot( p3, col=rgb(1,1,0,1/4), add=T)  
plot( p4, col=rgb(1,1,1,1/4), add=T)


p1=hist(SIRT1dat$B1,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
p2=hist(SIRT1dat[SIRT1dat$pLU<0.05,]$B1,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
plot(p1)  
plot( p2, col=2, add=T) 
abline(v=mean(SIRT1dat$B1),lwd=4)
abline(v=mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B1),lwd=4,col=2)

par(mfrow=c(1,1))

p1=hist(SIRT1dat$B3,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
p2=hist(SIRT1dat[SIRT1dat$pLU<0.05,]$B3,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10,add=T)
plot(p1)  
plot( p2, col=2, add=T) 
abline(v=mean(SIRT1dat$B3),lwd=4)
abline(v=mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B3),lwd=4,col=2)

p1=hist(SIRT1dat$B4,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
p2=hist(SIRT1dat[SIRT1dat$pLU<0.05,]$B4,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
plot(p1)  
plot( p2, col=2, add=T) 
abline(v=mean(SIRT1dat$B4),lwd=4)
abline(v=mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B4),lwd=4,col=2)

p1=hist(SIRT1dat$B10,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
p2=hist(SIRT1dat[SIRT1dat$pLU<0.05,]$B10,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
plot(p1)  
plot( p2, col=2, add=T) 
abline(v=mean(SIRT1dat$B10),lwd=4)
abline(v=mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B10),lwd=4,col=2)

p1=hist(SIRT1dat$B11,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
p2=hist(SIRT1dat[SIRT1dat$pLU<0.05,]$B11,xlab="CMC Beta",main="SIRT1",freq=T,col=2,breaks=10)
plot(p1)  
plot( p2, col=2, add=T) 
abline(v=mean(SIRT1dat$B11),lwd=4)
abline(v=mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B11),lwd=4,col=2)


abline(v=mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B3),lwd=4,col=rgb(0,0,1,1/4))
abline(v=mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B4),lwd=4,col=rgb(1,0,0,1/4))
abline(v=mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B10),lwd=4,col=rgb(1,0,0,1/4))
abline(v=mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B11),lwd=4,col=rgb(1,0,1,1/4))




length(SIRT1dat[SIRT1dat$pLU<0.05,]$CMCBeta);length(SIRT1dat[SIRT1dat$pLU<0.05,]$CMCBeta)/JJ
mean(SIRT1dat$CMCBeta);sd(SIRT1dat$CMCBeta)
mean(SIRT1dat[SIRT1dat$pLU<0.05,]$CMCBeta);sd(SIRT1dat[SIRT1dat$pLU<0.05,]$CMCBeta)
##Truth pLU
vec=vec1
vec[c(3,4,10,11)]
mean(vec[c(3,4,10,11)])
##Truth over replicates (JJ=2000)
mean(SIRT1dat$B3);mean(SIRT1dat$B4);mean(SIRT1dat$B10);mean(SIRT1dat$B11);
##Bias-Winner's Curse
mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B3)
mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B4)
mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B10)
mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B11)
##Absolute Bias
(mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B3)-vec[3])
(mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B4)-vec[4])
(mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B10)-vec[10])
(mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B11)-vec[11])
##Relative Bias
(mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B3)-vec[3])/vec[3]
(mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B4)-vec[4])/vec[4]
(mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B10)-vec[10])/vec[10]
(mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B11)-vec[11])/vec[11]



length(SIRT1dat[SIRT1dat$pLW<0.05,]$CMCBeta);length(SIRT1dat[SIRT1dat$pLW<0.05,]$CMCBeta)/JJ
mean(SIRT1dat$CMCBeta);sd(SIRT1dat$CMCBeta)
mean(SIRT1dat[SIRT1dat$pLW<0.05,]$CMCBeta);sd(SIRT1dat[SIRT1dat$pLW<0.05,]$CMCBeta)
##Truth  pLW
vec=vec1
vec[c(3,4,10,11)]
##Truth over replicates (JJ=2000)
mean(SIRT1dat$B3);mean(SIRT1dat$B4);mean(SIRT1dat$B10);mean(SIRT1dat$B11);
##Bias-Winner's Curse
mean(SIRT1dat[SIRT1dat$pLW<0.05,]$B3)
mean(SIRT1dat[SIRT1dat$pLW<0.05,]$B4)
mean(SIRT1dat[SIRT1dat$pLW<0.05,]$B10)
mean(SIRT1dat[SIRT1dat$pLW<0.05,]$B11)
##Absolute Bias
(mean(SIRT1dat[SIRT1dat$pLW<0.05,]$B3)-vec[3])
(mean(SIRT1dat[SIRT1dat$pLW<0.05,]$B4)-vec[4])
(mean(SIRT1dat[SIRT1dat$pLW<0.05,]$B10)-vec[10])
(mean(SIRT1dat[SIRT1dat$pLW<0.05,]$B11)-vec[11])
##Relative Bias
(mean(SIRT1dat[SIRT1dat$pLW<0.05,]$B3)-vec[3])/vec[3]
(mean(SIRT1dat[SIRT1dat$pLW<0.05,]$B4)-vec[4])/vec[4]
(mean(SIRT1dat[SIRT1dat$pLW<0.05,]$B10)-vec[10])/vec[10]
(mean(SIRT1dat[SIRT1dat$pLW<0.05,]$B11)-vec[11])/vec[11]


vec=c(0,0,0.83224,0.97060,0,0,0,0,0,0.93459,0.53073)
apply(Xsd,2,FUN=sum)[c(3,4,10,11)]/321

apply(Xsd,2,FUN=sum)
apply(Xsd,1,FUN=sum)
Xsd[266,c(4)]
Xsd[266,]


























############## Bootstrap replication under Hotelling T^2 test
###### SIRT1

which(SIRT1dat$pBeta<0.05&SIRT1dat$pQH<0.05)
which(SIRT1dat$pQH<0.05)
which(SIRT1dat$pBeta<0.05)

head(SIRT1dat)

SIRT1dat=newSIRT1.1
vec=c(0,0,0.83224,0.97060,0,0,0,0,0,0.93459,0.53073)

#SIRT1dat=newSIRT1.2
#SIRT1dat=newSIRT1.3

SIRT1dat$muBkmed=NA
SIRT1dat$muPkmed=NA
SIRT1dat$muBkmea=NA
SIRT1dat$muPkmea=NA


print(name)
DD = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/gene_info',header=T,sep=',')
ii = which(DD[,1]==name)
start = DD[ii,4]
last  = DD[ii,5]
chr = DD[ii,2]
name_to_read =  paste('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/c',chr,'_snps.unr',sep='')
A = read.csv(name_to_read,header=T,sep=',')
ID = A[,1]
subID = getPop(ID) 
B = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/snp_info',header=T,sep=',')
subSNP = subSet(B,names(A),start,last)
X = convertX(A[subID,subSNP])


MAF = getmaf(X)
SD = MAF<0.05 & MAF>0
Xsd  =(X[,SD])
if(name=='GCKR'){
Xsd  =matrix(X[,SD]) 
}
apply(Xsd,2,FUN=sum) 

#MAF = getmaf(X[BK,])
#SD = MAF<0.05 & MAF>0
#Xsd  =(X[,SD])
apply(Xsd[BK,],2,FUN=sum) 

XCMC=apply(Xsd,1,FUN=max)

dim(Xsd)
xx = names(A)
#xx[subSNP[MAF>0]]
#xx[subSNP[MAF<0.05]]

P = NULL

c=5
for(c in which(SIRT1dat$pQH<0.05)) {

print(paste(c,name))

BetaBk=NA
pBk=NA
BetaCk=NA
pCk=NA
pH=NULL
pt=NULL
sigK=0



for(k in 1:J){
	print(paste(k,"-",sigK))
	if(sigK==200) break
	BK=Bk[,k]
	CK=sam[-BK]	
	print(apply(Xsd[BK,],2,FUN=sum))
#	file = paste('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/unr_phen.',c,sep='')
#	Y = getY(ID[subID],file)
	Y = Phenos[,c]

	SQ = as.numeric(testQTHboot(Y[BK],Xsd[BK,],Xsd))
	perm=10^3
	storeQH = rep(0,perm)
	for (i in 1:perm){
		Ynew = sample(Y[BK])
		storeQH[i] = as.numeric(testQTHboot(Ynew,Xsd[BK,],Xsd))
	}
	pH[k] = sum(abs(storeQH)>=abs(SQ))/perm
	pH[k]
	
	
    

#	pt[k]=summary(lm(Y[BK]~XCMC[BK]))$coefficients[2,4]
#	pt[k]
	if(pH[k]<0.05){
		
		sigK=sigK+1
		
		pBk[k]=mean(XCMC[BK])
		if(pBk[k]>0){
		BetaBk[k]=summary(lm(Y[BK]~XCMC[BK]))$coefficients[2,1]
		}
		pCk[k]=mean(XCMC[CK])
		if(pCk[k]>0){
		BetaCk[k]=summary(lm(Y[CK]~XCMC[CK]))$coefficients[2,1]
		}
	}
	
}

#cbind(pH,pt,BetaBk,pBk,BetaCk,pCk)

muBkmed=median(BetaBk-BetaCk,na.rm = T)
muPkmed=median(pBk-pCk,na.rm = T)
SIRT1dat$muBkmed[c]=muBkmed
SIRT1dat$muPkmed[c]=muPkmed

muBkmea=mean(BetaBk-BetaCk,na.rm = T)
muPkmea=mean(pBk-pCk,na.rm = T)
SIRT1dat$muBkmea[c]=muBkmea
SIRT1dat$muPkmea[c]=muPkmea

}


# write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot',name,sep=''))
# bootSIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot',name,sep=''))

# write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot1',name,sep=''))
# boot1SIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot1',name,sep=''))

# write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot2',name,sep=''))
# boot2SIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot2',name,sep=''))

head(boot1SIRT1,100)


###Apply Ghosh et. al. to the significant replicates
SIRT1dat$BT1=NA
SIRT1dat$BT2=NA
SIRT1dat$BT3=NA

source("/Users/senadkokic/Desktop/NSERC/R code/GW-functions.R")# YOU READ IN THE FUNCTIONS FOR OUR METHOD


cee=abs(qnorm(.5*0.05)) # Bonferroni threshold for achieving study-wide significance = 0.1
                          # we use "cee" so R does not get confused with the function 'c'
cc=1
for(cc in which(SIRT1dat$pQH<0.05)) {
print(cc)

betahat=SIRT1dat[cc,]$CMCBeta # Reported OR = 1.54 
z=sign(betahat)*abs(qnorm(0.5*SIRT1dat[cc,]$pQH)) # Reported p-value = 5.7e-4, which we convert to a z-value


###################################################
#              THE PROPOSED APPROACH              #
###################################################

se=betahat/z # standard error of betahat

mutilde1=optimize(f=conditional.like,c(-20,20),maximum=T,z=z,cee=cee)$maximum # the conditional mle
betatilde1=mutilde1*se # the conditional mle on beta scale
SIRT1dat[cc,]$BT1=betatilde1# the condition mle OR estimate
SIRT1dat[cc,]$BT1 # print out OR_TILDE_1

a=sum(apply(as.matrix(seq(-50,50,by=0.01)),1,conditional.like.z,cee=cee,z=z))*0.01 # numeric integration
b=sum(apply(as.matrix(seq(-50,50,by=0.01)),1,conditional.like,cee=cee,z=z))*0.01
betatilde2=(a/b)*se
SIRT1dat[cc,]$BT2=betatilde2
SIRT1dat[cc,]$BT2 # print out OR_TILDE_2

betatilde3=(betatilde1+betatilde2)/2
SIRT1dat[cc,]$BT3=betatilde3
SIRT1dat[cc,]$BT3 # OR_TILDE_3

cimu=getCI(z=z,cee=cee,alpha=0.05)
cibeta=c(cimu$lowermu*(betahat/z),cimu$uppermu*(betahat/z))
cibeta # BIAS-CORRECTED 95% CI FOR OR
}

write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/bootQ',name,sep=''))
bootSIRT1Q=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/bootQ',name,sep=''))

write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot1Q',name,sep=''))
boot1SIRT1Q=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot1Q',name,sep=''))
head(boot1SIRT1Q,100)

write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot2Q',name,sep=''))
boot2SIRT1Q=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot2Q',name,sep=''))

























######## Now with dual criteria Hotelling <0.05 and CMCBeta<0.05
SIRT1dat$muBkmed2=NA
SIRT1dat$muPkmed2=NA
SIRT1dat$muBkmea2=NA
SIRT1dat$muPkmea2=NA


print(name)
DD = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/gene_info',header=T,sep=',')
ii = which(DD[,1]==name)
start = DD[ii,4]
last  = DD[ii,5]
chr = DD[ii,2]
name_to_read =  paste('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/c',chr,'_snps.unr',sep='')
A = read.csv(name_to_read,header=T,sep=',')
ID = A[,1]
subID = getPop(ID) 
B = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/snp_info',header=T,sep=',')
subSNP = subSet(B,names(A),start,last)
X = convertX(A[subID,subSNP])


MAF = getmaf(X)
SD = MAF<0.05 & MAF>0
Xsd  =(X[,SD])
if(name=='GCKR'){
Xsd  =matrix(X[,SD]) 
}
apply(Xsd,2,FUN=sum) 

#MAF = getmaf(X[BK,])
#SD = MAF<0.05 & MAF>0
#Xsd  =(X[,SD])
apply(Xsd[BK,],2,FUN=sum) 

XCMC=apply(Xsd,1,FUN=max)

dim(Xsd)
xx = names(A)
#xx[subSNP[MAF>0]]
#xx[subSNP[MAF<0.05]]

P = NULL


for(c in which(SIRT1dat$pBeta<0.05&SIRT1dat$pQH<0.05)) {

print(paste(c,name))

BetaBk=NA
pBk=NA
BetaCk=NA
pCk=NA
pH=NULL
pt=NULL
sigK=0


for(k in 1:J){
	print(paste(k,"-",sigK))
	if(sigK==100) break
	BK=Bk[,k]
	CK=sam[-BK]	
	print(apply(Xsd[BK,],2,FUN=sum))
#	file = paste('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/unr_phen.',c,sep='')
#	Y = getY(ID[subID],file)
	Y = Phenos[,c]
	
	SQ = as.numeric(testQTHboot(Y[BK],Xsd[BK,],Xsd))
	perm=10^3
	storeQH = rep(0,perm)
	for (i in 1:perm){
		Ynew = sample(Y[BK])
		storeQH[i] = as.numeric(testQTHboot(Ynew,Xsd[BK,],Xsd))
	}
	pH[k] = sum(abs(storeQH)>=abs(SQ))/perm
	pH[k]
	
	
 
	
	
	pt[k]=summary(lm(Y[BK]~XCMC[BK]))$coefficients[2,4]
	pt[k]
	if(pt[k]<0.05&pH[k]<0.05){
		
		sigK=sigK+1
		
		pBk[k]=mean(XCMC[BK])
		if(pBk[k]>0){
		BetaBk[k]=summary(lm(Y[BK]~XCMC[BK]))$coefficients[2,1]
		}
		pCk[k]=mean(XCMC[CK])
		if(pCk[k]>0){
		BetaCk[k]=summary(lm(Y[CK]~XCMC[CK]))$coefficients[2,1]
		}
	}
	
}


cbind(pH,pt,BetaBk,pBk,BetaCk,pCk)

muBkmed2=median(BetaBk-BetaCk,na.rm = T)
muPkmed2=median(pBk-pCk,na.rm = T)
SIRT1dat$muBkmed2[c]=muBkmed2
SIRT1dat$muPkmed2[c]=muPkmed2

muBkmea2=mean(BetaBk-BetaCk,na.rm = T)
muPkmea2=mean(pBk-pCk,na.rm = T)
SIRT1dat$muBkmea2[c]=muBkmea2
SIRT1dat$muPkmea2[c]=muPkmea2

}





























mean(SIRT1dat[SIRT1dat$pBeta<0.05&SIRT1dat$pQH<0.05,]$CMCBeta) # both Hotelling T2 and Beta (AGE) are significant
sd(SIRT1dat[SIRT1dat$pBeta<0.05&SIRT1dat$pQH<0.05,]$CMCBeta)

mean(with(SIRT1dat,CMCBeta-muBk),na.rm=T)    # CMC AGE bootstrap adjusted (Leal)
sd(with(SIRT1dat,CMCBeta-muBk),na.rm=T)      # CMC AGE bootstrap adjusted (Leal) - SD is slightly inflated!

sum(is.na(SIRT1dat$muBk)==FALSE)
length(SIRT1dat[SIRT1dat$pBeta<0.05&SIRT1dat$pQH<0.05,]$CMCBeta)

2*(1-pt(1.711,310))
qt(0.025,310)

hist(SIRT1dat$CMCBeta,xlab="CMC Beta",main=NULL,freq=T,breaks=20)
hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=2,breaks=5)


hist(SIRT1dat[SIRT1dat$pBeta<0.05&SIRT1dat$pQH<0.05,]$CMCBeta,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=2,breaks=10)
#hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=2,breaks=10)
hist(with(SIRT1dat,CMCBeta-muBk),xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=3,breaks=10)
abline(v=mean(SIRT1dat$CMCBeta),lwd=4,col=1)
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05&SIRT1dat$pQH<0.05,]$CMCBeta),lwd=4,col=2)
abline(v=mean(with(SIRT1dat,CMCBeta-muBk),na.rm=T),lwd=4,col=3)

title('SIRT1')









































j=29
XCMC=apply(Xsd,1,FUN=max)
XCMC[XCMC>0]=XCMC[XCMC>0]^0
qj[go]==mean(XCMC)
Betaj=NULL
pvalj=NULL
for (j in 1:200){
	file = paste('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/unr_phen.',j,sep='')
	Y = getY(ID[subID],file)
	Betaj[j]=summary(lm(Y~XCMC))$coefficients[2,1]
	pvalj[j]=summary(lm(Y~XCMC))$coefficients[2,4]
	print(j)

}
Bj=cbind(Bj,Betaj)
Pj=cbind(Pj,pvalj)










##### The true Effects!
go=0
qj=NULL
Betaj=NULL
pvalj=NULL
Bj=NULL
Pj=NULL



for (name in c('SIRT1','BCHE','PDGFD','SREBF1','GCKR','VLDLR','PLAT','RARB','INSIG1','VNN3','LPL','VWF','VNN1')){

print(name)
go=go+1
DD = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/gene_info',header=T,sep=',')
ii = which(DD[,1]==name)
start = DD[ii,4]
last  = DD[ii,5]
chr = DD[ii,2]
#
#   Read data from the ch3
name_to_read =  paste('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/c',chr,'_snps.unr',sep='')
A = read.csv(name_to_read,header=T,sep=',')
# helper functions
# helper functions

source('/Users/senadkokic/Desktop/NSERC/R code/gawhelp.r')
source('/Users/senadkokic/Desktop/NSERC/R code/gawpermM.r')
source('/Users/senadkokic/Desktop/NSERC/R code/methodsX-winners.r')



ID = A[,1]
subID = getPop(ID) 

B = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/snp_info',header=T,sep=',')
subSNP = subSet(B,names(A),start,last)
X = convertX(A[subID,subSNP])


# Get MAF

MAF = getmaf(X)
SD = MAF<0.05 & MAF>0
Xsd  =(X[,SD]) 
if(name=='GCKR'){
Xsd  =matrix(X[,SD]) 
}
dim(Xsd)
xx = names(A)
xx[subSNP[MAF>0]]
xx[subSNP[MAF<0.05]]
P = NULL

XCMC=apply(Xsd,1,FUN=max)
XCMC[XCMC>0]=XCMC[XCMC>0]^0
qj[go]==mean(XCMC)
Betaj=NULL
pvalj=NULL
for (j in 1:200){
	file = paste('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/unr_phen.',j,sep='')
	Y = getY(ID[subID],file)
	Betaj[j]=summary(lm(Y~XCMC))$coefficients[2,1]
	pvalj[j]=summary(lm(Y~XCMC))$coefficients[2,4]
	print(j)

}
Bj=cbind(Bj,Betaj)
Pj=cbind(Pj,pvalj)
}



head(Pj)
head(Bj)
apply(Bj,2,mean)
mean(Betaj)
colnames(Bj)=c('SIRT1','BCHE','PDGFD','SREBF1','GCKR' ,'VLDLR','PLAT','RARB','INSIG1','VNN3','LPL','VWF','VNN1')
colnames(Pj)=c('SIRT1','BCHE','PDGFD','SREBF1','GCKR' ,'VLDLR','PLAT','RARB','INSIG1','VNN3','LPL','VWF','VNN1')

name
Bj[,name]
Pj[,name]




####### These are the results (from the linear/quadratic/minP/Fisher tests plus the  results from the CMC[collapsing] paramtere estimation procedure)
name='GCKR'
#c('SIRT1','BCHE','PDGFD','SREBF1','GCKR','VLDLR','PLAT','RARB','INSIG1','VNN3','LPL','VWF','VNN1')
for (name in c('GCKR','VLDLR','PLAT','RARB','INSIG1','VNN3','LPL','VWF','VNN1')){
print(name)

DD = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/gene_info',header=T,sep=',')
ii = which(DD[,1]==name)
start = DD[ii,4]
last  = DD[ii,5]
chr = DD[ii,2]
#
#   Read data from the ch3
name_to_read =  paste('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/c',chr,'_snps.unr',sep='')
A = read.csv(name_to_read,header=T,sep=',')
# helper functions
# helper functions

source('/Users/senadkokic/Desktop/NSERC/R code/gawhelp.r')
source('/Users/senadkokic/Desktop/NSERC/R code/gawpermM.r')
source('/Users/senadkokic/Desktop/NSERC/R code/methodsX-winners.r')



ID = A[,1]
subID = getPop(ID) 

B = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/snp_info',header=T,sep=',')
subSNP = subSet(B,names(A),start,last)
X = convertX(A[subID,subSNP])


# Get MAF

MAF = getmaf(X)
SD = MAF<0.05 & MAF>0
Xsd  =(X[,SD]) 
if(name=='GCKR'){
Xsd  =matrix(X[,SD]) 
}

xx = names(A)
xx[subSNP[MAF>0]]
xx[subSNP[MAF<0.05]]
P = NULL


SW= NULL
SL= NULL
SWwin= NULL
SLwin= NULL
SC= NULL
SQ= NULL
### Demonstrate the winner's curse

for (j in 1:200){
file = paste('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/unr_phen.',j,sep='')
Y = getY(ID[subID],file)
X=Xsd

L = length(X[,1])
J = length(X[1,])

maf = colSums(X)/L
SWwin[j] = as.numeric(testQTLwin(Y,X,1/sqrt(maf*(1-maf))))
SW[j] =    as.numeric(testQTL(Y,X,1/sqrt(maf*(1-maf))))
SL[j] = as.numeric(testQTL(Y,X,rep(1,J)))
SLwin[j] = as.numeric(testQTLwin(Y,X,rep(1,J)))
SC[j] = as.numeric(testQTC(Y,X))
SQ[j] = as.numeric(testQTH(Y,X))
print(j)
}


for (j in 1:200){
# get Y

file = paste('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/unr_phen.',j,sep='')
Y = getY(ID[subID],file)
p = pvalCalc(Y,Xsd,10^4)
P = cbind(P,p)
cat(p,'\n')
print(j)
}

frame=cbind(t(P),SW,SL,SWwin,SLwin,SC,SQ,Bj[,name],Pj[,name])
colnames(frame)[c(1,2,3,4,5,6,13,14)]=c('pLW','pLU','pQU','pQH','pMin','pFisher','CMCBeta','pBeta')
write.csv(frame,paste('/Users/senadkokic/Desktop/NSERC/GAW17/',name,sep=''))


}






###################################################################################################
# Obtain the bias for each of the signif replicates
###### Bootstrapping Algorithm (1)
n=length(Y)
J=1000
Bk=matrix(rep(0,n*J),ncol=J)
dim(Bk)
for(k in 1:J){
sam=1:321
Bk[,k]=sample(sam,321,replace=T)
}
head(Bk)
#############################


###### SIRT1
which(SIRT1dat$pBeta<0.05&SIRT1dat$pQH<0.05)
j=NULL

SIRT1dat$muBk=NA
SIRT1dat$muPk=NA
head(SIRT1dat)


print(name)
DD = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/gene_info',header=T,sep=',')
ii = which(DD[,1]==name)
start = DD[ii,4]
last  = DD[ii,5]
chr = DD[ii,2]
name_to_read =  paste('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/c',chr,'_snps.unr',sep='')
A = read.csv(name_to_read,header=T,sep=',')
ID = A[,1]
subID = getPop(ID) 
B = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/snp_info',header=T,sep=',')
subSNP = subSet(B,names(A),start,last)
X = convertX(A[subID,subSNP])


MAF = getmaf(X)
SD = MAF<0.05 & MAF>0
Xsd  =(X[,SD])
if(name=='GCKR'){
Xsd  =matrix(X[,SD]) 
}
apply(Xsd,2,FUN=sum) 

#MAF = getmaf(X[BK,])
#SD = MAF<0.05 & MAF>0
#Xsd  =(X[,SD])
apply(Xsd[BK,],2,FUN=sum) 

dim(Xsd)
xx = names(A)
#xx[subSNP[MAF>0]]
#xx[subSNP[MAF<0.05]]
#P = NULL

which(SIRT1dat$pBeta<0.05&SIRT1dat$pQH<0.05)

for(c in which(SIRT1dat$pBeta<0.05&SIRT1dat$pQH<0.05)) {

print(paste(c,name))

BetaBk=NA
pBk=NA
BetaCk=NA
pCk=NA
pH=NULL
pt=NULL
sigK=0


for(k in 1:J){
	print(paste(k,"-",sigK))
	if(sigK==100) break
	BK=Bk[,k]
	CK=sam[-BK]	
	print(apply(Xsd[BK,],2,FUN=sum))
	file = paste('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/unr_phen.',c,sep='')
	Y = getY(ID[subID],file)

	SQ = as.numeric(testQTHboot(Y[BK],Xsd[BK,],Xsd))
	perm=10^3
	storeQH = rep(0,perm)
	for (i in 1:perm){
		Ynew = sample(Y[BK])
		storeQH[i] = as.numeric(testQTHboot(Ynew,Xsd[BK,],Xsd))
	}
	pH[k] = sum(abs(storeQH)>=abs(SQ))/perm
	pH[k]
	
	
    XCMC=apply(Xsd,1,FUN=max)
	
	
	pt[k]=summary(lm(Y[BK]~XCMC[BK]))$coefficients[2,4]
	pt[k]
	if(pt[k]<0.05&pH[k]<0.05){
		
		sigK=sigK+1
		
		pBk[k]=mean(XCMC[BK])
		if(pBk[k]>0){
		BetaBk[k]=summary(lm(Y[BK]~XCMC[BK]))$coefficients[2,1]
		}
		pCk[k]=mean(XCMC[CK])
		if(pCk[k]>0){
		BetaCk[k]=summary(lm(Y[CK]~XCMC[CK]))$coefficients[2,1]
		}
	}
	
}

cbind(pH,pt,BetaBk,pBk,BetaCk,pCk)

muBk=median(BetaBk-BetaCk,na.rm = T)
muPk=median(pBk-pCk,na.rm = T)
SIRT1dat$muBk[c]=muBk
SIRT1dat$muPk[c]=muPk

}
SIRT1dat$muBk
SIRT1dat$muPk
head(SIRT1dat,110)


write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot',name,sep=''))
bootSIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot',name,sep=''))
head(bootSIRT1)
plot(bootSIRT)





name='GCKR'
c('SIRT1','BCHE','PDGFD','SREBF1','GCKR','VLDLR','PLAT','RARB','INSIG1','VNN3','LPL','VWF','VNN1')
for(name in c('VLDLR','PLAT','RARB','INSIG1','VNN3','LPL','VWF','VNN1')){
fdat=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/',name,sep=''))
fdat
###### BCHE
which(fdat$pBeta<0.05& fdat$pQH<0.05)
j=NULL

fdat $muBk=NA
fdat $muPk=NA
#head(fdat)


print(name)
DD = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/gene_info',header=T,sep=',')
ii = which(DD[,1]==name)
start = DD[ii,4]
last  = DD[ii,5]
chr = DD[ii,2]
name_to_read =  paste('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/c',chr,'_snps.unr',sep='')
A = read.csv(name_to_read,header=T,sep=',')
ID = A[,1]
subID = getPop(ID) 
B = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/snp_info',header=T,sep=',')
subSNP = subSet(B,names(A),start,last)
X = convertX(A[subID,subSNP])


MAF = getmaf(X)
SD = MAF<0.05 & MAF>0
Xsd  =(X[,SD])
if(name=='GCKR'){
Xsd  =matrix(X[,SD]) 
}
apply(Xsd,2,FUN=sum) 

#MAF = getmaf(X[BK,])
#SD = MAF<0.05 & MAF>0
#Xsd  =(X[,SD])
#apply(Xsd[BK,],2,FUN=sum) 

dim(Xsd)
xx = names(A)
#xx[subSNP[MAF>0]]
#xx[subSNP[MAF<0.05]]
#P = NULL

c=4
for(c in which(fdat $pBeta<0.05& fdat $pQH<0.05)) {

print(paste(c,name))

BetaBk=NA
pBk=NA
BetaCk=NA
pCk=NA
pH=NULL
pt=NULL
sigK=0

k=1
for(k in 1:J){
	print(paste(k,"-",sigK))
	if(sigK==100) break
	BK=Bk[,k]
	CK=sam[-BK]	
	print(apply(Xsd[BK,],2,FUN=sum))
	#if (sum(apply(Xsd[BK,],2,FUN=sum)>0) {}
	file = paste('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/unr_phen.',c,sep='')
	Y = getY(ID[subID],file)

	SQ = as.numeric(testQTHboot(Y[BK],Xsd[BK,],Xsd))
	perm=10^3
	storeQH = rep(0,perm)
	for (i in 1:perm){
		Ynew = sample(Y[BK])
		storeQH[i] = as.numeric(testQTHboot(Ynew,Xsd[BK,],Xsd))
	}
	pH[k] = sum(abs(storeQH)>=abs(SQ))/perm
	pH[k]
	
	
    XCMC=apply(Xsd,1,FUN=max)
	
	
	pt[k]=summary(lm(Y[BK]~XCMC[BK]))$coefficients[2,4]
	pt[k]
	if(pt[k]<0.05&pH[k]<0.05){
		
		sigK=sigK+1
		
		pBk[k]=mean(XCMC[BK])
		if(pBk[k]>0){
		BetaBk[k]=summary(lm(Y[BK]~XCMC[BK]))$coefficients[2,1]
		}
		pCk[k]=mean(XCMC[CK])
		if(pCk[k]>0){
		BetaCk[k]=summary(lm(Y[CK]~XCMC[CK]))$coefficients[2,1]
		}
	}
	
}

#cbind(pH,pt,BetaBk,pBk,BetaCk,pCk)

muBk=median(BetaBk-BetaCk,na.rm = T)
muPk=median(pBk-pCk,na.rm = T)
fdat $muBk[c]=muBk
fdat $muPk[c]=muPk

}
#fdat $muBk
#fdat $muPk
#head(fdat,110)


write.csv(fdat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot',name,sep=''))

}





##########################################################################################################################################################################################################################################################################################################################################################################################
name='GCKR'
fdat=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/',name,sep=''))
fdat
###### BCHE
which(fdat$pBeta<0.05& fdat$pQH<0.05)
j=NULL

fdat $muBk=NA
fdat $muPk=NA
#head(fdat)


print(name)
DD = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/gene_info',header=T,sep=',')
ii = which(DD[,1]==name)
start = DD[ii,4]
last  = DD[ii,5]
chr = DD[ii,2]
name_to_read =  paste('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/c',chr,'_snps.unr',sep='')
A = read.csv(name_to_read,header=T,sep=',')
ID = A[,1]
subID = getPop(ID) 
B = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/snp_info',header=T,sep=',')
subSNP = subSet(B,names(A),start,last)
X = convertX(A[subID,subSNP])


MAF = getmaf(X)
SD = MAF<0.05 & MAF>0
Xsd  =(X[,SD])
if(name=='GCKR'){
Xsd  =matrix(X[,SD]) 
}
apply(Xsd,2,FUN=sum) /321

#MAF = getmaf(X[BK,])
#SD = MAF<0.05 & MAF>0
#Xsd  =(X[,SD])
#apply(Xsd[BK,],2,FUN=sum) 

dim(Xsd)
xx = names(A)
#xx[subSNP[MAF>0]]
#xx[subSNP[MAF<0.05]]
#P = NULL

c=4
for(c in which(fdat $pBeta<0.05& fdat $pQH<0.05)) {

print(paste(c,name))

BetaBk=NA
pBk=NA
BetaCk=NA
pCk=NA
pH=NULL
pt=NULL
sigK=0

k=1
for(k in 1:J){
	print(paste(k,"-",sigK))
	if(sigK==100) break
	BK=Bk[,k]
	CK=sam[-BK]	
	print(apply(matrix(Xsd[BK,]),2,FUN=sum))
	#if (sum(apply(Xsd[BK,],2,FUN=sum)>0) {}
	file = paste('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/unr_phen.',c,sep='')
	Y = getY(ID[subID],file)

	SQ = as.numeric(testQTHboot(Y[BK],matrix(Xsd[BK,]),Xsd))
	perm=10^3
	storeQH = rep(0,perm)
	for (i in 1:perm){
		Ynew = sample(Y[BK])
		storeQH[i] = as.numeric(testQTHboot(Ynew,matrix(Xsd[BK,]),Xsd))
	}
	pH[k] = sum(abs(storeQH)>=abs(SQ))/perm
	pH[k]
	
	
    XCMC=apply(Xsd,1,FUN=max)
	
	
	pt[k]=summary(lm(Y[BK]~XCMC[BK]))$coefficients[2,4]
	pt[k]
	if(pt[k]<0.05&pH[k]<0.05){
		
		sigK=sigK+1
		
		pBk[k]=mean(XCMC[BK])
		if(pBk[k]>0){
		BetaBk[k]=summary(lm(Y[BK]~XCMC[BK]))$coefficients[2,1]
		}
		pCk[k]=mean(XCMC[CK])
		if(pCk[k]>0){
		BetaCk[k]=summary(lm(Y[CK]~XCMC[CK]))$coefficients[2,1]
		}
	}
	
}

#cbind(pH,pt,BetaBk,pBk,BetaCk,pCk)

muBk=median(BetaBk-BetaCk,na.rm = T)
muPk=median(pBk-pCk,na.rm = T)
fdat $muBk[c]=muBk
fdat $muPk[c]=muPk

}
#fdat $muBk
#fdat $muPk
#head(fdat,110)


write.csv(fdat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot',name,sep=''))
###################################################################################################################################################################################################################################################








bootBCHE=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot',name,sep=''))
head(bootBCHE)
plot(bootBCHE)






bootSIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/bootSIRT1',sep=''))
bootBCHE=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/bootBCHE',sep=''))


SIRT1dat=read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/bootSIRT1')
SIRT1dat=read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/bootBCHE')
SIRT1dat=read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/bootPDGFD')
SIRT1dat=read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/bootSREBF1')
#SIRT1dat=read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/bootGCKR')
SIRT1dat=read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/bootRARB')
SIRT1dat=read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/bootPLAT')
SIRT1dat=read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/bootVLDLR')
SIRT1dat=read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/bootVNN3')
SIRT1dat=read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/bootINSIG1')
SIRT1dat=read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/bootLPL')
SIRT1dat=read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/bootVWF')



########### Now the analysis etc.
#AGE
#mean(SIRT1dat[SIRT1dat$pQH<0.05,]$CMCBeta)  # where Hotelling T2 is sig
#sd(SIRT1dat[SIRT1dat$pQH<0.05,]$CMCBeta)
#mean(SIRT1dat[SIRT1dat$pFisher<0.05,]$CMCBeta) # where Fisher is sig
#sd(SIRT1dat[SIRT1dat$pFisher<0.05,]$CMCBeta) 
#mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta)  # where Beta (AGE) is sig
#sd(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta) 

sum(SIRT1dat$pQH<0.05)/200
sum(SIRT1dat$pLW<0.05)/200
sum(SIRT1dat$pBeta<0.05)/200

mean(SIRT1dat$CMCBeta) # est of True AGE over the RV (causal and non causal)
sd(SIRT1dat$CMCBeta)

mean(SIRT1dat[SIRT1dat$pBeta<0.05&SIRT1dat$pQH<0.05,]$CMCBeta) # both Hotelling T2 and Beta (AGE) are significant
sd(SIRT1dat[SIRT1dat$pBeta<0.05&SIRT1dat$pQH<0.05,]$CMCBeta)

mean(with(SIRT1dat,CMCBeta-muBk),na.rm=T)    # CMC AGE bootstrap adjusted (Leal)
sd(with(SIRT1dat,CMCBeta-muBk),na.rm=T)      # CMC AGE bootstrap adjusted (Leal) - SD is slightly inflated!

sum(is.na(SIRT1dat$muBk)==FALSE)
length(SIRT1dat[SIRT1dat$pBeta<0.05&SIRT1dat$pQH<0.05,]$CMCBeta)


hist(SIRT1dat$CMCBeta,xlab="CMC Beta",main="SIRT1",freq=FALSE)
hist(SIRT1dat[SIRT1dat$pBeta<0.05&SIRT1dat$pQH<0.05,]$CMCBeta,xlab="CMC Beta",main="SIRT1",freq=FALSE,add=T)
hist(SIRT1dat[SIRT1dat$pBeta<0.05&SIRT1dat$pQH<0.05,]$CMCBeta,xlab="CMC Beta",main="SIRT1",freq=FALSE)

hist(SIRT1dat$CMCBeta,xlab="CMC Beta",main=NULL,freq=T,breaks=10)
hist(SIRT1dat[SIRT1dat$pBeta<0.05&SIRT1dat$pQH<0.05,]$CMCBeta,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=2,breaks=10)
#hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=2,breaks=10)
hist(with(SIRT1dat,CMCBeta-muBk),xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=3,breaks=10)
abline(v=mean(SIRT1dat$CMCBeta),lwd=4,col=1)
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05&SIRT1dat$pQH<0.05,]$CMCBeta),lwd=4,col=2)
abline(v=mean(with(SIRT1dat,CMCBeta-muBk),na.rm=T),lwd=4,col=3)

title('SIRT1')
title('BCHE')
title('PDGFD')
title('SREBF1')
title('RARB')
title('PLAT')
title('VLDLR')
title('VNN3')
title('INSIG1')
title('LPL')
title('VWF')







histogram(SIRT1dat$CMCBeta|SIRT1dat[SIRT1dat$pBeta<0.05&SIRT1dat$pQH<0.05,]$CMCBeta)
abline(v=mean(SIRT1dat$CMCBeta),lwd=4,col=4)
mean(SIRT1dat$CMCBeta)
text(19,55,"-12.6",col=4)
arrows(6,55,-12,55,length=0.1,col=4,lwd=2)

abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05&SIRT1dat$pQH<0.05,]$CMCBeta),lwd=4,col=2)
mean(SIRT1dat[SIRT1dat$pBeta<0.05&SIRT1dat$pQH<0.05,]$CMCBeta)
text(-65,55,"-39.0",col=2)
arrows(-55,55,-40,55,length=0.1,col=2,lwd=2)












































VNN1dat=read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/VNN1')
head(VNN1dat)

hist(VNN1dat$SWwin,xlab="W(L) Statistic",main="VNN1")
abline(v=mean(VNN1dat$SWwin),lwd=4,col=4)
mean(VNN1dat$SWwin)
text(19,55,"-12.6",col=4)
arrows(6,55,-12,55,length=0.1,col=4,lwd=2)

abline(v=mean(VNN1dat[VNN1dat$X.6<0.05,]$SWwin),lwd=4,col=2)
mean(VNN1dat[VNN1dat$X.6<0.05,]$SWwin)
text(-65,55,"-39.0",col=2)
arrows(-55,55,-40,55,length=0.1,col=2,lwd=2)



BCHEdat=read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/BCHE')

hist(BCHEdat$SWwin,xlab="W(L) Statistic",main="BCHE")
abline(v=mean(BCHEdat$SWwin),lwd=4,col=4)
text(10,55,"101",col=4)
mean(BCHEdat$SWwin)
arrows(30,55,99,55,length=0.1,col=4,lwd=2)

abline(v=mean(BCHEdat[BCHEdat$X.6<0.05,]$SWwin),lwd=4,col=2)
text(220,50,"148",col=2)
mean(BCHEdat[BCHEdat$X.6<0.05,]$SWwin)
arrows(190,50,149,50,length=0.1,col=2,lwd=2)


hist(BCHEdat$SW,xlab="W(L) Statistic",main="VNN1")
abline(v=mean(BCHEdat$SW),lwd=4,col=4)
text(10,55,"1.34",col=4)
arrows(1.5,240,1.36,240,length=0.1,col=4,lwd=2)

abline(v=mean(BCHEdat[BCHEdat$X.6<0.05,]$SW),lwd=4,col=2)
text(2.1,220,"1.81",col=2)
arrows(1.95,220,1.82,220,length=0.1,col=2,lwd=2)


SIRT1dat=read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/SIRT1')
head(SIRT1dat)
plot(SIRT1dat$SW,SIRT1dat$SL)
plot(SIRT1dat$SW,SIRT1dat$SWwin)
plot(SIRT1dat$SC,SIRT1dat$SQ)

# Master run for all genes
rm(list=ls(all=TRUE))

name = 'SIRT1'

DD = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/gene_info',header=T,sep=',')
ii = which(DD[,1]==name)
start = DD[ii,4]
last  = DD[ii,5]
chr = DD[ii,2]
#
#   Read data from the ch3
name_to_read =  paste('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/c',chr,'_snps.unr',sep='')
A = read.csv(name_to_read,header=T,sep=',')
# helper functions
# helper functions

source('/Users/senadkokic/Desktop/NSERC/R code/gawhelp.r')
source('/Users/senadkokic/Desktop/NSERC/R code/gawpermM.r')
source('/Users/senadkokic/Desktop/NSERC/R code/methodsX.r')



ID = A[,1]
subID = getPop(ID) 

B = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/snp_info',header=T,sep=',')
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


j=7

file = paste('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/unr_phen.',j,sep='')
Y = getY(ID[subID],file)
p = pvalCalc(Y,Xsd,10^4)
P = cbind(P,p)
cat(p,'\n')
print(j)

### AGE (Beta)
XCMC=apply(Xsd,1,FUN=sum)
XCMC
summary(lm(Y~XCMC))
summary(lm(Y~XCMC))$coefficients[2,1]












?hist
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



SW=0
SL=0
SWwin=0
SLwin=0
SC=0
SQ=0
### Demonstrate the winner's curse

for (j in 1:200){
file = paste('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/unr_phen.',j,sep='')
Y = getY(ID[subID],file)
Y
X=Xsd

L = length(X[,1])
J = length(X[1,])
maf = colSums(X)/L
SWwin[j] = as.numeric(testQTLwin(Y,X,1/sqrt(maf*(1-maf))))
SW[j] = as.numeric(testQTL(Y,X,1/sqrt(maf*(1-maf))))
SL[j] = as.numeric(testQTL(Y,X,rep(1,J)))
SLwin[j] = as.numeric(testQTLwin(Y,X,rep(1,J)))
SC[j] = as.numeric(testQTC(Y,X))
SQ[j] = as.numeric(testQTH(Y,X))
print(j)
}
SWwin
SLwin
SW
SL
SC
SQ

hist(SWwin)
hist(SLwin)
hist(SW)
hist(SL)
hist(SC)
hist(SQ)

mean(SL)
mean(SQ)




OR=0
pval=0
n=250
b0=log(.2)
b1=log(1.3)
for(i in 1:1000){

flog=function(x){
	rbinom(1,1,exp(b0+b1*x)/(1+exp(b0+b1*x)))
}
x=data.frame(c(rep(0,n),rep(1,n)))
names(x)=c("x")
x
x$y=apply(x,1,FUN=flog)
x
g1=glm(x~y,family=binomial,data=x)
OR[i]=exp(g1$coef[[2]])
pval[i]=summary(g1)$coefficients[2,4]
print(i)
}

hist(OR)
abline(v=mean(OR),lwd=4,col=4)
text(1.65,240,"1.34",col=4)
arrows(1.5,240,1.36,240,length=0.1,col=4,lwd=2)
abline(v=mean(OR[pval<0.05]),lwd=4,col=2)
text(2.1,220,"1.81",col=2)
arrows(1.95,220,1.82,220,length=0.1,col=2,lwd=2)

?abline
mean(OR)
mean(OR)
hist(OR[pval<0.05])
mean(OR[pval<0.05])
