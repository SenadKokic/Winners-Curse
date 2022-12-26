library(lattice)

# Master run for all genes
# rm(list=ls(all=TRUE))
# This script will anlyze GAW17 data.
#  gene under consideration
# 

##### The true Effects!
go=0
qj=NULL
Betaj=NULL
pvalj=NULL
Bj=NULL
Pj=NULL

for (name in c('SIRT1','BCHE','PDGFD','SREBF1','GCKR','VLDLR','PLAT','RARB','INSIG1','VNN3','LPL','VWF','VNN1')){

  # print the name of the gene
  print(name)
  # iterate counter
  go=go+1
  # read gene info file 
  DD = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/gene_info',header=T,sep=',')
  # find which index in the gene info matches the gene
  ii = which(DD[,1]==name)
  # get the start position, end position, and chromosome associated with gene
  start = DD[ii,4]
  last  = DD[ii,5]
  chr = DD[ii,2]
  #
  # Read data from the chromosome snps file
  name_to_read =  paste('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/c',chr,'_snps.unr',sep='')
  A = read.csv(name_to_read,header=T,sep=',')
  # helper functions
  # helper functions
  names(A)
  source('/Users/senadkokic/Desktop/NSERC/R code/gawhelp.r')
  source('/Users/senadkokic/Desktop/NSERC/R code/gawpermM.r')
  source('/Users/senadkokic/Desktop/NSERC/R code/methodsX-winners.r')
  # I believe these have to do with the names (?) of the SNPs
  ID = A[,1]
  # get the Han Chinese date
  subID = getPop(ID) 
  # read the file containing SNP  info
  B = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/snp_info',header=T,sep=',')
  # generates a subset from the SNP file
  subSNP = subSet(B,names(A),start,last)
  # convert the call matrix to a numeric matrix
  X = convertX(A[subID,subSNP])
  # Get MAF
  MAF = getmaf(X)
  # List indicating which MAF fit the criteria (rare variants)
  SD = MAF<0.05 & MAF>0
  # generate matrix based on above criteria
  Xsd  =(X[,SD]) 
  Xsd
  if(name=='GCKR'){
    # convert to a matrix
    Xsd  =matrix(X[,SD]) 
  }
  # of rows and columns in Xsd
  dim(Xsd)
  xx = names(A)
  # gets SNPs with desired MAF criteria
  xx[subSNP[MAF>0]]
  xx[subSNP[MAF<0.05]]
  P = NULL
  
  # applies maximum function to Xsd matrix
  XCMC=apply(Xsd,1,FUN=max)
  XCMC[XCMC>0]=XCMC[XCMC>0]^0
  XCMC[XCMC>0]
  # entry is the mean of the XCMC
  qj[go]==mean(XCMC)
  Betaj=NULL
  pvalj=NULL
  for (j in 1:200){
  	file = paste('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/unr_phen.',j,sep='')
  	# see gawhelp.r for function info
  	Y = getY(ID[subID],file)
  	# enter coefficients from linear model into beta and pval matrix
  	Betaj[j]=summary(lm(Y~XCMC))$coefficients[2,1]
  	pvalj[j]=summary(lm(Y~XCMC))$coefficients[2,4]
  	#print(j)
  
  }
  # combine matrices
  Bj=cbind(Bj,Betaj)
  Pj=cbind(Pj,pvalj)
}



head(Pj)
head(Bj)
# consider the means of the columns of Bj matrix
apply(Bj,2,mean)
# take the mean of Betaj
mean(Betaj)
# create column names for Bj and Pj matrices
colnames(Bj)=c('SIRT1','BCHE','PDGFD','SREBF1','GCKR' ,'VLDLR','PLAT','RARB','INSIG1','VNN3','LPL','VWF','VNN1')
colnames(Pj)=c('SIRT1','BCHE','PDGFD','SREBF1','GCKR' ,'VLDLR','PLAT','RARB','INSIG1','VNN3','LPL','VWF','VNN1')

# extract column from Bj and Pj matrix 
name
Bj[,name]
Pj[,name]




####### These are the results (from the linear/quadratic/minP/Fisher tests plus the  results from the CMC[collapsing] paramtere estimation procedure)
name='GCKR'
#c('SIRT1','BCHE','PDGFD','SREBF1','GCKR','VLDLR','PLAT','RARB','INSIG1','VNN3','LPL','VWF','VNN1')
for (name in c('GCKR','VLDLR','PLAT','RARB','INSIG1','VNN3','LPL','VWF','VNN1')){
  print(name)
  # read the gene info file
  DD = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/gene_info',header=T,sep=',')
  # take index in which the gene name appears
  ii = which(DD[,1]==name)
  # take the starting and ending position of the gene, and the chromosome
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
  
  
  # get the first column of the chromosome file
  ID = A[,1]
  # get the Han Chinese data based on the column
  subID = getPop(ID) 
  # read the SNP info
  B = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/snp_info',header=T,sep=',')
  # generate a subset of the SNPs 
  subSNP = subSet(B,names(A),start,last)
  # convert A[subID, subSNP] to a numeric matrix
  X = convertX(A[subID,subSNP])
  # Get MAF
  MAF = getmaf(X)
  # list indicating which elements of X are rare variants
  SD = MAF<0.05 & MAF>0
  # matrix based on SD aobve
  Xsd  =(X[,SD]) 
  # if the name matches the gene, create matrix 
  if(name=='GCKR'){
  Xsd  =matrix(X[,SD]) 
  }
  
  # take the names of the chromosome file
  xx = names(A)
  # see which subSNPs have a positive MAF, as well as those constituting as rare variants
  xx[subSNP[MAF>0]]
  xx[subSNP[MAF<0.05]]
  P = NULL
  
  # initialize matrices 
  SW= NULL
  SL= NULL
  SWwin= NULL
  SLwin= NULL
  SC= NULL
  SQ= NULL
  ### Demonstrate the winner's curse
  
  for (j in 1:200){
  # see gawhelp.r 
  file = paste('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/unr_phen.',j,sep='')
  Y = getY(ID[subID],file)
  X=Xsd
  # of rows and columns respectively
  L = length(X[,1])
  J = length(X[1,])
  # calculate the maf 
  maf = colSums(X)/L
  # calculate statistics, place in respective matrices
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
  # combine the matrices
  frame=cbind(t(P),SW,SL,SWwin,SLwin,SC,SQ,Bj[,name],Pj[,name])
  # delegate names to matrix columns
  colnames(frame)[c(1,2,3,4,5,6,13,14)]=c('pLW','pLU','pQU','pQH','pMin','pFisher','CMCBeta','pBeta')
  write.csv(frame,paste('/Users/senadkokic/Desktop/NSERC/GAW17/',name,sep=''))
}
###################################################################################################
# Obtain the bias for each of the signif replicates
###### Bootstrapping Algorithm (1)
# length of the Q2 values
n=length(Y)
J=1000
# create matrix 
Bk=matrix(rep(0,n*J),ncol=J)
dim(Bk)
# sample, insert into Bk matrix
for(k in 1:J){
sam=1:321
Bk[,k]=sample(sam,321,replace=T)
}
head(Bk)
#############################

###### SIRT1
head(SIRT1dat)
# determine indices in which pBeta and pQH match criteria
which(SIRT1dat$pBeta<0.05&SIRT1dat$pQH<0.05)
j=NULL
SIRT1dat$muBk=NA
SIRT1dat$muPk=NA
head(SIRT1dat)
# essentially the same as above
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
# apply sum to the columns of Xsd matrix
apply(Xsd,2,FUN=sum) 

#MAF = getmaf(X[BK,])
#SD = MAF<0.05 & MAF>0
#Xsd  =(X[,SD])

# BK is from the simulate file (line 333 in the for loop)
apply(Xsd[BK,],2,FUN=sum) 

dim(Xsd)
xx = names(A)
#xx[subSNP[MAF>0]]
#xx[subSNP[MAF<0.05]]
#P = NULL
# get indices of which SIRT1 elements meet criteria
which(SIRT1dat$pBeta<0.05&SIRT1dat$pQH<0.05)

# iterate through indices that meet criteria
for(c in which(SIRT1dat$pBeta<0.05&SIRT1dat$pQH<0.05)) {

print(paste(c,name))

# bootstrap sample coefficient
BetaBk=NA
# p value of bootstrap sample
pBk=NA
# coefficient for those left out of the sample
BetaCk=NA
# p value of those left out of bootstrap 
pCk=NA
# permutation p-value 
pH=NULL
# Wald p value
pt=NULL
# of significant reps
sigK=0

# there's a subscript out of bounds error that appears here... 
# should check that out...
for(k in 1:J){
  # print count, how many sig reps have been encountered
	print(paste(k,"-",sigK))
  # if 100 are counted stop the loop
	if(sigK==100) break
  # bootstrap 
	BK=Bk[,k]
	# left out of bootstrap
	CK=sam[-BK]	
  # print the sum of the covariates starting at BK
	print(apply(Xsd[BK,],2,FUN=sum))
	file = paste('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/unr_phen.',c,sep='')
	# get q2 values (see gawhelp.r)
	Y = getY(ID[subID],file)
  # convert QTH bootstrap tests character vector into numeric vector
	SQ = as.numeric(testQTHboot(Y[BK],Xsd[BK,],Xsd))
	# set the number of permutations
	perm=10^3
	# creates storage for bootstrap tests
	storeQH = rep(0,perm)
	for (i in 1:perm){
	  # sample from the bootstrap
		Ynew = sample(Y[BK])
		# store the computed statistic in the QH list
		storeQH[i] = as.numeric(testQTHboot(Ynew,Xsd[BK,],Xsd))
	}
	# calculate the permutation p-value
	pH[k] = sum(abs(storeQH)>=abs(SQ))/perm
	pH[k]
	
	# max of covariate matrix
  XCMC=apply(Xsd,1,FUN=max)
	
	# another p-value (Wald)
	pt[k]=summary(lm(Y[BK]~XCMC[BK]))$coefficients[2,4]
	pt[k]
	# if both the Wald and Permutation result are < 0.05, we have a 
	#   significant replicate 
	if(pt[k]<0.05&pH[k]<0.05){
		# increment the count of significant results
		sigK=sigK+1
	  # entry of the p value bootstrap matrix at k will be the mean of the 
		#   maximum covariate matrix BK'th row
		pBk[k]=mean(XCMC[BK])
		# if the above entry is > 0
		if(pBk[k]>0){
		  # coefficient drawn from the bootstrap sample
		  BetaBk[k]=summary(lm(Y[BK]~XCMC[BK]))$coefficients[2,1]
		}
		# entry of the p value bootstrap-excluded matrix at k will be the mean of 
		#   the maximum covariate matrix CK'th row
		pCk[k]=mean(XCMC[CK])
		if(pCk[k]>0){
		  # coefficient drawn from data left out of bootstrap
		  BetaCk[k]=summary(lm(Y[CK]~XCMC[CK]))$coefficients[2,1]
		}
	}
	
}
# there's a problem here too
# Warning message:
# In cbind(pH, pt, BetaBk, pBk, BetaCk, pCk) :
#  number of rows of result is not a multiple of vector length (arg 2)

# combine the vectors in order to make a table
cbind(pH,pt,BetaBk,pBk,BetaCk,pCk)

# take the medians of the estimates and p values 
#   (recall that there is no significant difference between the mean and median
#    for this, median more computationally efficient)
muBk=median(BetaBk-BetaCk,na.rm = T)
muPk=median(pBk-pCk,na.rm = T)

# insert as SIRT1dat attributes
# these will be inserted where the pBeta and pQH < 0.05
SIRT1dat$muBk[c]=muBk
SIRT1dat$muPk[c]=muPk

}
# display lists
SIRT1dat$muBk
SIRT1dat$muPk
# take 110 entries of SIRT1dat
head(SIRT1dat,110)

# write file
  # *****there's an error here... file doesn't exist
write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot',name,sep=''))
# read the file written above
bootSIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/boot',name,sep=''))
head(bootSIRT1)
plot(bootSIRT)





name='GCKR'
c('SIRT1','BCHE','PDGFD','SREBF1','GCKR','VLDLR','PLAT','RARB','INSIG1','VNN3','LPL','VWF','VNN1')
for(name in c('VLDLR','PLAT','RARB','INSIG1','VNN3','LPL','VWF','VNN1')){
  # read the csv that matches the gene name
  fdat=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/',name,sep=''))
  fdat
  ###### BCHE
  # determine the criteria that match 
  which(fdat$pBeta<0.05& fdat$pQH<0.05)
  j=NULL
  # set initial attributes to be null
  fdat $muBk=NA
  fdat $muPk=NA
  #head(fdat)
  # print the name of the gene
  print(name)
  # same as the beginning of the file 
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

# write bootstrap sample file associated with the gene
write.csv(fdat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot',name,sep=''))
###################################################################################################################################################################################################################################################
# read boot file for BCHE
bootBCHE=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot',name,sep=''))
# display data from bootBCHE
head(bootBCHE)
# plot the data 
plot(bootBCHE)

# read SIRT1 and BCHE data
bootSIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/bootSIRT1',sep=''))
bootBCHE=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/bootBCHE',sep=''))

# I'm assuming these are to be commented out depending on which gene is being examined
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
# sum up attributes of SIRT1dat 
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
