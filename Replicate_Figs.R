#########################################################################
# Clear global environment
rm(list = ls(all = TRUE))
# Load the lattice library for data viz
library(lattice)
# Load required 
source('/Users/senadkokic/Desktop/NSERC/R code/gawhelp.r')
source('/Users/senadkokic/Desktop/NSERC/R code/gawpermM.r')
source('/Users/senadkokic/Desktop/NSERC/R code/methodsX-winners.r')
# YOU READ IN THE FUNCTIONS FOR OUR METHOD
source("/Users/senadkokic/Desktop/NSERC/R code/GW-functions.R")
# Select the SIRT1 gene
name = 'SIRT1'
# Read the gene_info file in order to get gene information
geneInfo = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/gene_info',
                     header = TRUE, sep = ',')
# Use which in order to determine the column that matches SIRT1
matchingIndex = which(geneInfo[, 1] == name)
# Get starting/ending position of gene, and chromosome number
startingPos = geneInfo[matchingIndex, 4]
endingPos = geneInfo[matchingIndex, 5]
chrNum = geneInfo[matchingIndex, 2]
# Chromosome file name
chrFileName = paste('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/c',
                    chrNum,'_snps.unr',sep='')
# Read the file 
chrData = read.csv(chrFileName, header = TRUE, sep = ',')
# Take the headers (Chromosome IDs)
chrID = chrData[, 1]
# Call function to get Han Chinese data from the file
subID = getPop(chrID)
# Read snp_info file to get SNP data
snpData = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/snp_info',
                   header=T,sep=',')
# Call function to get subset of SNPs, based on positions of gene
subSNP = subSet(snpData, names(chrData), startingPos, endingPos)
# Convert call matrix into numeric matrix
numericMatrix = convertX(chrData[subID, subSNP])
# Determine the MAF
MAF = getmaf(numericMatrix)
# Generates a matrix of Boolean values based on "Rare Variants"
SD = MAF < 0.05 & MAF > 0
# Generate covariate matrix
covariateMatrix = (numericMatrix[, SD])
# Take the sum of columns of covariate matrix
apply(covariateMatrix, 2, FUN = sum)
apply(covariateMatrix, 2, FUN = sum)/321
chrNames = names(chrData)
P = NULL
# Apply CMC to the covariate matrix
covariateMatrixCMC = apply(covariateMatrix, 1, FUN = max)
# To make elements of above matrix all 0s and 1s
covariateMatrixCMC[covariateMatrixCMC > 0] =
  covariateMatrixCMC[covariateMatrixCMC > 0]^0
# Set parameters
JJ = 2000
Betaj = NULL
pvalj = NULL
qj = NULL
tBetaj = NULL
Bj = NULL
Tj = NULL
Pj = NULL
SW = NULL
SL = NULL
SWwin = NULL
SLwin = NULL
SC = NULL
SQ = NULL
NBetas = 11
IBetas=matrix(rep(NA,JJ*NBetas),ncol= NBetas)
TIBetas=matrix(rep(NA,JJ*NBetas),ncol= NBetas)
PIBetas=matrix(rep(NA,JJ*NBetas),ncol= NBetas)
Phenos=matrix(rep(NA,JJ*321),ncol= JJ)
# Set seed
set.seed(1)
# Start loop
for (j in 1:JJ) {
  # print(j)
  vec=c(0,0,0.83224,0.97060,0,0,0,0,0,0.93459,0.53073)
  # Create Design
  Design = apply(sweep(covariateMatrix, MARGIN = 2, vec, `*`), 1, FUN = sum)
  # Add normal simulations to the Design matrix
  Y = Design + rnorm(length(Design), 0, 1)
  # Insert elements of Y into Phenos matrix
  Phenos[, j] = Y
  # Assigns covariate CMC of estimate
  Betaj[j] = summary(lm(Y~covariateMatrixCMC))$coefficients[2, 1]
  # Covariate CMC of t value
  tBetaj[j] = summary(lm(Y~covariateMatrixCMC))$coefficients[2, 3]
  # Covariate CMC of Pr(>|t|)
  pvalj[j] = summary(lm(Y~covariateMatrixCMC))$coefficients[2, 4]
  # for # of variants
  for (i in 1:11) {
    # estimate
    IBetas[j,i]=summary(lm(Y~covariateMatrix[,i]))$coeff[2,1]
    # t value
    TIBetas[j,i]=summary(lm(Y~covariateMatrix[,i]))$coeff[2,3]
    # Pr(>|t|)
    PIBetas[j,i]=summary(lm(Y~covariateMatrix[,i]))$coeff[2,4]
  }
  numericMatrix = covariateMatrix
  # Get the number of rows and columns 
  numRows = length(numericMatrix[, 1])
  numCols = length(numericMatrix[1, ])
  # Taking the MAF as the # of occurrences divided by the sample size
  maf = colSums(numericMatrix)/numRows
  # Calculate statistics
  SWwin[j] = as.numeric(testQTLwin(Y,numericMatrix,1/sqrt(maf*(1-maf))))
  SW[j] =    as.numeric(testQTL(Y,numericMatrix,1/sqrt(maf*(1-maf))))
  SL[j] = as.numeric(testQTL(Y,numericMatrix,rep(1,numCols)))
  SLwin[j] = as.numeric(testQTLwin(Y,numericMatrix,rep(1,numCols))) 
  SC[j] = as.numeric(testQTC(Y,numericMatrix))
  SQ[j] = as.numeric(testQTH(Y,numericMatrix))
  # Calculate p values
  p = pvalCalc(Y, covariateMatrix, 10^3)
  # Add p value calculated above to P
  P = cbind(P, p)
  # Print p value
  # cat(p, '\n)
}
# Combine vectors/matrices by the columns
Bj=cbind(Bj,Betaj)
Tj=cbind(Tj,tBetaj)
Pj=cbind(Pj,pvalj)
# Data frame to combine the various estimates
frame=cbind(t(P),SW,SL,SWwin,SLwin,SC,SQ,Bj,Tj,Pj,IBetas)
# Name columns of data frame accordingly
colnames(frame)[c(1,2,3,4,5,6,13:26)]=c('pLW','pLU','pQU','pQH','pMin','pFisher',
                                        'CMCBeta','TCMC','pBeta','B1','B2','B3',
                                        'B4','B5','B6','B7','B8','B9','B10','B11')
# Create SIRT1 data based on the data frame
write.csv(frame, '/Users/senadkokic/Desktop/NSERC/GAW17/newsimSIRT1')
# Read the above data as SIRT1 data
newSIRT1 = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/newsimSIRT1')
# Change column names of Betas to match variants
colnames(IBetas)=c('B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11')
colnames(TIBetas)=c('B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11')
colnames(PIBetas)=c('B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11')
# Make PI and TI Betas into data frames
PIBetas = data.frame(PIBetas)
TIBetas = data.frame(TIBetas)
# Calculates correlation between TCMC and stats
with(newSIRT1,cor(TCMC,SW))
with(newSIRT1,cor(TCMC,SL))
with(newSIRT1,cor(TCMC,SC))
with(newSIRT1,cor(TCMC,SQ))

# In accordance with simulation type 1, create type 1 variables
IBetas1 = IBetas
TIBetas1 = TIBetas
PIBetas1 = PIBetas
newSIRT1.1 = newSIRT1
SIRT1dat = newSIRT1.1
################################################################################
# Obtain bias for all of the significant replicates
## Bootstrapping Algorithm (1)
n = length(Y)
numCols = 10000
# Initialize matrix for sample
Bk = matrix(rep(0, n*numCols), ncol = numCols)
# Delegate elements of matrix to be a sample
set.seed(2)
for (k in 1:numCols) {
  sam = 1:321
  Bk[, k] = sample(sam, 321, replace = TRUE)
}
################################################################################
# Bootstrap replication under the Wald Beta test
### SIRT1
SIRT1dat=newSIRT1.1
# Initialize elements to be NA
SIRT1dat$muBkmed=NA
SIRT1dat$muPkmed=NA
SIRT1dat$muBkmea=NA
SIRT1dat$muPkmea=NA

numericMatrix = convertX(chrData[subID, subSNP])
MAF = getmaf(numericMatrix)
SD = MAF < 0.05 & MAF > 0
covariateMatrix = (numericMatrix[, SD])
apply(covariateMatrix, 2, FUN = sum)
chrNames = names(chrData)
P = NULL
# Apply CMC to the covariate matrix
covariateMatrixCMC = apply(covariateMatrix, 1, FUN = max)
# Loop for elements that match criteria wrt pBeta
for (c in which(SIRT1dat$pBeta < 0.05)) {
  # print the index #
  # print(c)
  # Set initial values
  BetaBk=NA
  pBk=NA
  BetaCk=NA
  pCk=NA
  pH=NULL
  pt=NULL
  sigK=0
  for (k in 1:numCols) {
    # Print iteration as well as how many significant results have been seen
    # NOTE: this prints REALLY fast
    # print(paste(k, ' - ', sigK))
    # Break when 200 are observed (? not sure why)
    if (sigK == 200) break
    # Bootstrap sample
    BK = Bk[, k]
    # Remaining elements (left out of bootstrap)
    CK = sam[-BK]
    # Set Y = Phenos element based on index
    Y = Phenos[, c]
    # Compute linear model coefficient for Wald Statistic
    pt[k] = summary(lm(Y[BK]~covariateMatrixCMC[BK]))$coefficients[2,4]
    pt[k]
    # Significant result if the Wald stat < 0.05
    if (pt[k] < 0.05) {
      sigK = sigK + 1
      # Take the means of CMC matrix based on those sampled and unsampled
      pBk[k] = mean(covariateMatrixCMC[BK])
      if (pBk[k] > 0) {
        BetaBk[k]=summary(lm(Y[BK]~covariateMatrixCMC[BK]))$coefficients[2,1]
      }
      pCk[k] = mean(covariateMatrixCMC[CK])
      if (pCk[k] > 0) {
        BetaCk[k]=summary(lm(Y[CK]~covariateMatrixCMC[CK]))$coefficients[2,1]
      }
    }
  }
  # Insert differences between estimates into SIRT1dat
  # Median method
  muBkmed=median(BetaBk-BetaCk,na.rm = T)
  muPkmed=median(pBk-pCk,na.rm = T)
  SIRT1dat$muBkmed[c]=muBkmed
  SIRT1dat$muPkmed[c]=muPkmed
  # Mean method
  muBkmea=mean(BetaBk-BetaCk,na.rm = T)
  muPkmea=mean(pBk-pCk,na.rm = T)
  SIRT1dat$muBkmea[c]=muBkmea
  SIRT1dat$muPkmea[c]=muPkmea
}
write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot1',name,sep=''))
boot1SIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot1',name,sep=''))
### Apply Ghosh et. al. to the significant replicates
SIRT1dat$BT1=NA
SIRT1dat$BT2=NA
SIRT1dat$BT3=NA
# Bonferroni threshold for achieving study-wide significance = 0.1
# Named cee in order to avoid confusion b/w function c()
cee = abs(qnorm(0.5*0.05))
# Begin loop for pBeta 
for (cc in which(SIRT1dat$pBeta < 0.05)) {
  # Print the index
  # print(cc)
  # Proposed odds ratio = 1.54
  betahat = SIRT1dat[cc, ]$CMCBeta
  # Reported p-value = 5.7e-4, which we convert to a z-value
  z=sign(betahat)*abs(qnorm(0.5*SIRT1dat[cc,]$pBeta))
  ###################################################
  #              THE PROPOSED APPROACH              #
  ###################################################
  # Compute the standard error of betahat
  se = betahat/z
  # Compute the conditional MLE
  mutilde1 = optimize(f=conditional.like,
                      c(-20,20),maximum=T,z=z,cee=cee)$maximum
  # Compute the conditional MLE on the beta scale
  betatilde1 = mutilde1*se
  # Assign the conditional MLE of the odds ratio estimate 
  SIRT1dat[cc, ]$BT1 = betatilde1
  # Print out the estimate 
  SIRT1dat[cc, ]$BT1
  # Numeric integration
  a=sum(apply(as.matrix(seq(-50,50,by=0.01)),1,
              conditional.like.z,cee=cee,z=z))*0.01
  b=sum(apply(as.matrix(seq(-50,50,by=0.01)),1,
              conditional.like,cee=cee,z=z))*0.01
  # Alternate computation of odds ratio
  betatilde2 = (a/b)*se
  # Assign the estimate, print it
  SIRT1dat[cc,]$BT2=betatilde2
  # Compute the average of the previous ORs as the third estimate
  betatilde3 = (betatilde1 + betatilde2)/2
  # Assign the estimate, print it
  SIRT1dat[cc,]$BT3=betatilde3

  # Create confidence intervals for the true OR (bias-corrected, alpha = 0.05)
  cimu = getCI(z = z, cee = cee, alpha = 0.05)
  cibeta = c(cimu$lowermu*(betahat/z),cimu$uppermu*(betahat/z))
  # Print the bias corrected interval
  cibeta 
}
write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot',
                         name,sep=''))
bootSIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot',
                         name,sep=''))
write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot1',
                         name,sep=''))
boot1SIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot1',
                          name,sep=''))
write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot2',
                         name,sep=''))
boot2SIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot2',
                          name,sep=''))
##### Histogram comparisons:
SIRT1dat = bootSIRT1
IBetas=IBetas1
TIBetas=TIBetas1
PIBetas=PIBetas1

head(SIRT1dat)
tail(SIRT1dat)
# Set destination pdf for plots
#pdf('/Users/senadkokic/Desktop/NSERC/Results Files/fig1a.pdf')
hist(SIRT1dat$CMCBeta,xlab= expression(beta[CMC]),main=NULL,freq=T, col = "white", breaks=10)
hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta,xlab="CMC Beta",main="SIRT1",
     freq=T,add=T,col=2,breaks=10)
hist(with(SIRT1dat,CMCBeta-muBkmed),xlab="CMC Beta",main="SIRT1",
     freq=T,add=T,col=3,breaks=5)
abline(v=mean(SIRT1dat$CMCBeta),lwd=4,col=1)
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta),lwd=4,col=2)

abline(v=mean(with(SIRT1dat,CMCBeta-muBkmed),na.rm=T),lwd=4,col=3)
title('SIRT1 - Bootstrap correction')
#dev.off()

#pdf('/Users/senadkokic/Desktop/NSERC/Results Files/fig1c.pdf')
hist(SIRT1dat$CMCBeta,xlab=expression(beta[CMC]),main=NULL,freq=T,col = "white", breaks=10)
hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=2,breaks=10)
hist(with(SIRT1dat,BT1),xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=3,breaks=5)
abline(v=mean(SIRT1dat$CMCBeta),lwd=4,col=1)
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta),lwd=4,col=2)
abline(v=mean(with(SIRT1dat,BT1),na.rm=T),lwd=4,col=3)
title('SIRT1 - Likelihood Correction 1')
#dev.off()
###############################################################################
# For the Wald test (Table 1)
##Truth
sink('/Users/senadkokic/Desktop/NSERC/Results Files/Values_CMC_Uni')
round(vec[c(3,4,10,11)], 2)

##Truth over replicates (JJ=2000)
# These give B_k values for table
round(mean(SIRT1dat$B3), 2)
round(mean(SIRT1dat$B4), 2)
round(mean(SIRT1dat$B10), 2)
round(mean(SIRT1dat$B11), 2)
##Bias-Winner's Curse
# These give B_naive values for table
round(mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B3), 2)
round(mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B4), 2)
round(mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B10), 2)
round(mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B11), 2)
##Absolute Bias
round((mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B3)-vec[3]), 2)
round((mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B4)-vec[4]), 2)
round((mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B10)-vec[10]), 2)
round((mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B11)-vec[11]), 2)
##Relative Bias
round((mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B3)-vec[3])/vec[3], 2)
round((mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B4)-vec[4])/vec[4], 2)
round((mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B10)-vec[10])/vec[10], 2)
round((mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B11)-vec[11])/vec[11], 2)

round(sum(MAF[SD]*abs(vec))/sum(MAF[SD]), 2) # Calculation for true AGE
round(mean(SIRT1dat$CMCBeta), 2) # Estimate of True AGE over the RV (causal and non causal)
round(sd(SIRT1dat$CMCBeta), 2)
round(mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta), 2) # Naive estimate
round(sd(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta), 2)
round(mean(with(SIRT1dat,CMCBeta-muBkmed), na.rm = T), 2) # Bootstrap estimate
round(sd(with(SIRT1dat,CMCBeta-muBkmed), na.rm = T), 2)
round(mean(with(SIRT1dat,BT1),na.rm=T), 2) # Likelihood estimate
round(sd(with(SIRT1dat,BT1),na.rm=T), 2)
# Get the MAF
round(MAF[SD], 3)
# Power calculation
powerCalc = length(which(SIRT1dat$pBeta < 0.05))
paste("Power ", powerCalc, "/", JJ, " = ", round(powerCalc/JJ, 2), sep = '')
sink()
#########################################################################


#### STOP HERE
########################
# Bidirectional CMC Beta
########################
# Select the SIRT1 gene
name = 'SIRT1'
# Read the gene_info file in order to get gene information
geneInfo = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/gene_info',
                    header = TRUE, sep = ',')
# Use which in order to determine the column that matches SIRT1
matchingIndex = which(geneInfo[, 1] == name)
# Get starting/ending position of gene, and chromosome number
startingPos = geneInfo[matchingIndex, 4]
endingPos = geneInfo[matchingIndex, 5]
chrNum = geneInfo[matchingIndex, 2]
# Chromosome file name
chrFileName = paste('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/c',
                    chrNum,'_snps.unr',sep='')
# Read the file 
chrData = read.csv(chrFileName, header = TRUE, sep = ',')
# Take the headers (Chromosome IDs)
chrID = chrData[, 1]
# Call function to get Han Chinese data from the file
subID = getPop(chrID)
# Read snp_info file to get SNP data
snpData = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/snp_info',
                   header=T,sep=',')
# Call function to get subset of SNPs, based on positions of gene
subSNP = subSet(snpData, names(chrData), startingPos, endingPos)
# Convert call matrix into numeric matrix
numericMatrix = convertX(chrData[subID, subSNP])
# Determine the MAF
MAF = getmaf(numericMatrix)
# Generates a matrix of Boolean values based on "Rare Variants"
SD = MAF < 0.05 & MAF > 0
# Generate covariate matrix
covariateMatrix = (numericMatrix[, SD])
# Take the sum of columns of covariate matrix
apply(covariateMatrix, 2, FUN = sum)
apply(covariateMatrix, 2, FUN = sum)/321
chrNames = names(chrData)
P = NULL
# Apply CMC to the covariate matrix
covariateMatrixCMC = apply(covariateMatrix, 1, FUN = max)
# To make elements of above matrix all 0s and 1s
covariateMatrixCMC[covariateMatrixCMC > 0] =
  covariateMatrixCMC[covariateMatrixCMC > 0]^0
# Set parameters
JJ = 2000
Betaj = NULL
pvalj = NULL
qj = NULL
tBetaj = NULL
Bj = NULL
Tj = NULL
Pj = NULL
SW = NULL
SL = NULL
SWwin = NULL
SLwin = NULL
SC = NULL
SQ = NULL
NBetas = 11
IBetas=matrix(rep(NA,JJ*NBetas),ncol= NBetas)
TIBetas=matrix(rep(NA,JJ*NBetas),ncol= NBetas)
PIBetas=matrix(rep(NA,JJ*NBetas),ncol= NBetas)
Phenos=matrix(rep(NA,JJ*321),ncol= JJ)
# Start loop
set.seed(3)
for (j in 1:JJ) {
  # print(j)
  vec=c(0,0,0.83224,0.97060,0,0,0,0,0,-0.93459,0.53073)
  # Create Design
  Design = apply(sweep(covariateMatrix, MARGIN = 2, vec, `*`), 1, FUN = sum)
  # Add normal simulations to the Design matrix
  Y = Design + rnorm(length(Design), 0, 1)
  # Insert elements of Y into Phenos matrix
  Phenos[, j] = Y
  # Assigns covariate CMC of estimate
  Betaj[j] = summary(lm(Y~covariateMatrixCMC))$coefficients[2, 1]
  # Covariate CMC of t value
  tBetaj[j] = summary(lm(Y~covariateMatrixCMC))$coefficients[2, 3]
  # Covariate CMC of Pr(>|t|)
  pvalj[j] = summary(lm(Y~covariateMatrixCMC))$coefficients[2, 4]
  # for # of variants
  for (i in 1:11) {
    # estimate
    IBetas[j,i]=summary(lm(Y~covariateMatrix[,i]))$coeff[2,1]
    # t value
    TIBetas[j,i]=summary(lm(Y~covariateMatrix[,i]))$coeff[2,3]
    # Pr(>|t|)
    PIBetas[j,i]=summary(lm(Y~covariateMatrix[,i]))$coeff[2,4]
  }
  numericMatrix = covariateMatrix
  # Get the number of rows and columns 
  numRows = length(numericMatrix[, 1])
  numCols = length(numericMatrix[1, ])
  # Taking the MAF as the # of occurrences divided by the sample size
  maf = colSums(numericMatrix)/numRows
  # Calculate statistics
  SWwin[j] = as.numeric(testQTLwin(Y,numericMatrix,1/sqrt(maf*(1-maf))))
  SW[j] =    as.numeric(testQTL(Y,numericMatrix,1/sqrt(maf*(1-maf))))
  SL[j] = as.numeric(testQTL(Y,numericMatrix,rep(1,numCols)))
  SLwin[j] = as.numeric(testQTLwin(Y,numericMatrix,rep(1,numCols))) 
  SC[j] = as.numeric(testQTC(Y,numericMatrix))
  SQ[j] = as.numeric(testQTH(Y,numericMatrix))
  # Calculate p values
  p = pvalCalc(Y, covariateMatrix, 10^3)
  # Add p value calculated above to P
  P = cbind(P, p)
  # Print p value
  # cat(p, '\n)
}
# Combine vectors/matrices by the columns
Bj=cbind(Bj,Betaj)
Tj=cbind(Tj,tBetaj)
Pj=cbind(Pj,pvalj)
# Data frame to combine the various estimates
frame=cbind(t(P),SW,SL,SWwin,SLwin,SC,SQ,Bj,Tj,Pj,IBetas)
# Name columns of data frame accordingly
colnames(frame)[c(1,2,3,4,5,6,13:26)]=c('pLW','pLU','pQU','pQH','pMin','pFisher',
                                        'CMCBeta','TCMC','pBeta','B1','B2','B3',
                                        'B4','B5','B6','B7','B8','B9','B10','B11')
# Partially displays elements of above data frame
head(frame)
# Displays # of rows and columns of above data frame
dim(frame)
# Create SIRT1 data based on the data frame
write.csv(frame, '/Users/senadkokic/Desktop/NSERC/GAW17/newsimSIRT1')
# Read the above data as SIRT1 data
newSIRT1 = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/newsimSIRT1')
# Partially display elements of above data (data read correctly)
head(newSIRT1)
# Change column names of Betas to match variants
colnames(IBetas)=c('B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11')
colnames(TIBetas)=c('B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11')
colnames(PIBetas)=c('B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11')
# Make PI and TI Betas into data frames
PIBetas = data.frame(PIBetas)
TIBetas = data.frame(TIBetas)
# In accordance with simulation type 1, create type 1 variables
IBetas1 = IBetas
TIBetas1 = TIBetas
PIBetas1 = PIBetas
newSIRT1.1 = newSIRT1

################ Begin analysis
SIRT1dat = newSIRT1.1
################################################################################
# Obtain bias for all of the significant replicates
## Bootstrapping Algorithm (1)
n = length(Y)
numCols = 10000
# Initialize matrix for sample
Bk = matrix(rep(0, n*numCols), ncol = numCols)
# Delegate elements of matrix to be a sample
set.seed(4)
for (k in 1:numCols) {
  sam = 1:321
  Bk[, k] = sample(sam, 321, replace = TRUE)
}
################################################################################
# Bootstrap replication under the Wald Beta test
### SIRT1
# See which indices match criteria
which(SIRT1dat$pBeta<0.05&SIRT1dat$pQH<0.05)
which(SIRT1dat$pQH<0.05)
which(SIRT1dat$pBeta<0.05)
SIRT1dat=newSIRT1.1
#SIRT1dat=newSIRT1.2
#SIRT1dat=newSIRT1.3
# Initialize elements to be NA
SIRT1dat$muBkmed=NA
SIRT1dat$muPkmed=NA
SIRT1dat$muBkmea=NA
SIRT1dat$muPkmea=NA

numericMatrix = convertX(chrData[subID, subSNP])
MAF = getmaf(numericMatrix)
MAF
SD = MAF < 0.05 & MAF > 0
SD
covariateMatrix = (numericMatrix[, SD])
covariateMatrix
apply(covariateMatrix, 2, FUN = sum)
chrNames = names(chrData)
P = NULL
# Apply CMC to the covariate matrix
covariateMatrixCMC = apply(covariateMatrix, 1, FUN = max)
# Loop for elements that match criteria wrt pBeta
for (c in which(SIRT1dat$pBeta < 0.05)) {
  # print the index #
  # print(c)
  # Set initial values
  BetaBk=NA
  pBk=NA
  BetaCk=NA
  pCk=NA
  pH=NULL
  pt=NULL
  sigK=0
  for (k in 1:numCols) {
    # Print iteration as well as how many significant results have been seen
    # NOTE: this prints REALLY fast
    # print(paste(k, ' - ', sigK))
    # Break when 200 are observed (? not sure why)
    if (sigK == 200) break
    # Bootstrap sample
    BK = Bk[, k]
    # Remaining elements (left out of bootstrap)
    CK = sam[-BK]
    # Set Y = Phenos element based on index
    Y = Phenos[, c]
    # Compute linear model coefficient for Wald Statistic
    pt[k] = summary(lm(Y[BK]~covariateMatrixCMC[BK]))$coefficients[2,4]
    pt[k]
    
    # Significant result if the Wald stat < 0.05
    if (pt[k] < 0.05) {
      sigK = sigK + 1
      # Take the means of CMC matrix based on those sampled and unsampled
      pBk[k] = mean(covariateMatrixCMC[BK])
      if (pBk[k] > 0) {
        BetaBk[k]=summary(lm(Y[BK]~covariateMatrixCMC[BK]))$coefficients[2,1]
      }
      pCk[k] = mean(covariateMatrixCMC[CK])
      if (pCk[k] > 0) {
        BetaCk[k]=summary(lm(Y[CK]~covariateMatrixCMC[CK]))$coefficients[2,1]
      }
    }
  }
  # Insert differences between estimates into SIRT1dat
  # Median method
  muBkmed=median(BetaBk-BetaCk,na.rm = T)
  muPkmed=median(pBk-pCk,na.rm = T)
  SIRT1dat$muBkmed[c]=muBkmed
  SIRT1dat$muPkmed[c]=muPkmed
  # Mean method
  muBkmea=mean(BetaBk-BetaCk,na.rm = T)
  muPkmea=mean(pBk-pCk,na.rm = T)
  SIRT1dat$muBkmea[c]=muBkmea
  SIRT1dat$muPkmea[c]=muPkmea
}
write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot1',name,sep=''))
boot1SIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot1',name,sep=''))
### Apply Ghosh et. al. to the significant replicates
SIRT1dat$BT1=NA
SIRT1dat$BT2=NA
SIRT1dat$BT3=NA
# Bonferroni threshold for achieving study-wide significance = 0.1
# Named cee in order to avoid confusion b/w function c()
cee = abs(qnorm(0.5*0.05))
# Begin loop for pBeta 
for (cc in which(SIRT1dat$pBeta < 0.05)) {
  # Print the index
  # print(cc)
  # Proposed odds ratio = 1.54
  betahat = SIRT1dat[cc, ]$CMCBeta
  # Reported p-value = 5.7e-4, which we convert to a z-value
  z=sign(betahat)*abs(qnorm(0.5*SIRT1dat[cc,]$pBeta))
  ###################################################
  #              THE PROPOSED APPROACH              #
  ###################################################
  # Compute the standard error of betahat
  se = betahat/z
  # Compute the conditional MLE
  mutilde1 = optimize(f=conditional.like,
                      c(-20,20),maximum=T,z=z,cee=cee)$maximum
  # Compute the conditional MLE on the beta scale
  betatilde1 = mutilde1*se
  # Assign the conditional MLE of the odds ratio estimate 
  SIRT1dat[cc, ]$BT1 = betatilde1
  # Print out the estimate 
  SIRT1dat[cc, ]$BT1
  # Numeric integration
  a=sum(apply(as.matrix(seq(-50,50,by=0.01)),1,
              conditional.like.z,cee=cee,z=z))*0.01
  b=sum(apply(as.matrix(seq(-50,50,by=0.01)),1,
              conditional.like,cee=cee,z=z))*0.01
  # Alternate computation of odds ratio
  betatilde2 = (a/b)*se
  # Assign the estimate, print it
  SIRT1dat[cc,]$BT2=betatilde2
  # Compute the average of the previous ORs as the third estimate
  betatilde3 = (betatilde1 + betatilde2)/2
  # Assign the estimate, print it
  SIRT1dat[cc,]$BT3=betatilde3
  
  # Create confidence intervals for the true OR (bias-corrected, alpha = 0.05)
  cimu = getCI(z = z, cee = cee, alpha = 0.05)
  cibeta = c(cimu$lowermu*(betahat/z),cimu$uppermu*(betahat/z))
  # Print the bias corrected interval
  cibeta 
}
write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot',
                         name,sep=''))
bootSIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot',
                         name,sep=''))

write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot1',
                         name,sep=''))
boot1SIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot1',
                          name,sep=''))
head(boot1SIRT1,100)

write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot2',
                         name,sep=''))
boot2SIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot2',
                          name,sep=''))
##### Histogram comparisons:
SIRT1dat = bootSIRT1
IBetas=IBetas1
TIBetas=TIBetas1
PIBetas=PIBetas1
###############################################################################
#pdf('/Users/senadkokic/Desktop/NSERC/Results Files/fig1b.pdf')
hist(SIRT1dat$CMCBeta,xlab= expression(beta[CMC]),main=NULL,freq=T, col = "white",breaks=10)
hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta,xlab="CMC Beta",main="SIRT1",
     freq=T,add=T,col=2,breaks=10) 
hist(with(SIRT1dat,CMCBeta-muBkmed),xlab="CMC Beta",main="SIRT1",
     freq=T,add=T,col=3,breaks=5)
abline(v=mean(SIRT1dat$CMCBeta),lwd=4,col=1)
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta),lwd=4,col=2)

abline(v=mean(with(SIRT1dat,CMCBeta-muBkmed),na.rm=T),lwd=4,col=3)
title('SIRT1 - Bootstrap correction')
#dev.off()

#pdf('/Users/senadkokic/Desktop/NSERC/Results Files/fig1d.pdf')
hist(SIRT1dat$CMCBeta,xlab=expression(beta[CMC]),main=NULL,freq=T, col = "white", breaks=10)
hist(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta,xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=2,breaks=10)
hist(with(SIRT1dat,BT1),xlab="CMC Beta",main="SIRT1",freq=T,add=T,col=3,breaks=5)
abline(v=mean(SIRT1dat$CMCBeta),lwd=4,col=1)
abline(v=mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta),lwd=4,col=2)
abline(v=mean(with(SIRT1dat,BT1),na.rm=T),lwd=4,col=3)
title('SIRT1 - Likelihood Correction 1')
#dev.off()
###############################################################################
sink('/Users/senadkokic/Desktop/NSERC/Results Files/Values_CMC_Bi')
# For the Wald test (Table 1)
##Truth
round(vec[c(3,4,10,11)], 2)

##Truth over replicates (JJ=2000)
# These give B_k values for table
round(mean(SIRT1dat$B3), 2)
round(mean(SIRT1dat$B4), 2)
round(mean(SIRT1dat$B10), 2)
round(mean(SIRT1dat$B11), 2)
##Bias-Winner's Curse
# These give B_naive values for table
round(mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B3), 2)
round(mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B4), 2)
round(mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B10), 2)
round(mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B11), 2)
##Absolute Bias
round((mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B3)-vec[3]), 2)
round((mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B4)-vec[4]), 2)
round((mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B10)-vec[10]), 2)
round((mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B11)-vec[11]), 2)
##Relative Bias
round((mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B3)-vec[3])/vec[3], 2)
round((mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B4)-vec[4])/vec[4], 2)
round((mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B10)-vec[10])/vec[10], 2)
round((mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$B11)-vec[11])/vec[11], 2)

round(sum(MAF[SD]*vec)/sum(MAF[SD]), 2) # Calculation for true AGE
round(mean(SIRT1dat$CMCBeta), 2) # Estimate of True AGE over the RV (causal and non causal)
round(sd(SIRT1dat$CMCBeta), 2)
round(mean(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta), 2) # Naive estimate
round(sd(SIRT1dat[SIRT1dat$pBeta<0.05,]$CMCBeta), 2)
round(mean(with(SIRT1dat,CMCBeta-muBkmed), na.rm = T), 2) # Bootstrap estimate
round(sd(with(SIRT1dat,CMCBeta-muBkmed), na.rm = T), 2)
round(mean(with(SIRT1dat,BT1),na.rm=T), 2) # Likelihood estimate
round(sd(with(SIRT1dat,BT1),na.rm=T), 2)
# Get the MAF
round(MAF[SD], 3)
# Power calculation
powerCalc = length(which(SIRT1dat$pBeta < 0.05))
paste("Power ", powerCalc, "/", JJ, " = ", round(powerCalc/JJ, 2), sep = '')
sink()
###### END HERE (Resimulate Graphs)
#######################
# Unidirectional CAST (pLU values)
#######################
# Select the SIRT1 gene
name = 'SIRT1'
# Read the gene_info file in order to get gene information
geneInfo = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/gene_info',
                    header = TRUE, sep = ',')
# Use which in order to determine the column that matches SIRT1
matchingIndex = which(geneInfo[, 1] == name)
# Get starting/ending position of gene, and chromosome number
startingPos = geneInfo[matchingIndex, 4]
endingPos = geneInfo[matchingIndex, 5]
chrNum = geneInfo[matchingIndex, 2]
# Chromosome file name
chrFileName = paste('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/c',
                    chrNum,'_snps.unr',sep='')
# Read the file 
chrData = read.csv(chrFileName, header = TRUE, sep = ',')
# Take the headers (Chromosome IDs)
chrID = chrData[, 1]
# Call function to get Han Chinese data from the file
subID = getPop(chrID)
# Read snp_info file to get SNP data
snpData = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/snp_info',
                   header=T,sep=',')
# Call function to get subset of SNPs, based on positions of gene
subSNP = subSet(snpData, names(chrData), startingPos, endingPos)
# Convert call matrix into numeric matrix
numericMatrix = convertX(chrData[subID, subSNP])
# Determine the MAF
MAF = getmaf(numericMatrix)
# Display MAF matrix
MAF
# Generates a matrix of Boolean values based on "Rare Variants"
SD = MAF < 0.05 & MAF > 0
SD
# Generate covariate matrix
covariateMatrix = (numericMatrix[, SD])
covariateMatrix
# Take the sum of columns of covariate matrix
apply(covariateMatrix, 2, FUN = sum)
apply(covariateMatrix, 2, FUN = sum)/321
chrNames = names(chrData)
P = NULL
# Apply CMC to the covariate matrix
covariateMatrixCMC = apply(covariateMatrix, 1, FUN = max)
# To make elements of above matrix all 0s and 1s
covariateMatrixCMC[covariateMatrixCMC > 0] =
  covariateMatrixCMC[covariateMatrixCMC > 0]^0
# Set parameters
JJ = 2000
Betaj = NULL
pvalj = NULL
qj = NULL
tBetaj = NULL
Bj = NULL
Tj = NULL
Pj = NULL
SW = NULL
SL = NULL
SWwin = NULL
SLwin = NULL
SC = NULL
SQ = NULL
NBetas = 11
IBetas=matrix(rep(NA,JJ*NBetas),ncol= NBetas)
TIBetas=matrix(rep(NA,JJ*NBetas),ncol= NBetas)
PIBetas=matrix(rep(NA,JJ*NBetas),ncol= NBetas)
Phenos=matrix(rep(NA,JJ*321),ncol= JJ)
# Start loop
set.seed(5)
for (j in 1:JJ) {
  # print(j)
  vec=c(0,0,0.83224,0.97060,0,0,0,0,0,0.93459,0.53073)
  # Create Design
  Design = apply(sweep(covariateMatrix, MARGIN = 2, vec, `*`), 1, FUN = sum)
  # Add normal simulations to the Design matrix
  Y = Design + rnorm(length(Design), 0, 1)
  # Insert elements of Y into Phenos matrix
  Phenos[, j] = Y
  # Assigns covariate CMC of estimate
  Betaj[j] = summary(lm(Y~covariateMatrixCMC))$coefficients[2, 1]
  # Covariate CMC of t value
  tBetaj[j] = summary(lm(Y~covariateMatrixCMC))$coefficients[2, 3]
  # Covariate CMC of Pr(>|t|)
  pvalj[j] = summary(lm(Y~covariateMatrixCMC))$coefficients[2, 4]
  # for # of variants
  for (i in 1:11) {
    # estimate
    IBetas[j,i]=summary(lm(Y~covariateMatrix[,i]))$coeff[2,1]
    # t value
    TIBetas[j,i]=summary(lm(Y~covariateMatrix[,i]))$coeff[2,3]
    # Pr(>|t|)
    PIBetas[j,i]=summary(lm(Y~covariateMatrix[,i]))$coeff[2,4]
  }
  numericMatrix = covariateMatrix
  # Get the number of rows and columns 
  numRows = length(numericMatrix[, 1])
  numCols = length(numericMatrix[1, ])
  # Taking the MAF as the # of occurrences divided by the sample size
  maf = colSums(numericMatrix)/numRows
  # Calculate statistics
  SWwin[j] = as.numeric(testQTLwin(Y,numericMatrix,1/sqrt(maf*(1-maf))))
  SW[j] =    as.numeric(testQTL(Y,numericMatrix,1/sqrt(maf*(1-maf))))
  SL[j] = as.numeric(testQTL(Y,numericMatrix,rep(1,numCols)))
  SLwin[j] = as.numeric(testQTLwin(Y,numericMatrix,rep(1,numCols))) 
  SC[j] = as.numeric(testQTC(Y,numericMatrix))
  SQ[j] = as.numeric(testQTH(Y,numericMatrix))
  # Calculate p values
  p = pvalCalc(Y, covariateMatrix, 10^3)
  # Add p value calculated above to P
  P = cbind(P, p)
  # Print p value
  # cat(p, '\n)
}
# Combine vectors/matrices by the columns
Bj=cbind(Bj,Betaj)
Tj=cbind(Tj,tBetaj)
Pj=cbind(Pj,pvalj)
# Data frame to combine the various estimates
frame=cbind(t(P),SW,SL,SWwin,SLwin,SC,SQ,Bj,Tj,Pj,IBetas)
# Name columns of data frame accordingly
colnames(frame)[c(1,2,3,4,5,6,13:26)]=c('pLW','pLU','pQU','pQH','pMin','pFisher',
                                        'CMCBeta','TCMC','pBeta','B1','B2','B3',
                                        'B4','B5','B6','B7','B8','B9','B10','B11')
# Partially displays elements of above data frame
head(frame)
# Displays # of rows and columns of above data frame
dim(frame)
# Create SIRT1 data based on the data frame
write.csv(frame, '/Users/senadkokic/Desktop/NSERC/GAW17/newsimSIRT1')
# Read the above data as SIRT1 data
newSIRT1 = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/newsimSIRT1')
# Partially display elements of above data (data read correctly)
head(newSIRT1)
# Change column names of Betas to match variants
colnames(IBetas)=c('B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11')
colnames(TIBetas)=c('B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11')
colnames(PIBetas)=c('B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11')
# Make PI and TI Betas into data frames
PIBetas = data.frame(PIBetas)
TIBetas = data.frame(TIBetas)
# Calculates correlation between TCMC and stats
with(newSIRT1,cor(TCMC,SW))
with(newSIRT1,cor(TCMC,SL))
with(newSIRT1,cor(TCMC,SC))
with(newSIRT1,cor(TCMC,SQ))

# In accordance with simulation type 1, create type 1 variables
IBetas1 = IBetas
TIBetas1 = TIBetas
PIBetas1 = PIBetas
newSIRT1.1 = newSIRT1
##### Type 2
# IBetas2 = IBetas
# TIBetas2 = TIBetas
# PIBetas2 = PIBetas
# newSIRT1.2 = newSIRT1
##### Type 3
##### vec=c(0,0,0.83224,0.97060,0,0,0,0,0,-0.93459,0.53073)
# IBetas3=IBetas
# TIBetas3=TIBetas
# PIBetas3=PIBetas
# newSIRT1.3=newSIRT1

################ Begin analysis
# Declare SIRT1dat based on type
SIRT1dat = newSIRT1.1
# SIRT1dat = newSIRT1.2
# SIRT1dat = newSIRT1.3

# Calculate averages of values
sum(SIRT1dat$pQH<0.05)/JJ
sum(SIRT1dat$pQU<0.05)/JJ
sum(SIRT1dat$pLW<0.05)/JJ
sum(SIRT1dat$pBeta<0.05)/JJ
# Averages based on genetic effect
sum(PIBetas$B1<0.05)/JJ
sum(PIBetas$B3<0.05)/JJ
sum(PIBetas$B4<0.05)/JJ
sum(PIBetas$B10<0.05)/JJ
sum(PIBetas$B11<0.05)/JJ
# Calculate estimate of true AGE over both causal and non-causal RVs
mean(SIRT1dat$CMCBeta)
# Mean individual effects
mean(SIRT1dat$B1)
mean(SIRT1dat$B3)
mean(SIRT1dat$B4)
mean(SIRT1dat$B10)
mean(SIRT1dat$B11)
# Calculate estimate of  SD wrt true AGE over both causal and non-causal RVs
sd(SIRT1dat$CMCBeta)
# SD of indivdual effects
sd(SIRT1dat$B1)
sd(SIRT1dat$B3)
sd(SIRT1dat$B4)
sd(SIRT1dat$B10)
sd(SIRT1dat$B11)

################################################################################
# Obtain bias for all of the significant replicates
## Bootstrapping Algorithm (1)
n = length(Y)
numCols = 10000
# Initialize matrix for sample
Bk = matrix(rep(0, n*numCols), ncol = numCols)
# Delegate elements of matrix to be a sample
set.seed(6)
for (k in 1:numCols) {
  sam = 1:321
  Bk[, k] = sample(sam, 321, replace = TRUE)
}
################################################################################
# Bootstrap replication under the Wald Beta test
### SIRT1
# See which indices match criteria
which(SIRT1dat$pBeta<0.05&SIRT1dat$pQH<0.05)
which(SIRT1dat$pQH<0.05)
which(SIRT1dat$pBeta<0.05)
SIRT1dat=newSIRT1.1
#SIRT1dat=newSIRT1.2
#SIRT1dat=newSIRT1.3
# Initialize elements to be NA
SIRT1dat$muBkmed=NA
SIRT1dat$muPkmed=NA
SIRT1dat$muBkmea=NA
SIRT1dat$muPkmea=NA

numericMatrix = convertX(chrData[subID, subSNP])
MAF = getmaf(numericMatrix)
MAF
SD = MAF < 0.05 & MAF > 0
SD
covariateMatrix = (numericMatrix[, SD])
covariateMatrix
apply(covariateMatrix, 2, FUN = sum)
chrNames = names(chrData)
P = NULL
# Apply CMC to the covariate matrix
covariateMatrixCMC = apply(covariateMatrix, 1, FUN = max)
# Loop for elements that match criteria wrt pBeta
for (c in which(SIRT1dat$pLU < 0.05)) {
  # print the index #
  # print(c)
  # Set initial values
  BetaBk=NA
  pBk=NA
  BetaCk=NA
  pCk=NA
  pH=NULL
  pt=NULL
  sigK=0
  for (k in 1:numCols) {
    # Print iteration as well as how many significant results have been seen
    # NOTE: this prints REALLY fast
    # print(paste(k, ' - ', sigK))
    # Break when 200 are observed (? not sure why)
    if (sigK == 200) break
    # Bootstrap sample
    BK = Bk[, k]
    # Remaining elements (left out of bootstrap)
    CK = sam[-BK]
    # Set Y = Phenos element based on index
    Y = Phenos[, c]
    # Compute linear model coefficient for Wald Statistic
    pt[k] = summary(lm(Y[BK]~covariateMatrixCMC[BK]))$coefficients[2,4]
    pt[k]
    
    # Significant result if the Wald stat < 0.05
    if (pt[k] < 0.05) {
      sigK = sigK + 1
      # Take the means of CMC matrix based on those sampled and unsampled
      pBk[k] = mean(covariateMatrixCMC[BK])
      if (pBk[k] > 0) {
        BetaBk[k]=summary(lm(Y[BK]~covariateMatrixCMC[BK]))$coefficients[2,1]
      }
      pCk[k] = mean(covariateMatrixCMC[CK])
      if (pCk[k] > 0) {
        BetaCk[k]=summary(lm(Y[CK]~covariateMatrixCMC[CK]))$coefficients[2,1]
      }
    }
  }
  # Insert differences between estimates into SIRT1dat
  # Median method
  muBkmed=median(BetaBk-BetaCk,na.rm = T)
  muPkmed=median(pBk-pCk,na.rm = T)
  SIRT1dat$muBkmed[c]=muBkmed
  SIRT1dat$muPkmed[c]=muPkmed
  # Mean method
  muBkmea=mean(BetaBk-BetaCk,na.rm = T)
  muPkmea=mean(pBk-pCk,na.rm = T)
  SIRT1dat$muBkmea[c]=muBkmea
  SIRT1dat$muPkmea[c]=muPkmea
}
# Writing files
# write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot',name,sep=''))
# bootSIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot',name,sep=''))

write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot1',name,sep=''))
boot1SIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot1',name,sep=''))

# write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot2',name,sep=''))
# boot2SIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot2',name,sep=''))
head(boot1SIRT1, 100)
### Apply Ghosh et. al. to the significant replicates
SIRT1dat$BT1=NA
SIRT1dat$BT2=NA
SIRT1dat$BT3=NA
# Bonferroni threshold for achieving study-wide significance = 0.1
# Named cee in order to avoid confusion b/w function c()
cee = abs(qnorm(0.5*0.05))
# Begin loop for pBeta 
for (cc in which(SIRT1dat$pLU < 0.05)) {
  # Print the index
  # print(cc)
  # Proposed odds ratio = 1.54
  betahat = SIRT1dat[cc, ]$CMCBeta
  # Reported p-value = 5.7e-4, which we convert to a z-value
  z=sign(betahat)*abs(qnorm(0.5*SIRT1dat[cc,]$pLU))
  ###################################################
  #              THE PROPOSED APPROACH              #
  ###################################################
  # Compute the standard error of betahat
  se = betahat/z
  # Compute the conditional MLE
  mutilde1 = optimize(f=conditional.like,
                      c(-20,20),maximum=T,z=z,cee=cee)$maximum
  # Compute the conditional MLE on the beta scale
  betatilde1 = mutilde1*se
  # Assign the conditional MLE of the odds ratio estimate 
  SIRT1dat[cc, ]$BT1 = betatilde1
  # Print out the estimate 
  SIRT1dat[cc, ]$BT1
  # Numeric integration
  a=sum(apply(as.matrix(seq(-50,50,by=0.01)),1,
              conditional.like.z,cee=cee,z=z))*0.01
  b=sum(apply(as.matrix(seq(-50,50,by=0.01)),1,
              conditional.like,cee=cee,z=z))*0.01
  # Alternate computation of odds ratio
  betatilde2 = (a/b)*se
  # Assign the estimate, print it
  SIRT1dat[cc,]$BT2=betatilde2
  # Compute the average of the previous ORs as the third estimate
  betatilde3 = (betatilde1 + betatilde2)/2
  # Assign the estimate, print it
  SIRT1dat[cc,]$BT3=betatilde3
}
write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot',
                         name,sep=''))
bootSIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot',
                         name,sep=''))

write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot1',
                         name,sep=''))
boot1SIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot1',
                          name,sep=''))
head(boot1SIRT1,100)

write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot2',
                         name,sep=''))
boot2SIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot2',
                          name,sep=''))
##### Histogram comparisons:
SIRT1dat = bootSIRT1
IBetas=IBetas1
TIBetas=TIBetas1
PIBetas=PIBetas1
###############################################################################
sink('/Users/senadkokic/Desktop/NSERC/Results Files/Values_CAST_Uni')
# For the Wald test (Table 1)
##Truth
round(vec[c(3,4,10,11)], 2)

##Truth over replicates (JJ=2000)
# These give B_k values for table
round(mean(SIRT1dat$B3), 2)
round(mean(SIRT1dat$B4), 2)
round(mean(SIRT1dat$B10), 2)
round(mean(SIRT1dat$B11), 2)
##Bias-Winner's Curse
# These give B_naive values for table
round(mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B3), 2)
round(mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B4), 2)
round(mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B10), 2)
round(mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B11), 2)
##Absolute Bias
round((mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B3)-vec[3]), 2)
round((mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B4)-vec[4]), 2)
round((mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B10)-vec[10]), 2)
round((mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B11)-vec[11]), 2)
##Relative Bias
round((mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B3)-vec[3])/vec[3], 2)
round((mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B4)-vec[4])/vec[4], 2)
round((mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B10)-vec[10])/vec[10], 2)
round((mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B11)-vec[11])/vec[11], 2)
# Get the MAF
round(MAF[SD], 3)
# Power calculation
powerCalc = length(which(SIRT1dat$pLU < 0.05))
paste("Power ", powerCalc, "/", JJ, " = ", round(powerCalc/JJ, 2), sep = '')
sink()
#######################
# Bidirectional CAST
#######################
# Select the SIRT1 gene
name = 'SIRT1'
# Read the gene_info file in order to get gene information
geneInfo = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/gene_info',
                    header = TRUE, sep = ',')
# Use which in order to determine the column that matches SIRT1
matchingIndex = which(geneInfo[, 1] == name)
# Get starting/ending position of gene, and chromosome number
startingPos = geneInfo[matchingIndex, 4]
endingPos = geneInfo[matchingIndex, 5]
chrNum = geneInfo[matchingIndex, 2]
# Chromosome file name
chrFileName = paste('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/c',
                    chrNum,'_snps.unr',sep='')
# Read the file 
chrData = read.csv(chrFileName, header = TRUE, sep = ',')
# Take the headers (Chromosome IDs)
chrID = chrData[, 1]
# Call function to get Han Chinese data from the file
subID = getPop(chrID)
# Read snp_info file to get SNP data
snpData = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/snp_info',
                   header=T,sep=',')
# Call function to get subset of SNPs, based on positions of gene
subSNP = subSet(snpData, names(chrData), startingPos, endingPos)
# Convert call matrix into numeric matrix
numericMatrix = convertX(chrData[subID, subSNP])
# Determine the MAF
MAF = getmaf(numericMatrix)
# Display MAF matrix
MAF
# Generates a matrix of Boolean values based on "Rare Variants"
SD = MAF < 0.05 & MAF > 0
SD
# Generate covariate matrix
covariateMatrix = (numericMatrix[, SD])
covariateMatrix
# Take the sum of columns of covariate matrix
apply(covariateMatrix, 2, FUN = sum)
apply(covariateMatrix, 2, FUN = sum)/321
chrNames = names(chrData)
P = NULL
# Apply CMC to the covariate matrix
covariateMatrixCMC = apply(covariateMatrix, 1, FUN = max)
# To make elements of above matrix all 0s and 1s
covariateMatrixCMC[covariateMatrixCMC > 0] =
  covariateMatrixCMC[covariateMatrixCMC > 0]^0
# Set parameters
JJ = 2000
Betaj = NULL
pvalj = NULL
qj = NULL
tBetaj = NULL
Bj = NULL
Tj = NULL
Pj = NULL
SW = NULL
SL = NULL
SWwin = NULL
SLwin = NULL
SC = NULL
SQ = NULL
NBetas = 11
IBetas=matrix(rep(NA,JJ*NBetas),ncol= NBetas)
TIBetas=matrix(rep(NA,JJ*NBetas),ncol= NBetas)
PIBetas=matrix(rep(NA,JJ*NBetas),ncol= NBetas)
Phenos=matrix(rep(NA,JJ*321),ncol= JJ)
# Start loop
set.seed(7)
for (j in 1:JJ) {
  # print(j)
  vec=c(0,0,0.83224,0.97060,0,0,0,0,0,-0.93459,0.53073)
  # Create Design
  Design = apply(sweep(covariateMatrix, MARGIN = 2, vec, `*`), 1, FUN = sum)
  # Add normal simulations to the Design matrix
  Y = Design + rnorm(length(Design), 0, 1)
  # Insert elements of Y into Phenos matrix
  Phenos[, j] = Y
  # Assigns covariate CMC of estimate
  Betaj[j] = summary(lm(Y~covariateMatrixCMC))$coefficients[2, 1]
  # Covariate CMC of t value
  tBetaj[j] = summary(lm(Y~covariateMatrixCMC))$coefficients[2, 3]
  # Covariate CMC of Pr(>|t|)
  pvalj[j] = summary(lm(Y~covariateMatrixCMC))$coefficients[2, 4]
  # for # of variants
  for (i in 1:11) {
    # estimate
    IBetas[j,i]=summary(lm(Y~covariateMatrix[,i]))$coeff[2,1]
    # t value
    TIBetas[j,i]=summary(lm(Y~covariateMatrix[,i]))$coeff[2,3]
    # Pr(>|t|)
    PIBetas[j,i]=summary(lm(Y~covariateMatrix[,i]))$coeff[2,4]
  }
  numericMatrix = covariateMatrix
  # Get the number of rows and columns 
  numRows = length(numericMatrix[, 1])
  numCols = length(numericMatrix[1, ])
  # Taking the MAF as the # of occurrences divided by the sample size
  maf = colSums(numericMatrix)/numRows
  # Calculate statistics
  SWwin[j] = as.numeric(testQTLwin(Y,numericMatrix,1/sqrt(maf*(1-maf))))
  SW[j] =    as.numeric(testQTL(Y,numericMatrix,1/sqrt(maf*(1-maf))))
  SL[j] = as.numeric(testQTL(Y,numericMatrix,rep(1,numCols)))
  SLwin[j] = as.numeric(testQTLwin(Y,numericMatrix,rep(1,numCols))) 
  SC[j] = as.numeric(testQTC(Y,numericMatrix))
  SQ[j] = as.numeric(testQTH(Y,numericMatrix))
  # Calculate p values
  p = pvalCalc(Y, covariateMatrix, 10^3)
  # Add p value calculated above to P
  P = cbind(P, p)
  # Print p value
  # cat(p, '\n)
}
# Combine vectors/matrices by the columns
Bj=cbind(Bj,Betaj)
Tj=cbind(Tj,tBetaj)
Pj=cbind(Pj,pvalj)
# Data frame to combine the various estimates
frame=cbind(t(P),SW,SL,SWwin,SLwin,SC,SQ,Bj,Tj,Pj,IBetas)
# Name columns of data frame accordingly
colnames(frame)[c(1,2,3,4,5,6,13:26)]=c('pLW','pLU','pQU','pQH','pMin','pFisher',
                                        'CMCBeta','TCMC','pBeta','B1','B2','B3',
                                        'B4','B5','B6','B7','B8','B9','B10','B11')
# Partially displays elements of above data frame
head(frame)
# Displays # of rows and columns of above data frame
dim(frame)
# Create SIRT1 data based on the data frame
write.csv(frame, '/Users/senadkokic/Desktop/NSERC/GAW17/newsimSIRT1')
# Read the above data as SIRT1 data
newSIRT1 = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/newsimSIRT1')
# Partially display elements of above data (data read correctly)
head(newSIRT1)
# Change column names of Betas to match variants
colnames(IBetas)=c('B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11')
colnames(TIBetas)=c('B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11')
colnames(PIBetas)=c('B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11')
# Make PI and TI Betas into data frames
PIBetas = data.frame(PIBetas)
TIBetas = data.frame(TIBetas)
# Calculates correlation between TCMC and stats
with(newSIRT1,cor(TCMC,SW))
with(newSIRT1,cor(TCMC,SL))
with(newSIRT1,cor(TCMC,SC))
with(newSIRT1,cor(TCMC,SQ))

# In accordance with simulation type 1, create type 1 variables
IBetas1 = IBetas
TIBetas1 = TIBetas
PIBetas1 = PIBetas
newSIRT1.1 = newSIRT1
##### Type 2
# IBetas2 = IBetas
# TIBetas2 = TIBetas
# PIBetas2 = PIBetas
# newSIRT1.2 = newSIRT1
##### Type 3
##### vec=c(0,0,0.83224,0.97060,0,0,0,0,0,-0.93459,0.53073)
# IBetas3=IBetas
# TIBetas3=TIBetas
# PIBetas3=PIBetas
# newSIRT1.3=newSIRT1

################ Begin analysis
# Declare SIRT1dat based on type
SIRT1dat = newSIRT1.1
################################################################################
# Obtain bias for all of the significant replicates
## Bootstrapping Algorithm (1)
n = length(Y)
numCols = 10000
# Initialize matrix for sample
Bk = matrix(rep(0, n*numCols), ncol = numCols)
# Delegate elements of matrix to be a sample
set.seed(8)
for (k in 1:numCols) {
  sam = 1:321
  Bk[, k] = sample(sam, 321, replace = TRUE)
}
################################################################################
# Bootstrap replication under the Wald Beta test
### SIRT1
# See which indices match criteria
which(SIRT1dat$pBeta<0.05&SIRT1dat$pQH<0.05)
which(SIRT1dat$pQH<0.05)
which(SIRT1dat$pBeta<0.05)
SIRT1dat=newSIRT1.1
#SIRT1dat=newSIRT1.2
#SIRT1dat=newSIRT1.3
# Initialize elements to be NA
SIRT1dat$muBkmed=NA
SIRT1dat$muPkmed=NA
SIRT1dat$muBkmea=NA
SIRT1dat$muPkmea=NA

numericMatrix = convertX(chrData[subID, subSNP])
MAF = getmaf(numericMatrix)
MAF
SD = MAF < 0.05 & MAF > 0
SD
covariateMatrix = (numericMatrix[, SD])
covariateMatrix
apply(covariateMatrix, 2, FUN = sum)
chrNames = names(chrData)
P = NULL
# Apply CMC to the covariate matrix
covariateMatrixCMC = apply(covariateMatrix, 1, FUN = max)
# Loop for elements that match criteria wrt pBeta
for (c in which(SIRT1dat$pLU < 0.05)) {
  # print the index #
  # print(c)
  # Set initial values
  BetaBk=NA
  pBk=NA
  BetaCk=NA
  pCk=NA
  pH=NULL
  pt=NULL
  sigK=0
  for (k in 1:numCols) {
    # Print iteration as well as how many significant results have been seen
    # NOTE: this prints REALLY fast
    # print(paste(k, ' - ', sigK))
    # Break when 200 are observed (? not sure why)
    if (sigK == 200) break
    # Bootstrap sample
    BK = Bk[, k]
    # Remaining elements (left out of bootstrap)
    CK = sam[-BK]
    # Set Y = Phenos element based on index
    Y = Phenos[, c]
    # Compute linear model coefficient for Wald Statistic
    pt[k] = summary(lm(Y[BK]~covariateMatrixCMC[BK]))$coefficients[2,4]
    pt[k]
    
    # Significant result if the Wald stat < 0.05
    if (pt[k] < 0.05) {
      sigK = sigK + 1
      # Take the means of CMC matrix based on those sampled and unsampled
      pBk[k] = mean(covariateMatrixCMC[BK])
      if (pBk[k] > 0) {
        BetaBk[k]=summary(lm(Y[BK]~covariateMatrixCMC[BK]))$coefficients[2,1]
      }
      pCk[k] = mean(covariateMatrixCMC[CK])
      if (pCk[k] > 0) {
        BetaCk[k]=summary(lm(Y[CK]~covariateMatrixCMC[CK]))$coefficients[2,1]
      }
    }
  }
  # Insert differences between estimates into SIRT1dat
  # Median method
  muBkmed=median(BetaBk-BetaCk,na.rm = T)
  muPkmed=median(pBk-pCk,na.rm = T)
  SIRT1dat$muBkmed[c]=muBkmed
  SIRT1dat$muPkmed[c]=muPkmed
  # Mean method
  muBkmea=mean(BetaBk-BetaCk,na.rm = T)
  muPkmea=mean(pBk-pCk,na.rm = T)
  SIRT1dat$muBkmea[c]=muBkmea
  SIRT1dat$muPkmea[c]=muPkmea
}
# Writing files
# write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot',name,sep=''))
# bootSIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot',name,sep=''))

write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot1',name,sep=''))
boot1SIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot1',name,sep=''))

# write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot2',name,sep=''))
# boot2SIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot2',name,sep=''))
head(boot1SIRT1, 100)
### Apply Ghosh et. al. to the significant replicates
SIRT1dat$BT1=NA
SIRT1dat$BT2=NA
SIRT1dat$BT3=NA
# Bonferroni threshold for achieving study-wide significance = 0.1
# Named cee in order to avoid confusion b/w function c()
cee = abs(qnorm(0.5*0.05))
# Begin loop for pBeta 
for (cc in which(SIRT1dat$pLU < 0.05)) {
  # Print the index
  # print(cc)
  # Proposed odds ratio = 1.54
  betahat = SIRT1dat[cc, ]$CMCBeta
  # Reported p-value = 5.7e-4, which we convert to a z-value
  z=sign(betahat)*abs(qnorm(0.5*SIRT1dat[cc,]$pLU))
  ###################################################
  #              THE PROPOSED APPROACH              #
  ###################################################
  # Compute the standard error of betahat
  se = betahat/z
  # Compute the conditional MLE
  mutilde1 = optimize(f=conditional.like,
                      c(-20,20),maximum=T,z=z,cee=cee)$maximum
  # Compute the conditional MLE on the beta scale
  betatilde1 = mutilde1*se
  # Assign the conditional MLE of the odds ratio estimate 
  SIRT1dat[cc, ]$BT1 = betatilde1
  # Print out the estimate 
  SIRT1dat[cc, ]$BT1
  # Numeric integration
  a=sum(apply(as.matrix(seq(-50,50,by=0.01)),1,
              conditional.like.z,cee=cee,z=z))*0.01
  b=sum(apply(as.matrix(seq(-50,50,by=0.01)),1,
              conditional.like,cee=cee,z=z))*0.01
  # Alternate computation of odds ratio
  betatilde2 = (a/b)*se
  # Assign the estimate, print it
  SIRT1dat[cc,]$BT2=betatilde2
  # Compute the average of the previous ORs as the third estimate
  betatilde3 = (betatilde1 + betatilde2)/2
  # Assign the estimate, print it
  SIRT1dat[cc,]$BT3=betatilde3
}
write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot',
                         name,sep=''))
bootSIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot',
                         name,sep=''))

write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot1',
                         name,sep=''))
boot1SIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot1',
                          name,sep=''))
head(boot1SIRT1,100)

write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot2',
                         name,sep=''))
boot2SIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot2',
                          name,sep=''))
##### Histogram comparisons:
SIRT1dat = bootSIRT1
IBetas=IBetas1
TIBetas=TIBetas1
PIBetas=PIBetas1
###############################################################################
sink('/Users/senadkokic/Desktop/NSERC/Results Files/Values_CAST_Bi')
# For the Wald test (Table 1)
##Truth
round(vec[c(3,4,10,11)], 2)

##Truth over replicates (JJ=2000)
# These give B_k values for table
round(mean(SIRT1dat$B3), 2)
round(mean(SIRT1dat$B4), 2)
round(mean(SIRT1dat$B10), 2)
round(mean(SIRT1dat$B11), 2)
##Bias-Winner's Curse
# These give B_naive values for table
round(mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B3), 2)
round(mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B4), 2)
round(mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B10), 2)
round(mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B11), 2)
##Absolute Bias
round((mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B3)-vec[3]), 2)
round((mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B4)-vec[4]), 2)
round((mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B10)-vec[10]), 2)
round((mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B11)-vec[11]), 2)
##Relative Bias
round((mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B3)-vec[3])/vec[3], 2)
round((mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B4)-vec[4])/vec[4], 2)
round((mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B10)-vec[10])/vec[10], 2)
round((mean(SIRT1dat[SIRT1dat$pLU<0.05,]$B11)-vec[11])/vec[11], 2)
# Get the MAF
round(MAF[SD], 3)
# Power calculation
powerCalc = length(which(SIRT1dat$pLU < 0.05))
paste("Power ", powerCalc, "/", JJ, " = ", round(powerCalc/JJ, 2), sep = '')
sink()
#######################
# Unidirectional C-Alpha (pQU values)
#######################
# Select the SIRT1 gene
name = 'SIRT1'
# Read the gene_info file in order to get gene information
geneInfo = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/gene_info',
                    header = TRUE, sep = ',')
# Use which in order to determine the column that matches SIRT1
matchingIndex = which(geneInfo[, 1] == name)
# Get starting/ending position of gene, and chromosome number
startingPos = geneInfo[matchingIndex, 4]
endingPos = geneInfo[matchingIndex, 5]
chrNum = geneInfo[matchingIndex, 2]
# Chromosome file name
chrFileName = paste('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/c',
                    chrNum,'_snps.unr',sep='')
# Read the file 
chrData = read.csv(chrFileName, header = TRUE, sep = ',')
# Take the headers (Chromosome IDs)
chrID = chrData[, 1]
# Call function to get Han Chinese data from the file
subID = getPop(chrID)
# Read snp_info file to get SNP data
snpData = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/snp_info',
                   header=T,sep=',')
# Call function to get subset of SNPs, based on positions of gene
subSNP = subSet(snpData, names(chrData), startingPos, endingPos)
# Convert call matrix into numeric matrix
numericMatrix = convertX(chrData[subID, subSNP])
# Determine the MAF
MAF = getmaf(numericMatrix)
# Display MAF matrix
MAF
# Generates a matrix of Boolean values based on "Rare Variants"
SD = MAF < 0.05 & MAF > 0
SD
# Generate covariate matrix
covariateMatrix = (numericMatrix[, SD])
covariateMatrix
# Take the sum of columns of covariate matrix
apply(covariateMatrix, 2, FUN = sum)
apply(covariateMatrix, 2, FUN = sum)/321
chrNames = names(chrData)
P = NULL
# Apply CMC to the covariate matrix
covariateMatrixCMC = apply(covariateMatrix, 1, FUN = max)
# To make elements of above matrix all 0s and 1s
covariateMatrixCMC[covariateMatrixCMC > 0] =
  covariateMatrixCMC[covariateMatrixCMC > 0]^0
# Set parameters
JJ = 2000
Betaj = NULL
pvalj = NULL
qj = NULL
tBetaj = NULL
Bj = NULL
Tj = NULL
Pj = NULL
SW = NULL
SL = NULL
SWwin = NULL
SLwin = NULL
SC = NULL
SQ = NULL
NBetas = 11
IBetas=matrix(rep(NA,JJ*NBetas),ncol= NBetas)
TIBetas=matrix(rep(NA,JJ*NBetas),ncol= NBetas)
PIBetas=matrix(rep(NA,JJ*NBetas),ncol= NBetas)
Phenos=matrix(rep(NA,JJ*321),ncol= JJ)
# Start loop
set.seed(9)
for (j in 1:JJ) {
  # print(j)
  vec=c(0,0,0.83224,0.97060,0,0,0,0,0,0.93459,0.53073)
  # Create Design
  Design = apply(sweep(covariateMatrix, MARGIN = 2, vec, `*`), 1, FUN = sum)
  # Add normal simulations to the Design matrix
  Y = Design + rnorm(length(Design), 0, 1)
  # Insert elements of Y into Phenos matrix
  Phenos[, j] = Y
  # Assigns covariate CMC of estimate
  Betaj[j] = summary(lm(Y~covariateMatrixCMC))$coefficients[2, 1]
  # Covariate CMC of t value
  tBetaj[j] = summary(lm(Y~covariateMatrixCMC))$coefficients[2, 3]
  # Covariate CMC of Pr(>|t|)
  pvalj[j] = summary(lm(Y~covariateMatrixCMC))$coefficients[2, 4]
  # for # of variants
  for (i in 1:11) {
    # estimate
    IBetas[j,i]=summary(lm(Y~covariateMatrix[,i]))$coeff[2,1]
    # t value
    TIBetas[j,i]=summary(lm(Y~covariateMatrix[,i]))$coeff[2,3]
    # Pr(>|t|)
    PIBetas[j,i]=summary(lm(Y~covariateMatrix[,i]))$coeff[2,4]
  }
  numericMatrix = covariateMatrix
  # Get the number of rows and columns 
  numRows = length(numericMatrix[, 1])
  numCols = length(numericMatrix[1, ])
  # Taking the MAF as the # of occurrences divided by the sample size
  maf = colSums(numericMatrix)/numRows
  # Calculate statistics
  SWwin[j] = as.numeric(testQTLwin(Y,numericMatrix,1/sqrt(maf*(1-maf))))
  SW[j] =    as.numeric(testQTL(Y,numericMatrix,1/sqrt(maf*(1-maf))))
  SL[j] = as.numeric(testQTL(Y,numericMatrix,rep(1,numCols)))
  SLwin[j] = as.numeric(testQTLwin(Y,numericMatrix,rep(1,numCols))) 
  SC[j] = as.numeric(testQTC(Y,numericMatrix))
  SQ[j] = as.numeric(testQTH(Y,numericMatrix))
  # Calculate p values
  p = pvalCalc(Y, covariateMatrix, 10^3)
  # Add p value calculated above to P
  P = cbind(P, p)
  # Print p value
  # cat(p, '\n)
}
# Combine vectors/matrices by the columns
Bj=cbind(Bj,Betaj)
Tj=cbind(Tj,tBetaj)
Pj=cbind(Pj,pvalj)
# Data frame to combine the various estimates
frame=cbind(t(P),SW,SL,SWwin,SLwin,SC,SQ,Bj,Tj,Pj,IBetas)
# Name columns of data frame accordingly
colnames(frame)[c(1,2,3,4,5,6,13:26)]=c('pLW','pLU','pQU','pQH','pMin','pFisher',
                                        'CMCBeta','TCMC','pBeta','B1','B2','B3',
                                        'B4','B5','B6','B7','B8','B9','B10','B11')
# Partially displays elements of above data frame
head(frame)
# Displays # of rows and columns of above data frame
dim(frame)
# Create SIRT1 data based on the data frame
write.csv(frame, '/Users/senadkokic/Desktop/NSERC/GAW17/newsimSIRT1')
# Read the above data as SIRT1 data
newSIRT1 = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/newsimSIRT1')
# Partially display elements of above data (data read correctly)
head(newSIRT1)
# Change column names of Betas to match variants
colnames(IBetas)=c('B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11')
colnames(TIBetas)=c('B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11')
colnames(PIBetas)=c('B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11')
# Make PI and TI Betas into data frames
PIBetas = data.frame(PIBetas)
TIBetas = data.frame(TIBetas)
# Calculates correlation between TCMC and stats
with(newSIRT1,cor(TCMC,SW))
with(newSIRT1,cor(TCMC,SL))
with(newSIRT1,cor(TCMC,SC))
with(newSIRT1,cor(TCMC,SQ))

# In accordance with simulation type 1, create type 1 variables
IBetas1 = IBetas
TIBetas1 = TIBetas
PIBetas1 = PIBetas
newSIRT1.1 = newSIRT1
##### Type 2
# IBetas2 = IBetas
# TIBetas2 = TIBetas
# PIBetas2 = PIBetas
# newSIRT1.2 = newSIRT1
##### Type 3
##### vec=c(0,0,0.83224,0.97060,0,0,0,0,0,-0.93459,0.53073)
# IBetas3=IBetas
# TIBetas3=TIBetas
# PIBetas3=PIBetas
# newSIRT1.3=newSIRT1

################ Begin analysis
# Declare SIRT1dat based on type
SIRT1dat = newSIRT1.1

################################################################################
# Obtain bias for all of the significant replicates
## Bootstrapping Algorithm (1)
n = length(Y)
numCols = 10000
# Initialize matrix for sample
Bk = matrix(rep(0, n*numCols), ncol = numCols)
# Delegate elements of matrix to be a sample
set.seed(10)
for (k in 1:numCols) {
  sam = 1:321
  Bk[, k] = sample(sam, 321, replace = TRUE)
}
################################################################################
# Bootstrap replication under the Wald Beta test
### SIRT1
# See which indices match criteria
which(SIRT1dat$pBeta<0.05&SIRT1dat$pQH<0.05)
which(SIRT1dat$pQH<0.05)
which(SIRT1dat$pBeta<0.05)
SIRT1dat=newSIRT1.1
#SIRT1dat=newSIRT1.2
#SIRT1dat=newSIRT1.3
# Initialize elements to be NA
SIRT1dat$muBkmed=NA
SIRT1dat$muPkmed=NA
SIRT1dat$muBkmea=NA
SIRT1dat$muPkmea=NA

numericMatrix = convertX(chrData[subID, subSNP])
MAF = getmaf(numericMatrix)
MAF
SD = MAF < 0.05 & MAF > 0
SD
covariateMatrix = (numericMatrix[, SD])
covariateMatrix
apply(covariateMatrix, 2, FUN = sum)
chrNames = names(chrData)
P = NULL
# Apply CMC to the covariate matrix
covariateMatrixCMC = apply(covariateMatrix, 1, FUN = max)
# Loop for elements that match criteria wrt pBeta
for (c in which(SIRT1dat$pQU < 0.05)) {
  # print the index #
  # print(c)
  # Set initial values
  BetaBk=NA
  pBk=NA
  BetaCk=NA
  pCk=NA
  pH=NULL
  pt=NULL
  sigK=0
  for (k in 1:numCols) {
    # Print iteration as well as how many significant results have been seen
    # NOTE: this prints REALLY fast
    # print(paste(k, ' - ', sigK))
    # Break when 200 are observed (? not sure why)
    if (sigK == 200) break
    # Bootstrap sample
    BK = Bk[, k]
    # Remaining elements (left out of bootstrap)
    CK = sam[-BK]
    # Set Y = Phenos element based on index
    Y = Phenos[, c]
    # Compute linear model coefficient for Wald Statistic
    pt[k] = summary(lm(Y[BK]~covariateMatrixCMC[BK]))$coefficients[2,4]
    pt[k]
    
    # Significant result if the Wald stat < 0.05
    if (pt[k] < 0.05) {
      sigK = sigK + 1
      # Take the means of CMC matrix based on those sampled and unsampled
      pBk[k] = mean(covariateMatrixCMC[BK])
      if (pBk[k] > 0) {
        BetaBk[k]=summary(lm(Y[BK]~covariateMatrixCMC[BK]))$coefficients[2,1]
      }
      pCk[k] = mean(covariateMatrixCMC[CK])
      if (pCk[k] > 0) {
        BetaCk[k]=summary(lm(Y[CK]~covariateMatrixCMC[CK]))$coefficients[2,1]
      }
    }
  }
  # Insert differences between estimates into SIRT1dat
  # Median method
  muBkmed=median(BetaBk-BetaCk,na.rm = T)
  muPkmed=median(pBk-pCk,na.rm = T)
  SIRT1dat$muBkmed[c]=muBkmed
  SIRT1dat$muPkmed[c]=muPkmed
  # Mean method
  muBkmea=mean(BetaBk-BetaCk,na.rm = T)
  muPkmea=mean(pBk-pCk,na.rm = T)
  SIRT1dat$muBkmea[c]=muBkmea
  SIRT1dat$muPkmea[c]=muPkmea
}
# Writing files
# write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot',name,sep=''))
# bootSIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot',name,sep=''))

write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot1',name,sep=''))
boot1SIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot1',name,sep=''))

# write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot2',name,sep=''))
# boot2SIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot2',name,sep=''))
head(boot1SIRT1, 100)
### Apply Ghosh et. al. to the significant replicates
SIRT1dat$BT1=NA
SIRT1dat$BT2=NA
SIRT1dat$BT3=NA
# Bonferroni threshold for achieving study-wide significance = 0.1
# Named cee in order to avoid confusion b/w function c()
cee = abs(qnorm(0.5*0.05))
# Begin loop for pBeta 
for (cc in which(SIRT1dat$pQU < 0.05)) {
  # Print the index
  # print(cc)
  # Proposed odds ratio = 1.54
  betahat = SIRT1dat[cc, ]$CMCBeta
  # Reported p-value = 5.7e-4, which we convert to a z-value
  z=sign(betahat)*abs(qnorm(0.5*SIRT1dat[cc,]$pQU))
  ###################################################
  #              THE PROPOSED APPROACH              #
  ###################################################
  # Compute the standard error of betahat
  se = betahat/z
  # Compute the conditional MLE
  mutilde1 = optimize(f=conditional.like,
                      c(-20,20),maximum=T,z=z,cee=cee)$maximum
  # Compute the conditional MLE on the beta scale
  betatilde1 = mutilde1*se
  # Assign the conditional MLE of the odds ratio estimate 
  SIRT1dat[cc, ]$BT1 = betatilde1
  # Print out the estimate 
  SIRT1dat[cc, ]$BT1
  # Numeric integration
  a=sum(apply(as.matrix(seq(-50,50,by=0.01)),1,
              conditional.like.z,cee=cee,z=z))*0.01
  b=sum(apply(as.matrix(seq(-50,50,by=0.01)),1,
              conditional.like,cee=cee,z=z))*0.01
  # Alternate computation of odds ratio
  betatilde2 = (a/b)*se
  # Assign the estimate, print it
  SIRT1dat[cc,]$BT2=betatilde2
  # Compute the average of the previous ORs as the third estimate
  betatilde3 = (betatilde1 + betatilde2)/2
  # Assign the estimate, print it
  SIRT1dat[cc,]$BT3=betatilde3
}
write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot',
                         name,sep=''))
bootSIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot',
                         name,sep=''))

write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot1',
                         name,sep=''))
boot1SIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot1',
                          name,sep=''))
head(boot1SIRT1,100)

write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot2',
                         name,sep=''))
boot2SIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot2',
                          name,sep=''))
##### Histogram comparisons:
SIRT1dat = bootSIRT1
IBetas=IBetas1
TIBetas=TIBetas1
PIBetas=PIBetas1
###############################################################################
sink('/Users/senadkokic/Desktop/NSERC/Results Files/Values_CAlpha_Uni')
# For the Wald test (Table 1)
##Truth
round(vec[c(3,4,10,11)], 2)

##Truth over replicates (JJ=2000)
# These give B_k values for table
round(mean(SIRT1dat$B3), 2)
round(mean(SIRT1dat$B4), 2)
round(mean(SIRT1dat$B10), 2)
round(mean(SIRT1dat$B11), 2)
##Bias-Winner's Curse
# These give B_naive values for table
round(mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B3), 2)
round(mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B4), 2)
round(mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B10), 2)
round(mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B11), 2)
##Absolute Bias
round((mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B3)-vec[3]), 2)
round((mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B4)-vec[4]), 2)
round((mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B10)-vec[10]), 2)
round((mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B11)-vec[11]), 2)
##Relative Bias
round((mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B3)-vec[3])/vec[3], 2)
round((mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B4)-vec[4])/vec[4], 2)
round((mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B10)-vec[10])/vec[10], 2)
round((mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B11)-vec[11])/vec[11], 2)
# Get the MAF
round(MAF[SD], 3)
# Power calculation
powerCalc = length(which(SIRT1dat$pQU < 0.05))
paste("Power ", powerCalc, "/", JJ, " = ", round(powerCalc/JJ, 2), sep = '')
sink()
#######################
# Bidirectional C-Alpha
#######################
# Select the SIRT1 gene
name = 'SIRT1'
# Read the gene_info file in order to get gene information
geneInfo = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/gene_info',
                    header = TRUE, sep = ',')
# Use which in order to determine the column that matches SIRT1
matchingIndex = which(geneInfo[, 1] == name)
# Get starting/ending position of gene, and chromosome number
startingPos = geneInfo[matchingIndex, 4]
endingPos = geneInfo[matchingIndex, 5]
chrNum = geneInfo[matchingIndex, 2]
# Chromosome file name
chrFileName = paste('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/c',
                    chrNum,'_snps.unr',sep='')
# Read the file 
chrData = read.csv(chrFileName, header = TRUE, sep = ',')
# Take the headers (Chromosome IDs)
chrID = chrData[, 1]
# Call function to get Han Chinese data from the file
subID = getPop(chrID)
# Read snp_info file to get SNP data
snpData = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/snp_info',
                   header=T,sep=',')
# Call function to get subset of SNPs, based on positions of gene
subSNP = subSet(snpData, names(chrData), startingPos, endingPos)
# Convert call matrix into numeric matrix
numericMatrix = convertX(chrData[subID, subSNP])
# Determine the MAF
MAF = getmaf(numericMatrix)
# Display MAF matrix
MAF
# Generates a matrix of Boolean values based on "Rare Variants"
SD = MAF < 0.05 & MAF > 0
SD
# Generate covariate matrix
covariateMatrix = (numericMatrix[, SD])
covariateMatrix
# Take the sum of columns of covariate matrix
apply(covariateMatrix, 2, FUN = sum)
apply(covariateMatrix, 2, FUN = sum)/321
chrNames = names(chrData)
P = NULL
# Apply CMC to the covariate matrix
covariateMatrixCMC = apply(covariateMatrix, 1, FUN = max)
# To make elements of above matrix all 0s and 1s
covariateMatrixCMC[covariateMatrixCMC > 0] =
  covariateMatrixCMC[covariateMatrixCMC > 0]^0
# Set parameters
JJ = 2000
Betaj = NULL
pvalj = NULL
qj = NULL
tBetaj = NULL
Bj = NULL
Tj = NULL
Pj = NULL
SW = NULL
SL = NULL
SWwin = NULL
SLwin = NULL
SC = NULL
SQ = NULL
NBetas = 11
IBetas=matrix(rep(NA,JJ*NBetas),ncol= NBetas)
TIBetas=matrix(rep(NA,JJ*NBetas),ncol= NBetas)
PIBetas=matrix(rep(NA,JJ*NBetas),ncol= NBetas)
Phenos=matrix(rep(NA,JJ*321),ncol= JJ)
# Start loop
set.seed(11)
for (j in 1:JJ) {
  # print(j)
  vec=c(0,0,0.83224,0.97060,0,0,0,0,0,-0.93459,0.53073)
  # Create Design
  Design = apply(sweep(covariateMatrix, MARGIN = 2, vec, `*`), 1, FUN = sum)
  # Add normal simulations to the Design matrix
  Y = Design + rnorm(length(Design), 0, 1)
  # Insert elements of Y into Phenos matrix
  Phenos[, j] = Y
  # Assigns covariate CMC of estimate
  Betaj[j] = summary(lm(Y~covariateMatrixCMC))$coefficients[2, 1]
  # Covariate CMC of t value
  tBetaj[j] = summary(lm(Y~covariateMatrixCMC))$coefficients[2, 3]
  # Covariate CMC of Pr(>|t|)
  pvalj[j] = summary(lm(Y~covariateMatrixCMC))$coefficients[2, 4]
  # for # of variants
  for (i in 1:11) {
    # estimate
    IBetas[j,i]=summary(lm(Y~covariateMatrix[,i]))$coeff[2,1]
    # t value
    TIBetas[j,i]=summary(lm(Y~covariateMatrix[,i]))$coeff[2,3]
    # Pr(>|t|)
    PIBetas[j,i]=summary(lm(Y~covariateMatrix[,i]))$coeff[2,4]
  }
  numericMatrix = covariateMatrix
  # Get the number of rows and columns 
  numRows = length(numericMatrix[, 1])
  numCols = length(numericMatrix[1, ])
  # Taking the MAF as the # of occurrences divided by the sample size
  maf = colSums(numericMatrix)/numRows
  # Calculate statistics
  SWwin[j] = as.numeric(testQTLwin(Y,numericMatrix,1/sqrt(maf*(1-maf))))
  SW[j] =    as.numeric(testQTL(Y,numericMatrix,1/sqrt(maf*(1-maf))))
  SL[j] = as.numeric(testQTL(Y,numericMatrix,rep(1,numCols)))
  SLwin[j] = as.numeric(testQTLwin(Y,numericMatrix,rep(1,numCols))) 
  SC[j] = as.numeric(testQTC(Y,numericMatrix))
  SQ[j] = as.numeric(testQTH(Y,numericMatrix))
  # Calculate p values
  p = pvalCalc(Y, covariateMatrix, 10^3)
  # Add p value calculated above to P
  P = cbind(P, p)
  # Print p value
  # cat(p, '\n)
}
# Combine vectors/matrices by the columns
Bj=cbind(Bj,Betaj)
Tj=cbind(Tj,tBetaj)
Pj=cbind(Pj,pvalj)
# Data frame to combine the various estimates
frame=cbind(t(P),SW,SL,SWwin,SLwin,SC,SQ,Bj,Tj,Pj,IBetas)
# Name columns of data frame accordingly
colnames(frame)[c(1,2,3,4,5,6,13:26)]=c('pLW','pLU','pQU','pQH','pMin','pFisher',
                                        'CMCBeta','TCMC','pBeta','B1','B2','B3',
                                        'B4','B5','B6','B7','B8','B9','B10','B11')
# Partially displays elements of above data frame
head(frame)
# Displays # of rows and columns of above data frame
dim(frame)
# Create SIRT1 data based on the data frame
write.csv(frame, '/Users/senadkokic/Desktop/NSERC/GAW17/newsimSIRT1')
# Read the above data as SIRT1 data
newSIRT1 = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/newsimSIRT1')
# Partially display elements of above data (data read correctly)
head(newSIRT1)
# Change column names of Betas to match variants
colnames(IBetas)=c('B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11')
colnames(TIBetas)=c('B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11')
colnames(PIBetas)=c('B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11')
# Make PI and TI Betas into data frames
PIBetas = data.frame(PIBetas)
TIBetas = data.frame(TIBetas)

IBetas1 = IBetas
TIBetas1 = TIBetas
PIBetas1 = PIBetas
newSIRT1.1 = newSIRT1
################ Begin analysis
SIRT1dat = newSIRT1.1
################################################################################
# Obtain bias for all of the significant replicates
## Bootstrapping Algorithm (1)
n = length(Y)
numCols = 10000
# Initialize matrix for sample
Bk = matrix(rep(0, n*numCols), ncol = numCols)
# Delegate elements of matrix to be a sample
set.seed(12)
for (k in 1:numCols) {
  sam = 1:321
  Bk[, k] = sample(sam, 321, replace = TRUE)
}
################################################################################
# Bootstrap replication under the Wald Beta test
### SIRT1
SIRT1dat=newSIRT1.1
# Initialize elements to be NA
SIRT1dat$muBkmed=NA
SIRT1dat$muPkmed=NA
SIRT1dat$muBkmea=NA
SIRT1dat$muPkmea=NA

numericMatrix = convertX(chrData[subID, subSNP])
MAF = getmaf(numericMatrix)
SD = MAF < 0.05 & MAF > 0
covariateMatrix = (numericMatrix[, SD])
covariateMatrix
apply(covariateMatrix, 2, FUN = sum)
chrNames = names(chrData)
P = NULL
# Apply CMC to the covariate matrix
covariateMatrixCMC = apply(covariateMatrix, 1, FUN = max)
# Loop for elements that match criteria wrt pBeta
for (c in which(SIRT1dat$pQU < 0.05)) {
  # Set initial values
  BetaBk=NA
  pBk=NA
  BetaCk=NA
  pCk=NA
  pH=NULL
  pt=NULL
  sigK=0
  for (k in 1:numCols) {
    # Print iteration as well as how many significant results have been seen
    # NOTE: this prints REALLY fast
    # print(paste(k, ' - ', sigK))
    # Break when 200 are observed (? not sure why)
    if (sigK == 200) break
    # Bootstrap sample
    BK = Bk[, k]
    # Remaining elements (left out of bootstrap)
    CK = sam[-BK]
    # Set Y = Phenos element based on index
    Y = Phenos[, c]
    # Compute linear model coefficient for Wald Statistic
    pt[k] = summary(lm(Y[BK]~covariateMatrixCMC[BK]))$coefficients[2,4]
    pt[k]
    
    # Significant result if the Wald stat < 0.05
    if (pt[k] < 0.05) {
      sigK = sigK + 1
      # Take the means of CMC matrix based on those sampled and unsampled
      pBk[k] = mean(covariateMatrixCMC[BK])
      if (pBk[k] > 0) {
        BetaBk[k]=summary(lm(Y[BK]~covariateMatrixCMC[BK]))$coefficients[2,1]
      }
      pCk[k] = mean(covariateMatrixCMC[CK])
      if (pCk[k] > 0) {
        BetaCk[k]=summary(lm(Y[CK]~covariateMatrixCMC[CK]))$coefficients[2,1]
      }
    }
  }
  # Insert differences between estimates into SIRT1dat
  # Median method
  muBkmed=median(BetaBk-BetaCk,na.rm = T)
  muPkmed=median(pBk-pCk,na.rm = T)
  SIRT1dat$muBkmed[c]=muBkmed
  SIRT1dat$muPkmed[c]=muPkmed
  # Mean method
  muBkmea=mean(BetaBk-BetaCk,na.rm = T)
  muPkmea=mean(pBk-pCk,na.rm = T)
  SIRT1dat$muBkmea[c]=muBkmea
  SIRT1dat$muPkmea[c]=muPkmea
}
# Writing files
# write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot',name,sep=''))
# bootSIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot',name,sep=''))

write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot1',name,sep=''))
boot1SIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot1',name,sep=''))

# write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot2',name,sep=''))
# boot2SIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot2',name,sep=''))
head(boot1SIRT1, 100)
### Apply Ghosh et. al. to the significant replicates
SIRT1dat$BT1=NA
SIRT1dat$BT2=NA
SIRT1dat$BT3=NA
# Bonferroni threshold for achieving study-wide significance = 0.1
# Named cee in order to avoid confusion b/w function c()
cee = abs(qnorm(0.5*0.05))
# Begin loop for pBeta 
for (cc in which(SIRT1dat$pQU < 0.05)) {
  # Print the index
  # print(cc)
  # Proposed odds ratio = 1.54
  betahat = SIRT1dat[cc, ]$CMCBeta
  # Reported p-value = 5.7e-4, which we convert to a z-value
  z=sign(betahat)*abs(qnorm(0.5*SIRT1dat[cc,]$pQU))
  ###################################################
  #              THE PROPOSED APPROACH              #
  ###################################################
  # Compute the standard error of betahat
  se = betahat/z
  # Compute the conditional MLE
  mutilde1 = optimize(f=conditional.like,
                      c(-20,20),maximum=T,z=z,cee=cee)$maximum
  # Compute the conditional MLE on the beta scale
  betatilde1 = mutilde1*se
  # Assign the conditional MLE of the odds ratio estimate 
  SIRT1dat[cc, ]$BT1 = betatilde1
  # Print out the estimate 
  SIRT1dat[cc, ]$BT1
  # Numeric integration
  a=sum(apply(as.matrix(seq(-50,50,by=0.01)),1,
              conditional.like.z,cee=cee,z=z))*0.01
  b=sum(apply(as.matrix(seq(-50,50,by=0.01)),1,
              conditional.like,cee=cee,z=z))*0.01
  # Alternate computation of odds ratio
  betatilde2 = (a/b)*se
  # Assign the estimate, print it
  SIRT1dat[cc,]$BT2=betatilde2
  # Compute the average of the previous ORs as the third estimate
  betatilde3 = (betatilde1 + betatilde2)/2
  # Assign the estimate, print it
  SIRT1dat[cc,]$BT3=betatilde3
}
write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot',
                         name,sep=''))
bootSIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot',
                         name,sep=''))

write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot1',
                         name,sep=''))
boot1SIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot1',
                          name,sep=''))
head(boot1SIRT1,100)

write.csv(SIRT1dat,paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot2',
                         name,sep=''))
boot2SIRT1=read.csv(paste('/Users/senadkokic/Desktop/NSERC/GAW17/boot2',
                          name,sep=''))
##### Histogram comparisons:
SIRT1dat = bootSIRT1
IBetas=IBetas1
TIBetas=TIBetas1
PIBetas=PIBetas1
###############################################################################
sink('/Users/senadkokic/Desktop/NSERC/Results Files/Values_CAlpha_Bi')
# For the Wald test (Table 1)
##Truth
round(vec[c(3,4,10,11)], 2)

##Truth over replicates (JJ=2000)
# These give B_k values for table
round(mean(SIRT1dat$B3), 2)
round(mean(SIRT1dat$B4), 2)
round(mean(SIRT1dat$B10), 2)
round(mean(SIRT1dat$B11), 2)
##Bias-Winner's Curse
# These give B_naive values for table
round(mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B3), 2)
round(mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B4), 2)
round(mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B10), 2)
round(mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B11), 2)
##Absolute Bias
round((mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B3)-vec[3]), 2)
round((mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B4)-vec[4]), 2)
round((mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B10)-vec[10]), 2)
round((mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B11)-vec[11]), 2)
##Relative Bias
round((mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B3)-vec[3])/vec[3], 2)
round((mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B4)-vec[4])/vec[4], 2)
round((mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B10)-vec[10])/vec[10], 2)
round((mean(SIRT1dat[SIRT1dat$pQU<0.05,]$B11)-vec[11])/vec[11], 2)
# Get the MAF
round(MAF[SD], 3)
# Power calculation
powerCalc = length(which(SIRT1dat$pQU < 0.05))
paste("Power ", powerCalc, "/", JJ, " = ", round(powerCalc/JJ, 2), sep = '')
sink()












