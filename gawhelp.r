# Help function for GAW17 data

# Han Chinese data
getPop = function(ID){
A = read.csv('/Users/senadkokic/Desktop/NSERC/GAW17/GAW17 CD/unrelateds.ped',header=T,sep=',')
J = length(A[,1])
ii = NULL
for (i in 1:J){
if (substring(A[i,4],1,1) == 'H'){ii = c(ii,which(ID ==A[i,1]))}
if (substring(A[i,4],1,1) == 'J'){ii = c(ii,which(ID ==A[i,1]))}
if (substring(A[i,4],1,1) == 'D'){ii = c(ii,which(ID ==A[i,1]))}
}
return(ii)
}

# get subset SNPs
subSet = function(B,x,start,last){
cord = B[,3]
scord = which(cord>=start & cord<=last)
names = B[scord,1]
jj = NULL
l = length(names)
for (i in 1:l){
jj = c(jj,which(x==names[i]))
}
return(jj)
} 



# Convert call matrix to the numerical matrix
convertX = function(X){
if (is.null(dim(X))){
h=length(X)
l = 1
X = matrix(X)
} else{
h = length(X[,1])
l = length(X[1,]) 
}
S = matrix(0,h,l)
refS = getRef(X[1,])

for (i in 1:h){
  for (j in 1:l){
c = 0
A1 = substring(lapply(X[i,j],as.character),1,1)
if (A1 == refS[j]){c = c +1}
A1 = substring(lapply(X[i,j],as.character),3,3)
if (A1 == refS[j]){c = c +1}
S[i,j] = c  
 }
}
##lapply(Data1[pp,nn],as.character)

MAF = getmaf(S)
J = length(MAF)
for (j in 1:J){
if (MAF[j]>0.5){
S[,j] = abs(S[,j]-2)
}
}
return (S)
}
# get reference 
getRef = function(x){
l = length(x)
ref = rep(0,l)
for (i in 1:l){
cc = lapply(x[i],as.character)
ref[i] =  substring(cc,1,1)
}
return (ref)
}

# Get MAF
getmaf = function(X){
L = length(X[,1])
maf = colSums(X)/(2*L)
return (maf)
}


# Get Q2 values for individuals with ID from file
getY = function(ID,file){
A = read.csv(file,header=T,sep=',')
#A = read.csv('/Users/davids/Documents/Work/Thesis/GAW17/GAW17 CD/unr_phen.1',header=T,sep=',')
Y = NULL
jj = length(ID)
for (i in 1:jj){
kk = which(A[,1]==ID[i])
Y = c(Y,A[kk,6])
}
return(Y)
}
