# This is to creat functions
# It based on score and Section 4

# Calculte Score

Sstat = function(Y,X){
Ybar =mean(Y)
Yn  = Y - Ybar
S = t(X)%*%Yn
return (S)
}

mmatrix = function(X){
l = length(X[1,])
mX = NULL

for (i in 1:l){
mm = mean(X[,i]) 
mv = X[,i]-mm
mX = cbind(mX,mv)
}
return (mX)
}

testQTL = function(Y,X,w){
S = Sstat(Y,X)
sigma = var(Y)
linear  =  t(w)%*%S
return (abs(linear))
}

testQTLwin = function(Y,X,w){
S = Sstat(Y,X)
sigma = var(Y)
linear  =  t(w)%*%S
return ((linear))
}

testQTC = function(Y,X,w){
S = Sstat(Y,X)
sigma = var(Y)
quad  =  t(S)%*%S
return (quad)
}

testQTH = function(Y,X,w){
S = Sstat(Y,X)
sigma = var(Y)
quad  =  t(S)%*%S
mX = mmatrix(X)
varS = diag(diag((t(mX) %*% (mX))))
if(name=='GCKR'){
varS = t(mX) %*% (mX)	
}
quad  =  t(S)%*%solve(varS)%*%S
return (quad)
}

testQTHboot = function(Y,X,Xorig,w){
S = Sstat(Y,X)
sigma = var(Y)
quad  =  t(S)%*%S
mX = mmatrix(X)  # minor allele frequencies (I am forcing the use of MAF from original data to avoid MAF of 0 for some variants)
for(ll in which(mX[1,]==0)){
	mXorig=mmatrix(Xorig)
	mX[,ll]=mXorig[,ll]
}
varS = diag(diag((t(mX) %*% (mX))))
if(name=='GCKR'){
varS = t(mX) %*% (mX)	
}
quad  =  t(S)%*%solve(varS)%*%S
return (quad)
}

estBeta = function(Y,X,w){
	

}
