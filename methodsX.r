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
quad  =  t(S)%*%solve(varS)%*%S
return (quad)
}


