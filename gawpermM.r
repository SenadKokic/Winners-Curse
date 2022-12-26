
pvalCalc = function(Y,X,perm){
L = length(X[,1])
J = length(X[1,])
maf = colSums(X)/L
SW = as.numeric(testQTL(Y,X,1/sqrt(maf*(1-maf))))
SL = as.numeric(testQTL(Y,X,rep(1,J)))
SC = as.numeric(testQTC(Y,X))
SQ = as.numeric(testQTH(Y,X))

storeLW = rep(0,perm)
storeLS = rep(0,perm)
storeQC = rep(0,perm)
storeQH = rep(0,perm)

for (i in 1:perm){
Y = sample(Y)
storeLW[i] = as.numeric(testQTL(Y,X,1/sqrt(maf*(1-maf))))
storeLS[i] = as.numeric(testQTL(Y,X,rep(1,J)))
storeQC[i] = as.numeric(testQTC(Y,X))
storeQH[i] = as.numeric(testQTH(Y,X))
}

pW = sum(abs(storeLW)>=abs(SW))/perm
pL = sum(abs(storeLS)>=abs(SL))/perm
pC = sum(abs(storeQC)>=abs(SC))/perm
pH = sum(abs(storeQH)>=abs(SQ))/perm

LL = rep(0,perm)
LQ = rep(0,perm)
minP = rep(0,perm)
fP = rep(0,perm)

for (i in 1:perm){
LL[i] = sum(abs(storeLW)>=abs(storeLW[i]))/perm
LQ[i] = sum(abs(storeQH)>=abs(storeQH[i]))/perm
minP[i] = min(LL[i],LQ[i])
fP[i] = -2*log(LL[i]) - 2*log(LQ[i])
}
pFo = -2*log(pW) - 2*log(pH)
pMo = min(pW,pH)
pF = sum(fP>=pFo)/perm
pM = sum(minP<=pMo)/perm  
return (c(pW,pL,pC,pH,pM,pF))
}

singleSNP = function(Y,X){
Yb = mean(Y)

Sobs = abs(sum((Y-Yb)*X))
S = rep(0,10000)
for (i in 1:10000){
Y = sample(Y)
S[i] = abs(sum((Y-Yb)*X))
}
pval = sum(S>=Sobs)/10000
return (pval)
}


