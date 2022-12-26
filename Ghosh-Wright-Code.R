#########################################################
## THIS EXAMPLE R SCRIPT COMPUTES THE RESULTS 
## OF THE PROPOSED APPROACH REPORTED IN THE TOP LINE OF TABLE 1
## OF GHOSH ET AL.,
## FOR THE DATA USED BY YU ET AL. (2007) AM J HUM GENET 81:540-551.  
##
## WE ASSUME A BASIC FAMILIARITY WITH R,
## AND THAT THE FILE "functions.R" IS IN THE R WORKING DIRECTORY.
## 
##################################################

source("GW-functions.R")# YOU READ IN THE FUNCTIONS FOR OUR METHOD

##################################################
#              INPUT BY THE USER                 #
##################################################

cee=abs(qnorm(.5*0.1/48)) # Bonferroni threshold for achieving study-wide significance = 0.1
                          # we use "cee" so R does not get confused with the function 'c'
betahat=log(1.54) # Reported OR = 1.54 
z=sign(betahat)*abs(qnorm(0.5*5.7e-4)) # Reported p-value = 5.7e-4, which we convert to a z-value


###################################################
#              THE PROPOSED APPROACH              #
###################################################

se=betahat/z # standard error of betahat

mutilde1=optimize(f=conditional.like,c(-20,20),maximum=T,z=z,cee=cee)$maximum # the conditional mle
betatilde1=mutilde1*se # the conditional mle on beta scale
OR1=round(exp(betatilde1),2) # the condition mle OR estimate
OR1 # print out OR_TILDE_1

a=sum(apply(as.matrix(seq(-50,50,by=0.01)),1,conditional.like.z,cee=cee,z=z))*0.01 # numeric integration
b=sum(apply(as.matrix(seq(-50,50,by=0.01)),1,conditional.like,cee=cee,z=z))*0.01
betatilde2=(a/b)*se
OR2=round(exp(betatilde2),2)
OR2 # print out OR_TILDE_2

betatilde3=(betatilde1+betatilde2)/2
OR3=round(exp(betatilde3),2)
OR3 # OR_TILDE_3

cimu=getCI(z=z,cee=cee,alpha=0.05)
cibeta=round(exp(c(cimu$lowermu*(betahat/z),cimu$uppermu*(betahat/z))),2) 
cibeta # BIAS-CORRECTED 95% CI FOR OR