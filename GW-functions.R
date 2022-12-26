#########################################
## THE FUNCTIONS BELOW ARE USED TO OBTAIN THE 
## BIAS-CORRECTED ESTIMATES 
#########################################

conditional.like=function(mu,cee,z){
like=dnorm(z-mu)/(pnorm(mu-cee)+pnorm(-cee-mu))
return((abs(z)>cee)*like) }

conditional.like.z=function(mu,cee,z){
return(conditional.like(mu,cee,z)*mu)
}

#########################################
## THE FUNCTIONS BELOW ARE USED TO OBTAIN THE 
## BIAS-CORRECTED CONFIDENCE INTERVAL 
#########################################

ptruncnorm.lower=function(z,mu,cee,alpha){
 A=pnorm(-cee+mu)+pnorm(-cee-mu)
 term1=pnorm(z-mu)
 term2=pnorm(-cee-mu)
 term3=pnorm(-cee-mu)+pnorm(z-mu)-pnorm(cee-mu)
 result=(1/A)*(term1*(z<= -cee)+term2*(abs(z)<cee)+term3*(z>=cee))
 return(result-(alpha/2))
 }

ptruncnorm.upper=function(z,mu,cee,alpha){
 A=pnorm(-cee+mu)+pnorm(-cee-mu)
 term1=pnorm(z-mu)
 term2=pnorm(-cee-mu)
 term3=pnorm(-cee-mu)+pnorm(z-mu)-pnorm(cee-mu)
 result=(1/A)*(term1*(z<= -cee)+term2*(abs(z)<cee)+term3*(z>=cee))
 return(result-(1-alpha/2))
 }

find.lowerz=function(mu,z,cee,alpha){
lowerz=uniroot(ptruncnorm.lower,lower=-20,upper=20,mu=mu,cee=cee,alpha=alpha)$root
return(lowerz-z)
 }

find.upperz=function(mu,z,cee,alpha){
upperz=uniroot(ptruncnorm.upper,lower=-20,upper=20,mu=mu,cee=cee,alpha=alpha)$root
return(upperz-z)
 }

getCI=function(z,cee,alpha){
 uppermu=uniroot(find.lowerz,interval=c(-15,15),cee=cee,z=z,alpha=alpha)$root
 lowermu=uniroot(find.upperz,interval=c(-15,15),cee=cee,z=z,alpha=alpha)$root
 out=list(lowermu,uppermu)
 names(out)=c("lowermu","uppermu")
 return(out)
 }