####################################
## This program computes the point estimates and standard errors from the naive method 
##(ignoring exposure misclassification) using or not using the propensity scores, 
##and those from the estimation equation and jointly likelihood based methods. 
## Note that this is for external validation study, 
## and it takes care of left truncation.
## 
## Functions:
## LL_ex_ee: computes the likelihood for estimating equation method.
## hessian_ex_ee_noL3: computes the hessian matrix for estimating equation method.
## LL_ex_mle_noL3: compute the joint likelihood for joint likelihood method.

# Returned value into outfile: [uncorrected point estimate, 
#                  s.e. of uncorrected point estimate, 
#                  uncorrected point estimate with propensity score, 
#                  s.e. of uncorrected point estimate with propensity score,
#                  corrected point estimate using estimating equation method,
#                  s.e. of corrected point estimate using estimating equation #method,
#                  corrected point estimate using joint likelihood method,
#                  s.e. of corrected point estimate using joint likelihood #method]
#

########### Data Requirement##########
##
###validation study###
# data_V=data.frame(X,Z,W)
# X: true exposure  
# Z: surrogate exposure
# W: baseline covariates 
#!!!(varying names for W, but X,Z,T,D names should be fixed)

###main study###
# data_M=data.frame(Z,T,T0,D,W)
# W: baseline covariates
# Z: surrogate exposure
# T: End time
# T0: Start Time
# D: indicator of events
#######################################
library(survival)

LL_ex_ee <- function(b,PSXC, PXZ,data_M){
  
  D=data_M$D
  T=data_M$T
  T0=data_M$T0
  n=length(T)
  
  E=exp(b[1]+b[2]*PSXC)*PXZ+exp(b[2]*PSXC)*(1-PXZ)
  F=vector()
  for(k in 1:n){
    if(D[k]==1){
      Num=E[k]
      Den=sum(ifelse(T>T[k] & T0 < T[k],E,0))
      F[k]=Num/Den
    }
    else{F[k]=1}
  }
  LL=sum(log(F))  #prod(F)
  return(-LL)
}

hessian_ex_ee_noL3 <- function(b,gamma,PSXC,PXZ,PSZ,D,T,T0,Z,Z_V){
  e_it=exp(b[1]+b[2]*PSXC)*PXZ+exp(b[2]*PSXC)*(1-PXZ)
  b_it=exp(b[1]+b[2]*PSXC)*PXZ
  c_it=PSXC*e_it
  
  n_M=length(PSZ)
  n_V=length(Z_V)
  
  h22=0
  
  
  Ftemp22=vector()
  for(k in 1:n_M){
    if(D[k]==1){
      
      Num22_1=sum(ifelse(T>T[k] & T0 < T[k],c_it*PSXC,0))
      Den22_1=sum(ifelse(T>T[k] & T0 < T[k],e_it,0))
      
      Num22_2=(sum(ifelse(T>T[k] & T0 < T[k],c_it,0)))^2
      Den22_2=(sum(ifelse(T>T[k] & T0 < T[k],e_it,0)))^2
      
      Ftemp22[k]=-Num22_1/Den22_1+Num22_2/Den22_2
    }
    else{Ftemp22[k]=0}
  }
  
  
  h22=h22+sum(Ftemp22,na.rm =TRUE)
  
  h12=0
  
  Ftemp12=vector()
  for(k in 1:n_M){
    if(D[k]==1){
      
      Num12_1=sum(ifelse(T>T[k] & T0 < T[k],b_it*PSXC,0))
      Den12_1=sum(ifelse(T>T[k] & T0 < T[k],e_it,0))
      
      Num12_2=sum(ifelse(T>T[k] & T0 < T[k],b_it,0))*sum(ifelse(T>T[k],c_it,0))
      Den12_2=sum(ifelse(T>T[k] & T0 < T[k],e_it,0))^2
      
      Ftemp12[k]=-Num12_1/Den12_1+Num12_2/Den12_2
    }
    else{Ftemp12[k]=0}
  }
  
  
  h12=h12+sum(Ftemp12,na.rm =TRUE)
  
  h11=0
  
  Ftemp11=vector()
  for(k in 1:n_M){
    if(D[k]==1){
      h11=h11+b_it[k]/e_it[k]-b_it[k]^2/e_it[k]^2
      
      Num11_1=sum(ifelse(T>T[k] & T0 < T[k],b_it,0))
      Den11_1=sum(ifelse(T>T[k] & T0 < T[k],e_it,0))
      
      Num11_2=sum(ifelse(T>T[k] & T0 < T[k],b_it,0))^2
      Den11_2=sum(ifelse(T>T[k] & T0 < T[k],e_it,0))^2
      
      Ftemp11[k]=-Num11_1/Den11_1+Num11_2/Den11_2
    }
    else{Ftemp11[k]=0}
  }
  
  
  h11=h11+sum(Ftemp11,na.rm =TRUE)
  
  h33=0
  h34=0
  h44=0
  
  for(k in 1:n_V){
    h33=h33-exp(gamma[1]+gamma[2]*Z_V[k])/(1+exp(gamma[1]+gamma[2]*Z_V[k]))^2
    h34=h34-exp(gamma[1]+gamma[2]*Z_V[k])*Z_V[k]/(1+exp(gamma[1]+gamma[2]*Z_V[k]))^2
    h44=h44-exp(gamma[1]+gamma[2]*Z_V[k])*Z_V[k]^2/(1+exp(gamma[1]+gamma[2]*Z_V[k]))+exp(2*(gamma[1]+gamma[2]*Z_V[k]))*Z_V[k]^2/(1+exp(gamma[1]+gamma[2]*Z_V[k]))^2
  }
  
  ##########
  px1z1=exp(gamma[1]+gamma[2])/(1+exp(gamma[1]+gamma[2]))
  px0z1=1-px1z1
  px1z0=exp(gamma[1])/(1+exp(gamma[1]))
  px0z0=1-px1z0
  
  dpsxcdgamma0=px1z1*px0z1*PSZ+px1z0*px0z0*(1-PSZ)
  dpsxcdgamma1=px1z1*px0z1*PSZ
  
  px1z=exp(gamma[1]+gamma[2]*Z)/(1+exp(gamma[1]+gamma[2]*Z))
  px0z=1-px1z
  
  dedgamma0=exp(b[1]+b[2]*PSXC)*px1z*(dpsxcdgamma0*b[2]+px0z)+exp(b[2]*PSXC)*px0z*(dpsxcdgamma0*b[2]-px1z)
  dedgamma1=exp(b[1]+b[2]*PSXC)*px1z*(dpsxcdgamma1*b[2]+Z*px0z)+exp(b[2]*PSXC)*px0z*(dpsxcdgamma1*b[2]-Z*px1z)
  
  dcdgamma0=dpsxcdgamma0*e_it+PSXC*dedgamma0
  dcdgamma1=dpsxcdgamma1*e_it+PSXC*dedgamma1
  
  dbdgamma0=exp(b[1]+b[2]*PSXC)*px1z*(dpsxcdgamma0*b[2]+px0z)
  dbdgamma1=exp(b[1]+b[2]*PSXC)*px1z*(dpsxcdgamma1*b[2]+Z*px0z)
  
  h23=0
  Ftemp23=vector()
  for(k in 1:n_M){
    if(D[k]==1){
      h23=h23+dcdgamma0[k]/e_it[k]-c_it[k]*dedgamma0[k]/e_it[k]^2
      
      Num23_1=sum(ifelse(T>T[k] & T0 < T[k],dcdgamma0,0))
      Den23_1=sum(ifelse(T>T[k] & T0 < T[k],e_it,0))
      
      Num23_2=sum(ifelse(T>T[k] & T0 < T[k],dedgamma0,0))*sum(ifelse(T>T[k],c_it,0))
      Den23_2=sum(ifelse(T>T[k] & T0 < T[k],e_it,0))^2
      
      Ftemp23[k]=-Num23_1/Den23_1+Num23_2/Den23_2
    }
    else{Ftemp23[k]=0}
  }
  
  
  h23=h23+sum(Ftemp23,na.rm =TRUE)
  
  h24=0
  Ftemp24=vector()
  for(k in 1:n_M){
    if(D[k]==1){
      h24=h24+dcdgamma1[k]/e_it[k]-c_it[k]*dedgamma1[k]/e_it[k]^2
      
      Num24_1=sum(ifelse(T>T[k] & T0 < T[k],dcdgamma1,0))
      Den24_1=sum(ifelse(T>T[k] & T0 < T[k],e_it,0))
      
      Num24_2=sum(ifelse(T>T[k] & T0 < T[k],dedgamma1,0))*sum(ifelse(T>T[k],c_it,0))
      Den24_2=sum(ifelse(T>T[k] & T0 < T[k],e_it,0))^2
      
      Ftemp24[k]=-Num24_1/Den24_1+Num24_2/Den24_2
    }
    else{Ftemp24[k]=0}
  }
  
  
  h24=h24+sum(Ftemp24,na.rm =TRUE)
  
  h13=0
  Ftemp13=vector()
  for(k in 1:n_M){
    if(D[k]==1){
      h13=h13+dbdgamma0[k]/e_it[k]-b_it[k]*dedgamma0[k]/e_it[k]^2
      
      Num13_1=sum(ifelse(T>T[k] & T0 < T[k],dbdgamma0,0))
      Den13_1=sum(ifelse(T>T[k] & T0 < T[k],e_it,0))
      
      Num13_2=sum(ifelse(T>T[k] & T0 < T[k],dedgamma0,0))*sum(ifelse(T>T[k],b_it,0))
      Den13_2=sum(ifelse(T>T[k] & T0 < T[k],e_it,0))^2
      
      Ftemp13[k]=-Num13_1/Den13_1+Num13_2/Den13_2
    }
    else{Ftemp13[k]=0}
  }
  
  
  h13=h13+sum(Ftemp13,na.rm =TRUE)
  
  h14=0
  Ftemp14=vector()
  for(k in 1:n_M){
    if(D[k]==1){
      h14=h14+dbdgamma1[k]/e_it[k]-b_it[k]*dedgamma1[k]/e_it[k]^2
      
      Num14_1=sum(ifelse(T>T[k] & T0 < T[k],dbdgamma1,0))
      Den14_1=sum(ifelse(T>T[k] & T0 < T[k],e_it,0))
      
      Num14_2=sum(ifelse(T>T[k] & T0 < T[k],dedgamma1,0))*sum(ifelse(T>T[k],b_it,0))
      Den14_2=sum(ifelse(T>T[k] & T0 < T[k],e_it,0))^2
      
      Ftemp14[k]=-Num14_1/Den14_1+Num14_2/Den14_2
    }
    else{Ftemp14[k]=0}
  }
  
  
  h14=h14+sum(Ftemp14,na.rm =TRUE)
  
  hessian_result=matrix(data=c(h11,h12,h13,h14,h12,h22,h23,h24,0,0,h33,h34,0,0,h34,h44),nrow=4,byrow=TRUE)
  
  return(hessian_result)
}

LL_ex_mle_noL3 <- function(b,data_M,data_V,PSZ){
  n_M=dim(data_M)[1]
  n_V=dim(data_V)[1]
  D=data_M$D
  T=data_M$T
  T0=data_M$T0
  Z=c(data_M[,1], data_V[,2])
  X=data_V[,1]
  PXZ = exp(b[3]+Z*b[4])/(1+exp(b[3]+Z*b[4]))
  PSXC =  exp(b[3]+b[4])/(1+exp(b[3]+ b[4]))* PSZ + exp(b[3])/(1+exp(b[3]))*(1-PSZ)
  
  E=exp(b[1]+b[2]*PSXC)*PXZ+exp(b[2]*PSXC)*(1-PXZ)
  
  L1=vector()
  L2=vector()
  for(k in 1:n_M){
    if(D[k]==1){
      Num=E[k]
      Den=sum(ifelse(T>T[k] & T0 < T[k],E,0))
      L1[k]=Num/Den
    }
    else{L1[k]=1}
  }
  
  logL1=sum(log(L1))  #prod(F)2
  
  L2=ifelse(X[1:n_V]==1,exp(b[3]+b[4]*Z[(n_M+1):(n_M+n_V)])/(1+exp(b[3]+b[4]*Z[(n_M+1):(n_M+n_V)])),1/(1+exp(b[3]+b[4]*Z[(n_M+1):(n_M+n_V)])))
  logL2=sum(log(L2))
  
  LL=logL1+logL2
  
  return(-LL)
}

########################
# 1. read the files in #
########################

data_M =read.csv("data_main.csv")
data_V = read.csv("data_validation.csv")
outfile = "out.csv"

#####################################################
#First, get the relationship in the validation study#
#####################################################

##put the vector of coefficients as gamma##
XWZ_val = glm(X ~ Z, family=binomial(logit), data=data_V)
gamma=c(XWZ_val$coefficients[1],XWZ_val$coefficients[2])

########################################################################
#second, get the crude PSZ in the main study and PSXC in the main study#
########################################################################

###This piece can be directly get from SAS logistic regression(main study), ###
###only need the predicted value P(Z=1|W)##
##########put the predicted value into PSZ##########################################################################
W_M = data_M[,5:ncol(data_M)]
ZWformula = paste(names(W_M),collapse="+")
ZWformula1 = paste("Z ~ ",ZWformula,sep="")

ZW = glm(ZWformula1, family=binomial(logit), data=data_M)
PSZ=fitted.values(ZW)

# uncorrected models
r3<-coxph(as.formula(paste("Surv(time=T0, time2=T, event=D) ~ Z+",ZWformula,sep="")), data=data_M) 
r4<-coxph(Surv(time=T0, time2=T, event=D) ~ Z+PSZ, data_M) 

PXZ=exp(gamma[1] + gamma[2]*data_M$Z)/(1+exp(gamma[1]+gamma[2]*data_M$Z))
PSXC=exp(gamma[1]+gamma[2])/(1+exp(gamma[1]+gamma[2]))*PSZ+1/(1+exp(gamma[1]))*(1-PSZ)


# EE
b0_ee=c(0,0)

res2_ex_ee<-optim(b0_ee,LL_ex_ee,PSXC=PSXC, PXZ=PXZ,data_M=data_M,method="BFGS",hessian=TRUE)
b2_ex_ee=res2_ex_ee$par


hes1=hessian_ex_ee_noL3(b2_ex_ee,gamma,PSXC=PSXC,PXZ,PSZ,D=data_M$D,T=data_M$T,T0=data_M$T0,Z=data_M$Z, Z_V=data_V$Z)

varg=matrix(data=c(res2_ex_ee$hessian[1,1],res2_ex_ee$hessian[1,2],0,0,res2_ex_ee$hessian[2,1],res2_ex_ee$hessian[2,2],0,0,
                   0,0,-hes1[3,3],-hes1[3,4],0,0,-hes1[4,3],-hes1[4,4]),nrow=4,byrow=TRUE)


b2_ex_sd_ee_noL3=try(sqrt((solve(hes1)%*%varg%*%solve(hes1))[1,1]),silent=TRUE)
if ('try-error' %in% class(b2_ex_sd_ee_noL3)) b2_ex_sd_ee_noL3=NaN

# JL
b0_mle_noL3=c(0,0,0,0)


res2_ex_mle<-optim(b0_mle_noL3,LL_ex_mle_noL3,data_M=data_M,data_V=data_V,PSZ=PSZ,method="BFGS",hessian = TRUE)
b2_ex_mle=res2_ex_mle$par
b2_ex_sd_mle=try(sqrt(solve(res2_ex_mle$hessian)[1,1]),silent=TRUE)
if ('try-error' %in% class(b2_ex_sd_mle)) b2_ex_sd_mle=NaN  

# returned results
returnValue = c(r3$coefficients[1],sqrt(r3$var[1,1]),r4$coefficients[1],sqrt(r4$var[1,1]),b2_ex_ee[1],b2_ex_sd_ee_noL3,b2_ex_mle[1],b2_ex_sd_mle)

write.csv(returnValue,outfile)




