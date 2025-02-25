
## this code assumes the existence of 2 files:
## (1) the main study CSV file
## (2) the vali study CSV file
## the MS file has columns containing Y, Z1, Z2, C, all binary variables
## Z1 and Z2 are potentially misclassified variables. C is a correctly classified variable, e.g. age group
## Y is the outcome variable
## the VS file has columns containing Y, Z1, Z2, C, X1, X2
## the subjects in both files do not overlap

##############################################################################################################

# Steps for correcting for misclassification bias
# 1. read the main and validation study files in
# 2. define the functions for MS/IVS design under double and single transportability assumptions
# 3. use 'optim' to find the MLEs for the parameters, inputting initial values for the paramters into 'optim'
# 4. calculate the corrected pPAR and its CI under the 2 assumptions


##########################################################################################
################## 1. read the files in ##################################################
##########################################################################################

mainstudy = read.csv("main.csv")
ym = mainstudy[,1]
zm1 = mainstudy[,2]
zm2 = mainstudy[,3]
cm = mainstudy[,4]
zm = cbind(zm1,zm2)

valistudy = read.csv("vali.csv")
yv = valistudy[,1]
zv1 = valistudy[,2]
zv2 = valistudy[,3]
cv = valistudy[,4]
xv1 = valistudy[,5]
xv2 = valistudy[,6]
zv = cbind(zv1,zv2)
xv = cbind(xv1,xv2)

### if functions have not been defined, go to bottom of document and run the code under 'define the functions'

ppar = matrix(NA,3,1)	##### create arrays for storing estimates
CI_ppar = matrix(NA,3,2)

fit1 = glm(c(ym,yv)~c(zm1,zv1)+c(zm2,zv2)+c(cm1,cv1),family="binomial")

	# =============================================================================
	# ==== Calculate the Naive pPAR ===============================================
	# =============================================================================
	beta = coef(fit1)
	pi_z = table(4*cm1+(zm1+2*zm2))
	pi_z = pi_z/sum(pi_z)
	pi = pi_z
	RR0 = 1
	RR1 = exp(beta[1]+beta[2])/(1+exp(beta[1]+beta[2])) / ( exp(beta[1])/(1+exp(beta[1])) )
	RR2 = exp(beta[1]+beta[3])/(1+exp(beta[1]+beta[3])) / ( exp(beta[1])/(1+exp(beta[1])) )
	RR3 = exp(beta[1]+beta[2]+beta[3])/(1+exp(beta[1]+beta[2]+beta[3])) / ( exp(beta[1])/(1+exp(beta[1])) )
	RR0a = exp(beta[1]+beta[4])/(1+exp(beta[1]+beta[4])) / ( exp(beta[1])/(1+exp(beta[1])) )
	RR1a = exp(beta[1]+beta[2]+beta[4])/(1+exp(beta[1]+beta[2]+beta[4])) / ( exp(beta[1])/(1+exp(beta[1])) )
	RR2a = exp(beta[1]+beta[3]+beta[4])/(1+exp(beta[1]+beta[3]+beta[4])) / ( exp(beta[1])/(1+exp(beta[1])) )
	RR3a = exp(beta[1]+beta[2]+beta[3]+beta[4])/(1+exp(beta[1]+beta[2]+beta[3]+beta[4])) / ( exp(beta[1])/(1+exp(beta[1])) )
 	RR_vec = c(1,RR1,RR2,RR3,RR0a,RR1a,RR2a,RR3a)	;	RR_num = c(1,1,RR2,RR2,RR0a,RR0a,RR2a,RR2a)
	ppar[1,1] =	1 - sum(pi*RR_num)/sum(pi*RR_vec)	

### if desired, obtain the naive CI under "Calculate the Uncorrected CI", then come back to obtain corrected point estimates
	
	# ====================================================================================================
	# ========== Evaluate the IVS2-corrected pPAR=========================================================
	# ====================================================================================================
	gamma = log(pi_z[2:8]/pi_z[1])	
	initial.value = c(beta,rep(0,12),gamma)		## naive values used as initial estimates for beta/gamma
	print(ivs.doubletransportability(initial.value))
	ivs2output = optim(initial.value,ivs.doubletransportability,method="BFGS",control=list(abstol=0.0000000000000001),hessian=T)
	print(ivs2output)
	betatheta = ivs2output$par
	vcov1 = solve(ivs2output$hessian)
	beta = betatheta[1:4]
	xeta = betatheta[17:23]	
	pi = c(0,xeta)
	pi = exp(pi)
	sum_exp_pi = sum(pi)
	pi = pi/sum_exp_pi
	RR0 = 1
	RR1 = exp(beta[1]+beta[2])/(1+exp(beta[1]+beta[2])) / ( exp(beta[1])/(1+exp(beta[1])) )
	RR2 = exp(beta[1]+beta[3])/(1+exp(beta[1]+beta[3])) / ( exp(beta[1])/(1+exp(beta[1])) )
	RR3 = exp(beta[1]+beta[2]+beta[3])/(1+exp(beta[1]+beta[2]+beta[3])) / ( exp(beta[1])/(1+exp(beta[1])) )
	RR0a = exp(beta[1]+beta[4])/(1+exp(beta[1]+beta[4])) / ( exp(beta[1])/(1+exp(beta[1])) )
	RR1a = exp(beta[1]+beta[2]+beta[4])/(1+exp(beta[1]+beta[2]+beta[4])) / ( exp(beta[1])/(1+exp(beta[1])) )
	RR2a = exp(beta[1]+beta[3]+beta[4])/(1+exp(beta[1]+beta[3]+beta[4])) / ( exp(beta[1])/(1+exp(beta[1])) )
	RR3a = exp(beta[1]+beta[2]+beta[3]+beta[4])/(1+exp(beta[1]+beta[2]+beta[3]+beta[4])) / ( exp(beta[1])/(1+exp(beta[1])) )
	 RR_vec = c(1,RR1,RR2,RR3,RR0a,RR1a,RR2a,RR3a)
	 RR_num = c(1,1,RR2,RR2,RR0a,RR0a,RR2a,RR2a)
	 ppar[2,1] = 1 - sum(pi*RR_num)/sum(pi*RR_vec)

	# ====================================================================================================
	# ====Evaluate the CI ================================================================================
	# ====================================================================================================

	df_dpi = matrix(0,8,1)
	 RR_vec = c(1,RR1,RR2,RR3,RR0a,RR1a,RR2a,RR3a)
	 RR_num = c(1,1,RR2,RR2,RR0a,RR0a,RR2a,RR2a)
	 a = sum(pi*RR_num)
	 b = sum(pi*RR_vec)
	 b2 = b^2
	for(i in 1:8){
		df_dpi[i,1] = (b*RR_num[i] - a*RR_vec[i])/(b2)
	}

	dpi_dxeta = matrix(0,7,7)
	for (i in 1:6){
		dpi_dxeta[i,i] = exp(xeta[i])*(sum_exp_pi-exp(xeta[i]))/sum_exp_pi/sum_exp_pi
  		for (j in (i+1):7){
    			dpi_dxeta[i,j] = -exp(xeta[i]+xeta[j])/sum_exp_pi/sum_exp_pi
    			dpi_dxeta[j,i] = dpi_dxeta[i,j]
		}
	}
	i = 7
	dpi_dxeta[i,i] = exp(xeta[i])*(sum_exp_pi-exp(xeta[i]))/sum_exp_pi/sum_exp_pi
 	 		
	top_row = rep(NA,7)
	for (i in 1:7){
  		top_row[i] = -exp(xeta[i])/sum_exp_pi/sum_exp_pi
	}
	dpi_dxeta = cbind(top_row, dpi_dxeta)

	dRR_dbeta = matrix(0,4,7)
	# dRR1/dbeta0
	dRR_dbeta[1,1] = (1-exp(beta[2]))*(exp(beta[1]+beta[2]))/(1+exp(beta[1]+beta[2]))/(1+exp(beta[1]+beta[2]))
	# dRR2/dbeta0
	dRR_dbeta[1,2] = (1-exp(beta[3]))*(exp(beta[1]+beta[3]))/(1+exp(beta[1]+beta[3]))/(1+exp(beta[1]+beta[3]))
	# dRR3/dbeta0
	dRR_dbeta[1,3] = (1-exp(beta[2]+beta[3]))*(exp(beta[1]+beta[2]+beta[3]))/(1+exp(beta[1]+beta[2]+beta[3]))/(1+exp(beta[1]+beta[2]+beta[3]))
	# dRR0a/dbeta0
	dRR_dbeta[1,4] = (1-exp(beta[4]))*(exp(beta[1]+beta[4]))/(1+exp(beta[1]+beta[4]))/(1+exp(beta[1]+beta[4]))
	# dRR1a/dbeta0
	dRR_dbeta[1,5] = (1-exp(beta[2]+beta[4]))*(exp(beta[1]+beta[2]+beta[4]))/(1+exp(beta[1]+beta[2]+beta[4]))/(1+exp(beta[1]+beta[2]+beta[4]))
	# dRR2a/dbeta0
	dRR_dbeta[1,6] = (1-exp(beta[3]+beta[4]))*(exp(beta[1]+beta[3]+beta[4]))/(1+exp(beta[1]+beta[3]+beta[4]))/(1+exp(beta[1]+beta[3]+beta[4]))
	# dRR3a/dbeta0
	dRR_dbeta[1,7] = (1-exp(beta[2]+beta[3]+beta[4]))*(exp(beta[1]+beta[2]+beta[3]+beta[4]))/(1+exp(beta[1]+beta[2]+beta[3]+beta[4]))/(1+exp(beta[1]+beta[2]+beta[3]+beta[4]))
	
	# dRR1/dbeta1
	dRR_dbeta[2,1] = exp(beta[2])*(1+exp(beta[1]))/(1+exp(beta[1]+beta[2]))/(1+exp(beta[1]+beta[2]))  
	
	# dRR2/dbeta2
	dRR_dbeta[3,2] = exp(beta[3])*(1+exp(beta[1]))/(1+exp(beta[1]+beta[3]))/(1+exp(beta[1]+beta[3]))  

	# dRR3/dbeta1
	dRR_dbeta[2,3] = exp(beta[2]+beta[3])*(1+exp(beta[1]))/(1+exp(beta[1]+beta[2]+beta[3]))/(1+exp(beta[1]+beta[2]+beta[3]))  
	# dRR3/dbeta2
	dRR_dbeta[3,3] = dRR_dbeta[2,3]

	# dRR0a/dbeta_age
	dRR_dbeta[4,4] = exp(beta[4])*(1+exp(beta[1]))/(1+exp(beta[1]+beta[4]))/(1+exp(beta[1]+beta[4]))  

	# dRR1a/dbeta1
	dRR_dbeta[2,5] = exp(beta[2]+beta[4])*(1+exp(beta[1]))/(1+exp(beta[1]+beta[2]+beta[4]))/(1+exp(beta[1]+beta[2]+beta[4]))  
	# dRR1a/dbeta_age
	dRR_dbeta[4,5] = dRR_dbeta[2,5]

	# dRR2a/dbeta2
	dRR_dbeta[3,6] = exp(beta[3]+beta[4])*(1+exp(beta[1]))/(1+exp(beta[1]+beta[3]+beta[4]))/(1+exp(beta[1]+beta[3]+beta[4]))  
	# dRR2a/dbeta_age
	dRR_dbeta[4,6] = dRR_dbeta[3,6]

	# dRR3a/dbeta1
	dRR_dbeta[2,7] = exp(beta[2]+beta[3]+beta[4])*(1+exp(beta[1]))/(1+exp(beta[1]+beta[2]+beta[3]+beta[4]))/(1+exp(beta[1]+beta[2]+beta[3]+beta[4]))  
	# dRR3a/dbeta2
	dRR_dbeta[3,7] = dRR_dbeta[2,7]
	# dRR3a/dbeta_age
	dRR_dbeta[4,7] = dRR_dbeta[2,7]

	df_dRR = matrix(0,7,1)	
	for (i in c(2,4,6,8)){
		df_dRR[i-1,1] = pi[i]*a/b2
	}
	for (i in c(3,5,7)){
		df_dRR[i-1,1] = (pi[i]*a-b*(pi[i]+pi[i+1]))/b2

	}

	var_ppar = t(df_dRR) %*% t(dRR_dbeta) %*% vcov1[1:4,1:4] %*% dRR_dbeta %*% df_dRR +
			t(df_dpi) %*% t(dpi_dxeta) %*% vcov1[17:23,17:23] %*% dpi_dxeta %*% df_dpi +
			2* ( t(df_dpi) %*% t(dpi_dxeta) %*% vcov1[17:23,1:4] %*% dRR_dbeta %*% df_dRR )
	
	CI_ppar[2,1] = ppar[2,1] - qnorm(0.975)*sqrt(var_ppar)
	CI_ppar[2,2] = ppar[2,1] + qnorm(0.975)*sqrt(var_ppar)



	# ====================================================================================================
	# ========== Evaluate the IVS1-corrected pPAR=========================================================
	# ====================================================================================================
	
initial.value = c(beta,rep(0,12),gamma)
print(ivs.singletransportability(initial.value))
ivs1output = optim(initial.value,ivs.singletransportability,method="BFGS",control=list(abstol=0.0000000000000001),hessian=T)
print(ivs1output)
	betatheta = ivs1output$par
	beta = betatheta[1:4]
	xeta = betatheta[17:23]	
	pi = c(0,xeta)
	pi = exp(pi)
	sum_exp_pi = sum(pi)
	pi = pi/sum_exp_pi
	RR0 = 1
	RR1 = exp(beta[1]+beta[2])/(1+exp(beta[1]+beta[2])) / ( exp(beta[1])/(1+exp(beta[1])) )
	RR2 = exp(beta[1]+beta[3])/(1+exp(beta[1]+beta[3])) / ( exp(beta[1])/(1+exp(beta[1])) )
	RR3 = exp(beta[1]+beta[2]+beta[3])/(1+exp(beta[1]+beta[2]+beta[3])) / ( exp(beta[1])/(1+exp(beta[1])) )
	RR0a = exp(beta[1]+beta[4])/(1+exp(beta[1]+beta[4])) / ( exp(beta[1])/(1+exp(beta[1])) )
	RR1a = exp(beta[1]+beta[2]+beta[4])/(1+exp(beta[1]+beta[2]+beta[4])) / ( exp(beta[1])/(1+exp(beta[1])) )
	RR2a = exp(beta[1]+beta[3]+beta[4])/(1+exp(beta[1]+beta[3]+beta[4])) / ( exp(beta[1])/(1+exp(beta[1])) )
	RR3a = exp(beta[1]+beta[2]+beta[3]+beta[4])/(1+exp(beta[1]+beta[2]+beta[3]+beta[4])) / ( exp(beta[1])/(1+exp(beta[1])) )
	 RR_vec = c(1,RR1,RR2,RR3,RR0a,RR1a,RR2a,RR3a)
	 RR_num = c(1,1,RR2,RR2,RR0a,RR0a,RR2a,RR2a)
	 ppar[3,1] = 1 - sum(pi*RR_num)/sum(pi*RR_vec)

	# ====================================================================================================
	# ====Evaluate the CI ================================================================================
	# ====================================================================================================

	df_dpi = matrix(0,8,1)
	 RR_vec = c(1,RR1,RR2,RR3,RR0a,RR1a,RR2a,RR3a)
	 RR_num = c(1,1,RR2,RR2,RR0a,RR0a,RR2a,RR2a)
	 a = sum(pi*RR_num)
	 b = sum(pi*RR_vec)
	 b2 = b^2
	for(i in 1:8){
		df_dpi[i,1] = (b*RR_num[i] - a*RR_vec[i])/(b2)
	}

	dpi_dxeta = matrix(0,7,7)
	for (i in 1:6){
		dpi_dxeta[i,i] = exp(xeta[i])*(sum_exp_pi-exp(xeta[i]))/sum_exp_pi/sum_exp_pi
  		for (j in (i+1):7){
    			dpi_dxeta[i,j] = -exp(xeta[i]+xeta[j])/sum_exp_pi/sum_exp_pi
    			dpi_dxeta[j,i] = dpi_dxeta[i,j]
		}
	}
	i = 7
	dpi_dxeta[i,i] = exp(xeta[i])*(sum_exp_pi-exp(xeta[i]))/sum_exp_pi/sum_exp_pi
 	 		
	top_row = rep(NA,7)
	for (i in 1:7){
  		top_row[i] = -exp(xeta[i])/sum_exp_pi/sum_exp_pi
	}
	dpi_dxeta = cbind(top_row, dpi_dxeta)

	dRR_dbeta = matrix(0,4,7)
	# dRR1/dbeta0
	dRR_dbeta[1,1] = (1-exp(beta[2]))*(exp(beta[1]+beta[2]))/(1+exp(beta[1]+beta[2]))/(1+exp(beta[1]+beta[2]))
	# dRR2/dbeta0
	dRR_dbeta[1,2] = (1-exp(beta[3]))*(exp(beta[1]+beta[3]))/(1+exp(beta[1]+beta[3]))/(1+exp(beta[1]+beta[3]))
	# dRR3/dbeta0
	dRR_dbeta[1,3] = (1-exp(beta[2]+beta[3]))*(exp(beta[1]+beta[2]+beta[3]))/(1+exp(beta[1]+beta[2]+beta[3]))/(1+exp(beta[1]+beta[2]+beta[3]))
	# dRR0a/dbeta0
	dRR_dbeta[1,4] = (1-exp(beta[4]))*(exp(beta[1]+beta[4]))/(1+exp(beta[1]+beta[4]))/(1+exp(beta[1]+beta[4]))
	# dRR1a/dbeta0
	dRR_dbeta[1,5] = (1-exp(beta[2]+beta[4]))*(exp(beta[1]+beta[2]+beta[4]))/(1+exp(beta[1]+beta[2]+beta[4]))/(1+exp(beta[1]+beta[2]+beta[4]))
	# dRR2a/dbeta0
	dRR_dbeta[1,6] = (1-exp(beta[3]+beta[4]))*(exp(beta[1]+beta[3]+beta[4]))/(1+exp(beta[1]+beta[3]+beta[4]))/(1+exp(beta[1]+beta[3]+beta[4]))
	# dRR3a/dbeta0
	dRR_dbeta[1,7] = (1-exp(beta[2]+beta[3]+beta[4]))*(exp(beta[1]+beta[2]+beta[3]+beta[4]))/(1+exp(beta[1]+beta[2]+beta[3]+beta[4]))/(1+exp(beta[1]+beta[2]+beta[3]+beta[4]))
	
	# dRR1/dbeta1
	dRR_dbeta[2,1] = exp(beta[2])*(1+exp(beta[1]))/(1+exp(beta[1]+beta[2]))/(1+exp(beta[1]+beta[2]))  
	
	# dRR2/dbeta2
	dRR_dbeta[3,2] = exp(beta[3])*(1+exp(beta[1]))/(1+exp(beta[1]+beta[3]))/(1+exp(beta[1]+beta[3]))  

	# dRR3/dbeta1
	dRR_dbeta[2,3] = exp(beta[2]+beta[3])*(1+exp(beta[1]))/(1+exp(beta[1]+beta[2]+beta[3]))/(1+exp(beta[1]+beta[2]+beta[3]))  
	# dRR3/dbeta2
	dRR_dbeta[3,3] = dRR_dbeta[2,3]

	# dRR0a/dbeta_age
	dRR_dbeta[4,4] = exp(beta[4])*(1+exp(beta[1]))/(1+exp(beta[1]+beta[4]))/(1+exp(beta[1]+beta[4]))  

	# dRR1a/dbeta1
	dRR_dbeta[2,5] = exp(beta[2]+beta[4])*(1+exp(beta[1]))/(1+exp(beta[1]+beta[2]+beta[4]))/(1+exp(beta[1]+beta[2]+beta[4]))  
	# dRR1a/dbeta_age
	dRR_dbeta[4,5] = dRR_dbeta[2,5]

	# dRR2a/dbeta2
	dRR_dbeta[3,6] = exp(beta[3]+beta[4])*(1+exp(beta[1]))/(1+exp(beta[1]+beta[3]+beta[4]))/(1+exp(beta[1]+beta[3]+beta[4]))  
	# dRR2a/dbeta_age
	dRR_dbeta[4,6] = dRR_dbeta[3,6]

	# dRR3a/dbeta1
	dRR_dbeta[2,7] = exp(beta[2]+beta[3]+beta[4])*(1+exp(beta[1]))/(1+exp(beta[1]+beta[2]+beta[3]+beta[4]))/(1+exp(beta[1]+beta[2]+beta[3]+beta[4]))  
	# dRR3a/dbeta2
	dRR_dbeta[3,7] = dRR_dbeta[2,7]
	# dRR3a/dbeta_age
	dRR_dbeta[4,7] = dRR_dbeta[2,7]

	df_dRR = matrix(0,7,1)	
	for (i in c(2,4,6,8)){
		df_dRR[i-1,1] = pi[i]*a/b2
	}
	for (i in c(3,5,7)){
		df_dRR[i-1,1] = (pi[i]*a-b*(pi[i]+pi[i+1]))/b2

	}

	var_ppar = t(df_dRR) %*% t(dRR_dbeta) %*% vcov1[1:4,1:4] %*% dRR_dbeta %*% df_dRR +
			t(df_dpi) %*% t(dpi_dxeta) %*% vcov1[17:23,17:23] %*% dpi_dxeta %*% df_dpi +
			2* ( t(df_dpi) %*% t(dpi_dxeta) %*% vcov1[17:23,1:4] %*% dRR_dbeta %*% df_dRR )
	

	CI_ppar[3,1] = ppar[3,1] - qnorm(0.975)*sqrt(var_ppar)
	CI_ppar[3,2] = ppar[3,1] + qnorm(0.975)*sqrt(var_ppar)


print(ppar)
print(CI_ppar)


write.csv(cbind(ppar,CI_ppar),"pPAR_and_CI.csv")




####################################################################################################################################################################################
################## define the functions ############################################################################################################################################
####################################################################################################################################################################################

ivs.doubletransportability = function(betatheta){
	nz = ncol(zm) 			# nz is the number of columns of the Z matrix in the Main study, i.e. the number of mismeasured variables
	nc = length(table(cm1)) 	# nc is the number of categories of AGE
	n = 1+nz+nc-1				# n is the total number of parameters in f1
	beta0 = betatheta[1]
	beta1 = betatheta[2]
	beta2 = betatheta[3]
	betac1 = betatheta[4]
	gamma = betatheta[(n+1):(n+19)]
	 gamma10 = gamma[1]
	 gamma20 = gamma[2]
	 gamma30 = gamma[3]
	 gamma11 = gamma[4]
	 gamma21 = gamma[5]	
	 gamma31 = gamma[6]
	 gamma12 = gamma[7]
	 gamma22 = gamma[8]
	 gamma32 = gamma[9]
	 gamma13 = gamma[10]
	 gamma23 = gamma[11]
	 gamma33 = gamma[12]	 
	xeta1 = gamma[13]
	xeta2 = gamma[14]
	xeta3 = gamma[15]
	xeta4 = gamma[16]
	xeta5 = gamma[17]
	xeta6 = gamma[18]
	xeta7 = gamma[19]
	sum_exp_xeta = sum(exp(gamma[13:19]))+1
	sumlogf1vs = sum( yv * (cbind(1,xv,cv1)%*%c(beta0,beta1,beta2,betac1)   ) ) - sum(log(1+exp(cbind(1,xv,cv1)%*%c(beta0,beta1,beta2,betac1) )))

	zv1 = zv[,1]
	zv2 = zv[,2]
	xv1 = xv[,1]
	xv2 = xv[,2]
	f2vs =   (           1/(1+exp(gamma10)+exp(gamma20)+exp(gamma30))*as.numeric(xv1==0)*as.numeric(xv2==0) +            1/(1+exp(gamma11)+exp(gamma21)+exp(gamma31))*as.numeric(xv1==1)*as.numeric(xv2==0)                                                + 1/(1+exp(gamma12)+exp(gamma22)+exp(gamma32))*as.numeric(xv1==0)*as.numeric(xv2==1)                                               +1/(1+exp(gamma13)+exp(gamma23)+exp(gamma33))*as.numeric(xv1==1)*as.numeric(xv2==1))*as.numeric(zv1==0)*as.numeric(zv2==0)+
		 (exp(gamma10)/(1+exp(gamma10)+exp(gamma20)+exp(gamma30))*as.numeric(xv1==0)*as.numeric(xv2==0) + exp(gamma11)/(1+exp(gamma11)+exp(gamma21)+exp(gamma31))*as.numeric(xv1==1)*as.numeric(xv2==0)                                     + exp(gamma12)/(1+exp(gamma12)+exp(gamma22)+exp(gamma32))*as.numeric(xv1==0)*as.numeric(xv2==1)                                    +exp(gamma13)/(1+exp(gamma13)+exp(gamma23)+exp(gamma33))*as.numeric(xv1==1)*as.numeric(xv2==1))*as.numeric(zv1==1)*as.numeric(zv2==0)+
		 (exp(gamma20)/(1+exp(gamma10)+exp(gamma20)+exp(gamma30))*as.numeric(xv1==0)*as.numeric(xv2==0) + exp(gamma21)/(1+exp(gamma11)+exp(gamma21)+exp(gamma31))*as.numeric(xv1==1)*as.numeric(xv2==0)                                     + exp(gamma22)/(1+exp(gamma12)+exp(gamma22)+exp(gamma32))*as.numeric(xv1==0)*as.numeric(xv2==1)                                    +exp(gamma23)/(1+exp(gamma13)+exp(gamma23)+exp(gamma33))*as.numeric(xv1==1)*as.numeric(xv2==1))*as.numeric(zv1==0)*as.numeric(zv2==1)+
		 (exp(gamma30)/(1+exp(gamma10)+exp(gamma20)+exp(gamma30))*as.numeric(xv1==0)*as.numeric(xv2==0) + exp(gamma31)/(1+exp(gamma11)+exp(gamma21)+exp(gamma31))*as.numeric(xv1==1)*as.numeric(xv2==0)                                     + exp(gamma32)/(1+exp(gamma12)+exp(gamma22)+exp(gamma32))*as.numeric(xv1==0)*as.numeric(xv2==1)                                    +exp(gamma33)/(1+exp(gamma13)+exp(gamma23)+exp(gamma33))*as.numeric(xv1==1)*as.numeric(xv2==1))*as.numeric(zv1==1)*as.numeric(zv2==1)
      sumlogf2vs = sum(log(f2vs))			# add up the log-likelihoods over all Z
	
	
	f6vs = 		 1/(sum_exp_xeta) * as.numeric(xv1==0)*as.numeric(xv2==0)*as.numeric(cv1==0) +
		exp(xeta1)/(sum_exp_xeta) * as.numeric(xv1==1)*as.numeric(xv2==0)*as.numeric(cv1==0) +
		exp(xeta2)/(sum_exp_xeta) * as.numeric(xv1==0)*as.numeric(xv2==1)*as.numeric(cv1==0) +
		exp(xeta3)/(sum_exp_xeta) * as.numeric(xv1==1)*as.numeric(xv2==1)*as.numeric(cv1==0) + 
		exp(xeta4)/(sum_exp_xeta) * as.numeric(xv1==0)*as.numeric(xv2==0)*as.numeric(cv1==1) +
		exp(xeta5)/(sum_exp_xeta) * as.numeric(xv1==1)*as.numeric(xv2==0)*as.numeric(cv1==1) +
		exp(xeta6)/(sum_exp_xeta) * as.numeric(xv1==0)*as.numeric(xv2==1)*as.numeric(cv1==1) +
		exp(xeta7)/(sum_exp_xeta) * as.numeric(xv1==1)*as.numeric(xv2==1)*as.numeric(cv1==1)  
	sumlogf6vs = sum(log(f6vs))			# add up the log-likelihoods over all Z

	age0 = cm1==0
	ym_0 = ym[age0]; zm_0 = zm[age0,]; cm_0 = cm1[age0]
	zm1 = zm_0[,1]
	zm2 = zm_0[,2]

	f3 = matrix(0,length(ym_0),(2^nz))
	

	x1 = 0; x2 = 0
	f3[,1] = exp(ym_0*cbind(1,x1,x2,cm_0)%*%c(beta0,beta1,beta2,0))/(1+exp(cbind(1,x1,x2,cm_0)%*%c(beta0,beta1,beta2,0))) * (   	           1/(1+exp(gamma10)+exp(gamma20)+exp(gamma30))*as.numeric(zm1==0)*as.numeric(zm2==0)	+
												exp(gamma10)/(1+exp(gamma10)+exp(gamma20)+exp(gamma30))*as.numeric(zm1==1)*as.numeric(zm2==0)	+
												exp(gamma20)/(1+exp(gamma10)+exp(gamma20)+exp(gamma30))*as.numeric(zm1==0)*as.numeric(zm2==1)	+
												exp(gamma30)/(1+exp(gamma10)+exp(gamma20)+exp(gamma30))*as.numeric(zm1==1)*as.numeric(zm2==1)	) * 1/(sum_exp_xeta)
	x1 = 1; x2 = 0
	f3[,2] = exp(ym_0*cbind(1,x1,x2,cm_0)%*%c(beta0,beta1,beta2,0))/(1+exp(cbind(1,x1,x2,cm_0)%*%c(beta0,beta1,beta2,0))) * (   	           1/(1+exp(gamma11)+exp(gamma21)+exp(gamma31))*as.numeric(zm1==0)*as.numeric(zm2==0)	+
												exp(gamma11)/(1+exp(gamma11)+exp(gamma21)+exp(gamma31))*as.numeric(zm1==1)*as.numeric(zm2==0)	+
												exp(gamma21)/(1+exp(gamma11)+exp(gamma21)+exp(gamma31))*as.numeric(zm1==0)*as.numeric(zm2==1)	+
												exp(gamma31)/(1+exp(gamma11)+exp(gamma21)+exp(gamma31))*as.numeric(zm1==1)*as.numeric(zm2==1)	) * exp(xeta1)/(sum_exp_xeta)
	x1 = 0; x2 = 1
	f3[,3] = exp(ym_0*cbind(1,x1,x2,cm_0)%*%c(beta0,beta1,beta2,0))/(1+exp(cbind(1,x1,x2,cm_0)%*%c(beta0,beta1,beta2,0))) * (   	           1/(1+exp(gamma12)+exp(gamma22)+exp(gamma32))*as.numeric(zm1==0)*as.numeric(zm2==0)	+
												exp(gamma12)/(1+exp(gamma12)+exp(gamma22)+exp(gamma32))*as.numeric(zm1==1)*as.numeric(zm2==0)	+
												exp(gamma22)/(1+exp(gamma12)+exp(gamma22)+exp(gamma32))*as.numeric(zm1==0)*as.numeric(zm2==1)	+
												exp(gamma32)/(1+exp(gamma12)+exp(gamma22)+exp(gamma32))*as.numeric(zm1==1)*as.numeric(zm2==1)	) * exp(xeta2)/(sum_exp_xeta)
	x1 = 1; x2 = 1
	f3[,4] = exp(ym_0*cbind(1,x1,x2,cm_0)%*%c(beta0,beta1,beta2,0))/(1+exp(cbind(1,x1,x2,cm_0)%*%c(beta0,beta1,beta2,0))) * (   	           1/(1+exp(gamma13)+exp(gamma23)+exp(gamma33))*as.numeric(zm1==0)*as.numeric(zm2==0)	+
												exp(gamma13)/(1+exp(gamma13)+exp(gamma23)+exp(gamma33))*as.numeric(zm1==1)*as.numeric(zm2==0)	+
												exp(gamma23)/(1+exp(gamma13)+exp(gamma23)+exp(gamma33))*as.numeric(zm1==0)*as.numeric(zm2==1)	+
												exp(gamma33)/(1+exp(gamma13)+exp(gamma23)+exp(gamma33))*as.numeric(zm1==1)*as.numeric(zm2==1)	) * exp(xeta3)/(sum_exp_xeta)
	logf3 = log(apply(f3,1,sum))
	sumlogf3a = sum(logf3)	

	age1 = cm1==1
	ym_1 = ym[age1]; zm_1 = zm[age1,]; cm_1 = cm1[age1]
	zm1 = zm_1[,1]
	zm2 = zm_1[,2]
	f3 = matrix(0,length(ym_1),(2^nz))

	x1 = 0; x2 = 0
	f3[,1] = exp(ym_1*cbind(1,x1,x2,cm_1)%*%c(beta0,beta1,beta2,betac1))/(1+exp(cbind(1,x1,x2,cm_1)%*%c(beta0,beta1,beta2,betac1))) * (   	           1/(1+exp(gamma10)+exp(gamma20)+exp(gamma30))*as.numeric(zm1==0)*as.numeric(zm2==0)	+
												exp(gamma10)/(1+exp(gamma10)+exp(gamma20)+exp(gamma30))*as.numeric(zm1==1)*as.numeric(zm2==0)	+
												exp(gamma20)/(1+exp(gamma10)+exp(gamma20)+exp(gamma30))*as.numeric(zm1==0)*as.numeric(zm2==1)	+
												exp(gamma30)/(1+exp(gamma10)+exp(gamma20)+exp(gamma30))*as.numeric(zm1==1)*as.numeric(zm2==1)	) * exp(xeta4)/(sum_exp_xeta)
	x1 = 1; x2 = 0
	f3[,2] = exp(ym_1*cbind(1,x1,x2,cm_1)%*%c(beta0,beta1,beta2,betac1))/(1+exp(cbind(1,x1,x2,cm_1)%*%c(beta0,beta1,beta2,betac1))) * (   	           1/(1+exp(gamma11)+exp(gamma21)+exp(gamma31))*as.numeric(zm1==0)*as.numeric(zm2==0)	+
												exp(gamma11)/(1+exp(gamma11)+exp(gamma21)+exp(gamma31))*as.numeric(zm1==1)*as.numeric(zm2==0)	+
												exp(gamma21)/(1+exp(gamma11)+exp(gamma21)+exp(gamma31))*as.numeric(zm1==0)*as.numeric(zm2==1)	+
												exp(gamma31)/(1+exp(gamma11)+exp(gamma21)+exp(gamma31))*as.numeric(zm1==1)*as.numeric(zm2==1)	) * exp(xeta5)/(sum_exp_xeta)
	x1 = 0; x2 = 1
	f3[,3] = exp(ym_1*cbind(1,x1,x2,cm_1)%*%c(beta0,beta1,beta2,betac1))/(1+exp(cbind(1,x1,x2,cm_1)%*%c(beta0,beta1,beta2,betac1))) * (   	           1/(1+exp(gamma12)+exp(gamma22)+exp(gamma32))*as.numeric(zm1==0)*as.numeric(zm2==0)	+
												exp(gamma12)/(1+exp(gamma12)+exp(gamma22)+exp(gamma32))*as.numeric(zm1==1)*as.numeric(zm2==0)	+
												exp(gamma22)/(1+exp(gamma12)+exp(gamma22)+exp(gamma32))*as.numeric(zm1==0)*as.numeric(zm2==1)	+
												exp(gamma32)/(1+exp(gamma12)+exp(gamma22)+exp(gamma32))*as.numeric(zm1==1)*as.numeric(zm2==1)	) * exp(xeta6)/(sum_exp_xeta)
	x1 = 1; x2 = 1
	f3[,4] = exp(ym_1*cbind(1,x1,x2,cm_1)%*%c(beta0,beta1,beta2,betac1))/(1+exp(cbind(1,x1,x2,cm_1)%*%c(beta0,beta1,beta2,betac1))) * (   	           1/(1+exp(gamma13)+exp(gamma23)+exp(gamma33))*as.numeric(zm1==0)*as.numeric(zm2==0)	+
												exp(gamma13)/(1+exp(gamma13)+exp(gamma23)+exp(gamma33))*as.numeric(zm1==1)*as.numeric(zm2==0)	+
												exp(gamma23)/(1+exp(gamma13)+exp(gamma23)+exp(gamma33))*as.numeric(zm1==0)*as.numeric(zm2==1)	+
												exp(gamma33)/(1+exp(gamma13)+exp(gamma23)+exp(gamma33))*as.numeric(zm1==1)*as.numeric(zm2==1)	) * exp(xeta7)/(sum_exp_xeta)
	logf3 = log(apply(f3,1,sum))
	sumlogf3 = sum(logf3)	+ sumlogf3a

	-(sumlogf1vs+sumlogf2vs+sumlogf6vs+sumlogf3)			# return the negative of the sum of the log-likelihoods
}




ivs.singletransportability = function(betatheta){
	nz = ncol(zm) 			# nz is the number of columns of the Z matrix in the Main study, i.e. the number of mismeasured variables
	nc = length(table(cm1)) 	# nc is the number of categories of AGE
	n = 1+nz+nc-1				# n is the total number of parameters in f1
	beta0 = betatheta[1]
	beta1 = betatheta[2]
	beta2 = betatheta[3]
	betac1 = betatheta[4]
	gamma = betatheta[(n+1):(n+19)]
	 gamma10 = gamma[1]
	 gamma20 = gamma[2]
	 gamma30 = gamma[3]
	 gamma11 = gamma[4]
	 gamma21 = gamma[5]	
	 gamma31 = gamma[6]
	 gamma12 = gamma[7]
	 gamma22 = gamma[8]
	 gamma32 = gamma[9]
	 gamma13 = gamma[10]
	 gamma23 = gamma[11]
	 gamma33 = gamma[12]	 
	xeta1 = gamma[13]
	xeta2 = gamma[14]
	xeta3 = gamma[15]
	xeta4 = gamma[16]
	xeta5 = gamma[17]
	xeta6 = gamma[18]
	xeta7 = gamma[19]
	sum_exp_xeta = sum(exp(gamma[13:19]))+1
	sumlogf1vs = sum( yv * (cbind(1,xv,cv1)%*%c(beta0,beta1,beta2,betac1)   ) ) - sum(log(1+exp(cbind(1,xv,cv1)%*%c(beta0,beta1,beta2,betac1) )))

	zv1 = zv[,1]
	zv2 = zv[,2]
	xv1 = xv[,1]
	xv2 = xv[,2]
	f2vs =   (           1/(1+exp(gamma10)+exp(gamma20)+exp(gamma30))*as.numeric(xv1==0)*as.numeric(xv2==0) +            1/(1+exp(gamma11)+exp(gamma21)+exp(gamma31))*as.numeric(xv1==1)*as.numeric(xv2==0)                                                + 1/(1+exp(gamma12)+exp(gamma22)+exp(gamma32))*as.numeric(xv1==0)*as.numeric(xv2==1)                                               +1/(1+exp(gamma13)+exp(gamma23)+exp(gamma33))*as.numeric(xv1==1)*as.numeric(xv2==1))*as.numeric(zv1==0)*as.numeric(zv2==0)+
		 (exp(gamma10)/(1+exp(gamma10)+exp(gamma20)+exp(gamma30))*as.numeric(xv1==0)*as.numeric(xv2==0) + exp(gamma11)/(1+exp(gamma11)+exp(gamma21)+exp(gamma31))*as.numeric(xv1==1)*as.numeric(xv2==0)                                     + exp(gamma12)/(1+exp(gamma12)+exp(gamma22)+exp(gamma32))*as.numeric(xv1==0)*as.numeric(xv2==1)                                    +exp(gamma13)/(1+exp(gamma13)+exp(gamma23)+exp(gamma33))*as.numeric(xv1==1)*as.numeric(xv2==1))*as.numeric(zv1==1)*as.numeric(zv2==0)+
		 (exp(gamma20)/(1+exp(gamma10)+exp(gamma20)+exp(gamma30))*as.numeric(xv1==0)*as.numeric(xv2==0) + exp(gamma21)/(1+exp(gamma11)+exp(gamma21)+exp(gamma31))*as.numeric(xv1==1)*as.numeric(xv2==0)                                     + exp(gamma22)/(1+exp(gamma12)+exp(gamma22)+exp(gamma32))*as.numeric(xv1==0)*as.numeric(xv2==1)                                    +exp(gamma23)/(1+exp(gamma13)+exp(gamma23)+exp(gamma33))*as.numeric(xv1==1)*as.numeric(xv2==1))*as.numeric(zv1==0)*as.numeric(zv2==1)+
		 (exp(gamma30)/(1+exp(gamma10)+exp(gamma20)+exp(gamma30))*as.numeric(xv1==0)*as.numeric(xv2==0) + exp(gamma31)/(1+exp(gamma11)+exp(gamma21)+exp(gamma31))*as.numeric(xv1==1)*as.numeric(xv2==0)                                     + exp(gamma32)/(1+exp(gamma12)+exp(gamma22)+exp(gamma32))*as.numeric(xv1==0)*as.numeric(xv2==1)                                    +exp(gamma33)/(1+exp(gamma13)+exp(gamma23)+exp(gamma33))*as.numeric(xv1==1)*as.numeric(xv2==1))*as.numeric(zv1==1)*as.numeric(zv2==1)
      sumlogf2vs = sum(log(f2vs))			# add up the log-likelihoods over all Z
	
	
	age0 = cm1==0
	ym_0 = ym[age0]; zm_0 = zm[age0,]; cm_0 = cm1[age0]
	zm1 = zm_0[,1]
	zm2 = zm_0[,2]

	f3 = matrix(0,length(ym_0),(2^nz))
	

	x1 = 0; x2 = 0
	f3[,1] = exp(ym_0*cbind(1,x1,x2,cm_0)%*%c(beta0,beta1,beta2,0))/(1+exp(cbind(1,x1,x2,cm_0)%*%c(beta0,beta1,beta2,0))) * (   	           1/(1+exp(gamma10)+exp(gamma20)+exp(gamma30))*as.numeric(zm1==0)*as.numeric(zm2==0)	+
												exp(gamma10)/(1+exp(gamma10)+exp(gamma20)+exp(gamma30))*as.numeric(zm1==1)*as.numeric(zm2==0)	+
												exp(gamma20)/(1+exp(gamma10)+exp(gamma20)+exp(gamma30))*as.numeric(zm1==0)*as.numeric(zm2==1)	+
												exp(gamma30)/(1+exp(gamma10)+exp(gamma20)+exp(gamma30))*as.numeric(zm1==1)*as.numeric(zm2==1)	) * 1/(sum_exp_xeta)
	x1 = 1; x2 = 0
	f3[,2] = exp(ym_0*cbind(1,x1,x2,cm_0)%*%c(beta0,beta1,beta2,0))/(1+exp(cbind(1,x1,x2,cm_0)%*%c(beta0,beta1,beta2,0))) * (   	           1/(1+exp(gamma11)+exp(gamma21)+exp(gamma31))*as.numeric(zm1==0)*as.numeric(zm2==0)	+
												exp(gamma11)/(1+exp(gamma11)+exp(gamma21)+exp(gamma31))*as.numeric(zm1==1)*as.numeric(zm2==0)	+
												exp(gamma21)/(1+exp(gamma11)+exp(gamma21)+exp(gamma31))*as.numeric(zm1==0)*as.numeric(zm2==1)	+
												exp(gamma31)/(1+exp(gamma11)+exp(gamma21)+exp(gamma31))*as.numeric(zm1==1)*as.numeric(zm2==1)	) * exp(xeta1)/(sum_exp_xeta)
	x1 = 0; x2 = 1
	f3[,3] = exp(ym_0*cbind(1,x1,x2,cm_0)%*%c(beta0,beta1,beta2,0))/(1+exp(cbind(1,x1,x2,cm_0)%*%c(beta0,beta1,beta2,0))) * (   	           1/(1+exp(gamma12)+exp(gamma22)+exp(gamma32))*as.numeric(zm1==0)*as.numeric(zm2==0)	+
												exp(gamma12)/(1+exp(gamma12)+exp(gamma22)+exp(gamma32))*as.numeric(zm1==1)*as.numeric(zm2==0)	+
												exp(gamma22)/(1+exp(gamma12)+exp(gamma22)+exp(gamma32))*as.numeric(zm1==0)*as.numeric(zm2==1)	+
												exp(gamma32)/(1+exp(gamma12)+exp(gamma22)+exp(gamma32))*as.numeric(zm1==1)*as.numeric(zm2==1)	) * exp(xeta2)/(sum_exp_xeta)
	x1 = 1; x2 = 1
	f3[,4] = exp(ym_0*cbind(1,x1,x2,cm_0)%*%c(beta0,beta1,beta2,0))/(1+exp(cbind(1,x1,x2,cm_0)%*%c(beta0,beta1,beta2,0))) * (   	           1/(1+exp(gamma13)+exp(gamma23)+exp(gamma33))*as.numeric(zm1==0)*as.numeric(zm2==0)	+
												exp(gamma13)/(1+exp(gamma13)+exp(gamma23)+exp(gamma33))*as.numeric(zm1==1)*as.numeric(zm2==0)	+
												exp(gamma23)/(1+exp(gamma13)+exp(gamma23)+exp(gamma33))*as.numeric(zm1==0)*as.numeric(zm2==1)	+
												exp(gamma33)/(1+exp(gamma13)+exp(gamma23)+exp(gamma33))*as.numeric(zm1==1)*as.numeric(zm2==1)	) * exp(xeta3)/(sum_exp_xeta)
	logf3 = log(apply(f3,1,sum))
	sumlogf3a = sum(logf3)	

	age1 = cm1==1
	ym_1 = ym[age1]; zm_1 = zm[age1,]; cm_1 = cm1[age1]
	zm1 = zm_1[,1]
	zm2 = zm_1[,2]
	f3 = matrix(0,length(ym_1),(2^nz))

	x1 = 0; x2 = 0
	f3[,1] = exp(ym_1*cbind(1,x1,x2,cm_1)%*%c(beta0,beta1,beta2,betac1))/(1+exp(cbind(1,x1,x2,cm_1)%*%c(beta0,beta1,beta2,betac1))) * (   	           1/(1+exp(gamma10)+exp(gamma20)+exp(gamma30))*as.numeric(zm1==0)*as.numeric(zm2==0)	+
												exp(gamma10)/(1+exp(gamma10)+exp(gamma20)+exp(gamma30))*as.numeric(zm1==1)*as.numeric(zm2==0)	+
												exp(gamma20)/(1+exp(gamma10)+exp(gamma20)+exp(gamma30))*as.numeric(zm1==0)*as.numeric(zm2==1)	+
												exp(gamma30)/(1+exp(gamma10)+exp(gamma20)+exp(gamma30))*as.numeric(zm1==1)*as.numeric(zm2==1)	) * exp(xeta4)/(sum_exp_xeta)
	x1 = 1; x2 = 0
	f3[,2] = exp(ym_1*cbind(1,x1,x2,cm_1)%*%c(beta0,beta1,beta2,betac1))/(1+exp(cbind(1,x1,x2,cm_1)%*%c(beta0,beta1,beta2,betac1))) * (   	           1/(1+exp(gamma11)+exp(gamma21)+exp(gamma31))*as.numeric(zm1==0)*as.numeric(zm2==0)	+
												exp(gamma11)/(1+exp(gamma11)+exp(gamma21)+exp(gamma31))*as.numeric(zm1==1)*as.numeric(zm2==0)	+
												exp(gamma21)/(1+exp(gamma11)+exp(gamma21)+exp(gamma31))*as.numeric(zm1==0)*as.numeric(zm2==1)	+
												exp(gamma31)/(1+exp(gamma11)+exp(gamma21)+exp(gamma31))*as.numeric(zm1==1)*as.numeric(zm2==1)	) * exp(xeta5)/(sum_exp_xeta)
	x1 = 0; x2 = 1
	f3[,3] = exp(ym_1*cbind(1,x1,x2,cm_1)%*%c(beta0,beta1,beta2,betac1))/(1+exp(cbind(1,x1,x2,cm_1)%*%c(beta0,beta1,beta2,betac1))) * (   	           1/(1+exp(gamma12)+exp(gamma22)+exp(gamma32))*as.numeric(zm1==0)*as.numeric(zm2==0)	+
												exp(gamma12)/(1+exp(gamma12)+exp(gamma22)+exp(gamma32))*as.numeric(zm1==1)*as.numeric(zm2==0)	+
												exp(gamma22)/(1+exp(gamma12)+exp(gamma22)+exp(gamma32))*as.numeric(zm1==0)*as.numeric(zm2==1)	+
												exp(gamma32)/(1+exp(gamma12)+exp(gamma22)+exp(gamma32))*as.numeric(zm1==1)*as.numeric(zm2==1)	) * exp(xeta6)/(sum_exp_xeta)
	x1 = 1; x2 = 1
	f3[,4] = exp(ym_1*cbind(1,x1,x2,cm_1)%*%c(beta0,beta1,beta2,betac1))/(1+exp(cbind(1,x1,x2,cm_1)%*%c(beta0,beta1,beta2,betac1))) * (   	           1/(1+exp(gamma13)+exp(gamma23)+exp(gamma33))*as.numeric(zm1==0)*as.numeric(zm2==0)	+
												exp(gamma13)/(1+exp(gamma13)+exp(gamma23)+exp(gamma33))*as.numeric(zm1==1)*as.numeric(zm2==0)	+
												exp(gamma23)/(1+exp(gamma13)+exp(gamma23)+exp(gamma33))*as.numeric(zm1==0)*as.numeric(zm2==1)	+
												exp(gamma33)/(1+exp(gamma13)+exp(gamma23)+exp(gamma33))*as.numeric(zm1==1)*as.numeric(zm2==1)	) * exp(xeta7)/(sum_exp_xeta)
	logf3 = log(apply(f3,1,sum))
	sumlogf3 = sum(logf3)	+ sumlogf3a
	-(sumlogf1vs+sumlogf2vs+sumlogf3)			# return the negative of the sum of the log-likelihoods
}



	# ====================================================================================================================================================================
	# ==== Calculate the Uncorrected CI ==================================================================================================================================
	# ====================================================================================================================================================================
	## calculate the naive joint prevalences
	cov_pz = matrix(NA,8,8)
	nm = nrow(zm) + nrow(zv)
	for(i in 1:8){
		cov_pz[i,i] = pi[i]*(1-pi[i])/nm
	}
	for(i in 1:7){
		for(j in (i+1):8){
			cov_pz[i,j] = - pi[i]*pi[j]/nm
			cov_pz[j,i] = cov_pz[i,j]
		}
	}

	df_dpi = matrix(0,8,1)
	 RR_vec = c(1,RR1,RR2,RR3,RR0a,RR1a,RR2a,RR3a)
	 RR_num = c(1,1,RR2,RR2,RR0a,RR0a,RR2a,RR2a)
	 a = sum(pi*RR_num)
	 b = sum(pi*RR_vec)
	 b2 = b^2
	for(i in 1:8){
		df_dpi[i,1] = (b*RR_num[i] - a*RR_vec[i])/(b2)
	}

	### dRR_dbeta = matrix(0,3,2)
	# dRR1/dbeta0
	# dRR_dbeta[1,1] = (1-exp(beta[2]))*(exp(beta[1]+beta[2]))/(1+exp(beta[1]+beta[2]))/(1+exp(beta[1]+beta[2]))
	# dRR2/dbeta0
	# dRR_dbeta[1,2] = (1-exp(beta[3]))*(exp(beta[1]+beta[3]))/(1+exp(beta[1]+beta[3]))/(1+exp(beta[1]+beta[3]))
	# dRR1/dbeta1
	# dRR_dbeta[2,1] = exp(beta[2])*(1+exp(beta[1]))/(1+exp(beta[1]+beta[2]))/(1+exp(beta[1]+beta[2]))  
	# dRR2/dbeta2
	# dRR_dbeta[3,2] = exp(beta[3])*(1+exp(beta[1]))/(1+exp(beta[1]+beta[3]))/(1+exp(beta[1]+beta[3]))  

	dRR_dbeta = matrix(0,4,7)
	# dRR1/dbeta0
	dRR_dbeta[1,1] = (1-exp(beta[2]))*(exp(beta[1]+beta[2]))/(1+exp(beta[1]+beta[2]))/(1+exp(beta[1]+beta[2]))
	# dRR2/dbeta0
	dRR_dbeta[1,2] = (1-exp(beta[3]))*(exp(beta[1]+beta[3]))/(1+exp(beta[1]+beta[3]))/(1+exp(beta[1]+beta[3]))
	# dRR3/dbeta0
	dRR_dbeta[1,3] = (1-exp(beta[2]+beta[3]))*(exp(beta[1]+beta[2]+beta[3]))/(1+exp(beta[1]+beta[2]+beta[3]))/(1+exp(beta[1]+beta[2]+beta[3]))
	# dRR0a/dbeta0
	dRR_dbeta[1,4] = (1-exp(beta[4]))*(exp(beta[1]+beta[4]))/(1+exp(beta[1]+beta[4]))/(1+exp(beta[1]+beta[4]))
	# dRR1a/dbeta0
	dRR_dbeta[1,5] = (1-exp(beta[2]+beta[4]))*(exp(beta[1]+beta[2]+beta[4]))/(1+exp(beta[1]+beta[2]+beta[4]))/(1+exp(beta[1]+beta[2]+beta[4]))
	# dRR2a/dbeta0
	dRR_dbeta[1,6] = (1-exp(beta[3]+beta[4]))*(exp(beta[1]+beta[3]+beta[4]))/(1+exp(beta[1]+beta[3]+beta[4]))/(1+exp(beta[1]+beta[3]+beta[4]))
	# dRR3a/dbeta0
	dRR_dbeta[1,7] = (1-exp(beta[2]+beta[3]+beta[4]))*(exp(beta[1]+beta[2]+beta[3]+beta[4]))/(1+exp(beta[1]+beta[2]+beta[3]+beta[4]))/(1+exp(beta[1]+beta[2]+beta[3]+beta[4]))
	
	# dRR1/dbeta1
	dRR_dbeta[2,1] = exp(beta[2])*(1+exp(beta[1]))/(1+exp(beta[1]+beta[2]))/(1+exp(beta[1]+beta[2]))  
	
	# dRR2/dbeta2
	dRR_dbeta[3,2] = exp(beta[3])*(1+exp(beta[1]))/(1+exp(beta[1]+beta[3]))/(1+exp(beta[1]+beta[3]))  

	# dRR3/dbeta1
	dRR_dbeta[2,3] = exp(beta[2]+beta[3])*(1+exp(beta[1]))/(1+exp(beta[1]+beta[2]+beta[3]))/(1+exp(beta[1]+beta[2]+beta[3]))  
	# dRR3/dbeta2
	dRR_dbeta[3,3] = dRR_dbeta[2,3]

	# dRR0a/dbeta_age
	dRR_dbeta[4,4] = exp(beta[4])*(1+exp(beta[1]))/(1+exp(beta[1]+beta[4]))/(1+exp(beta[1]+beta[4]))  

	# dRR1a/dbeta1
	dRR_dbeta[2,5] = exp(beta[2]+beta[4])*(1+exp(beta[1]))/(1+exp(beta[1]+beta[2]+beta[4]))/(1+exp(beta[1]+beta[2]+beta[4]))  
	# dRR1a/dbeta_age
	dRR_dbeta[4,5] = dRR_dbeta[2,5]

	# dRR2a/dbeta2
	dRR_dbeta[3,6] = exp(beta[3]+beta[4])*(1+exp(beta[1]))/(1+exp(beta[1]+beta[3]+beta[4]))/(1+exp(beta[1]+beta[3]+beta[4]))  
	# dRR2a/dbeta_age
	dRR_dbeta[4,6] = dRR_dbeta[3,6]

	# dRR3a/dbeta1
	dRR_dbeta[2,7] = exp(beta[2]+beta[3]+beta[4])*(1+exp(beta[1]))/(1+exp(beta[1]+beta[2]+beta[3]+beta[4]))/(1+exp(beta[1]+beta[2]+beta[3]+beta[4]))  
	# dRR3a/dbeta2
	dRR_dbeta[3,7] = dRR_dbeta[2,7]
	# dRR3a/dbeta_age
	dRR_dbeta[4,7] = dRR_dbeta[2,7]

	df_dRR = matrix(0,7,1)	
	for (i in c(2,4,6,8)){
		df_dRR[i-1,1] = pi[i]*a/b2
	}
	for (i in c(3,5,7)){
		df_dRR[i-1,1] = (pi[i]*a-b*(pi[i]+pi[i+1]))/b2

	}

	var_ppar = t(df_dRR) %*% t(dRR_dbeta) %*% vcov(fit1) %*% dRR_dbeta %*% df_dRR +
			t(df_dpi) %*%  cov_pz %*%  df_dpi 

	
	CI_ppar[1,1] = ppar[1,1] - qnorm(0.975)*sqrt(var_ppar)
	CI_ppar[1,2] = ppar[1,1] + qnorm(0.975)*sqrt(var_ppar)
