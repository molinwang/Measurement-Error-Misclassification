########################
# regCalibCRS function #
########################

# authors: Wenze Tang and Molin Wang
# version 1: Jan 5, 2022

# Abstract:
# This function corrects for measurement error in exposure or exposure and covariates and gives corrected coefficients associated standard errors, p values
# as well as variance-covariance matrix of corrected coefficients. Linear model for calibration regression and generalized linear model for outcome regression
# are currently supported under external validation study design. Standard errors are obtained by bootstrapping. Non-linear terms such as interaction terms 
# need to be included as a permanent covariate in input datasets as opposed to be included in the regression formula. 
# A validation study is required to empirically characterize the measurement error calibration model. Options are given for main study/external validation study 
# design (Carroll, Ruppert, Stefanski and Cainiceanu; 2004).

# Reference:
# Carroll, Ruppert, Stefanski and Cainiceanu. "Measurement error in nonlinear models". Chapman & Hall/CRC. 2006.

#' Find corrected parameters (both single or multiple) and their standard deviations and p values for 
#' mismeasured continuous exposure
#' @param ms The input main study data set as dataframe. This dataframe should minimally include variables specified in 
#' `sur`, `outcome` and `covOutcome` (if any).
#' @param vs The input internal/external validation data set as dataframe. This dataframe should minimally include 
#' variables indicated in `exp`, `sur` and `covCalib` (if any).
#' @param sur character vector of mismeasured exposure and covariates (i.e. surrogates) in the main study dataset.
#' @param exp character vector of correctly measured exposure and covariates that has a one-to-one correspondence to those 
#' specified in `sur` in validation dataset. Must have same length as `sur`.
#' @param vsIndicator character of the indicator variable in main study indicating the subject having a validation record, ==1 if having validation record. 
#' Only required when external==FALSE. 
#' @param covCalib=NULL character vector of names of correctly measured covariates to adjust for in calibration model, included non-linear terms.
#' @param covOutcome=NULL character vector of names of correctly measured covariates to adjust for in outcome model with corrected exposure, 
#' included non-linear terms. 
#' @param outcome Outcome variable.
#' @param method="lm" Methods for outcome modeling, currently only `lm` or `glm` methods are available. Required.    
#' @param external=TRUE Indicates whether `vs` is an external validation set. If external=FALSE, then vs must contain
#'  variable specified in `outcome`. If external=FALSE, then user cannot supply uncorrected estimates, i.e. supplyEstimates=FALSE.
#' @param family Supply family parameter to pass to glm function. Not a character. Required if method="glm".
#' @param link Supply link parameter to pass to glm function. Should be character. Required if method="glm". 
#' @return printable dataframe from standard regression results (when supplyEstimates==FALSE) as well as corrected results


regCalibCRS<-function(ms,vs, 
                      sur, exp, vsIndicator,covCalib=NULL, covOutcome=NULL, outcome=NA, 
                      method="lm",family=NA,link=NA,external=TRUE){
  require("stats")
  require("earth")
  require("matrixcalc")
  require("Matrix")
  require("dplyr")
  require("R.utils")
  require("glue")
  
  # show intermediate quantities controller
  show_detail = FALSE
  
  ###################  
  # check arguments #
  ###################
  ## data related warnings
  if(external==TRUE){
    if(missing(vs)){
      stop("Input data vs not supplied.")
    }else if(class(vs)!="data.frame"){
      stop("Input data vs must be of data.frame class.")
    }
  }

  
  if(missing(ms)){
    stop("Input data ms not supplied.")
  }else if(class(ms)!="data.frame"){
    stop("Input data ms must be of data.frame class.")
  }
  
  if(missing(sur)|missing(exp)){
    stop("Missing exposure variable.")
  }else if(class(sur)!="character"|class(exp)!="character"){
    stop("mExp or exp is not supplied with character vector.")
  }else if(length(sur)!=length(exp)){
    stop("Length of correctly measured variables differs from length of mismeasured variables.")
  }
  
  if(missing(covCalib)){
    warning("No covariates supplied.")
  }else if(length(covCalib)!=0&class(covCalib)!="character"|(length(covOutcome)!=0&class(covOutcome)!="character")){
    stop("covCalib or covOutcome is not supplied with character vector.")
  }
  
  if(missing(outcome)){
    stop("Outcome is missing.")
  }else if(class(outcome)!="character"|outcome==""|outcome==" "){
    stop("outcome is not supplied with appropriate character.")
  }
  
  # if(method!="lm"){
  #   if(is.na(family)|is.na(link)){
  #     stop("You must supply `family` and `link` parameters if you do not wish to use least square linear outcome model.")
  #   }
  # }
  
  ## 1. check if MS contains data indicated by sur and covOutcome
  
  # if(length(covOutcome)>0&length(covCalib)>0){
  #   MSVars_spec<-(c(sur,covCalib,covOutcome))
  # }else if(length(covCalib)==0&length(covOutcome)==0){
  #   MSVars_spec<-(c(sur))
  # }
  MSVars_spec<-(c(sur,covOutcome))
  MSVars<-colnames(ms)
  inMSVars<-(MSVars_spec%in%MSVars)
  if(sum(inMSVars)!=length(MSVars_spec)){
    stop("Main study dataset does not contain all the necessary variables specified in one of the following parameter: id, sur, covCalib and covOutcome.")
  } 
  
  
  ## 2. check if EVS contains data indicate by sur, exp and covCalib; for IVS, check additionally for outcome
  if(external==TRUE){
    EVSVars<-colnames(vs)
    EVSVars_spec<-(c(sur,exp,covCalib))
    inMSVars<-(EVSVars_spec%in%EVSVars)
    if(sum(inMSVars)!=length(EVSVars_spec)){
      stop("Validation study dataset does not contain all the necessary variables specified in one of the following parameter: id, sur, exp and covCalib.")
    }
  }

  
  ######################
  # Embedded functions #
  ######################
  #### create a function similar to the design function in SAS code
  # design<-function(TMatrix,PVector){
  #   zeroMatrix<-matrix(rep(0,(PVector^2)^2),nrow=PVector^2)
  #   # replace 0 in the zeroMatrix's positions indicated in TMatrix as 1
  #   for(i in 1:length(TMatrix)){
  #     tPos<-TMatrix[i]
  #     zeroMatrix[i,tPos]=1
  #   }
  #   return(zeroMatrix)
  # }
  
  # #  for testing purpose
  # supplyEstimates = FALSE
  # ms=ds1
  # vs=ds1
  # outcome="cvd_bi"
  # exp=c("dr_fiber")
  # sur=c("aofib86av")
  # covCalib=covariateV3V4
  # covOutcome=V1List
  # linear=FALSE
  # external=TRUE
  # family=binomial
  # link="logit"
  
  # vs=EVS
  # ms=MS
  # outcome="Y_bi"
  # exp=c("X1")
  # sur=c("Z1")
  # covCalib=c("V2","V4")
  # covOutcome=c("V4","R")
  # method="glm"
  # family=binomial
  # link="logit"
  # external=TRUE
  
  # ms=df1
  # vsIndicator="indicatorVS"
  # outcome="Y"
  # exp=c("X1","X2")
  # sur=c("Z1","Z2")
  # covCalib=c("V2","V4","V5")
  # covOutcome=c("V2","V4","R")
  # method="lm"
  # external=FALSE
  
  # ms=ds1
  # vs=ds2 
  # outcome="cvd_bi"
  # exp=c("dr_fiber")
  # sur=c("aofib90av")
  # covCalib=covariateV2V3V4
  # covOutcome=covariateV2V3V4
  # method="glm"
  # external=TRUE
  # family=binomial
  # link="logit"
  # 
  # load("C:/Users/wenze/Downloads/main.Rdata")
  # load("C:/Users/wenze/Downloads/valid.Rdata")
  # outcome<-"newcase"
  # expList<-c("drtotinc","drcalinc","dralcinc")
  # surList<-c("fqtotinc","fqcalinc","fqalcinc")
  # adjList<-c("agec")
  # ms=main
  # vs=valid 
  # sur=surList
  # exp=expList
  # covCalib=adjList
  # covOutcome=NULL
  # outcome=outcome 
  # method="glm"
  # family=binomial
  # link="logit"
  # external=TRUE
  # 
  # vs=EVS
  # ms=MS
  # outcome="Y"
  # exp=c("X1","X2")
  # sur=c("Z1","Z2")
  # covCalib=c("V2","V3","V4")
  # covOutcome=c("V3","V4")
  # method="lm"
  # external=TRUE
  # 
  # vs=EVS
  # ms=MS
  # outcome="Y"
  # exp=c("X")
  # sur=c("Z")
  # covCalib=c()
  # covOutcome=c()
  # method="lm"
  # external=TRUE
  ######################
  # Computation starts #
  ######################
  #####################
  # 1. Point estimate #
  #####################
  ## obtain complete-case data for validation and main dataset - external study
  commonCovariates <- intersect(covCalib,covOutcome)
  covCalibOnly<-setdiff(covCalib,commonCovariates)
  covOutcomeOnly<-setdiff(covOutcome,commonCovariates)
  covAll<-c(commonCovariates,covCalibOnly,covOutcomeOnly)
  
  if(external==TRUE){
    allVars_vs<-c(exp,sur,covCalib)
    vs_complete<-vs%>%dplyr::select(all_of(allVars_vs))%>%na.omit()
    
    allVars_ms<-c(sur,covAll,outcome)
    ms_complete<-ms%>%dplyr::select(all_of(allVars_ms))%>%na.omit() 
  }else if(external==FALSE){
    ## obtain complete-case data for validation and main dataset - internal study
    ### if internal validation then draw from main study
    allVars_ms<-c(sur,exp,covAll,outcome,vsIndicator)
    allVars_ms_miss<-c(sur,covAll,outcome,vsIndicator)
    ms_complete<-ms[complete.cases(ms[,allVars_ms_miss]),]
    ms_complete<-ms_complete%>%dplyr::select(all_of(allVars_ms))
    
    allVars_vs<-c(exp,sur,covCalib,outcome)
    vs_complete<-ms_complete%>%dplyr::select(all_of(allVars_vs))%>%na.omit()
  }
  
  ## create design matrix of calibration model
  if(length(covCalib)==0){
    exposureFormulaX<-paste0("~",paste0(sur,collapse="+"))
    exposureFormulaY<-paste0("~",paste0(exp,collapse="+"))
  }else{
    exposureFormulaX<-paste0("~",paste0(sur,collapse="+"),"+",paste0(covCalib,collapse="+"))
    exposureFormulaY<-paste0("~",paste0(exp,collapse="+"),"+",paste0(covCalib,collapse="+")) #Y here represents true exposure      
  }
  
  X_VS <- model.matrix(object= as.formula(exposureFormulaX) ,data=vs_complete)
  X_MS <- model.matrix(object= as.formula(exposureFormulaX) ,data=ms_complete) # for prediction purpose
  Y_VS <- model.matrix(object= as.formula(exposureFormulaY) ,data=vs_complete)[,-1] #remove intercept
  exposureModelVarNames<-colnames(X_VS)  
  
  #step 1: calibration model  
  ### exposure measurement error model, instead of using modeling, use matrix operations in a least square linear regression context.
  #### reminder: X, Z can be k x 1 dimension 
  
  X=as.matrix(X_VS)  #### design matrix, including intercept
  X.MS=as.matrix(X_MS)  #### design matrix, including intercept, will be used for prediction
  Y=as.matrix(Y_VS)  # no intercep
  
  #### Let B = B*%*%inv(GAMMAs) is the corrected parameter vector
  #### Variance(B) = B* SIGMA(inv(GAMMAs)) t(B*) + t(inv(GAMMAs)) SIGMA(B*) inv(GAMMAs), 
  ####        where SIGMA(inv(GAMMAs)) = t(d inv(GAMMAs)/d GAMMAs) Cov(GAMMAs) (d inv(GAMMAs)/d GAMMAs)
  ####        Cov(GAMMAs)=Var(G) is just VG in the following with dimension p^2 x p^2  
  #### We have so far B* = Bstar, SIGMA(B*)=VBstar, 
  #### in the following we derive first inv(GAMMAs) and then SIGMA(inv(GAMMAs)), where GAMMAs will exclude intercept 
  
  n<-nrow(X)   # number of obs in X
  pMeModel<-ncol(X)-1 #the original number of parameters in the calibration model except intercept. 
  p<-ncol(X)-1 #the original number of parameters in the calibration model except intercept. 
  
  F=solve(t(X)%*%X) #shorthand
  GWI=F%*%t(X)%*%Y  #estimates of GAMMAs, the parameter matrix in exp (vector) ~ sur (vector) + covCalib (vector), including intercept, dim: (p+1) x p (p equations embedded in X ~ Z + covCalib)
  if(length(covCalib)==0){
    GEV=t(GWI[sur,])
  }else{
    GEV=t(GWI[,exp])
  }
  ERR=(Y-X%*%GWI) # dimension n x p, doe snot depend on 
  
  # print(GWI)
  #### remove objects that are no longer used (free up memory)
  #remove(list=(c("X","Y")))
  
  #### calculate regular covariance matrix for the validation regression coefficients Gammas, which has dimension (p+1) x p 
  # S = (t(ERR)%*%ERR)/(n-pMeModel-1) # MSE, dimension p x p, does not depend on whether length(covOutcome)>0
  # if(length(covCalib)==0){
  #   colnames(S)=exp
  #   rownames(S)=exp
  #   VEV=S
  # }else{
  #   VEV = S[exp,exp] #variance-covariance matrix of the error variables Var(error) Z in exp ~ sur + covCalib) + error
  # }
  # G = t(GWI[2:(p+1),]) # dimension: p x p; This is subset of all GAMMAs <p+1 x p> (excluding the intercepts for each exp ~ sur + covCalib)
  # VG = matrixcalc::direct.prod(x=S,y=F[2:(p+1),2:(p+1)])# This is similar to the sigma^2 (X'X)^-1 in the case of unidimensional vector 
  # 
  # IGT = t(solve(t(G)))
  
  ## prediction
  if(ncol(X.MS%*%GWI)==1){
    X.hat_MS = as.matrix((X.MS%*%GWI)) # in main study
    X.hat_VS = as.matrix((X%*%GWI)) # in validation study      
  }else{
    X.hat_MS = as.matrix((X.MS%*%GWI)[,exp]) # in main study
    X.hat_VS = as.matrix((X%*%GWI)[,exp]) # in validation study    
  }

  colnames(X.hat_MS) <- paste0(exp,".hat") # as.matrix is used for cases where X.hat_MS is a vector
  colnames(X.hat_VS) <- paste0(exp,".hat")
  
  MS_new<-ms_complete
  nrow_ms_complete<-nrow(ms_complete)
  # now create (predicted) true variables
  ## if external validation study, then directly create new X using predicted X
  for(i in 1:length(exp)){
    exp_i <- exp[i]
    exphat_i <- paste0(exp,".hat")[i]
    if(external==TRUE){
      MS_new[[exp_i]]<-X.hat_MS[,exphat_i]
    }else if(external==FALSE){
      MS_new[MS_new[vsIndicator]==0,][[exp_i]]<-X.hat_MS[MS_new[,vsIndicator]==0,exphat_i]
      # only overwrite X that are missing in main study
    }
  }

  # remove(list=(c("X","Y")))

  #step 2: outcome model
  ## create formula for outcome model
  if (length(covOutcome)==0){
    outcomeFormula<-paste0("~",paste0(exp,collapse="+"))
  }else{
    outcomeFormula<-paste0("~",paste0(exp,collapse="+"),"+",paste0(covOutcome,collapse="+"))
  }
  
  ## CRS outcome model modeling (additive)
  if(method=="lm"){
    outModel<-lm(formula=as.formula(paste0(outcome,outcomeFormula)),data=MS_new)
  }else{
    outModel<-glm(formula=as.formula(paste0(outcome,outcomeFormula)),
                  data=MS_new,family=family)
  }
  
  # prediction
  if(method=="lm"){
    Y.hat_MS<-predict(outModel)
  }else if(method=="glm"){
    Y.hat_MS<-predict.glm(outModel,type="response")
  }
  
  Y_MS<-MS_new[,outcome] # outcome in main study
  # if internal validation study additionally calculate
    if(external==FALSE){
      Y.hat_VS<-Y.hat_MS[MS_new[,vsIndicator]==1]
      Y_VS<-MS_new[MS_new[,vsIndicator]==1,outcome]
    }
  
  # point estimate
  pointEstimate=coef(outModel)
  
  if(show_detail == TRUE){
    print("point estimate =:")
    print(pointEstimate)
  }
  #################################
  # 2. 95% CI and Standard Errors #
  #################################
  # validation sample size
  n=nrow(vs_complete)
  # main sample size
  m=nrow(MS_new)

  
  ms_complete_VS<-ms_complete
  
  ## prep work: create design matrix
  if(length(covCalib)==0){
    exposureFormulaX<-paste0("~",paste0(sur,collapse="+"))
  }else{
    exposureFormulaX<-paste0("~",paste0(sur,collapse="+"),"+",paste0(covCalib,collapse="+"))
  }
     
  if(length(covOutcome)==0){
    exposureFormulaY<-paste0("~",paste0(exp,collapse="+")) # the Y here now is the right hand side of Y ~ X + V
  }else{
    exposureFormulaY<-paste0("~",paste0(exp,collapse="+"),"+",paste0(covOutcome,collapse="+"))   # the Y here now is the right hand side of Y ~ X + V    
  }     
  
  Z_EVS <- model.matrix(object= as.formula(exposureFormulaX) ,data=vs_complete)%>%as.data.frame()
  Z_MS  <- model.matrix(object= as.formula(exposureFormulaX) ,data=ms_complete_VS)%>%as.data.frame()
  X_MS  <- model.matrix(object= as.formula(exposureFormulaY) ,data=MS_new)%>%as.data.frame()
  
  if(external==FALSE){
    EVS_new<-MS_new[MS_new[,vsIndicator]==1,]
    X_EVS <- model.matrix(object= as.formula(exposureFormulaY) ,data=EVS_new)%>%as.data.frame()
  }
  
  ## Overwrite all covariate list to include possible expansion in name of categorical variable(s)
  covOutcome     <- setdiff(colnames(X_MS)[-1],exp)
  covCalib       <- setdiff(colnames(Z_EVS)[-1],sur)
  covBoth        <- intersect(covOutcome,covCalib)
  covOutcomeOnly <- setdiff(covOutcome,covBoth)
  covCalibOnly   <- setdiff(covCalib,covBoth)
  covAll         <- union(covCalib,covOutcome)
  
  #step 1: force zeros for 
  ## (1) covariates in covCalib but not in MS and
  ## (2) covariates in covOutcome but not in EVS
  
  ## (1) covariates in covCalib but not in MS and
  if(length(covCalibOnly)>0){
    for(i in 1:length(covCalibOnly)){
      name <- covCalibOnly[i]
      assign(name,rep(0,m))
      X_MS[,name]<-get(name)
      if(external==FALSE){
        assign(name,rep(0,n))
        X_EVS[,name]<-get(name)       
      }
    }
  }
  
  ## (2) covariates in covOutcome but not in EVS
  if(length(covOutcomeOnly)>0){
    for(i in 1:length(covOutcomeOnly)){
      name <- covOutcomeOnly[i]
      assign(name,rep(0,n))
      Z_EVS[,name]<-get(name)
      assign(name,rep(0,m))
      Z_MS[,name]<-get(name)
    }
  }

  # reorder all variable names as exp/sur + covAll
  varOrderZ <-c("(Intercept)",sur,covAll) 
  varOrderX <-c("(Intercept)",exp,covAll) 
  Z_EVS   <- Z_EVS[,varOrderZ]%>%as.matrix()
  Z_MS    <- Z_MS[,varOrderZ]%>%as.matrix()
  X_MS    <- X_MS[,varOrderX]%>%as.matrix()
  if(external==FALSE){
    X_EVS   <- X_EVS[,varOrderX]%>%as.matrix()
  }
  
  #step 2. embedded function, obtain Pm Z_i X_i where Z and X has the same dimension and Z_i and X_i are the row vector of the matrix, Yhat is the predicted outcome for GLM models
  meanOfRowVectorInMatrix1<-function(Vhat,Z,X){
    # Z=Z_EVS
    # X=Z_EVS
    # Vhat=rep(1,nrow(Z_EVS))
    matrixRowValue<-matrix(NA,nrow=nrow(Z),ncol=ncol(Z)*ncol(X))
    # obtain the matrix where each row contains the value of the product of Z_i and X_i
    for(i in 1:nrow(Z)){
      Vhat_i <- as.matrix(Vhat[i])
      Z_i <- as.matrix(Z[i,])
      X_i <- as.matrix(X[i,])
      ZX_i <- (Z_i)%*%Vhat_i%*%t(X_i)
      matrixRowValue[i,]<-ZX_i
    }
    # obtain mean for each column in the matrix
    meanMatrixRowValue<-apply(X=matrixRowValue,MARGIN=2,FUN=mean)
    outputMatrix<-matrix(meanMatrixRowValue,nrow=ncol(Z),ncol=ncol(Z))
    colnames(outputMatrix)<-c("Intercept",exp,covAll)
    rownames(outputMatrix)<-c("Intercept",exp,covAll)
    return(outputMatrix)
  }
  # 
  meanOfRowVectorInMatrix2<-function(DhatZ,DhatX,Z,X){
    # DhatZ<-diff_X_ZV[,"X1"]
    # DhatX<-diff_X_ZV[,"X2"]
    # Z<-Z_EVS
    # X<-Z_EVS
    matrixRowValue<-matrix(NA,nrow=nrow(Z),ncol=ncol(Z)*ncol(X))
    # obtain the matrix where each row contains the value of the product of Z_i and X_i
    for(i in 1:nrow(Z)){
      DhatZ_i <- DhatZ[i]
      DhatX_i <- DhatX[i]
      Z_i <- as.matrix(Z[i,])
      X_i <- as.matrix(X[i,])

      ZX_i <- (Z_i)%*%DhatZ_i%*%t((X_i)%*%DhatX_i)
      matrixRowValue[i,]<-ZX_i
    }
    # obtain mean for each column in the matrix
    meanMatrixRowValue<-apply(X=matrixRowValue,MARGIN=2,FUN=mean)
    outputMatrix<-matrix(meanMatrixRowValue,nrow=ncol(Z),ncol=ncol(Z))
    return(outputMatrix)
  }
  
  #step 3: prepare items related to the situation where there is >1 mismeasured exposure
  mismeasureNumber <- length(exp)
  
  #step 4: PREPARE FOR GLM's COEFFICIENTs
  # 1. coefficient in matrix A that varies by family GLM and link function
  if(method=="lm"){
    V_hat = rep(1,length(Y.hat_MS))
  }else if(link=="logit"){
    V_hat = Y.hat_MS*(1-Y.hat_MS)
  }else if(link=="log"){
    V_hat = Y.hat_MS
  }
  

  # Now obtain the block matrices in the Carroll's method
  ## preparation
  ### obtain variance matrix: (X- Z alpha) (X - Z alpha)^T
  if(ncol(Y)==1){
    var_X_ZV  <- (Y - X.hat_VS)^2 # Y here represents the X in the X ~ Z + V model
    diff_X_ZV <- as.matrix(Y - X.hat_VS) # Y here represents the X in the X ~ Z + V model, the as.matrix is to make sure it is in matrix form even for a vector in the multivariate case.
    
  }else{
    var_X_ZV  <- (Y[,exp] - X.hat_VS)^2 # Y here represents the X in the X ~ Z + V model
    diff_X_ZV <- as.matrix(Y[,exp] - X.hat_VS) # Y here represents the X in the X ~ Z + V model, the as.matrix is to make sure it is in matrix form even for a vector in the multivariate case.
  }
  
  if(show_detail == TRUE){
    print("(X-X.hat)^2/n =:")
    print(mean(var_X_ZV))
    print(sqrt(mean(var_X_ZV)))
  }
  colnames(diff_X_ZV)<-exp
  ### obtain variance matrix: (Y- X beta) (Y- X beta)^T
  var_Y_XV  <- (MS_new[,outcome] - Y.hat_MS)^2 # Y here represents the Y in the Y ~ Xhat + V model
  # sqrt(mean(var_Y_XV))
  diff_Y_XV <- (MS_new[,outcome] - Y.hat_MS) # Y here represents the Y in the Y ~ Xhat + V model
  ### for internal validity study
  if(external==FALSE){
    var_Y_XV_VS  <- (MS_new[MS_new[,vsIndicator]==1,outcome] - Y.hat_VS)^2 # Y here represents the Y in the Y ~ Xhat + V model
    diff_Y_XV_VS <- (MS_new[MS_new[,vsIndicator]==1,outcome] - Y.hat_VS) # Y here represents the Y in the Y ~ Xhat + V model        
  }
  if(show_detail == TRUE){
    print("(Y-Y.hat)^2/m =:")
    print(mean(var_Y_XV))
    print(sqrt(mean(var_Y_XV)))
  }

  ## A MATRICES
  ### A_11
  # A_11 <- meanOfRowVectorInMatrix1(rep(1,nrow(Z_EVS)), Z_EVS,Z_EVS)
  
  A_11 <- t(Z_EVS)%*%(Z_EVS)/n
  A_11_i <- A_11
  nrow_A_11_i<- nrow(A_11_i)
  ncol_A_11_i<- ncol(A_11_i)
  if(mismeasureNumber>1){
    for(i in 2:(mismeasureNumber)){
      nrow_A_11<-nrow(A_11)
      ncol_A_11<-ncol(A_11)
      # recursively create the diagnal matrices made up of t(Z)Z
      A_11_j<-as.matrix(rbind(cbind(A_11,matrix(0,nrow=nrow_A_11,ncol=ncol_A_11_i)),cbind(matrix(0,nrow=nrow_A_11_i,ncol=ncol_A_11),A_11_i)))
      A_11<- A_11_j    
    }    
    rownames(A_11)<- rep(c("Intercept",exp,covAll),mismeasureNumber)
    colnames(A_11)<- rep(c("Intercept",exp,covAll),mismeasureNumber)
  }
  # remove outcome model only rows
  A_11 <- A_11[!row.names(A_11)%in%covOutcomeOnly,!colnames(A_11)%in%covOutcomeOnly]
  
  if(show_detail == TRUE){
    print("A_11 =:")
    print(A_11)
  }

  ### A_21
  # A_21 <- pointEstimate[exp][1]*meanOfRowVectorInMatrix1(V_hat, X_MS, Z_MS)
  
  # the following correction made by Jingyu Cui at Feb 20, 2026
  # running error: dimension does not comparable in calculating A_21_X_MS
  # minus_matri_4_X_MS = cbind(rep(0,length(diff_Y_XV)),
  #                            diff_Y_XV,
  #                            matrix(0, nrow=nrow(X_MS), ncol=(ncol(X_MS)- mismeasureNumber - 1)))
  
  minus_matri_4_X_MS = cbind(rep(0,length(diff_Y_XV)),
                            diff_Y_XV,
                             matrix(0, nrow=nrow(X_MS), ncol=(ncol(X_MS)- 2)))
  A_21_X_MS = pointEstimate[exp][1]*X_MS -minus_matri_4_X_MS
  A_21 <- (V_hat*t(A_21_X_MS))%*%Z_MS/m
  if(mismeasureNumber>1){
    # A_21_unit<-meanOfRowVectorInMatrix1(V_hat, X_MS, Z_MS)
    A_21_unit <- (V_hat*t(A_21_X_MS))%*%Z_MS/m
    A_21_i <- A_21
    for(i in 2:(mismeasureNumber)){
      pointEstaimte_i<-pointEstimate[exp][i]
      A_21_iplus1<- pointEstaimte_i*A_21_unit
      # recursively create the diagnal matrices made up of t(Z)Z
      A_21_i<-as.matrix(cbind(A_21_i,A_21_iplus1))
    }    
    A_21<- A_21_i    
  } 
  rownames(A_21) <- colnames(X_MS)
  # remove covCalib/covOutcome only that respectively have row/col value 0 (by design)
  A_21 <- A_21[!row.names(A_21)%in%covCalibOnly,!colnames(A_21)%in%covOutcomeOnly]
  if(show_detail == TRUE){
    print("A_21 =:")
    print(A_21)
  }
  
  A_12 = matrix(0, nrow = ncol(A_21), ncol = nrow(A_21))
  ### A_22
  # A_22 <- meanOfRowVectorInMatrix1(V_hat, X_MS, X_MS)
  A_22 <- (V_hat*t(X_MS))%*%X_MS/m
  A_22 <- A_22[!row.names(A_22)%in%covCalibOnly,!colnames(A_22)%in%covCalibOnly]
  if(show_detail == TRUE){
    print("A_22 =:")
    print(A_22)
  }
  
  A = cbind(rbind(A_11,A_21), rbind(A_12,A_22))
  invA <- solve(A)
  ## B MATRICES
  
  ### reminder for the following quantities
  # Z_EVS, Z_MS
  # X_EVS, X_MS
  # Y.hat_MS, Y.hat_VS <only available for IVS>
  # X.hat_MS, X.hat_VS
  
  
  B_11 <- matrix(NA,nrow=length(exp)*ncol(Z_EVS),ncol=length(exp)*ncol(Z_EVS))
  for(i in 1:(mismeasureNumber)){
    for(j in 1:(mismeasureNumber)){
      exp_i <- exp[i]
      exp_j <- exp[j]
      # B_11_ij<- meanOfRowVectorInMatrix2(diff_X_ZV[,exp_i],diff_X_ZV[,exp_j],Z_EVS,Z_EVS)/n
      #
      e_i_matrix = diff_X_ZV[,exp_i]
      e_j_matrix = diff_X_ZV[,exp_j]
      sigma2X_ZV_ij = sum(e_i_matrix*e_j_matrix)/(n - length(covCalib) - mismeasureNumber - 1) # dof for MEM 
      B_11_ij<- sigma2X_ZV_ij*t(Z_EVS)%*%Z_EVS/(n^2)
      B_11[((i-1)*(ncol(Z_EVS))+1):(i*(ncol(Z_EVS))),((j-1)*(ncol(Z_EVS))+1):(j*(ncol(Z_EVS)))]<-B_11_ij
    }
  }
  rownames(B_11)<- rep(c("Intercept",exp,covAll),mismeasureNumber)
  colnames(B_11)<- rep(c("Intercept",exp,covAll),mismeasureNumber)
  if(show_detail == TRUE){
    print("B_11 =:")
    print(B_11)
  }
  # remove covOutcome/covCalib only that respectively have row/col value 0 (by design)
  B_11 <- B_11[!row.names(B_11)%in%covOutcomeOnly,!colnames(B_11)%in%covOutcomeOnly]
  
  # B_22 <-  meanOfRowVectorInMatrix2(diff_Y_XV,diff_Y_XV, X_MS,X_MS)/m
  sigma2Y_XhatV = sum(var_Y_XV)/(m - length(covOutcome) - mismeasureNumber - 1) # dof for outcome Y model
  
  B_22 <-  sigma2Y_XhatV*(t(X_MS)%*%X_MS)/m^2
  B_22 <-  B_22[!row.names(B_22)%in%covCalibOnly,!colnames(B_22)%in%covCalibOnly]
  if(show_detail == TRUE){
    print("B_22 =:")
    print(B_22)
  }
  B_12 <- matrix(0,nrow=length(exp)*ncol(Z_EVS),ncol=ncol(Z_EVS))
  if(external==FALSE){
    for(j in 1:(mismeasureNumber)){
      exp_j <- exp[j]
      # B_12_j<- meanOfRowVectorInMatrix2(diff_X_ZV[,exp_j],diff_Y_XV_VS,Z_EVS,X_EVS)/m
      B_12_j<- (diff_X_ZV[,exp_j]*t(Z_EVS))%*%(diff_Y_XV_VS*X_EVS)/(m*n)
      B_12[((j-1)*(ncol(Z_EVS))+1):(j*(ncol(Z_EVS))),]<-B_12_j
    }
  }
  rownames(B_12)<- rep(c("Intercept",exp,covAll),mismeasureNumber)
  colnames(B_12)<- c("Intercept",exp,covAll)
  # remove covOutcome/covCalib only that respectively have row/col value 0 (by design)
  B_12 <- B_12[!row.names(B_12)%in%covOutcomeOnly,!colnames(B_12)%in%covCalibOnly]
  
  B = rbind(cbind(B_11,B_12), cbind(t(B_12),B_22))
  temp = invA %*% B %*% t(invA)
  Z_EVS_varNames = colnames(Z_EVS)
  covarianceEstimates <- temp[!row.names(temp)%in%Z_EVS_varNames,!colnames(temp)%in%Z_EVS_varNames]
  if(show_detail == TRUE){
    print("VAR_B =:")
    print(covarianceEstimates)
  }
  outcomeModelVars_all <- colnames(X_MS)[-1] # this includes all covariates (expanded)  (model matrix version)
  outcomeModel_only<-names(coef(outModel))[-1] # this includes only outcome variable specified in the outcome model (model matrix version)
  
  # the block below is revised by Jingyu Cui on Feb 20, 2026 since running error
  # # point estimate only contains outcome model variables (model matrix version)
  # Bhat  <- pointEstimate[outcomeModel_only]
  # if(!is.matrix(covarianceEstimates)){
  #   VB <- covarianceEstimates  
  #   SE_B<-sqrt(VB)
  # }else{
  #   VB <- covarianceEstimates[outcomeModel_only,outcomeModel_only]
  #   SE_B<-sqrt(diag(VB))
  # }
  
  #### Jingyu's correction ####
  
  kZ <- nrow(A_11)
  kT <- nrow(temp)
  
  covarianceEstimates <- temp[(kZ+1):kT, (kZ+1):kT, drop = FALSE]
  rownames(covarianceEstimates) <- colnames(covarianceEstimates) <- rownames(A_22)
  
  # ---- revised: point estimate + VCOV + SE for outcome-model coefficients ----
  
  # coefficients actually in the fitted outcome model (model-matrix / expanded factor names)
  outcomeModel_only <- names(coef(outModel))[-1]   # same as your existing line
  Bhat <- pointEstimate[outcomeModel_only]
  
  # covarianceEstimates should already be the OUTCOME block VCOV
  # (i.e., bottom-right block of temp) with dimnames matching the outcome block.
  # Ensure it's a matrix with dimnames:
  covarianceEstimates <- as.matrix(covarianceEstimates)
  
  # Subset to coefficients in outModel (drop=FALSE keeps matrix form even if 1 coef)
  missing <- setdiff(outcomeModel_only, rownames(covarianceEstimates))
  if (length(missing) > 0) {
    stop("Outcome coefficients missing from covarianceEstimates: ",
         paste(missing, collapse = ", "))
  }
  
  VB <- covarianceEstimates[outcomeModel_only, outcomeModel_only, drop = FALSE]
  SE_B <- sqrt(diag(VB))
  
  #################### end of Jingyu Cui revision ####
  
  # compose output tables

  zValue<-Bhat/SE_B
  pValue<-2*pnorm(-abs(as.numeric(zValue)))
  lcl<-Bhat - qnorm(0.975)*SE_B
  ucl<-Bhat + qnorm(0.975)*SE_B
  BSebP<-cbind(Bhat,SE_B,zValue,pValue,lcl,ucl)
  
  colnames(BSebP)<-c("Estimate","Std. Error","Z Value","Pr(>|Z|)","lower 95%CI","upper 95%CI")
  # print(BSebP)
  # output objects
  outputList<-list(BSebP,VB)
  names(outputList)<-c("correctedCoefTable","correctedVCOV")
  
  return(outputList)
}




