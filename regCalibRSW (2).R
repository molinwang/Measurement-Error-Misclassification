########################
# regCalib function #
########################

# authors: Wenze Tang and Molin Wang
# version 1: jul 29, 2021

# Abstract:
# This function corrects for measurement error in exposure or exposure and covariates and gives corrected coefficients associated standard errors, p values
# as well as variance-covariance matrix of corrected coefficients. Users can choose to either use built-in outcome model where only linear and generalized linear
# model is currently supported under both internal and external validation study design, or supply their own coefficients and variance-covariance matrix from 
# outcome models such as logistic model or Cox models under external validation study design. 
# A validation study is required to empirically characterize the measurement error calibration model. Options are given for main study/external validation study 
# design and main study/internal validation study design (not an option when users supply their own uncorrected estimates) (Spiegelman, Carrol, Kipnis; 2001).
#  More technical details are given in Rosner et al (1989), Rosner et al (1990) and Spiegelman et al (1997). 

# Reference:
# Rosner B, Willett WC, Spiegelman D "Correction of logistic relative risk estimates and confidence
# intervals for systematic withinperson measurement error". Statistics in Medicine 8: 1051-1069, 1989.
# Rosner B, Spiegelman D, Willett WC "Correction of logistic regression relative risk estimates and
# confidence intervals for measurement error: the case of multiple covariates measured with error".
# American Journal of Epidemiology 1990;132: 734-735.
# Spiegelman D, McDermott A, Rosner B "The many uses of the 'regression calibration' method
# for measurement error bias correction in nutritional epidemiology". American Journal of Clinical
# Nutrition, 1997; 65:1179S-1186S.
# Spiegelman D, Carroll RJ, Kipnis V "Efficient regression calibration for logistic regression in main
# study/internal validation study designs with an imperfect regerence instrument". Statistics in
# Medicine 2001;20:139-160.


#' Find corrected parameters (both single or multiple) and their standard deviations and p values for 
#' mismeasured continuous exposure
#' @param supplyEstimates=FALSE Indicates whether uncorrected estimates will be supplied by the user. 
#' If supplied by the user, then `ms` is optional. Standard regression results will not be returned if TRUE.     
#' @param ms The input main study data set as dataframe. This dataframe should minimally include variables specified in 
#' `sur`, `covCalib` and `covOutcome` (if any).
#' @param vs The input internal/external validation data set as dataframe. This dataframe should minimally include 
#' variables indicated in `exp`, `sur` and `covCalib`.
#' @param sur character vector of mismeasured exposure and covariates (i.e. surrogates) in the main study dataset.
#' @param exp character vector of correctly measured exposure and covariates that has a one-to-one correspondence to those 
#' specified in `sur` in validation dataset. Must have same length as `sur`.
#' @param covCalib=NULL character vector of names of correctly measured covariates to adjust for in calibration model and outcome model.
#' @param covOutcome=NULL character vector of names of correctly measured risk factors (that are not associated with exposure or surrogate) 
#' for outcome in the main study data to adjust for in outcome model in additional to those specified in covCalib. 
#' Should not include any variable from `covCalib`. Leave unspecified if no such risk factors. 
#' @param outcome Outcome variable, required if method is lm or glm.
#' @param event Required if method is cox. The event status indicator, normally 0=alive/censored, 1=dead/event. 
#'  For interval censored data, the status indicator is 0=right censored, 1= event at time, 2=left censored, 3=interval censored.
#' @param time Required if method is cox. For right censored data, this is the follow up time. For interval data, 
#' the first argument is the starting time for the interval.
#' @param method="lm" Methods for modeling, currently only `lm`, `glm` and `cox` methods are available. Required.    
#' @param family Supply family parameter to pass to glm function. Not a character. Required if method="glm".
#' @param link Supply link parameter to pass to glm function. Should be character. Required if method="glm". 
#' @param external=TRUE Indicates whether `vs` is an external validation set. If external=FALSE, then vs must contain
#'  variable specified in `outcome`. If external=FALSE, then user cannot supply uncorrected estimates, i.e. supplyEstimates=FALSE.
#' @param pointEstimates A numeric vector of point estimates from uncorrected standard regression analysis. 
#' Intercept estimate must be removed. Must include names for each point estimates corresponding to the (transformed) 
#' names from `covCalib` followed by `covOutcome`. Must be supplied if supplyEstimates = TRUE. 
#' @param vcovEstimates A p by p Variance-covariance matrix estimates from uncorrected standard regression analysis. 
#' Intercept estimates must be removed. Must include column names (excluding intercept) for the estimates corresponding 
#' to the (transformed) names from `covCalib` followed by `covOutcome`. Must be supplied if supplyEstimates = TRUE.  
#' @return printable dataframe from standard regression results (when supplyEstimates==FALSE) as well as corrected results

regCalibRSW<-function(supplyEstimates=FALSE, ms,vs, 
                   sur, exp, covCalib=NULL, covOutcome=NULL, outcome=NA, event, time, method="lm",family=NA,link=NA, external=TRUE
                   ,pointEstimates=NA, vcovEstimates=NA){
require("stats")
require("earth")
require("matrixcalc")
require("Matrix")
require("dplyr")
require("glm")
require("lm")
require("survival")
###################  
# check arguments #
###################
  ## data related warnings
  if(missing(vs)){
    stop("Input data vs not supplied.")
  }else if(class(vs)!="data.frame"){
    stop("Input data vs must be of data.frame class.")
  }
  
  if(supplyEstimates==FALSE){
    if(missing(ms)){
      stop("Input data ms not supplied.")
    }else if(class(ms)!="data.frame"){
      stop("Input data ms must be of data.frame class.")
    }
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
  }else if(length(base::intersect(covCalib,covOutcome))>0){
    stop("There should be no overlapping variables in `covCalib` and `covOutcome`.")
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
  
  ## 1. check if MS contains data indicated by sur, covCalib and covOutcome
  if(supplyEstimates==FALSE){
    # if(length(covOutcome)>0&length(covCalib)>0){
    #   MSVars_spec<-(c(sur,covCalib,covOutcome))
    # }else if(length(covCalib)==0&length(covOutcome)==0){
    #   MSVars_spec<-(c(sur))
    # }
    MSVars_spec<-(c(sur,covCalib,covOutcome))
    MSVars<-colnames(ms)
    inMSVars<-(MSVars_spec%in%MSVars)
    if(sum(inMSVars)!=length(MSVars_spec)){
      stop("Main study dataset does not contain all the necessary variables specified in one of the following parameter: id, sur, covCalib and covOutcome.")
    } 
  }

  ## 2. check if EVS contains data indicate by sur, exp and covCalib; for IVS, check additionally for outcome
  EVSVars<-colnames(vs)
  EVSVars_spec<-(c(sur,exp,covCalib))
  inMSVars<-(EVSVars_spec%in%EVSVars)
  if(sum(inMSVars)!=length(EVSVars_spec)){
    stop("Validation study dataset does not contain all the necessary variables specified in one of the following parameter: id, sur, exp and covCalib.")
  }
  if(external==FALSE){
    if(!(outcome%in%colnames(vs))){
      stop("Outcome variable is not available in the supplied internal validation data.")
    }
  }
  
  if(supplyEstimates==FALSE){
    ## check if variables in covOutcome are strongly associated with exposure or surrogates, by regress exposure on covOutcome
    if(length(covOutcome)>0&class(covOutcome)=="character"){
      trackerPvalueGE005=0
      for(i in 1:length(sur)){
        surUnivariate<-sur[i]
        checkCovOutcomeFormula<-paste0(surUnivariate,"~",paste0(covOutcome,collapse="+"))
        checkCovOutcomeModel<-lm(data=ms,formula=as.formula(checkCovOutcomeFormula))
        checkCovOutcomePValues<-summary(checkCovOutcomeModel)$coefficients[-1,4]
        if(sum(checkCovOutcomePValues<0.05)>0){
          trackerPvalueGE005=trackerPvalueGE005+1
        }
        # print(list(checkCovOutcomeFormula,summary(checkCovOutcomeModel),checkCovOutcomePValues,trackerPvalueGE005)) # testing purpose
      }
      
      if(trackerPvalueGE005>0){
        warning("At least one of the risk factors specified in covOutcome is strongly associated with exposure(s) or surrogate(s) and should be reconsidered to be specified in covCalib instead.")
      }  
    }   
  }

  if(supplyEstimates==TRUE){
  #   ## check if length of point estiamtes and vcov estimates supplied are the same as specified in covCalib and covOutcome combined. 
  #   if(length(covOutcome)==0){
  #     l_est<-length(covCalib) # needs to account for NA which has length 1
  #   }else{
  #     l_est<-length(covCalib) + length(covOutcome)
  #   }
  # 
  #   if(l_est!=length(pointEstimates)){
  #     stop("Length of point estimates from uncorrected regression analysis is not equal to the length of covCalib and covOutcome combined!")
  #   }else if(l_est!=dim(vcovEstimates)[1]|l_est!=dim(vcovEstimates)[2]){
  #     stop("Dimension of variance-covariance estimates from uncorrected regression analysis is not equal to the length of covCalib and covOutcome combined!")
  #   }
    
    ## check whether there are names for the point Estimates and column names for the vcov estimates
    if(length(names(pointEstimates))==0|length(colnames(vcovEstimates))==0){
      stop("There must be names for user supplied point estiamtes and column names for variance covariance matrix!")
    }
    
    ## check if vcov estiamtes are a square matrix
    if(dim(vcovEstimates)[1]!=dim(vcovEstimates)[2]){
      stop("Covariance matrix supplied is not a symmetric square matrix. Please check!")
    }
  }

  ######################
  # Embedded functions #
  ######################
  #### create a function similar to the design function in SAS code
  design<-function(TMatrix,PVector){
    zeroMatrix<-matrix(rep(0,(PVector^2)^2),nrow=PVector^2)
    # replace 0 in the zeroMatrix's positions indicated in TMatrix as 1
    for(i in 1:length(TMatrix)){
      tPos<-TMatrix[i]
      zeroMatrix[i,tPos]=1
    }
    return(zeroMatrix)
  }
  
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

  # supplyEstimates = TRUE
  # vs=EVS
  # outcome="Y"
  # exp=c("X1","X2")
  # sur=c("Z1","Z2")
  # covCalib=c("V1","V2","V3","V4")
  # covOutcome=NA
  # pointEstimates = pointEst
  # vcovEstimates = vcovEst
  # external=TRUE


  # supplyEstimates = FALSE
  # ms=MS
  # vs=EVS
  # outcome="Y"
  # exp=c("X1","X2")
  # sur=c("Z1","Z2")
  # covCalib=c("V1","V2","V3","V4")
  # covOutcome=c("R")
  # external=TRUE
  # method="lm"

  
  # supplyEstimates = FALSE
  # vs=valid
  # ms=main
  # outcome="case"
  # exp=c("X_cont")
  # sur=c("Z_cont")
  # external=TRUE
  # covCalib=NULL
  # pointEstimates = NULL
  # vcovEstimates = NULL
  # covOutcome=NULL
  
  # supplyEstimates = FALSE
  # ms=ds1
  # vs=ds1
  # outcome=c("cvd_bi")
  # exp=c("dr_fiber")
  # sur=c("fib")
  # covCalib=covariateV3V4
  # covOutcome=NULL
  # method="glm"
  # external=TRUE
  # family=binomial
  # link="logit"
  
  # supplyEstimates = FALSE
  # ms=MS
  # vs=EVS
  # event="Y_bi"
  # time ="Y_time"
  # outcome=NULL
  # exp=c("X1","X2")
  # sur=c("Z1","Z2")
  # covCalib=covCalib=c("V1","V2","V3","V4")
  # covOutcome=NULL
  # method="cox"
  # external=TRUE

  ######################
  # Computation starts #
  ######################
  #step 0: create design matrix
  ## create formula for outcome model
  if(supplyEstimates==FALSE){
    if (length(covOutcome)==0&length(covCalib)==0){
      outcomeFormula<-paste0("~",paste0(sur,collapse="+"))
      allVars_ms<-c(outcome,sur)      
    }else if(length(covOutcome)==0){
      outcomeFormula<-paste0("~",paste0(sur,collapse="+"),"+",paste0(covCalib,collapse="+"))
      allVars_ms<-c(outcome,sur,covCalib)
    }else{
      outcomeFormula<-paste0("~",paste0(sur,collapse="+"),"+",paste0(covCalib,collapse="+"),"+",paste0(covOutcome,collapse="+"))
      allVars_ms<-c(outcome,sur,covCalib, covOutcome)
    }
    
    # rewrite if method is cox
    if(method=="cox"){
      if (length(covOutcome)==0&length(covCalib)==0){
        outcomeFormula<-paste0("~",paste0(sur,collapse="+"))
        allVars_ms<-c(time,event,sur)      
      }else if(length(covOutcome)==0){
        outcomeFormula<-paste0("~",paste0(sur,collapse="+"),"+",paste0(covCalib,collapse="+"))
        allVars_ms<-c(time,event,sur,covCalib)
      }else{
        outcomeFormula<-paste0("~",paste0(sur,collapse="+"),"+",paste0(covCalib,collapse="+"),"+",paste0(covOutcome,collapse="+"))
        allVars_ms<-c(time,event,sur,covCalib, covOutcome)
      }
    }
  }


  ## identify complete cases
  if(external==TRUE){
    allVars_vs<-c(exp,sur,covCalib)
  }else if(external==FALSE){
    if (length(covOutcome)==0&length(covCalib)==0){
      allVars_vs<-c(exp,sur,outcome)
    }else if(length(covOutcome)==0){
      allVars_vs<-c(exp,sur,covCalib,outcome)
    }else{
      allVars_vs<-c(exp,sur,covCalib,covOutcome,outcome)
    }
  }
  
  if(length(covOutcome)==0&length(covCalib)==0){
    exposureFormulaX<-paste0("~",paste0(sur,collapse="+"))
    exposureFormulaY<-paste0("~",paste0(exp,collapse="+"))  
  }else{
    exposureFormulaX<-paste0("~",paste0(sur,collapse="+"),"+",paste0(covCalib,collapse="+"))
    exposureFormulaY<-paste0("~",paste0(exp,collapse="+"),"+",paste0(covCalib,collapse="+"))      
    
  }
  
  if(supplyEstimates==FALSE){
    ms_complete<-ms%>%dplyr::select(all_of(allVars_ms))%>%na.omit()
    X_MS <- model.matrix(object= as.formula(outcomeFormula),data=ms_complete)
    Y_MS<- ms_complete[,outcome]
    outcomeModelVarNames<-colnames(X_MS)
  }else if(supplyEstimates==TRUE){
    outcomeModelVarNames<-c("(Intercept)",names(pointEstimates))
  }

  vs_complete<-vs%>%dplyr::select(all_of(allVars_vs))%>%na.omit()
  
  X_VS <- model.matrix(object= as.formula(exposureFormulaX) ,data=vs_complete)
  Y_VS<- model.matrix(object= as.formula(exposureFormulaY) ,data=vs_complete)[,-1] #remove intercept
  exposureModelVarNames<-colnames(X_VS)
  
  ## calculate the number of design variables derived from covOutcome
  ## i.e. the variables that are in outcomeModelVarNames but not in exposureModelVarName

  riskFactorModelVarNames<-setdiff(outcomeModelVarNames,exposureModelVarNames) # in case intercept is not included as in outcomeModelVarNames<-names(pointEstimates)

  #step 1: outcome model
  ## RSW outcome model modeling (additive)
  if(supplyEstimates==FALSE){
    
    if(method=="lm"){
      outModel<-lm(formula=as.formula(paste0(outcome,outcomeFormula)),data=ms_complete)
      outcomeParam=coef(outModel)
      outcomeParamVCOV=vcov(outModel)
      outcomeModelResults<-(list(outcomeParam,outcomeParamVCOV))
    #}else if(survival==TRUE){
    }else if(method=="glm"){
      outModel<-glm(formula=as.formula(paste0(outcome,outcomeFormula)),
                    data=ms_complete,family=family(link))
      outcomeParam=coef(outModel)
      outcomeParamVCOV=vcov(outModel)
      outcomeModelResults<-(list(outcomeParam,outcomeParamVCOV))
    }else if(method=="cox"){
      outModel<-coxph(formula=as.formula(paste0("Surv(",time,",", event,")",outcomeFormula)),
                    data=ms_complete)
      outcomeParam=coef(outModel)
      outcomeParamVCOV=vcov(outModel)
      outcomeModelResults<-(list(outcomeParam,outcomeParamVCOV))      
    }
    
    ## before calculating Betas from outcome model
    ## check that the exposure model design matrix variables are part of the design matrix variables from
    ## the outcome model (sometimes for categorical variables exposure model covariates do not contain 
    ## all the values presented in outcome model --> dimension no long matches)

      if(length(outcomeModelVarNames)!=(length(c(exposureModelVarNames,covOutcome)))){
        stop("At least one categorical variable in main data set does not have the same length of values as in the validation data set. This violates the positivity required for 
            for the transportability of validation model. Consider data restriction or using continuous variable for extrapolation.")
      }
    
    if(method=="cox"){
      Bstar<-t(t(outcomeParam)) #p' x 1
      VBstar<-outcomeParamVCOV[1:length(outcomeParam),1:length(outcomeParam)] # p' x p'
      BstarSebstarP<-summary(outModel)$coef[,] # entire table   
    }else{
      Bstar<-t(t(outcomeParam[2:length(outcomeParam)])) #p' x 1
      VBstar<-outcomeParamVCOV[2:length(outcomeParam),2:length(outcomeParam)] # p' x p'
      BstarSebstarP<-summary(outModel)$coef[-1,] # entire table      
    }

  }else if(supplyEstimates==TRUE){
    Bstar<-as.matrix(pointEstimates)
    VBstar<-vcovEstimates
  }
  
  
  #step 2: calibration model  
  ### exposure measurement error model, instead of using modeling, use matrix operations in a least square linear regression context.
  #### reminder: X, Z can be k x 1 dimension 

  X=as.matrix(X_VS)  #### design matrix, including intercept
  Y=as.matrix(Y_VS) # no intercept

  #### Let B = B*%*%inv(GAMMAs) is the corrected parameter vector
  #### Variance(B) = B* SIGMA(inv(GAMMAs)) t(B*) + t(inv(GAMMAs)) SIGMA(B*) inv(GAMMAs), 
  ####        where SIGMA(inv(GAMMAs)) = t(d inv(GAMMAs)/d GAMMAs) Cov(GAMMAs) (d inv(GAMMAs)/d GAMMAs)
  ####        Cov(GAMMAs)=Var(G) is just VG in the following with dimension p^2 x p^2  
  #### We have so far B* = Bstar, SIGMA(B*)=VBstar, 
  #### in the following we derive first inv(GAMMAs) and then SIGMA(inv(GAMMAs)), where GAMMAs will exclude intercept 
  
  lCovOutcome<-length(riskFactorModelVarNames)
  n<-nrow(X)   # number of obs in X
  pMeModel<-ncol(X)-1 #the original number of parameters in the calibration model except intercept. 
  p<-ncol(X)-1+lCovOutcome
  
if(length(covOutcome)==0){
    F=solve(t(X)%*%X) #shorthand
    GWI=F%*%t(X)%*%Y  #estimates of GAMMAs, the parameter matrix in exp (vector) ~ sur (vector) + covCalib (vector), including intercept, dim: (p+1) x p (p equations embedded in X ~ Z + covCalib)
    if(length(covCalib)==0){
      GEV=t(GWI[sur,])
    }else{
      GEV=t(GWI[,exp])
    }
    ERR=(Y-X%*%GWI) # dimension n x p, doe snot depend on 
  }else if(length(covOutcome)>0){  #### rewrite GAMMAs matrix and dimension of GAMMAs
    # add an identity matrix with length of covOUtcome to the right lower end of GWI (GAMMAs) matrix
    I<-diag(lCovOutcome)
    rightLowerZeroMatrix<-matrix(0,nrow=lCovOutcome,ncol=lCovOutcome)
    #Now recalculate F, GWI, GEV and ERR
    # add zeros to X matrix to obtain as filler for covOutcome variables, order matters
    F=solve(t(X)%*%X) #shorthand (original F)
    # Expand F, GWI, i.e. GAMMAs
    GWI=as.matrix(Matrix::bdiag((F%*%t(X)%*%Y),I))
    F=as.matrix(Matrix::bdiag(F,rightLowerZeroMatrix))
    # GEV=t(GWI[,exp])
    # expand X, Y
    X=cbind(X,matrix(1,nrow=nrow(X),ncol=lCovOutcome))
    Y=cbind(Y,matrix(1,nrow=nrow(Y),ncol=lCovOutcome)) 
    ERR=(Y-X%*%GWI)
    
    colnames(GWI)<-colnames(Y)
    
    if(length(covCalib)==0){
      GEV=t(GWI[,sur])
    }else{
      GEV=t(GWI[,exp])
    }
  } 
  # print(GWI)
  #### remove objects that are no longer used (free up memory)
  #remove(list=(c("X","Y")))
  
  #### calculate regular covariance matrix for the validation regression coefficients Gammas, which has dimension (p+1) x p 
  S = (t(ERR)%*%ERR)/(n-pMeModel-1) # MSE, dimension p x p, does not depend on whether length(covOutcome)>0
  if(length(covCalib)==0){
    colnames(S)=exp
    rownames(S)=exp
    VEV=S
  }else{
    VEV = S[exp,exp] #variance-covariance matrix of the error variables Var(error) Z in exp ~ sur + covCalib) + error
  }
  G = t(GWI[2:(p+1),]) # dimension: p x p; This is subset of all GAMMAs <p+1 x p> (excluding the intercepts for each exp ~ sur + covCalib)
  VG = matrixcalc::direct.prod(x=S,y=F[2:(p+1),2:(p+1)])# This is similar to the sigma^2 (X'X)^-1 in the case of unidimensional vector 
                                                        # POTENTIAL ERROR HERE - direct product between matrix S (p x p) and F[2:p+1,2:p+1]                                                       # VG is p^2 x p^2 and is the estimate of Var(G) where G again is p x p
  # remove(list=(c("S","F","ERR")))
  IGT = t(solve(t(G)))
  #### create matrix(k) to reaarange elements of (direct.dot(IGT,t(IGT))) into the order needed; 
  t = as.vector(matrix(1:p^2,nrow=p,byrow=TRUE))
  # this creates a vector that corresponds to the positions indicator of a matrix
  # .e.g., if p=3, then k is {1,4,7,2,5,8,3,6,9}, corresponding to a 
  # matrix(c(1,2,3,4,5,6,7,8,9),nrow=3,byrow=TRUE)
  
  k = design(t,p) # k is the (d GAMMAs)/(d GAMMAs) that is why it is p^2 x p^2 matrix
  m_k = k%*%(direct.prod(IGT,t(IGT))) # this is d inv(GAMMAs)/d (GAMMAs) 
                                      # = t(inv(GAMMAs))%*%(d GAMMAs)/(d GAMMAs)%*%inv(GAMMAs)
  #### Finally calculate corrected betas and variables
  #### Let B = B*%*%inv(GAMMAs) is the corrected parameter vector
  # print(IGT)
  B =  t(IGT)%*%Bstar 
  
  #### Variance(B) = B* SIGMA(inv(GAMMAs)) t(B*) + t(inv(GAMMAs)) SIGMA(B*) inv(GAMMAs), 
  ####        where SIGMA(inv(GAMMAs)) = t(d inv(GAMMAs)/d GAMMAs) Cov(GAMMAs) (d inv(GAMMAs)/d GAMMAs)
  ####        Cov(GAMMAs)=Var(G) is just VG in the following with dimension p^2 x p^2  
  VB = direct.prod(diag(p),t(Bstar))%*%m_k%*%VG%*%t(m_k)%*%t(direct.prod(diag(p),t(Bstar))) + t(IGT)%*%VBstar%*%IGT
 
#   expFormula<-paste0(exp,"~",paste0(sur,collapse="+"),"+",paste0(covCalib,collapse="+"))
#   expModel<-lm(formula=as.formula(expFormula),data=EVS)
#   expParam=coef(expModel)
#   expParamVCOV=vcov(expModel)
#   exposureModelR+-esults<-(list(expParam,expParamVCOV))
  
  ##################################################
  # Additional steps for internal validation study #
  ##################################################
  if(external==FALSE){
    # additionally run outcome ~ exp + covCalib within validation study
 #   if(survival==FALSE){
      if(length(covCalib)==0&length(covOutcome)==0){
        expFormula<-paste0("~",paste0(exp,collapse="+"))
      }else if(length(covOutcome)==0){
        expFormula<-paste0("~",paste0(exp,collapse="+"),"+",paste0(covCalib,collapse="+"))
      }else{
        expFormula<-paste0("~",paste0(exp,collapse="+"),"+",paste0(covCalib,collapse="+"),"+",paste0(covOutcome,collapse="+"))
      }
#    }else{
      # for survival need to take into account additional variables, censoring etc. 
#    }
    # regression
    if(method=="lm"){
      expModel_internal<-lm(data=vs_complete,formula=as.formula(paste0(outcome,expFormula)))
#    }else if(survival==TRUE){
    # survival
    }else{
      # vs_complete<-vs_complete%>%select(-"cvd_bi")%>%
      #           mutate(cvd_bi=rbinom(n=n(),size=1,prob=0.15))
      expModel_internal<-glm(data=vs_complete,formula=as.formula(paste0(outcome,expFormula)),family=family(link))
    }
    BV  = expModel_internal$coef[-1] #remove intercept
    VBV = vcov(expModel_internal)[-1,-1]
    SEBV = sqrt(diag(VBV))
    
    # Combine internal validation study estimates with corrected main study estimates
    # (weighted average). Weighting done with inverse covariance of est. (properly scaled);
    EST1 = B
    VEST1 = VB
    
    VE = as.matrix(Matrix::bdiag(VEST1, VBV))
    XC = rbind(diag(p),diag(p))
    # The denominator of the weight is 
    # solve(solve(VEST1) + solve(VBC)), with 
    # solve(solve(VEST1) and solve(VBC)) being the weight
    M = solve(t(XC)%*%solve(VE)%*%XC)%*%t(XC)%*%solve(VE)
    # combined (weighted) point and cov estimates
    BRCI = M%*%c(EST1, BV)
    VBRCI = M%*%VE%*%t(M)

    # overwrite original B and VB
    B = BRCI
    VB = VBRCI
    
 #   remove(list=c("EST1","VEST1","VE","XC","M","BRCI","VBRCI"))
    
    }
 
  remove(list=(c("k","m_k","GWI","VG","G")))
  # dim of outcomeCov
  
  correctedVarNames<-c(exp,setdiff(outcomeModelVarNames[-1],sur))
  # compose output tables
  SE_B<-sqrt(diag(VB))
  zValue<-B/SE_B
  pValue<-2*pnorm(-abs(as.numeric(zValue)))
  lcl<-B - qnorm(0.975)*SE_B
  ucl<-B + qnorm(0.975)*SE_B
  BSebP<-cbind(B,SE_B,zValue,pValue,lcl,ucl)
  
  colnames(BSebP)<-c("Estimate","Std. Error","Z Value","Pr(>|Z|)","lower 95%CI","upper 95%CI")
  rownames(BSebP)<-names(SE_B)
  # output objects
  if(supplyEstimates==FALSE){
    outputList<-list(BSebP,VB,BstarSebstarP,VBstar,GEV,S)
    names(outputList)<-c("correctedCoefTable","correctedVCOV","standardCoefTable","standardVCOV","calibrationModelCoefTable","calibrationModelVCOV")
  }else if(supplyEstimates==TRUE){
    outputList<-list(BSebP,VB,GEV,S)
    names(outputList)<-c("correctedCoefTable","correctedVCOV","calibrationModelCoefTable","calibrationModelVCOV")
  }
  
  return(outputList)
  
}