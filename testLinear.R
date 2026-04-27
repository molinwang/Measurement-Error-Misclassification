#######################
# testLinear function #
#######################

# authors: Wenze Tang and Molin Wang
# version 1: Aug 3, 2021

# Abstract:
# This function fits restricted cubic splines to linear and generalized linear models to detect non-linear relationship between 
# outcome (dependent variable) and exposure, a single independent variable. It allows for controlling for covariates. The output includes (1) a set of p-values
# from the likelihood ratio tests for non-linearity (test 1), overall association (test 2) and linear relation between independent variable and outcome (test 3) 
# and (2) a graph of the predicted outcome values and pointwise confidence bands versus exposure of interest, together with a density plot of the exposure. 
# For binary outcome, the predicted outcome is odds ratio. 

#' Use restricted cubic spline model to test linear assumption on linear, log and logit scale and graph fitted curve.   
#' @param nknots=5 Number of knots. Default is 5. The minimum value is 3 (i.e. one inner knot two outer knots).Required.   
#' @param knotsvalues knot locations. If not given, knots will be estimated using default quantiles of `var`. For 3 knots, the outer quantiles used are 0.10 and 0.90. 
#' For 4-6 knots, the outer quantiles used are 0.05 and 0.95. For nk>6, the outer quantiles are 0.025 and 0.975. The knots are equally spaced between these 
#' on the quantile scale. For fewer than 100 non-missing values of `var`, the outer knots are the 5th smallest and largest `var`. Must be unique values. 
#' @param ds The input dataframe. This dataframe should minimally include variables indicated in `adj`, `var` and `outcome`.Required.   
#' @param var single character representing independent variable being evaluated for linear relationship with `outcome`. Required.   
#' @param outcome single character representing dependent variable being evaluated for linear relationship with `var`.Required.   
#' @param adj=NULL Character vector of covariates to be adjusted in the model. For graphing purpose, mean values will be used for numeric values and 
#' mode will be used for binary and character variables.   
#' @param method="lm" Methods for modeling, currently only `lm` or `glm` methods are available. Required.    
#' @param family Family parameter to pass to glm function. Not a character. Required if method="glm".
#' @param link Link parameter to pass to glm function. Should be character. Required if method="glm". 
#' @param ref Reference value to be used. If not specified, use minimum of exposure values. Optional. 
#' @param densplot Produce exposure density plot. Default is FALSE. 

#' @return Two-object list that includes a printable dataframe of non-linearity test results and a fitted curve graph.  
#

testLinear<-function(ds,var,outcome,adj=NULL,nknots=5,knotsvalues=NULL,method="lm", family,link,ref=NULL,densplot=FALSE){
    require(dplyr)
    require(lmtest)
    require(Hmisc)
    require(ggpubr)
    require(data.table)    
  
    #### data preparation
    getmode <- function(v) {
      uniqv <- unique(v)
      uniqv[which.max(tabulate(match(v, uniqv)))]
    }
    #### (1) assess fqalc
    ##### pre-specify knots (use quintiles)
        ##### for testing purpose
        # nknots=5
        # outcome="dralc"
        # var="fqalc"
        # adj=c("fqcal","fqtot","agec")
        # ds=valid
        # method="lm"
        # # family=binomial
        # # link="logit"
        # knotsvalues=NULL
        # densplot=TRUE
        # ref=5
        # knotsvalues=quantile((ds%>%dplyr::select(var))[,1],c(0.1,0.3,0.5,0.7,0.9))

    
    # test6<-testLinear(ds=main,outcome="newcase",var="fqalc",adj=c("fqtot","fqcal","agec"),knotsvalues=quantile(main$fqalc,c(0.1,0.3,0.5,0.7,0.9)),method="glm",family=binomial,link="logit")
    # 
    
    
    if(length(knotsvalues)==0){ #if we are going to use nkots only but not knotsvalues, 
      # then in case default percentile value is not unique, only use the unique cutoff values and update knot numbers
      knotsvalues=rcspline.eval(pull(ds%>%dplyr::select(var)),nk=nknots,knots.only=TRUE)
    }
        
    if(length(unique(knotsvalues))!=nknots){
      warning("Default knot values are not unique. Right now we only use unique knot values. You can supply your own knot values or reduce number of knots (>=3).")
      nknots=length(unique(knotsvalues))
    }
    rcsTerms<-paste0("rcs",seq(1:(nknots-2))) # cubic spline terms
    
    ###########################
    # Begin Test of Linearity #
    ###########################
    # data preparation
    # create restricted cubic spline data based on original data and knots
    if(length(adj)>0){
      dsRCspline<-cbind(ds[,c(outcome,adj)],rcspline.eval(pull(ds%>%dplyr::select(var)),nk=nknots,knots=knotsvalues,inclx=TRUE))%>%as.data.frame()
    }else{
      dsRCspline<-cbind(ds[,c(outcome)],rcspline.eval(pull(ds%>%dplyr::select(var)),nk=nknots,knots=knotsvalues,inclx=TRUE))%>%as.data.frame()
    }

    colnames(dsRCspline)<-c(outcome,adj,var,rcsTerms)
    covAllTerms<-c(adj,var,rcsTerms)
    covLinearTerms<-c(adj,var)
    # test: only keep first RCSpline term
    # covAllTerms<-c(adj,var,rcsTerms[1])
    # covLinearTerms<-c(adj,var)    

    
    # for(i in 1:(nknots-1)){
    #   nameRCS=paste0("rcs",i)
    #   ds<-ds%>%mutate(!!nameRCS:=ifelse(get(var)>=knotsValueToUse[i]&get(var)<knotsValueToUse[i+1],(get(var)-knotsValueToUse[i])^3,0)
    #                   )
    #   if(i==(nknots-1)){
    #     ds<-ds%>%mutate(lastLinearTerm=ifelse(get(var)>=knotsValueToUse[nknots],(get(var)-knotsValueToUse[nknots]),0)
    #                     )
    #   }
    # }
    
    if(method=="lm"){
      # 1. Test of non-linear terms: Use LRT to compare model with vs w/o non-linear terms
      formulaFull<-paste0(outcome,"~",paste0(covAllTerms,collapse="+"))
      formulaLinearOnly<-paste0(outcome,"~",paste0(covLinearTerms,collapse="+"))
      
      modelFull<-lm(data=dsRCspline, formula=as.formula(formulaFull))
      modelLinearOnly<-lm(data=dsRCspline, formula=as.formula(formulaLinearOnly))
      test1<-lrtest(modelLinearOnly,modelFull)
      test1.df<-test1$Df[2]
      test1.Chisq<-test1$Chisq[2]
      test1.p<-test1$`Pr(>Chisq)`[2]
      
      # 2. Test of overall curvature: Use LRT to compare full model with model with covariates only
      # formulaFull<-paste0(outcome,"~",var,"+lastLinearTerm +", paste0(rcsTerms,collapse="+"),"+",paste0(adj,collapse="+"))
      if(length(adj)>0){
        formulaAdjOnly<-paste0(outcome,"~",paste0(adj,collapse="+"))
      }else{formulaAdjOnly<-paste0(outcome,"~ 1")}

      # modelFull<-lm(data=dsRCspline, formula=as.formula(formulaFull))
      modelAdjOnly<-lm(data=dsRCspline, formula=as.formula(formulaAdjOnly))
      test2<-lrtest(modelAdjOnly,modelFull)
      test2.df<-test2$Df[2]
      test2.Chisq<-test2$Chisq[2]
      test2.p<-test2$`Pr(>Chisq)`[2]
      
      # 3. Test of linear term: Use LRT to compare model w/o non-linear terms with model with covariates only
      test3<-lrtest(modelAdjOnly,modelLinearOnly)
      test3.df<-test3$Df[2]
      test3.Chisq<-test3$Chisq[2]
      test3.p<-test3$`Pr(>Chisq)`[2]      
      
    }else{
      # 1. Test of non-linear terms: Use LRT to compare model with vs w/o non-linear terms
      formulaFull<-paste0(outcome,"~",paste0(covAllTerms,collapse="+"))
      formulaLinearOnly<-paste0(outcome,"~",paste0(covLinearTerms,collapse="+"))
      
      modelFull<-glm(data=dsRCspline, formula=as.formula(formulaFull),family=family(link=link))
      modelLinearOnly<-glm(data=dsRCspline, formula=as.formula(formulaLinearOnly),family=family(link=link))
      test1<-lrtest(modelLinearOnly,modelFull)
      test1.df<-test1$Df[2]
      test1.Chisq<-test1$Chisq[2]
      test1.p<-test1$`Pr(>Chisq)`[2]
      
      # 2. Test of overall curvature: Use LRT to compare full model with model with only covariates
      # formulaFull<-paste0(outcome,"~",var,"+lastLinearTerm +", paste0(rcsTerms,collapse="+"),"+",paste0(adj,collapse="+"))
      if(length(adj)>0){
        formulaAdjOnly<-paste0(outcome,"~",paste0(adj,collapse="+"))
      }else{formulaAdjOnly<-paste0(outcome,"~ 1")}
      modelFull<-glm(data=dsRCspline, formula=as.formula(formulaFull),family=family(link=link))
      modelAdjOnly<-glm(data=dsRCspline, formula=as.formula(formulaAdjOnly),family=family(link=link))
      test2<-lrtest(modelAdjOnly,modelFull)
      test2.df<-test2$Df[2]
      test2.Chisq<-test2$Chisq[2]
      test2.p<-test2$`Pr(>Chisq)`[2]
      
      # 3. Test of linear term: Use LRT to compare model w/o non-linear terms with model with covariates only
      test3<-lrtest(modelAdjOnly,modelLinearOnly)
      test3.df<-test3$Df[2]
      test3.Chisq<-test3$Chisq[2]
      test3.p<-test3$`Pr(>Chisq)`[2]      
    }
    
    # create an dataframe object with test two test results; 
    testResults<-rbind(cbind(test1.Chisq,test1.df,round(test1.p,4)),
                       cbind(test2.Chisq,test2.df,round(test2.p,4)),
                       cbind(test3.Chisq,test3.df,round(test3.p,4)))
    colnames(testResults)<-c("LRT-based chisq","df","p value")
    rownames(testResults)<-c("Test of non-linear cubic spline terms","Test of overall curvature (i.e. all cubic spline terms)","Test of linear term only.")
    
    ##################
    # Begin Graphing #
    ##################
    
    # prepare for graphing
    #### range of mismeasured exposure
    lims<-range(ds%>%dplyr::select(eval(var)))
    grid<-seq(from=lims[1], to = lims[2],length.out=150)
    gridlength<-length(grid)

    # nknots=5

    grid.CBspline<-rcspline.eval(grid,knots=knotsvalues,inclx=TRUE) # we use original knot Values


    ##### find mean and mode for covariates that are not being assessed for linearity (this will be used for lm model only). 
    if(length(adj)>0){
      adjNewData<-matrix(NA,ncol=length(adj),nrow=gridlength)%>%as.data.frame()
      
      for(i in 1:length(adj)){
        cov=adj[i]
        vectorOfInterest= (ds%>%dplyr::select(eval(cov)))[,1]
        if(class(vectorOfInterest)=="numeric"&length(unique(vectorOfInterest))!=2){
          cov.grid<-as.numeric(rep(times=gridlength,mean(vectorOfInterest)))
        }else if(length(unique(vectorOfInterest))==2|class(vectorOfInterest)!="numeric"){ # if binary (or 2 category), or non-numeric variables use mode
          cov.grid<-rep(times=gridlength,getmode(vectorOfInterest))
        }
        # print(class(cov.grid))
        adjNewData[,i]<-cov.grid
      }
      
      newData<-cbind(adjNewData,grid.CBspline)%>%as.data.frame()%>%as.data.frame()
      
    }else{newData<-grid.CBspline%>%as.data.frame()}

    # get prediction data ready 
    ## (1) newData will be used for lm prediction
    colnames(newData)<-c(adj,var,rcsTerms)
    covariateValues<-newData[1,adj]
    # test: only keep first RCSpline term
    # newData<-cbind(adjNewData,grid.CBspline[,1:2])%>%as.data.frame()
    # colnames(newData)<-c(adj,var,rcsTerms[1])
    
    ## (2) for glm prediction, we will calculate predicted value based on (self-selected) reference level
    refValue=lims[1]
    if(length(ref)!=0){
      refValue=ref
    }
    
    ### agumented reference value data 

    HrefValue<-rcspline.eval(refValue,knots=knotsvalues,inclx=TRUE)
    HrefValues<-matrix(rep(HrefValue,times=gridlength),nrow=gridlength,byrow=TRUE)%>%as.data.frame()

    # name spline terms
    RCSplineTerms<-c(var,rcsTerms)
    
    # Test only with first term
    # RCSplineTerms<-c(var,rcsTerms[1]) #note the var here is actully x-x', where x' is the min value of all x's
    
    # construct data
    reducedData<-newData[,RCSplineTerms]-HrefValues
    # extract betas and cov matrix
    pointEstimates=coef(modelFull)[RCSplineTerms]
    vcovEstimates=vcov(modelFull)[RCSplineTerms,RCSplineTerms]
    
    ##### predict
    if(method=="lm"){
      pd<-predict(modelFull,newdata =newData,interval="confidence",level=0.95)
      pd<-cbind(grid,pd)%>%as.data.frame 
      ylims<-c(min(pd$lwr),max(pd$upr))    
      
    }else if(method=="glm"&(link=="logit"|link=="log")){

      
      # predict 
      pd<-data.table(logOR=as.vector(as.matrix(reducedData)%*%pointEstimates),
                     varLogOR=diag(as.matrix(reducedData)%*%vcovEstimates%*%t(reducedData))
      )%>%mutate(
        fit=exp(logOR),
        lwr=exp(logOR - 1.96*sqrt(varLogOR)),
        upr=exp(logOR + 1.96*sqrt(varLogOR))        
      )%>%dplyr::select(fit,lwr,upr)%>%as.data.frame()
      
      ylims<-c(0,min(10,max(pd$upr)))    
      
    }else{

      pd<-predict.glm(modelFull,newdata =newData,type="response",se.fit=TRUE)%>%as.data.frame%>%
        mutate(lwr=fit - 1.96*se.fit,
               upr=fit + 1.96*se.fit)%>%
        dplyr::select(c("fit","lwr","upr"))
      pd<-cbind(grid,pd)%>%as.data.frame    
      ylims<-c(0,min(10,max(pd$upr)))    
    }

    ##### plotting the Regression Line to the scatterplot   
    
    if(method=="lm"){
      p1 <- ggplot(data=pd, aes(grid,fit)) +
        geom_point(size = 1) +
        geom_line(aes(y = pd[,"fit"]), size = 1.2) +
        geom_ribbon(aes(ymin = pd[,"lwr"], ymax = pd[,"upr"]), alpha = 0.2, fill="darkgreen") +
        ggtitle("Fitted Restricted Cubic Spline Line and Its Pointwise 95% Confidence  Interval") +
        coord_cartesian(xlim=lims,ylim=ylims) +
        geom_vline(xintercept =knotsvalues , linetype="twodash", 
                   color = "darkblue", size=0.5) +
        xlab(var) + ylab(outcome)   
    }else if(method!="lm"&length(link)!=0){
      if(link=="logit"){
        p1 <- ggplot(data=pd, aes(grid,fit)) +
          geom_point(size = 1) +
          geom_line(aes(y = pd[,"fit"]), size = 1.2) +
          geom_ribbon(aes(ymin = pd[,"lwr"], ymax = pd[,"upr"]), alpha = 0.2, fill="darkgreen") +
          ggtitle("Fitted Restricted Cubic Spline Line and Its Pointwise 95% Confidence Interval") +
          coord_cartesian(xlim=lims,ylim=ylims) +
          geom_vline(xintercept =knotsvalues , linetype="twodash", 
                     color = "darkblue", size=0.5) +
          geom_hline(yintercept = 1 , linetype="solid", 
                     color = "lightblue", size=0.7) +
          xlab(var) + ylab(paste0("Odds Ratio: ",outcome)) 
        
        
      }else if(link=="log"){
        p1 <- ggplot(data=pd, aes(grid,fit)) +
          geom_point(size = 1) +
          geom_line(aes(y = pd[,"fit"]), size = 1.2) +
          geom_ribbon(aes(ymin = pd[,"lwr"], ymax = pd[,"upr"]), alpha = 0.2, fill="darkgreen") +
          ggtitle("Fitted Restricted Cubic Spline Line and Its Pointwise 95% Confidence Interval") +
          coord_cartesian(xlim=lims,ylim=ylims) +
          geom_vline(xintercept =knotsvalues , linetype="twodash", 
                     color = "darkblue", size=0.5) +
          geom_hline(yintercept =1 , linetype="solid", 
                     color = "lightblue", size=0.7) +
          xlab(var) + ylab(paste0("Incidence Rate Ratio: ",outcome))      
      }
         
    }

    p2<-ggplot(data=ds,aes(get(var))) + ggtitle("Density of the Exposure") +
      geom_density(bw="ucv") + 
      geom_vline(xintercept =knotsvalues , linetype="twodash", 
                 color = "darkblue", size=0.5) +
      coord_cartesian(xlim=lims) +
      xlab(var) + ylab("density")
    
    graph<-p1
    
    if(densplot==TRUE){
      graph<-ggarrange(p1,p2, labels=c("A","B"),ncol=1,nrow=2,heights=c(2,1))
    }
    
    outputList<-list(testResults,graph,covariateValues,knotsvalues)
    names(outputList)<-c("testOfLinearity","fittedCBSplineCurve","covariateFixedValues","knotValues")
    
    return(outputList)

}

