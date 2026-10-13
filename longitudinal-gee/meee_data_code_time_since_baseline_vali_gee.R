########################################### load and prepare main study data #################################################
library(geepack)
main_study_long<- read.csv("main_study_long.csv",header = T)
class(main_study_long)# data.frame
dim(main_study_long)# 58245*31
length(unique(main_study_long$id)) # 19415 unique id (people)

# make sure at the same starting point as the original one: only delete people who don't have pm2.5 exposure; now there are 19409 people: each person has been ordered from t= 1 to 3
main_study_long_no_pm<- main_study_long[(is.na(main_study_long$pm25c14)==TRUE),] # no pm2.5 data at all 
main_study_long<- main_study_long[(is.na(main_study_long$pm25c14)==FALSE),] 

# create a variable named "age_new", indicating varying age at each visit; "intage" meaning "baseline age"
main_study_long$age_new<- 0

for (i in 1:nrow(main_study_long)){
  
  if (main_study_long$t[i]==1){
    main_study_long$age_new[i]<- main_study_long$intage[i]
  }
  if (main_study_long$t[i]==2){
    main_study_long$age_new[i]<- main_study_long$intage[i]
  }
  if (main_study_long$t[i]==3){
    main_study_long$age_new[i]<- (main_study_long$intage[i]+main_study_long$intage2[i])/2
  }
  
}

########################################### load the validation study data ###########################################
vali_study<- read.csv("vali_study.csv",header = T)
vali_study<- vali_study[is.na(vali_study$fine_jeff)==FALSE & is.na(vali_study$fine_ambor)==FALSE & is.na(vali_study$age)==FALSE,]

vali_study$fine_jeff<- (vali_study$fine_jeff-14)/10 # surrogate exposure
vali_study$fine_ambor<- (vali_study$fine_ambor-14)/10 # true exposure 
vali_study$inter_pm25_age<- vali_study$age*vali_study$fine_jeff

###  for people from four cities + seatle 
# check localized assumption 
num_rows<- dim(vali_study)[1]
vali_study$ave_exposure<- 0
vali_study$ave_exposure[1]<- vali_study[1,]$fine_jeff
for (i in 2:num_rows){
  id_current<- vali_study[i,]$index_id
  surrogate_exposure<- vali_study[i,]$fine_jeff
  j<- (i-1)
  while(vali_study[j,]$index_id==id_current){
    surrogate_exposure<- c(surrogate_exposure,vali_study[j,]$fine_jeff)
    j<- (j-1)
  }
  if (length(surrogate_exposure)>1){
    vali_study$ave_exposure[i]<- mean(surrogate_exposure[-1])
  } else{
    vali_study$ave_exposure[i]<- NA
  }
}

vali_study_check<- na.omit(vali_study)
vali_study_check_gee <- geeglm(fine_ambor ~ fine_jeff + age + inter_pm25_age +ave_exposure, data=vali_study_check, family= gaussian(link="identity"),id=index_id, corstr = "ind", std.err="san.se")
anova(vali_study_check_gee) # sequential anova test suggests that ave_exposure doesn't have significant impact 

# remove unnecessary variables and only keep index_id, fine_ambor, fine_jeff, age, inter_pm25_age
colnames(vali_study)
#vali_study<- vali_study[,c(1,4,6,8,11)]
vali_study_gee <- geeglm(fine_ambor ~ fine_jeff + age + inter_pm25_age , data=vali_study, family= gaussian(link="identity"),id=index_id, corstr = "ind", std.err="san.se")

# correlation between c and C
cor(vali_study$fine_jeff,vali_study$fine_ambor)

#time-varying function for main study (true x)
alpha0_estimated<- vali_study_gee$coefficients[1]
alpha1_estimated<- vali_study_gee$coefficients[2]
alpha2_estimated<- vali_study_gee$coefficients[3]
alpha3_estimated<- vali_study_gee$coefficients[4]
alpha_estimated<- c(alpha0_estimated,alpha1_estimated,alpha2_estimated,alpha3_estimated)

main_study_long$pm25c14_corrected<- ((alpha0_estimated+ alpha1_estimated*(main_study_long$pm25c14)+alpha2_estimated*main_study_long$age_new+alpha3_estimated*(main_study_long$pm25c14)*main_study_long$age_new))

###  for people from four cities only 
# select records with chosen cities
vali_study_four_city<- vali_study[vali_study$city %in% c(0,1,2,3),] # this further excludes people from seatle
table(vali_study_four_city$city)
length(unique(vali_study_four_city$index_id))

vali_study<- vali_study_four_city

# check localized assumption 
num_rows<- dim(vali_study)[1]
vali_study$ave_exposure<- 0
vali_study$ave_exposure[1]<- vali_study[1,]$fine_jeff
for (i in 2:num_rows){
  id_current<- vali_study[i,]$index_id
  surrogate_exposure<- vali_study[i,]$fine_jeff
  j<- (i-1)
  while(vali_study[j,]$index_id==id_current){
    surrogate_exposure<- c(surrogate_exposure,vali_study[j,]$fine_jeff)
    j<- (j-1)
  }
  if (length(surrogate_exposure)>1){
    vali_study$ave_exposure[i]<- mean(surrogate_exposure[-1])
  } else{
    vali_study$ave_exposure[i]<- NA
  }
}

vali_study_check<- na.omit(vali_study)
vali_study_check_gee <- geeglm(fine_ambor ~ fine_jeff + age + inter_pm25_age +ave_exposure, data=vali_study_check, family= gaussian(link="identity"),id=index_id, corstr = "ind", std.err="san.se")
anova(vali_study_check_gee) # sequential anova test suggests that ave_exposure doesn't have significant impact 

# remove unnecessary variables and only keep index_id, fine_ambor, fine_jeff, age, inter_pm25_age
colnames(vali_study)
vali_study_gee <- geeglm(fine_ambor ~ fine_jeff + age + inter_pm25_age , data=vali_study, family= gaussian(link="identity"),id=index_id, corstr = "ind", std.err="san.se")

# correlation between c and C
cor(vali_study$fine_jeff,vali_study$fine_ambor)

#time-varying function for main study (true x)
alpha0_estimated<- vali_study_gee$coefficients[1]
alpha1_estimated<- vali_study_gee$coefficients[2]
alpha2_estimated<- vali_study_gee$coefficients[3]
alpha3_estimated<- vali_study_gee$coefficients[4]
alpha_estimated<- c(alpha0_estimated,alpha1_estimated,alpha2_estimated,alpha3_estimated)

main_study_long$pm25c14_corrected<- ((alpha0_estimated+ alpha1_estimated*(main_study_long$pm25c14)+alpha2_estimated*main_study_long$age_new+alpha3_estimated*(main_study_long$pm25c14)*main_study_long$age_new))

############################################################## fit the GEE model ##########################################################################################################
### time since baseline as time scale: 2-year decline. We want to assess the relationship between global cognitive score and time passed since baseline, this can be affected by baseline age, edu/husband's edu, long-term physical activity and alcohol consumption
# uncorrected: almost the same result, same as original both for full and partial model
fit.uns_z_time_since_baseline_full <- geeglm(gs ~  pm25c14*time+  intage*time + edu2*time + edu3*time + husbem*time + husbe2*time + husbe3*time + husbe4*time + husbe5*time + phymvm*time + phymv1*time + phymv2*time +  phymv3*time +  alcmvm*time + alcmv2*time + alcmv3*time + alcmv4*time, data=main_study_long, family= gaussian(link="identity"),id=id, corstr = "uns", std.err="san.se")
fit.uns_z_time_since_baseline_partial<- geeglm(gs ~  pm25c14*time+  intage*time, data=main_study_long, family= gaussian(link="identity"),id=id, corstr = "uns", std.err="san.se")

# corrected: when calculating se, needs to consider extra interaction term... but they're actually just zero
fit.uns_x_time_since_baseline_full <- geeglm(gs ~  pm25c14_corrected*time+  intage*time + edu2*time + edu3*time + husbem*time + husbe2*time + husbe3*time + husbe4*time + husbe5*time + phymvm*time + phymv1*time + phymv2*time +  phymv3*time +  alcmvm*time + alcmv2*time + alcmv3*time + alcmv4*time, data=main_study_long, family= gaussian(link="identity"),id=id, corstr = "uns", std.err="san.se")
fit.uns_x_time_since_baseline_partial<- geeglm(gs ~  pm25c14_corrected*time+  intage*time, data=main_study_long, family= gaussian(link="identity"),id=id, corstr = "uns", std.err="san.se")

######################################################### corrected se of measurement-error corrected method #####################################################################################################
### for validation study: create a new id
n2<- length(unique(vali_study$index_id))#see how many nonmissing observations we have
id_num<- 1
vali_study$id_new<- id_num
for (i in 2:dim(vali_study)[1]){
  if (vali_study$index_id[i]==vali_study$index_id[(i-1)]){
    vali_study$id_new[i]<- vali_study$id_new[(i-1)]
  }
  
  else {
    id_num<- id_num + 1 
    vali_study$id_new[i]<- id_num
  }
  
}

# changes to A(\theta)
a_star<- list()
sum_matrix_alpha<- 0
for (i in 1:n2){
  a_star[[i]]<- (cbind(rep(1,sum(vali_study$id_new==i)),vali_study[(vali_study$id_new==i),]$fine_jeff,vali_study[(vali_study$id_new==i),]$age,vali_study[(vali_study$id_new==i),]$inter_pm25_age))
  sum_matrix_alpha<- sum_matrix_alpha+ (as.matrix(t(a_star[[i]]))%*%as.matrix(vali_study[(vali_study$id_new==i),]$fine_ambor-a_star[[i]]%*%alpha_estimated))%*%t((t(as.matrix(a_star[[i]]))%*%as.matrix(vali_study[(vali_study$id_new==i),]$fine_ambor-a_star[[i]]%*%alpha_estimated)))
  print(i)
}

### for main study: 
beta_estimated<- as.vector(fit.uns_x_time_since_baseline_partial$coefficients)

# need to ensure beta3 denotes interaction terms and beta4 dentoes the confounder!! otherwise the calculation will be wrong below
beta_estimated[4]<- as.vector(fit.uns_x_time_since_baseline_partial$coefficients)[5]
beta_estimated[5]<- as.vector(fit.uns_x_time_since_baseline_partial$coefficients)[4]
num_month<- 3
p_beta<- 6
p_alpha<- 4

# unstructured correlation
v_hat<- matrix(1,nrow=num_month,ncol=num_month)
l<- 0
for (i in 1:(num_month-1)){
  for (j in (i+1):num_month){
    l<- l+1
    v_hat[i,j]<- fit.uns_x_time_since_baseline_partial$geese$alpha[l]
    v_hat[j,i]<- v_hat[i,j]
  }
}

v_hat<- v_hat*fit.uns_x_time_since_baseline_partial$geese$gamma

### in the se correction process: delete missing records and create a new id from 1 to 19002

observed<- main_study_long
observed_no_gs<- observed[(is.na(main_study_long$gs)==TRUE),] 
observed<- observed[(is.na(main_study_long$gs)==FALSE),] # final people used for analysis 

# create a new id: finished
id_num<- 1
observed$id_new<- id_num
for (i in 2:dim(observed)[1]){
  if (observed$id[i]==observed$id[(i-1)]){
    observed$id_new[i]<- observed$id_new[(i-1)]
  }
  
  else {
    id_num<- id_num + 1 
    observed$id_new[i]<- id_num
  }
  
}

n1<- length(unique(observed$id_new))

# this is used to check the missing pattern - intermittent missing (38) & dropout missing of 'gs' (6661)
observed$individual_time_calculation<- 0 
for (i in 1:nrow(observed)){
  observed[i,]$individual_time_calculation<- sum(observed[observed$id_new==observed[i,]$id_new,]$t)
  print(i)
}

# split people into observed_full_three & observed_not_three
observed_full_three<- observed[observed$individual_time_calculation==6,]
observed_not_three<- observed[observed$individual_time_calculation!=6,]
length(unique(observed_full_three$id_new))
length(unique(observed_not_three$id_new))
  
colnames(observed)[2]<- "w"
colnames(observed)[33]<- "x"
colnames(observed)[29]<- "z"
colnames(observed)[30]<- "y"
observed$interaction_x<- observed$x*observed$time
observed$interaction_z<- observed$z*observed$time
observed$interaction_w<- observed$w*observed$time
sum_matrix<- 0

x_star<- list()
sigma_hat<- list()
sigma_hat_inv<- list()

for (i in 1:n1){
  x_star[[i]]<- (cbind(rep(1,sum(observed$id_new==i)),observed[(observed$id_new==i),]$x,observed[(observed$id_new==i),]$time,observed[(observed$id_new==i),]$interaction_x,observed[(observed$id_new==i),]$w,observed[(observed$id_new==i),]$interaction_w))
  sigma_hat_inv[[i]]<- solve(v_hat[observed[(observed$id_new==i),]$t,observed[(observed$id_new==i),]$t])
  sum_matrix<- sum_matrix+ (as.matrix(t(x_star[[i]]))%*%as.matrix(sigma_hat_inv[[i]])%*%as.matrix(observed[(observed$id_new==i),]$y-x_star[[i]]%*%beta_estimated))%*%t((t(as.matrix(x_star[[i]]))%*%as.matrix(sigma_hat_inv[[i]])%*%as.matrix(observed[(observed$id_new==i),]$y-x_star[[i]]%*%beta_estimated)))
  print(i)
}

# need to save this object because it takes a rather long time to run 
# save(x_star,sigma_hat_inv,sum_matrix,file = "x_star_sigma_hat_sum_matrix.RData")

a_hat<- matrix(0,nrow=(p_alpha+p_beta),ncol=(p_alpha+p_beta))
a_hat[(1:p_alpha),(1:p_alpha)]<- sum_matrix_alpha/(n2^2)
a_hat[(p_alpha+1):(p_alpha+p_beta),(p_alpha+1):(p_alpha+p_beta)]<- sum_matrix/(n1^2)

# changes to B(\theta)
b_hat<- matrix(0,nrow=(p_alpha+p_beta),ncol=(p_alpha+p_beta))
sum_matrix_2_alpha<- 0
for (i in 1:n2){
  sum_matrix_2_alpha<- sum_matrix_2_alpha- t(a_star[[i]])%*%a_star[[i]]
}

b_hat[(1:p_alpha),(1:p_alpha)]<- sum_matrix_2_alpha/n2

sum_matrix_2<- 0
for (i in 1:n1){
  sum_matrix_2<- sum_matrix_2- t(x_star[[i]])%*%sigma_hat_inv[[i]]%*%x_star[[i]]
}

b_hat[(p_alpha+1):(p_alpha+p_beta),(p_alpha+1):(p_alpha+p_beta)]<- sum_matrix_2/n1

x_star_der_alpha0<- list()
sum_first_alpha0<- 0
sum_second_alpha0<- 0
for (i in 1:n1){
  x_star_der_alpha0[[i]]<- matrix(0,nrow=sum(observed$id_new==i),ncol=p_beta)
  x_star_der_alpha0[[i]][,2]<- rep(1,times=(sum(observed$id_new==i)))
  x_star_der_alpha0[[i]][,4]<- observed[observed$id_new==i,]$time
  
  sum_first_alpha0<- sum_first_alpha0+ t(x_star_der_alpha0[[i]])%*%sigma_hat_inv[[i]]%*%observed[observed$id_new==i,]$y
  
  sum_second_alpha0<- sum_second_alpha0+ (t(x_star_der_alpha0[[i]])%*%sigma_hat_inv[[i]]%*% x_star[[i]]+t(x_star[[i]])%*%sigma_hat_inv[[i]]%*%x_star_der_alpha0[[i]])%*%beta_estimated
  print(i)
}

sum_matrix_alpha0<- sum_first_alpha0- sum_second_alpha0

b_hat[(p_alpha+1):(p_alpha+p_beta),1]<- sum_matrix_alpha0/n1

# need to save this object because it takes a rather long time to run 
# save(a_hat,b_hat,file = "a_hat_b_hat.RData")

x_star_der_alpha1<- list()
sum_first_alpha1<- 0
sum_second_alpha1<- 0
for (i in 1:n1){
  x_star_der_alpha1[[i]]<- matrix(0,nrow=sum(observed$id_new==i),ncol=p_beta)
  
  x_star_der_alpha1[[i]][,2]<- observed[observed$id_new==i,]$z
  x_star_der_alpha1[[i]][,4]<- observed[observed$id_new==i,]$interaction_z
  
  sum_first_alpha1<- sum_first_alpha1+ t(x_star_der_alpha1[[i]])%*%sigma_hat_inv[[i]]%*% observed[observed$id_new==i,]$y
  
  sum_second_alpha1<- sum_second_alpha1+ (t(x_star_der_alpha1[[i]])%*%sigma_hat_inv[[i]]%*% x_star[[i]]+t(x_star[[i]])%*%sigma_hat_inv[[i]]%*%x_star_der_alpha1[[i]])%*%beta_estimated
  print(i) 
}

sum_matrix_alpha1<- sum_first_alpha1- sum_second_alpha1

b_hat[(p_alpha+1):(p_alpha+p_beta),2]<- sum_matrix_alpha1/n1

# need to save this object because it takes a rather long time to run 
# save(a_hat,b_hat,file = "a_hat_b_hat.RData")

x_star_der_alpha2<- list()
sum_first_alpha2<- 0
sum_second_alpha2<- 0
for (i in 1:n1){
  x_star_der_alpha2[[i]]<- matrix(0,nrow=sum(observed$id_new==i),ncol=p_beta)
  
  x_star_der_alpha2[[i]][1,2]<- observed[observed$id_new==i,]$time[1]
  x_star_der_alpha2[[i]][1,4]<- x_star_der_alpha2[[i]][1,2]*observed[observed$id_new==i,]$time[1]
  
  if (sum(observed$id_new==i)>=2){
    for (j in 2:sum(observed$id_new==i)){
      x_star_der_alpha2[[i]][j,2]<- observed[observed$id_new==i,]$time[(1:j-1)]%*%(observed[observed$id_new==i,]$time[(2:j)]-observed[observed$id_new==i,]$time[(1:j-1)])/(observed[observed$id_new==i,]$time[j]-observed[observed$id_new==i,]$time[1])
      x_star_der_alpha2[[i]][j,4]<-  x_star_der_alpha2[[i]][j,2]*observed[observed$id_new==i,]$time[j]
    }
  }
  
  
  sum_first_alpha2<- sum_first_alpha2+ t(x_star_der_alpha2[[i]])%*%sigma_hat_inv[[i]]%*% observed[observed$id_new==i,]$y
  
  sum_second_alpha2<- sum_second_alpha2+ (t(x_star_der_alpha2[[i]])%*%sigma_hat_inv[[i]]%*% x_star[[i]]+t(x_star[[i]])%*%sigma_hat_inv[[i]]%*%x_star_der_alpha2[[i]])%*%beta_estimated
  
  print(i)
}

sum_matrix_alpha2<- sum_first_alpha2- sum_second_alpha2

b_hat[(p_alpha+1):(p_alpha+p_beta),3]<- sum_matrix_alpha2/n1

# need to save this object because it takes a rather long time to run 
# save(a_hat,b_hat,file = "a_hat_b_hat.RData")

x_star_der_alpha3<- list()
sum_first_alpha3<- 0
sum_second_alpha3<- 0
for (i in 1:n1){
  x_star_der_alpha3[[i]]<- matrix(0,nrow=sum(observed$id_new==i),ncol=p_beta)
  
  x_star_der_alpha3[[i]][1,2]<- observed[observed$id_new==i,]$time[1]*observed[observed$id_new==i,]$z[1]
  x_star_der_alpha3[[i]][1,4]<- x_star_der_alpha3[[i]][1,2]*observed[observed$id_new==i,]$time[1]
  
  if (sum(observed$id_new==i)>=2){
    for (j in 2:sum(observed$id_new==i)){
      x_star_der_alpha3[[i]][j,2]<- sum(observed[observed$id_new==i,]$z[(1:(j-1))]*observed[observed$id_new==i,]$time[(1:(j-1))]*(observed[observed$id_new==i,]$time[(2:j)]-observed[observed$id_new==i,]$time[(1:(j-1))]))/(observed[observed$id_new==i,]$time[j]-observed[observed$id_new==i,]$time[1])
      x_star_der_alpha3[[i]][j,4]<-  x_star_der_alpha3[[i]][j,2]*observed[observed$id_new==i,]$time[j]
    }
  }
  
  sum_first_alpha3<- sum_first_alpha3+ t(x_star_der_alpha3[[i]])%*%sigma_hat_inv[[i]]%*% observed[observed$id_new==i,]$y
  
  sum_second_alpha3<- sum_second_alpha3+ (t(x_star_der_alpha3[[i]])%*%sigma_hat_inv[[i]]%*% x_star[[i]]+t(x_star[[i]])%*%sigma_hat_inv[[i]]%*%x_star_der_alpha3[[i]])%*%beta_estimated
  
  print(i)
}

sum_matrix_alpha3<- sum_first_alpha3- sum_second_alpha3

b_hat[(p_alpha+1):(p_alpha+p_beta),4]<- sum_matrix_alpha3/n1

# need to save this object because it takes a rather long time to run 
# save(a_hat,b_hat,file = "a_hat_b_hat.RData")

#final 

b_hat_inv<- solve(b_hat)

corrected_estimator<- b_hat_inv%*%a_hat%*%t(b_hat_inv)
corrected_se <- sqrt(diag(corrected_estimator)[(p_alpha+1):(p_alpha+p_beta)])

#################################################### summarize the results into a table ############################################################################################### 
# for naive method
naive_results<- matrix(0,nrow=p_beta,ncol=4)
colnames(naive_results)<- c("est","se","lower ci","upper ci")
rownames(naive_results)<- c("intercept","pm25c14","timesincebaseline","interaction","intage","intage_interaction")
naive_results[(1:p_beta),1]<- fit.uns_z_time_since_baseline_partial$coef
naive_results[4,1]<- fit.uns_z_time_since_baseline_partial$coef[5]
naive_results[5,1]<- fit.uns_z_time_since_baseline_partial$coef[4]

naive_results[(1:p_beta),2]<- sqrt(diag(fit.uns_z_time_since_baseline_partial$geese$vbeta))
naive_results[4,2]<- sqrt(diag(fit.uns_z_time_since_baseline_partial$geese$vbeta))[5]
naive_results[5,2]<- sqrt(diag(fit.uns_z_time_since_baseline_partial$geese$vbeta))[4]

naive_results[(1:p_beta),3]<- naive_results[1:p_beta,1]-1.96* naive_results[1:p_beta,2]
naive_results[(1:p_beta),4]<- naive_results[1:p_beta,1]+1.96* naive_results[1:p_beta,2]

# for corrected method
correction_results<- matrix(0,nrow=p_beta,ncol=4)
colnames(correction_results)<- c("est","se","lower ci","upper ci")
rownames(correction_results)<- c("intercept","pm25c14","timesincebaseline","interaction","intage","intage_interaction")
correction_results[(1:p_beta),1]<- beta_estimated
correction_results[(1:p_beta),2]<- corrected_se

correction_results[(1:p_beta),3]<- correction_results[1:p_beta,1]-1.96* correction_results[1:p_beta,2]
correction_results[(1:p_beta),4]<- correction_results[1:p_beta,1]+1.96* correction_results[1:p_beta,2]

####################################### plot relationship between cognition and pm2.5 based on the estimated model ###############################################################
### summary of main study
# spatial-temporal predicted
summary(main_study_long$pm25c14*10+14)  
sqrt(var(main_study_long$pm25c14*10+14))

summary(main_study_long$timec)

# baseline age 
intage<- NULL
for (i in 1:n1){
  intage[i]<- unique(observed[observed$id_new==i,]$intage)
  print(i)
}

summary(intage)
sqrt(var(intage))

# age
summary(main_study_long$age_new) 

# cognitive scores
summary(main_study_long$gs)
sqrt(var(main_study_long$gs))

### summary of validation study 
# personal pm2.5 (true)
summary(vali_study$fine_ambor)*10+14
sqrt(var(vali_study$fine_ambor*10+14))

# spatial-temporal predicted
summary(vali_study$fine_jeff)*10+14
sqrt(var(vali_study$fine_jeff*10+14))

# age 
summary(vali_study$age)
sqrt(var(vali_study$age))

pm25c14min<- -1.214
pm25c14mean<- 0.016
pm25c14max<- 1.145

### calculation of slopes and intercepts for lines 
# uncorrected, pm25c14min
uncorrected_pm25c14min_intercept<- 3.262 + 0.003 * pm25c14min - 0.044 * 74.2
uncorrected_pm25c14min_slope<- ((0.712) + (-0.018 * pm25c14min)+ (-0.01*74.2))/2

# uncorrected, pm25c14mean
uncorrected_pm25c14mean_intercept<- 3.262 + 0.003 * pm25c14mean - 0.044 * 74.2
uncorrected_pm25c14mean_slope<- ((0.712) + (-0.018 * pm25c14mean)+ (-0.01*74.2))/2

# uncorrected, pm25c14max
uncorrected_pm25c14max_intercept<- ((0.712) + (-0.018 * pm25c14max)+ (-0.01*74.2))/2
uncorrected_pm25c14max_slope<- ((0.712) + (-0.018 * pm25c14max)+ (-0.01*74.2))/2

# corrected, pm25c14min
corrected_pm25c14min_intercept<- 3.265 + 0.006 * pm25c14min - 0.044 * 74.2
corrected_pm25c14min_slope<- ((0.697) + (-0.029 * pm25c14min) + (-0.01*74.2))/2

# corrected, pm25c14mean
corrected_pm25c14mean_intercept<-  3.265 + 0.006 * pm25c14mean - 0.044 * 74.2
corrected_pm25c14mean_slope<- ((0.697) + (-0.029 * pm25c14mean) + (-0.01*74.2))/2

# corrected, pm25c14max
corrected_pm25c14max_intercept<-  3.265 + 0.006 * pm25c14max - 0.044 * 74.2
corrected_pm25c14max_slope<- ((0.697) + (-0.029 * pm25c14max) + (-0.01*74.2))/2

### draw the plot

# corrected max 
time_vector<- seq(from=0,to=10,by=0.01)
gs_vector<- corrected_pm25c14max_intercept  + corrected_pm25c14max_slope * time_vector
plot(time_vector,gs_vector,col="blue",lty=1,type="n",xlab="Time since Baseline (in years)",ylab="Global Cognitive Score")
lines(time_vector,gs_vector,col="blue",lty=1)

# uncorrected max 
time_vector<- seq(from=0,to=10,by=0.01)
gs_vector<- uncorrected_pm25c14max_intercept  + uncorrected_pm25c14max_slope * time_vector
lines(time_vector,gs_vector,col="red",lty=1)

# corrected mean
time_vector<- seq(from=0,to=10,by=0.01)
gs_vector<- corrected_pm25c14mean_intercept  + corrected_pm25c14mean_slope * time_vector
plot(time_vector,gs_vector,col="blue",lty=1,type="n",xlab="Time since Baseline (in years)",ylab="Global Cognitive Score")
lines(time_vector,gs_vector,col="blue",lty=1)

# uncorrected mean
time_vector<- seq(from=0,to=10,by=0.01)
gs_vector<- uncorrected_pm25c14mean_intercept  + uncorrected_pm25c14mean_slope * time_vector
lines(time_vector,gs_vector,col="red",lty=1)

# corrected_min
time_vector<- seq(from=0,to=10,by=0.01)
gs_vector<- corrected_pm25c14min_intercept  + corrected_pm25c14min_slope * time_vector
lines(time_vector,gs_vector,col="blue",lty=4)

# uncorrected_min
time_vector<- seq(from=0,to=10,by=0.01)
gs_vector<- uncorrected_pm25c14min_intercept  + uncorrected_pm25c14min_slope * time_vector
lines(time_vector,gs_vector,col="red",lty=4)

# use solid/dotted lines to differentiate corrected and uncorrected lines on one plot
# corrected mean
time_vector<- seq(from=0,to=10,by=0.01)
gs_vector<- corrected_pm25c14mean_intercept  + corrected_pm25c14mean_slope * time_vector
plot(time_vector,gs_vector,col="blue",lty=1,
     type="n",
     xlab="Time since Baseline (in years)",ylab="Global Cognitive Score")
lines(time_vector,gs_vector,lty=1)

# uncorrected mean
time_vector<- seq(from=0,to=10,by=0.01)
gs_vector<- uncorrected_pm25c14mean_intercept  + uncorrected_pm25c14mean_slope * time_vector
lines(time_vector,gs_vector,lty=2)
legend("bottomleft",
       legend = c("Corrected - PM2.5 mean", "Uncorrected - PM2.5 mean"),
       lty = c(1,2),
       pt.cex = 1,
       bty = "n",
       text.col = "black",
       horiz = F
)

########################################### summary of key variables for all people (19002) ###########################################
### education
index<- unique(observed$id_new)
edu_2<- NULL
edu_3<- NULL
for (i in 1:n1){
  edu_2[i]<- unique(observed[(observed$id_new==index[i]),]$edu2)
  edu_3[i]<- unique(observed[(observed$id_new==index[i]),]$edu3)
  print(i)
}

sum(edu_2)
sum(edu_3)

### husband's education
husbe_2<- NULL
husbe_3<- NULL
husbe_4<- NULL
husbe_5<- NULL
husbe_m<- NULL
for (i in 1:n1){
  husbe_2[i]<- unique(observed[(observed$id_new==index[i]),]$husbe2)
  husbe_3[i]<- unique(observed[(observed$id_new==index[i]),]$husbe3)
  husbe_4[i]<- unique(observed[(observed$id_new==index[i]),]$husbe4)
  husbe_5[i]<- unique(observed[(observed$id_new==index[i]),]$husbe5)
  husbe_m[i]<- unique(observed[(observed$id_new==index[i]),]$husbem)
  print(i)
}

sum(husbe_2)
sum(husbe_3)
sum(husbe_4)
sum(husbe_5)
sum(husbe_m)

### alcohol consumption 
alc_2<- NULL
alc_3<- NULL
alc_4<- NULL
alc_m<- NULL
for (i in 1:n1){
  alc_2[i]<- unique(observed[(observed$id_new==index[i]),]$alcmv2)
  alc_3[i]<- unique(observed[(observed$id_new==index[i]),]$alcmv3)
  alc_4[i]<- unique(observed[(observed$id_new==index[i]),]$alcmv4)
  alc_m[i]<- unique(observed[(observed$id_new==index[i]),]$alcmvm)
  print(i)
}

sum(alc_2)
sum(alc_3)
sum(alc_4)
sum(alc_m)

### physical activity 
phy_1<- NULL
phy_2<- NULL
phy_3<- NULL
phy_m<- NULL
for (i in 1:n1){
  phy_1[i]<- unique(observed[(observed$id_new==index[i]),]$phymv1)
  phy_2[i]<- unique(observed[(observed$id_new==index[i]),]$phymv2)
  phy_3[i]<- unique(observed[(observed$id_new==index[i]),]$phymv3)
  phy_m[i]<- unique(observed[(observed$id_new==index[i]),]$phymvm)
  print(i)
}

sum(phy_1)
sum(phy_2)
sum(phy_3)
sum(phy_m)


########################################### summary of key variables for people with all three (12303) & not three cognitive assessments (6699) ###########################################
### summary of pm2.5
# all three assessments
summary(observed_full_three$pm25c14*10+14)  
sqrt(var(observed_full_three$pm25c14*10+14))

# not_three 
summary(observed_not_three$pm25c14*10+14)  
sqrt(var(observed_not_three$pm25c14*10+14))

### summary of baseline age 
# full_three 
intage_full_three<- NULL
n1_full_three<- length(unique(observed_full_three$id_new))
index_full_three<- unique(observed_full_three$id_new)

for (i in 1:n1_full_three){
  intage_full_three[i]<- unique(observed_full_three[(observed_full_three$id_new==index_full_three[i]),]$intage)
  print(i)
}

summary(intage_full_three)
sqrt(var(intage_full_three))

# not_three 
intage_not_three<- NULL
n1_not_three<- length(unique(observed_not_three$id_new))
index_not_three<- unique(observed_not_three$id_new)

for (i in 1:n1_not_three){
  intage_not_three[i]<- unique(observed_not_three[(observed_not_three$id_new==index_not_three[i]),]$intage)
  print(i)
}

summary(intage_not_three)
sqrt(var(intage_not_three))

### summary of global cognitive scores
# full_three 
summary(observed_full_three$gs)
sqrt(var(observed_full_three$gs))

# not_three 
summary(observed_not_three$gs)
sqrt(var(observed_not_three$gs))

### education
# full_three 
edu_full_three_2<- NULL
edu_full_three_3<- NULL
for (i in 1:n1_full_three){
  edu_full_three_2[i]<- unique(observed_full_three[(observed_full_three$id_new==index_full_three[i]),]$edu2)
  edu_full_three_3[i]<- unique(observed_full_three[(observed_full_three$id_new==index_full_three[i]),]$edu3)
  print(i)
}

sum(edu_full_three_2)
sum(edu_full_three_3)

# not_three 
edu_not_three_2<- NULL
edu_not_three_3<- NULL
for (i in 1:n1_not_three){
  edu_not_three_2[i]<- unique(observed_not_three[(observed_not_three$id_new==index_not_three[i]),]$edu2)
  edu_not_three_3[i]<- unique(observed_not_three[(observed_not_three$id_new==index_not_three[i]),]$edu3)
  print(i)
}

sum(edu_not_three_2)
sum(edu_not_three_3)

# 19002 - above = 14265

### husband's education

# full_three 
husbe_full_three_2<- NULL
husbe_full_three_3<- NULL
husbe_full_three_4<- NULL
husbe_full_three_5<- NULL
husbe_full_three_m<- NULL
for (i in 1:n1_full_three){
  husbe_full_three_2[i]<- unique(observed_full_three[(observed_full_three$id_new==index_full_three[i]),]$husbe2)
  husbe_full_three_3[i]<- unique(observed_full_three[(observed_full_three$id_new==index_full_three[i]),]$husbe3)
  husbe_full_three_4[i]<- unique(observed_full_three[(observed_full_three$id_new==index_full_three[i]),]$husbe4)
  husbe_full_three_5[i]<- unique(observed_full_three[(observed_full_three$id_new==index_full_three[i]),]$husbe5)
  husbe_full_three_m[i]<- unique(observed_full_three[(observed_full_three$id_new==index_full_three[i]),]$husbem)
  print(i)
}

sum(husbe_full_three_2)
sum(husbe_full_three_3)
sum(husbe_full_three_4)
sum(husbe_full_three_5)
sum(husbe_full_three_m)

# not_three
husbe_not_three_2<- NULL
husbe_not_three_3<- NULL
husbe_not_three_4<- NULL
husbe_not_three_5<- NULL
husbe_not_three_m<- NULL
for (i in 1:n1_not_three){
  husbe_not_three_2[i]<- unique(observed_not_three[(observed_not_three$id_new==index_not_three[i]),]$husbe2)
  husbe_not_three_3[i]<- unique(observed_not_three[(observed_not_three$id_new==index_not_three[i]),]$husbe3)
  husbe_not_three_4[i]<- unique(observed_not_three[(observed_not_three$id_new==index_not_three[i]),]$husbe4)
  husbe_not_three_5[i]<- unique(observed_not_three[(observed_not_three$id_new==index_not_three[i]),]$husbe5)
  husbe_not_three_m[i]<- unique(observed_not_three[(observed_not_three$id_new==index_not_three[i]),]$husbem)
  print(i)
}

sum(husbe_not_three_2)
sum(husbe_not_three_3)
sum(husbe_not_three_4)
sum(husbe_not_three_5)
sum(husbe_not_three_m)

### alcohol consumption 

# full three
alc_full_2<- NULL
alc_full_3<- NULL
alc_full_4<- NULL
alc_full_m<- NULL
for (i in 1:n1_full_three){
  alc_full_2[i]<- unique(observed_full_three[(observed_full_three$id_new==index_full_three[i]),]$alcmv2)
  alc_full_3[i]<- unique(observed_full_three[(observed_full_three$id_new==index_full_three[i]),]$alcmv3)
  alc_full_4[i]<- unique(observed_full_three[(observed_full_three$id_new==index_full_three[i]),]$alcmv4)
  alc_full_m[i]<- unique(observed_full_three[(observed_full_three$id_new==index_full_three[i]),]$alcmvm)
  print(i)
}

sum(alc_full_2)
sum(alc_full_3)
sum(alc_full_4)
sum(alc_full_m)

# not three
alc_not_2<- NULL
alc_not_3<- NULL
alc_not_4<- NULL
alc_not_m<- NULL
for (i in 1:n1_not_three){
  alc_not_2[i]<- unique(observed_not_three[(observed_not_three$id_new==index_not_three[i]),]$alcmv2)
  alc_not_3[i]<- unique(observed_not_three[(observed_not_three$id_new==index_not_three[i]),]$alcmv3)
  alc_not_4[i]<- unique(observed_not_three[(observed_not_three$id_new==index_not_three[i]),]$alcmv4)
  alc_not_m[i]<- unique(observed_not_three[(observed_not_three$id_new==index_not_three[i]),]$alcmvm)
  print(i)
}

sum(alc_not_2)
sum(alc_not_3)
sum(alc_not_4)
sum(alc_not_m)

### physical activity 

# full three
phy_full_1<- NULL
phy_full_2<- NULL
phy_full_3<- NULL
phy_full_m<- NULL
for (i in 1:n1_full_three){
  phy_full_1[i]<- unique(observed_full_three[(observed_full_three$id_new==index_full_three[i]),]$phymv1)
  phy_full_2[i]<- unique(observed_full_three[(observed_full_three$id_new==index_full_three[i]),]$phymv2)
  phy_full_3[i]<- unique(observed_full_three[(observed_full_three$id_new==index_full_three[i]),]$phymv3)
  phy_full_m[i]<- unique(observed_full_three[(observed_full_three$id_new==index_full_three[i]),]$phymvm)
  print(i)
}

sum(phy_full_1)
sum(phy_full_2)
sum(phy_full_3)
sum(phy_full_m)

# not three
phy_not_1<- NULL
phy_not_2<- NULL
phy_not_3<- NULL
phy_not_m<- NULL
for (i in 1:n1_not_three){
  phy_not_1[i]<- unique(observed_not_three[(observed_not_three$id_new==index_not_three[i]),]$phymv1)
  phy_not_2[i]<- unique(observed_not_three[(observed_not_three$id_new==index_not_three[i]),]$phymv2)
  phy_not_3[i]<- unique(observed_not_three[(observed_not_three$id_new==index_not_three[i]),]$phymv3)
  phy_not_m[i]<- unique(observed_not_three[(observed_not_three$id_new==index_not_three[i]),]$phymvm)
  print(i)
}

sum(phy_not_1)
sum(phy_not_2)
sum(phy_not_3)
sum(phy_not_m)
