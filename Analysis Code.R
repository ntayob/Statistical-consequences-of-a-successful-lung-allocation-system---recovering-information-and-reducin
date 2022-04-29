library(gee)
library(nleqslv)
library(gdata)
library(survival)

#############################################
# Load functions
#############################################
sum_function=function(X)
{
  temp=array(X,c(length(X),length(X)))
  lowerTriangle(temp,diag=F)=0
  apply(temp,2,sum,na.rm=TRUE)
}

#Function for Step 1a: Fit a Cox model for censored outcomes
get_beta_censor=function(X,delta,t,Z1,Z2,Z3,W)
{
  n=length(X)
  delta_censor=1-delta
  start_time=rep(t,rep(n,length(t)))
  temp=rep(X,length(t))-start_time
  temp2=c(rep(t[-1],rep(n,length(t)-1)),X)
  stop_time=apply(cbind(temp,temp2),1,min)
  start_time_2=start_time
  stop_time_2=stop_time
  start_time_2[start_time>=stop_time]=NA
  stop_time_2[start_time>=stop_time]=NA
  delta_time_dependent=rep(delta_censor,length(t))
  delta_time_dependent[stop_time_2==temp2]=0
  Z=c(Z1,Z2,Z3)
  W1=c(W,W,W)
  censor_model=coxph(Surv(start_time_2,stop_time_2,delta_time_dependent)~Z+W1)
  #print(censor_model)
  coef(censor_model)
}

#Function for Step 1b and c: Construct inverse weights and compute the IPCW survival estimate
get_IPCW_surv=function(X,delta,t,Z1,Z2,Z3,W,beta_censor)
{
  n=length(X)
  time=sort(unique(X))
  n_time=length(time)
  
  temp3=array(X,c(n,n_time))
  temp4=t(array(time,c(n_time,n)))
  delta_array=array(delta,c(n,n_time))
  dN_T=array(as.numeric(temp3==temp4 & delta_array==1),c(n,n_time))
  dN_Q=array(as.numeric(temp3==temp4 & delta_array==0),c(n,n_time))
  Y=array(as.numeric(temp3>=temp4),c(n,n_time))
  n_Z_1=sum(time<t[2])
  n_Z_2=sum(time<t[3] & time>t[2])
  n_Z_3=n_time-(n_Z_1+n_Z_2)
  Z_array=array(rep(Z1,n_Z_1),c(n,n_Z_1))
  Z_array=cbind(Z_array,array(rep(Z2,n_Z_2),c(n,n_Z_2)))
  Z_array=cbind(Z_array,array(rep(Z3,n_Z_3),c(n,n_Z_3)))
  
  W_array=array(W,c(n,n_time))
  
  denominator_weight=apply(Y*exp(beta_censor[1]*Z_array + beta_censor[2]*W_array),2,sum,na.rm=TRUE)
  numerator_weight=t(exp(beta_censor[1]*Z_array + beta_censor[2]*W_array))*apply(dN_Q,2,sum,na.rm=TRUE)
  weight=exp(apply(numerator_weight/denominator_weight,2,sum_function))
  
  denominator=apply(Y*t(weight),2,sum,na.rm=TRUE)
  temp5=apply(t(dN_T)*weight/denominator,1,sum,na.rm=TRUE)
  IPCW_CH=sum_function(temp5)
  IPCW_surv=exp(-IPCW_CH)
  list(time=time, IPCW_surv=IPCW_surv,weight=weight)
}

#Function for Step 2
get_IPCW_pseudo_obs=function(X,delta,t,Z1,Z2,Z3,W)
{
  n=length(X)
  
  beta_censor=get_beta_censor(X,delta,t,Z1,Z2,Z3,W)
  IPCW_surv_results=get_IPCW_surv(X,delta,t,Z1,Z2,Z3,W,beta_censor)
  
  time=IPCW_surv_results$time
  IPCW_surv=IPCW_surv_results$IPCW_surv
  
  surv_times=seq(from=0,to=t[length(t)]+Tau,by=0.001)
  a=seq(from=0.0001,to=Tau,by=0.001)
  K=length(a)-1
  d=a[2:(K+1)]-a[1:K]
  d=c(0,d,0)
  b=(d[1:(K+1)]+d[2:(K+2)])/2
  IPCW_survival=rep(1,length(surv_times))
  for(j in 2:length(surv_times))
  {
    if(time[1]<=surv_times[j])
    {
      IPCW_survival[j]=min(IPCW_surv[time<=surv_times[j]])
    }
  }
  dIPCW_survival=IPCW_survival[2:length(surv_times)]-IPCW_survival[1:(length(surv_times)-1)]
  #Get Pseudo obs
  delta_hat_1=-sum(b*log(a)*dIPCW_survival[(surv_times<Tau)]) + log(Tau)*min(IPCW_surv[time<=Tau])
  delta_hat_2=(-sum(b*log(a)*dIPCW_survival[(surv_times<(Tau+t[2]) & surv_times>=t[2])]) + log(Tau)*min(IPCW_surv[time<=(Tau+t[2])]))/min(IPCW_surv[time<=t[2]])
  delta_hat_3=(-sum(b*log(a)*dIPCW_survival[(surv_times<(Tau+t[3]) & surv_times>=t[3])]) + log(Tau)*min(IPCW_surv[time<=(Tau+t[3])]))/min(IPCW_surv[time<=t[3]])
  delta_hat=c(delta_hat_1,delta_hat_2,delta_hat_3)
  list(beta_censor=beta_censor,IPCW_survival=IPCW_survival,delta_hat=delta_hat)
}

#############################################
# Load Data
#############################################
input_data=read.csv("data.csv")
names(input_data)
n=nrow(input_data)

SID=input_data$SID
X=input_data$X
delta=input_data$delta
Z_0=input_data$Z_0
Z_6=input_data$Z_6
Z_12=input_data$Z_12
W=input_data$W

#Format data for GEE analysis (Figure 1 in manuscript)
Z1=c(Z_0,Z_6,Z_12)
Z2=c(W,W,W)
ID=c(SID,SID,SID)
wave=c(rep(1,n),rep(2,n),rep(3,n))
t=c(0,6,12)
Tau=12
  
##########################################################################################
# Step 2: Restricted Mean Model via Pseudo Observations Adjusted for Dependent Censoring
##########################################################################################
all_obs_results=get_IPCW_pseudo_obs(X,delta,t,Z_0,Z_6,Z_12,W)
beta_censor=all_obs_results$beta_censor
n_t=apply(array(X,c(n,length(t)))>=t(array(t,c(length(t),n))),2,sum)

pseudo_obs=array(NA,c(n,length(t)))
for(k in 1:n)
{
  remove_k_results=get_IPCW_pseudo_obs(X[-k],delta[-k],t,Z_0[-k],Z_6[-k],Z_12[-k],W[-k])
  pseudo_obs[k,]=n_t*all_obs_results$delta_hat-(n_t-1)*remove_k_results$delta_hat #Formula (2) in manuscript
  print(k)
}
Y=array(pseudo_obs,c(n*length(t),1))

#Format and fit GEE model to pseudo obs
data=data.frame(ID,Y,Z1,Z2,wave)
data.noNA=subset(data[order(ID),],is.na(Z1)==FALSE)
pseudo_obs_model=gee(Y~Z1+Z2,id=ID, corstr="unstructured",data=data.noNA)
beta_pseudo_obs=coef(pseudo_obs_model)
var_beta_pseudo_obs=diag(pseudo_obs_model$robust.variance)

beta_pseudo_obs
#(Intercept)          Z1          Z2 
# 1.73680853 -0.02466993  0.06389635 
lower_CI_beta=beta_pseudo_obs - 1.96*sqrt(var_beta_pseudo_obs)
#(Intercept)          Z1          Z2 
# 1.48682323 -0.04768082 -0.39432241 
upper_CI_beta=beta_pseudo_obs + 1.96*sqrt(var_beta_pseudo_obs)
# (Intercept)           Z1           Z2 
# 1.986793830 -0.001659039  0.522115097 

#Save values for use in MI
beta0_PO=beta_pseudo_obs[1]
beta1_PO=beta_pseudo_obs[2]
beta2_PO=beta_pseudo_obs[3]

#One window only
pseudo_obs_model_1window=lm(Y~Z1+Z2,data=data.noNA[data.noNA$wave==1,])
beta_pseudo_obs_1window=coef(pseudo_obs_model_1window)
#(Intercept)          Z1          Z2 
# 1.68720447 -0.01241842  0.14176693 
var_beta_pseudo_obs_1window=diag(vcov(pseudo_obs_model_1window))

lower_CI_beta=beta_pseudo_obs_1window - 1.96*sqrt(var_beta_pseudo_obs_1window)
#(Intercept)          Z1          Z2 
# 1.27152527 -0.06563672 -0.45639153
upper_CI_beta=beta_pseudo_obs_1window + 1.96*sqrt(var_beta_pseudo_obs_1window)
#(Intercept)          Z1          Z2 
# 2.10288366  0.04079989  0.73992539 

##########################################################################################
# Step 3: Multiple Imputation Method for dependently censored data
##########################################################################################
M=10
set.seed(1)

#Identify patients needing imputation
impute_subjects=SID[X<(Tau+t[length(t)]) & delta==0]
n_impute_subjects=length(impute_subjects)
imputed_data=array(apply(cbind(X,rep(Tau+t[length(t)],n)),1,min),c(n,M))
imputed_data[impute_subjects,]=NA

#Identify times when dN>0
time_X=sort(unique(X))
n_Z_1=sum(time_X<t[2])
n_Z_2=sum(time_X<t[3] & time_X>t[2])
n_Z_3=length(time_X)-n_Z_1-n_Z_2
#Construct matrices for time-dependent and time-independent variables at times when dN>0
Z_array=array(rep(Z_0,n_Z_1),c(n,n_Z_1))
Z_array=cbind(Z_array,array(rep(Z_6,n_Z_2),c(n,n_Z_2)))
Z_array=cbind(Z_array,array(rep(Z_12,n_Z_3),c(n,n_Z_3)))
W_array=array(W,c(n,n_Z_1+n_Z_2+n_Z_3))

#Get unconditional weight matrix
IPCW_surv_results=get_IPCW_surv(X,delta,t,Z_0,Z_6,Z_12,W,beta_censor)
weight=t(IPCW_surv_results$weight) #each row corresponds to a subject (K in step 3(b) of manuscript)

for(k in 1:n_impute_subjects)
{
  #Step 3(a): Identify risk set for patient k
  patient_id=impute_subjects[k]
  censoring_time_index=seq(1:length(time_X))[time_X==X[patient_id]]
  largest_possible_t=max(t[t<X[patient_id]]) #time of last window for patient k
  risk_set=rep(0,n)
  i=0
  epsilon=0.01
  while(sum(risk_set)<5)
  {
    epsilon=epsilon+i*0.001
    risk_set[X>X[patient_id] & abs(beta1_PO*Z_array[,censoring_time_index] + beta2_PO*W_array[,censoring_time_index] - beta1_PO*Z_array[patient_id,censoring_time_index] - beta2_PO*W_array[patient_id,censoring_time_index])<epsilon]=1
    i=i+1
  }
  
  X_riskset=X[risk_set==1]
  delta_riskset=delta[risk_set==1]
  W_riskset=W[risk_set==1]
  
  n_riskset=length(X_riskset)
  time_riskset=sort(unique(X_riskset))
  n_time_riskset=length(time_riskset)
  risk_set_index=array(seq(1:length(time_X)),c(length(time_X),n_time_riskset))[array(time_X,c(length(time_X),n_time_riskset))==t(array(time_riskset,c(n_time_riskset,length(time_X))))]
  
  #Step 3(b): Inverse Weighted Survival Estimation Within Risk Set
  weight_denominator=array(weight[risk_set==1,censoring_time_index],c(n_riskset,n_time_riskset))
  weight_riskset=weight[risk_set==1,risk_set_index]/weight_denominator #conditional weight matrix within risk set times (W in manuscript)
  
  #Estimate survival 
  temp3=array(X_riskset,c(n_riskset,n_time_riskset))
  temp4=t(array(time_riskset,c(n_time_riskset,n_riskset)))
  delta_array_riskset=array(delta_riskset,c(n_riskset,n_time_riskset))
  dN_T_riskset=array(as.numeric(temp3==temp4 & delta_array_riskset==1),c(n_riskset,n_time_riskset))
  Y_riskset=array(as.numeric(temp3>=temp4),c(n_riskset,n_time_riskset))
  
  denominator=apply(Y_riskset*weight_riskset,2,sum,na.rm=TRUE)
  temp5=apply(t(dN_T_riskset)*t(weight_riskset)/denominator,1,sum,na.rm=TRUE)
  IPCW_CH_Rk=sum_function(temp5)
  IPCW_surv_Rk=exp(-IPCW_CH_Rk)
  IPCW_time_Rk=time_riskset
  IPCW_surv_Rk[IPCW_time_Rk>(Tau+largest_possible_t)]=0
  #plot(IPCW_time_Rk,IPCW_surv_Rk,type="s",col="red")
  
  #Step 3(c): Sampling a valid impute
  valid_imputation=0
  while(valid_imputation==0)
  {
    U=runif(M)
    impute_index=(length(IPCW_time_Rk)-apply(array(IPCW_surv_Rk,c(length(IPCW_surv_Rk),M))<t(array(U,c(M,length(IPCW_surv_Rk)))),2,sum)) + 1
    imputation_times=IPCW_time_Rk[impute_index]
    
    temp=array(X_riskset,c(n_riskset,M))==t(array(IPCW_time_Rk[impute_index],c(M,n_riskset)))
    risk_set_imputes_subjects=array(SID[risk_set==1],c(length(SID[risk_set==1]),M))[temp]
    residuals=log(apply(cbind(X[risk_set_imputes_subjects]-largest_possible_t,rep(Tau,M)),1,min))-(beta0_PO + beta1_PO*Z_array[risk_set_imputes_subjects,censoring_time_index] + beta2_PO*W_array[risk_set_imputes_subjects,censoring_time_index])
    subject_imputes=exp(beta0_PO + beta1_PO*Z_array[patient_id,censoring_time_index] + beta2_PO*W_array[patient_id,censoring_time_index] + residuals) + largest_possible_t
    subject_imputes[imputation_times>(Tau+largest_possible_t)]=(Tau+largest_possible_t)
    if(sum(subject_imputes>X[patient_id])==M)
    {
      valid_imputation=1
    }
    print("imputing")
  }
  
  imputed_data[patient_id,]=subject_imputes
  print(k)
}

#Analysis of the M multiply imputed datasets
b1=array(NA,c(M,3))
var_b1=array(NA,c(M,3))
b1_1window=array(NA,c(M,3))
var_b1_1window=array(NA,c(M,3))
for(m in 1:M)
{
  #Transformed variables
  T_imputed=imputed_data[,m]
  T_t1=T_imputed-t[1]
  T_t1[T_imputed<=t[1]]=NA
  Y1=log(apply(cbind(T_t1,rep(Tau,n)),1,min))
  
  T_t2=T_imputed-t[2]
  T_t2[T_imputed<=t[2]]=NA
  Y2=log(apply(cbind(T_t2,rep(Tau,n)),1,min))
  
  T_t3=T_imputed-t[3]
  T_t3[T_imputed<=t[3]]=NA
  Y3=log(apply(cbind(T_t3,rep(Tau,n)),1,min))
  
  # Fit GEE model
  Y=c(Y1,Y2,Y3)
  data=data.frame(ID,Y,Z1,Z2,wave)
  data.noNA=subset(data[order(ID),],is.na(Y)==FALSE)
  model_MI=gee(Y~Z1+Z2,id=ID, corstr="unstructured",data=data.noNA)
  
  b1[m,]=model_MI$coefficients
  var_b1[m,]=diag(model_MI$robust.variance)
  
  # one window only
  model_MI_1window=lm(Y~Z1+Z2,data=data.noNA[data.noNA$wave==1,])
  b1_1window[m,]=model_MI_1window$coefficients
  var_b1_1window[m,]=diag(vcov(model_MI_1window))
}

beta_MI=apply(b1,2,mean)
# 2.086291240 -0.004914556  0.081668610
W=apply(var_b1,2,mean)
B=diag(cov(b1))
var_beta_MI=W + (1+1/M)*B

lower_CI=beta_MI-1.96*sqrt(var_beta_MI)
# 1.89283215 -0.02618138 -0.17967344
upper_CI=beta_MI+1.96*sqrt(var_beta_MI)
# 2.27975033 0.01635227 0.34301066

# one window only
beta_MI_1window=apply(b1_1window,2,mean)
# 2.161092575 -0.009589567  0.098010557
W=apply(var_b1_1window,2,mean)
B=diag(cov(b1_1window))
var_beta_MI_1window=W + (1+1/M)*B

lower_CI=beta_MI_1window-1.96*sqrt(var_beta_MI_1window)
# 1.92343563 -0.04006516 -0.24606131
upper_CI=beta_MI_1window+1.96*sqrt(var_beta_MI_1window)
# 2.39874952 0.02088603 0.44208242