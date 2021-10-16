#import data as obs, delta and covs
#requires ma_auxiliaryFunctions

create.imp.dr.l2 = function(obs,delta,covs,tau=NULL,quant=0.05,time.point){
  n = length(obs)
  
  parms = parameter_BrierTreeKM(obs=obs,delta=delta,covs=covs, tau=tau, quant=quant, time.point)
  #parms equals c(a0,a1,b0,b1,c0,c1,obs,delta,diag(m1),n)
  a1 <- parms[(n+1):(2*n)]
  b1 <- parms[(3*n+1):(4*n)]
  c1 <- parms[(5*n+1):(6*n)]
  y.imp.all <- a1 + b1 - c1
  
  return(y.imp.all)
}

# Create imputation for buckley James and l2 loss
create.imp.bj.l2 = function(obs,delta,covs,time.point){
  
  y.imp.all = parameter_BJ(obs, delta, covs, time.point)
  
  return(y.imp.all)
}

# Create imputation for dr and brier loss
create.imp.dr.brier = function(obs,delta,covs,tau=NULL,quant=0.05,time.point){
  n = length(obs)
  parms = parameter_BrierTreeKM(obs,delta,covs,tau,quant,time.point)
  a1 <- parms[(n+1):(2*n)]
  b1 <- parms[(3*n+1):(4*n)]
  c1 <- parms[(5*n+1):(6*n)]
  y.brier <- a1 + b1 - c1
  return(y.brier)
}

# Create imputation for dr and brier loss
create.imp.dr.brier.km = function(obs,delta,covs,tau=NULL, quant=0.05, time.point){
  n = length(obs)
  parms = parameter_BrierTreeKM(obs,delta,covs,tau,quant, time.point)
  a1 <- parms[(n+1):(2*n)]
  b1 <- parms[(3*n+1):(4*n)]
  c1 <- parms[(5*n+1):(6*n)]
  y.brier <- a1 + b1 - c1
  return(y.brier)
}

# Create imputation for dr and brier loss
create.imp.bj.brier = function(obs,delta,covs, time.point){
  y.brier = parameter_BJBrier(obs=obs,delta=delta,covs=covs, time.point=time.point)
  
  return(y.brier)
}

# Create imputation for dr and brier loss
create.ipcw.weights = function(obs,delta,covs,tau=NULL, quant=0.05, time.point){
  
  # Creating the new T(t) dataset
  data.tmp <- truncate_data(obs, delta, tau=time.point)
  obs.t = data.tmp$obs
  delta.t = data.tmp$delta
  n = length(obs)
  
  # Calculating the conditional censoring distribution.
  tem = SurvTreeEstimate_G(obs.t,delta.t,covs = covs, tau=tau,quant=quant)
  # Calculating the censoring distribution
  surC_rf = as.matrix(subset(tem, select=-(1:3)))#surC_rf=tem$cens.est
  
  # Calculating a0, a1, b0, b1, c0, c1
  ipcw.weights=delta.t/diag(surC_rf)
  data.ipcw = data.frame(obs.t, delta.t, ipcw.weights) 
  return(data.ipcw)
}

