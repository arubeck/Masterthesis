library(rpart)

#delta is event indicator (not censoring indicator), i.e. delta = 0 for censored outcomes
truncate_data <- function(obs, delta, tau = NULL, quant = 0.05){
  n = length(obs)
  #cut observations after specified time point tau, e.g. important for estimation of restricted mean survival
  if (!is.null(tau)){
    delta = delta * (obs <= tau) + (obs > tau)
    obs = pmin(obs, tau)
  }
  #dtype specifies how many % of data are truncated
  else {
    delta[order(obs)][floor(n*(1-quant)):n] = FALSE
    obs[order(obs)][floor(n*(1-quant)):n] = obs[order(obs)][floor(n*(1-quant))]
  }
  data_out <- data.frame(
    obs = obs,
    delta = delta
  )
  return(data_out)
}

#functions that estimates conditional censoring distribution
#first approach: Kaplan-Meier estimator
KMestimate_G <- function(obs, delta, tau, quant){
  #number of observations
  n = length(obs)
  #truncate data, if no truncation is desired, set quant=0
  trunc = truncate_data(obs, delta, tau, quant)
  obs = trunc$obs
  delta = trunc$delta
  
  #calculate KM estimator
  hazC=mapply(function(xx,dd){dd/sum((obs>=xx))},xx=obs,dd=delta)
  surC_km=mapply(function(xx){prod(1-hazC[obs<=xx])},xx=obs)
  est_G = data.frame(
    surC_km = surC_km,
    obs = obs,
    delta = delta
  )
  return(est_G)
}
  
SurvTreeEstimate_G <- function(obs, delta, covs, tau, quant){
  #number of observations
  n = length(obs)
  #change perspective: for G we need indicator to indicate censoring, not the event
  data.used <- data.frame(obs, 1 - delta, covs)
  names(data.used)[1:2] <- c("obs", "delta.cens")
  #Fit regression tree
  surv.tree = rpart(Surv(obs,delta.cens)~., data = data.used, minbucket = 30)
  #survival curve
  pred.surv.tree <- predict(surv.tree, proximity = FALSE)# Finding the terminal nodes
  sett=unique(surv.tree[['where']])
  nset=length(sett)
  
  cens.est = matrix(0, nrow = n, ncol = n)
  obs.used = rep(NA, n)
  delta.used = rep(NA, n)
  
  for (i in 1:nset){
    # Finding the subset corresponding the ith node of the tree
    subset=(1:n)[surv.tree[['where']]==sett[i]]
    nlen=length(subset)
    # Observed time within the node
    sobs=obs[subset]
    # Failure indicators within each node.
    sdelta=delta[subset]
    
    # Doing truncation within that subset
    # Changing the dataset to account for truncation. 
    aa=truncate_data(sobs,1-sdelta,tau = tau, quant = quant)
    # Observed time after truncation
    sobs=aa[[1]]
    # Failure indicator after truncation.
    sdelta = aa[[2]]
    
    obs.used[subset] = sobs
    delta.used[subset] = sdelta
    
    # Calculating the KM estimator censoring curve within a node
    # Calculating the jumps in the KM estimator
    hazC=mapply(function(xx,dd){dd/sum((sobs>=xx))},xx=sobs,dd=sdelta)
    surC_km=mapply(function(xx){prod(1-hazC[sobs<=xx])},xx=sobs)
    cens.est[subset, ] = matrix(surC_km,nrow=length(subset),ncol=length(surC_km),byrow=TRUE)
  }
  SurvTree.estimate = data.frame(
    cens.est = cens.est,
    obs.used = obs.used,
    delta.used = delta.used,
    term.nodes = surv.tree[['where']]
  )
  return(SurvTree.estimate)
}

#estimate conditional expectations using random forests
RandomForestEstimate_m <- function(obs, delta, covs){
  #number of observations
  n = length(obs)
  
  # Creating the data frame
  data.used <- data.frame(obs, delta, covs)
  names(data.used)[1:2] <- c("obs", "delta")
  
  #Fitting a randm forest
  rand.for = rfsrc(Surv(obs,delta)~ ., data = data.used, importance = "none")
  # Getting the Survival Curves.
  pred.rf <- predict(rand.for, proximity = FALSE)
  
  # Calculation of the survivor function
  m1=matrix(0,n, n)
  # Finding unique event times
  time.used <- pred.rf$time.interest
  # Finding the jumps in the estimator \hat P(T >t|W_i)
  surv.diff <- matrix(0, ncol = sum(delta), nrow = n)
  
  for(i in 1:n){
    # Calculating the jumps in the random forest model survival curve estimator
    surv.diff[i, ] <- c(1, pred.rf$survival[i, ][-length(pred.rf$survival[i, ])]) - pred.rf$survival[i, ]
  }
  for(j in 1:n){
    if(delta[j]==FALSE){
      for(i in 1:n){
        if(obs[j]<=obs[i]){
          if(sum(surv.diff[i, ][time.used > obs[j]]) != 0){
            # Calculating the conditional expectation
            m1[j,i]=  sum(log(time.used[time.used > obs[j]]) * surv.diff[i, ][time.used > obs[j]])/sum(surv.diff[i, ][time.used > obs[j]])
          }
        }
      }
      if (sum(surv.diff[i, ][time.used > obs[j]]) == 0){
        m1[j,]=log(obs[j])
      }
    }
  }
  return(m1)
}

BrierRandomForest_m <- function(obs, delta, covs, time.point){
  n = length(obs)
  
  # Creating the data frame
  data.used <- data.frame(obs, delta, covs)
  names(data.used)[1:2] <- c("obs", "delta")
  # Fitting the Tree
  rand.for = rfsrc(Surv(obs,delta)~ ., data = data.used, importance = "none")
  
  # Getting the Survival Curves. 
  pred.rf <- predict(rand.for, proximity = FALSE)
  
  # Calculating the cox model
  #m1[j, i] <- P(T >t|T > T_j, W_i)
  #P(T > tau|T >u,W) = P(T>tau|W)/P(T>u|W)
  
  # predsurvAFT compute P(T>t|W) tval is t, fit is the AFT fit and zval is the covariate values.
  predsurvRF = function(time, cov.index){
    time.point = sum(pred.rf$time.interest < time) + 1
    surv = c(1, pred.rf$survival[cov.index, ])[time.point]
    return(surv)
  }
  
  # Defininf the matrix where m1[j, i] = P(T>t|T >T_j, W_i)
  m1 <- matrix(0, ncol = n, nrow = n)
  
  # Calculating the cox model
  #m1[j, i] <- P(T >t|T > T_j, W_i)
  for(j in 1:n){
    if(delta[j] == 0 & obs[j] < time.point){
      for(i in 1:n){
        # Calculate the denominator and the numerator in the desired probability
        prob.est.den = predsurvRF(obs[j], i)
        prob.est.num = predsurvRF(time.point, i)
        if(prob.est.den != 0){
          m1[j, i] = prob.est.num/prob.est.den
        }
        if(prob.est.den == 0){
          m1[j, i] = 0.5
        }
      }
    }
  }
  
  # Return the probabilities
  return(m1)
}

#calculation of parameters needed for trafo
parameter_RegularTreeKM <- function(obs, delta, covs, tau, quant){
  n = length(obs)
  
  # Calculating the conditional expectation    
  m1 = RandomForestEstimate_m(obs,delta,covs)
  
  # Calculating the conditional censoring distribution.
  tem=SurvTreeEstimate_G(obs,delta,covs = covs, tau = tau, quant = quant)
  # Calculating the censoring distribution
  surC_rf=tem$cens.est
  # Observed event times for adjusted for truncation
  obs=tem$obs.used
  # Failure indicator adjusted for truncation
  delta=tem$delta.used
  # Finding which observations fall in which terminal node
  term.nodes = tem$term.nodes
  
  # Calculating a0, a1, b0, b1, c0, c1
  a0=delta/diag(surC_rf)
  a1=a0*log(obs)
  
  b0=(1-delta)/diag(surC_rf)
  b1=b0 * diag(m1)
  
  c0 <- rep(NA, n)
  c1 <- rep(NA, n)
  
  # Creating the ordered data
  ord.used = order(obs)
  obs.order = obs[ord.used]
  delta.order = delta[ord.used]
  
  # Finding the terminal nodes
  sett=unique(term.nodes)
  nset=length(sett)
  
  for (i in 1:nset){
    # Finding the subset corresponding the ith node of the tree
    subset=(1:n)[term.nodes==sett[i]]
    nlen=length(subset)
    # Observed time within the node
    sobs=obs[subset]
    # Failure indicators within each node.
    sdelta=delta[subset]
    
    kk=mapply(function(tt){sum((tt<=sobs))},tt=sobs)
    c0[subset]=mapply(function(tt){sum(b0[subset]*(sobs<=tt)/kk)},tt=sobs)
    c1[subset]=mapply(function(tt,i){sum(b0[subset]*(sobs<=tt)*m1[subset,i]/kk)},tt=sobs,i=1:nlen)
  }
  parms = c(a0,a1,b0,b1,c0,c1,obs,delta,diag(m1),n)
  return(parms)
}
  
parameter_BrierTreeKM <- function(obs, delta, covs, tau, quant, time.point){
  n = length(obs)
  
  # Creating the new T(t) dataset
  obs.t = pmin(obs, time.point)
  delta.t = delta * (obs <= time.point) + (obs > time.point)
  
  # Calculating the conditional expectation    
  m1=BrierRandomForest_m(obs, delta, covs, time.point)
  
  # Calculating the conditional censoring distribution.
  tem = SurvTreeEstimate_G(obs.t,delta.t,tau = tau, quant = quant, covs = covs)
  # Calculating the censoring distribution
  surC_rf=tem$cens.est
  # Finding which observations fall in which terminal node
  term.nodes = tem$term.nodes
  
  # Calculating a0, a1, b0, b1, c0, c1
  a0=delta.t/diag(surC_rf)
  a1 = a0 * (obs > time.point)
  
  b0=(1-delta.t)/diag(surC_rf)
  b1=b0 * diag(m1)
  
  c0 <- rep(NA, n)
  c1 <- rep(NA, n)
  
  # Finding the terminal nodes
  sett=unique(term.nodes)
  nset=length(sett)
  
  for (i in 1:nset){
    # Finding the subset corresponding the ith node of the tree
    subset=(1:n)[term.nodes==sett[i]]
    nlen=length(subset)
    # Observed time within the node
    sobs=obs.t[subset]
    # Failure indicators within each node.
    sdelta=delta.t[subset]
    
    kk=mapply(function(tt){sum((tt<=sobs))},tt=sobs)
    c0[subset]=mapply(function(tt){sum(b0[subset]*(sobs<=tt)/kk)},tt=sobs)
    c1[subset]=mapply(function(tt,i){sum(b0[subset]*(sobs<=tt)*m1[subset,i]/kk)},tt=sobs,i=1:nlen)
  }
  
  parms = c(a0,a1,b0,b1,c0,c1,obs,delta,diag(m1),n)
  return(parms)
}

parameter_BJBrier <- function(obs, delta, covs, tau, quant, time.point){
  obs.t = pmin(obs, time.point)
  delta.t = delta * (obs <= time.point) + (obs > time.point)
  
  # Calculating the conditional expectation    
  m1 = BrierRandomForest_m(obs,delta,covs, time.point)
  
  a1 = delta.t *  (obs > time.point) + (1 - delta.t) * diag(m1)

  return(a1)
}

parameter_BJ <- function(obs, delta, covs, tau, quant, time.point){
  # Calculating the conditional expectation    
  m1 = RandomForestEstimate_m(obs,delta,covs)
  a1 = delta *  log(obs) + (1 - delta) * diag(m1)
  
  return(a1)
}
  
