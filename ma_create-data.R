#Setting 1
##For censoring rate approx 0.31 choose cens_rate = 0.5
create_data_setting1 <- function(nobs, ncov, cens_rate=0.5){
  #Initialize Variables
  covariates <- matrix(data = 0,nrow = nobs, ncol = ncov)
  surv_time <- rep(0, nobs)
  cens_time <- rep(0, nobs)
  obs_time <- rep(0,nobs)
  delta <- rep(0, nobs)
  
  for (i in 1:nobs){
    #Calculate Covariance Matrix
    covariance_matrix <- matrix(data=NA, nrow=ncov, ncol=ncov)
    for (j in 1:ncov){
      for (k in 1:ncov){
        covariance_matrix[j,k] = 0.9^(abs(j-k))
      }
    }
    #Calculate Covariates from multivariate normal distribution with mean 0
    covariates[i,] <- mvrnorm(n=1, mu=rep(0,ncov), Sigma=covariance_matrix)
    
    #Calculate Survival Times from exponential distribution
    ##Calculate Mean
    mu_surv = exp(0.1*sum(covariates[i,c(9:18)]))
    surv_time[i] = rexp(1, 1/mu_surv)
    
    #Calculate Censoring Destribution
    cens_time[i] = rexp(n=1, rate=cens_rate)
    
    #Update observed time and indicator
    delta[i] = (surv_time[i]<=cens_time[i])
    obs_time[i] = min(surv_time[i], cens_time[i])
  }
  
  #Store simulated data as data.frame
  obs_tmp <- data.frame(
    obs = obs_time,
    delta = delta
  )
  covs_tmp <- data.frame(covariates)
  obs_data <- cbind(obs_tmp, covs_tmp)
  return(obs_data)
}

#Setting 2
create_data_setting2 <- function(nobs, ncov){
  #Initialize Variables
  covariates <- matrix(data = 0,nrow = nobs, ncol = ncov)
  surv_time <- rep(0, nobs)
  cens_time <- rep(0, nobs)
  obs_time <- rep(0,nobs)
  delta <- rep(0, nobs)
  
  for (i in 1:nobs){
    #Calculate Covariates from uniform random distribution on the interval [0,1]
    covariates[i,] = runif(ncov)
    
    #Calculate Survival Times from exponential distribution
    ##Calculate Mean
    mu_surv = sin(covariates[i,1]*pi)+2*abs(covariates[i,2]-0.5)+covariates[i,3]^3
    surv_time[i] = rexp(1, 1/mu_surv)
    
    #Calculate Censoring Destribution
    cens_time[i] = runif(1, min = 0, max = 6)
    
    #Update observed time and indicator
    delta[i] = (surv_time[i]<=cens_time[i])
    obs_time[i] = min(surv_time[i], cens_time[i])
  }
  
  #Store simulated data as data.frame
  obs_tmp <- data.frame(
    obs = obs_time,
    delta = delta
  )
  covs_tmp <- data.frame(covariates)
  obs_data <- cbind(obs_tmp, covs_tmp)
  return(obs_data)
}

#Setting 3
create_data_setting3 <- function(nobs, ncov){
  #Initialize Variables
  covariates <- matrix(data = 0,nrow = nobs, ncol = ncov)
  surv_time <- rep(0, nobs)
  cens_time <- rep(0, nobs)
  obs_time <- rep(0,nobs)
  delta <- rep(0, nobs)
  
  for (i in 1:nobs){
    #Calculate Covariance Matrix
    covariance_matrix <- matrix(data=NA, nrow=ncov, ncol=ncov)
    for (j in 1:ncov){
      for (k in 1:ncov){
        covariance_matrix[j,k] = 0.75^(abs(j-k))
      }
    }
    #Calculate Covariates from multivariate normal distribution with mean 0
    covariates[i,] <- mvrnorm(n=1, mu=rep(0,ncov), Sigma=covariance_matrix)
    
    #Calculate Survival Times from gamma distribution with scale parameter = 2
    ##Calculate Shape Parameter
    mu_surv = 0.5 + 0.3*abs(sum(covariates[c(11:15)]))
    surv_time[i] = rgamma(1, shape = 2/mu_surv, scale = 2)
    
    #Calculate Censoring Destribution
    cens_time[i] = runif(1, min=0, max=15)
    
    #Update observed time and indicator
    delta[i] = (surv_time[i]<=cens_time[i])
    obs_time[i] = min(surv_time[i], cens_time[i])
  }
  
  #Store simulated data as data.frame
  obs_tmp <- data.frame(
    obs = obs_time,
    delta = delta
  )
  covs_tmp <- data.frame(covariates)
  obs_data <- cbind(obs_tmp, covs_tmp)
  return(obs_data)
}

#Setting 4
create_data_setting4 <- function(nobs, ncov){
  #Initialize Variables
  covariates <- matrix(data = 0,nrow = nobs, ncol = ncov)
  surv_time <- rep(0, nobs)
  cens_time <- rep(0, nobs)
  obs_time <- rep(0,nobs)
  delta <- rep(0, nobs)
  
  for (i in 1:nobs){
    #Calculate Covariance Matrix
    covariance_matrix <- matrix(data=NA, nrow=ncov, ncol=ncov)
    for (j in 1:ncov){
      for (k in 1:ncov){
        covariance_matrix[j,k] = 0.75^(abs(j-k))
      }
    }
    #Calculate Covariates from multivariate normal distribution with mean 0
    covariates[i,] <- mvrnorm(n=1, mu=rep(0,ncov), Sigma=covariance_matrix)
    
    #Calculate Survival Times from log-normal distribution
    ##Calculate Mean
    mu_surv = 0.1*abs(sum(covariates[c(1:5)]))+0.1*abs(sum(covariates[c(16:20)]))
    surv_time[i] = rlnorm(1, meanlog=mu_surv)
    
    #Calculate Censoring Times from log normal distribution, scale parameter =1
    mu_cens = mu_surv+0.5
    cens_time[i] = rlnorm(1, meanlog = mu_cens)
    
    #Update observed time and indicator
    delta[i] = (surv_time[i]<=cens_time[i])
    obs_time[i] = min(surv_time[i], cens_time[i])
  }
  
  #Store simulated data as data.frame
  obs_tmp <- data.frame(
    obs = obs_time,
    delta = delta
  )
  covs_tmp <- data.frame(covariates)
  obs_data <- cbind(obs_tmp, covs_tmp)
  return(obs_data)
}


#data.used: dataset with first column named obs, second column named delta, rest covariates
#v gives number of cross validations
#time.point specifies parameter of interest for Brier loss
create.crossval.data = function(data.used, v, time.point){
  #change names as required
  names(data.used)[c(1,2)] = c("obs","delta")
  ## Creating the cross validated group. 
  get.sam <- rep(c(1:v),sum(data.used$delta==TRUE)/v)
  #rep rounds down, add "incomplete" enumeration
  if(length(get.sam)<sum(data.used$delta==TRUE)){
    get.sam <- c(get.sam,c(1:(sum(data.used$delta==TRUE)-length(get.sam))))
  }
  #draw sample of get.sam (vector of indices), with the size of #of events, without replacement
  delta.1 <- sample(get.sam,sum(data.used$delta==TRUE),replace=F)
  #draw sample of get.sam (with adjusted length) with the size of #of censorings, without replacement
  delta.0 <- sample(rep(c(1:v),length(data.used$delta)-length(delta.1)),sum(data.used$delta==FALSE),replace=F) 
  k <- 1
  l <- 1
  val.sample <- NULL
  
  for(m in 1:length(data.used$delta)){
    #if data is censored
    if(data.used$delta[m]==0){
      #add kth element of delta.0 (which contains 1:5 number-of-censored-observation times)
      val.sample[m] <- delta.0[k]
      k <- k+1
    }
    #if data is not censored
    if(data.used$delta[m]==1){
      #take element of delta.1, which contains replicats of 1:5 number-of-observed-events times
      val.sample[m] <- delta.1[l]
      l <- l+1
    }
  }
  #perform crossvalidation
  for(j in 1:v){
    #draw nobs samples from 1:v, use each index with val.sample != 1 as trainings data, rest as test data
    val.sample = sample(1:v, dim(data.used)[1], replace = TRUE)
    data.train = data.used[val.sample != 1, ]
    obs.train = data.train$obs
    delta.train = data.train$delta
    cov.train = data.train[, -c(1,2)]
    data.test = data.used[val.sample == 1, ]
    obs.test = data.test$obs
    delta.test = data.test$delta
    cov.test = data.test[, -c(1,2)]
    
    # Create imputation
    imp.val.dr = create.imp.dr.brier.km(obs = obs.train, delta = delta.train, covs = cov.train, quant=0.1, time.point = time.point)
    imp.val.bj = create.imp.bj.brier(obs = obs.train, delta = delta.train, covs = cov.train, time.point = time.point)
    data.ipcw = create.ipcw.weights(obs = data.used$time, delta = data.used$event, covs = data.used[, -c(1,2)], quant=0.1, time.point = time.point)
    data.ipcw = data.ipcw[val.sample == 1, 3]
    
    x_train = cov.train
    y_train.dr = imp.val.dr
    y_train.bj = imp.val.bj
    x_test = cov.test
    y_test = obs.test
    
    data.surv = data.frame(data.used$obs, data.used$delta, data.used[, -c(1,2)])[val.sample != 1, ]
    names(data.surv)[1:2] = c("obs", "delta")
  }
  tmp <- data.frame(
    x_train = x_train,
    y_train.dr = y_train.dr,
    y_train.bj = y_train.bj,
    x_test = x_test,
    y_test = y_test,
    data.durv = data.surv,
    data.ipcw = data.ipcw,
    time.point = time.point
  )
  return(tmp)
}


