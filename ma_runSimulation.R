library(randomForestSRC)
library(rpart)
library(survival)
library(penalized)
library(ranger)
library(survivalmodels)
library(rockchalk)
library(glmnet)
library(keras)
library(dplyr)
library(MASS)


#run code with simulated data
setwd("C:/Users/tobia/Desktop/Masterthesis-main/Masterthesis-main")
source("ma_auxillary-functions.R")
source("ma_cox-procedure.R")
source("ma_create-data.R")
source("ma_deeplearning.R")
source("ma_trafo.R")

path="C:/Users/tobia/Desktop/Masterthesis-main/Masterthesis-main/data/"

nobs = 300
ncov = 25

numb = 50

######Setting1
comp_cox_pen=rep(0,numb)
comp_cox=rep(0,numb)
comp_rand=rep(0,numb)

comp_bj_cox_pen=rep(0,numb)
comp_bj_cox=rep(0,numb)
comp_bj_rand=rep(0,numb)

comp_dr_cox_pen=rep(0,numb)
comp_dr_cox=rep(0,numb)
comp_dr_rand=rep(0,numb)

comp_dl_dr=rep(0,numb)
comp_dl_dr_crossval=rep(0,numb)
comp_dl_bj=rep(0,numb)
comp_dl_bj_crossval=rep(0,numb)

cens_rate = rep(0,numb)

for (i in 1:numb){
  tmp = create_data_setting1(nobs,ncov)
  obs.time = tmp$obs
  obs.delta = tmp$delta
  obs.covs = subset(tmp, select = -c(obs, delta))
  colnames(obs.covs) <- paste0("W", colnames(obs.covs), "_")
  
  cens_rate[i] = length(which(obs.delta==0))/length(obs.time)
  
  obs.data <- data.frame(
    obs = obs.time,
    delta = obs.delta,
    covs = obs.covs
  )
  
  time.point = 0.4189946
  tmp <- create.crossval.data(obs.data, v=5, time.point=time.point)
  
  x_train <- tmp$x_train
  y_train <- tmp$y_train
  y_train.dr <- tmp$y_train.dr
  y_train.bj <- tmp$y_train.bj
  x_test <- tmp$x_test
  y_test <- tmp$y_test
  data.surv <- tmp$data.surv
  ipcw_weights <- tmp$data.ipcw
  time.point <- tmp$time.point
  
  comp <- compare.methods.brier(x_train = x_train, y_train = y_train, x_test = x_test, y_test = y_test, data_surv = data.surv, ipcw_weights = ipcw_weights, time.point = time.point)
  comp_cox_pen[i] = comp$pred.err.cox.pen
  comp_cox[i] = comp$pred.err.cox
  comp_rand[i] = comp$pred.err.rand
  
  comp_bj <- compare.methods.brier(x_train = x_train, y_train = y_train.bj, x_test = x_test, y_test = y_test, data_surv = data.surv, ipcw_weights = ipcw_weights, time.point = time.point)
  comp_bj_cox_pen[i] = comp_bj$pred.err.cox.pen
  comp_bj_cox[i] = comp_bj$pred.err.cox
  comp_bj_rand[i] = comp_bj$pred.err.rand
  
  comp_dr <- compare.methods.brier(x_train = x_train, y_train = y_train.dr, x_test = x_test, y_test = y_test, data_surv = data.surv, ipcw_weights = ipcw_weights, time.point = time.point)
  comp_dr_cox_pen[i] = comp_bj$pred.err.cox.pen
  comp_dr_cox[i] = comp_bj$pred.err.cox
  comp_dr_rand[i] = comp_bj$pred.err.rand
  
  comp_1 <- one.sim.dl.dr(x_train = x_train, y_train.dr = y_train.dr, x_test = x_test, y_test = y_test, data_surv = data.surv, ipcw_weights = ipcw_weights, time.point = time.point)
  comp_dl_dr[i] = comp_1$pred.err.dl
  comp_dl_dr_crossval[i] = comp_1$pred.err.dl.cv
  
  comp_2 <- one.sim.dl.bj(x_train = x_train, y_train.bj = y_train.bj, x_test = x_test, y_test = y_test, data_surv = data.surv, ipcw_weights = ipcw_weights, time.point = time.point)
  comp_dl_bj[i] = comp_2$pred.err.dl
  comp_dl_bj_crossval[i] = comp_2$pred.err.dl.cv
  
  #store results
  df = bind_rows(data.frame(x_train=x_train, y_train.bj = y_train.bj, y_train.dr=y_train.dr),
                 data.frame(x_test=x_test,y_test=y_test), data.frame(data.surv=data.surv),
                 data.frame(ipcw_weights = ipcw_weights),data.frame(time.point=time.point),
                 data.frame(comp = comp, comp_bj=comp_bj, comp_dr=comp_dr, comp_dl_bj=comp_1, comp_dl_dr=comp_2))
  data_name = paste(path,"setting1-iter",i,"data",".csv",sep="")
  write.csv(df, data_name)
}

save.image(file="setting1.RData")

#########Setting2
comp_cox_pen=rep(0,numb)
comp_cox=rep(0,numb)
comp_rand=rep(0,numb)

comp_bj_cox_pen=rep(0,numb)
comp_bj_cox=rep(0,numb)
comp_bj_rand=rep(0,numb)

comp_dr_cox_pen=rep(0,numb)
comp_dr_cox=rep(0,numb)
comp_dr_rand=rep(0,numb)

comp_dl_dr=rep(0,numb)
comp_dl_dr_crossval=rep(0,numb)
comp_dl_bj=rep(0,numb)
comp_dl_bj_crossval=rep(0,numb)

cens_rate = rep(0,numb)


for (i in 1:numb){
  tmp = create_data_setting2(nobs,ncov)
  obs.time = tmp$obs
  obs.delta = tmp$delta
  obs.covs = subset(tmp, select = -c(obs, delta))
  colnames(obs.covs) <- paste0("W", colnames(obs.covs), "_")
  
  cens_rate[i] = length(which(obs.delta==0))/length(obs.time)
  
  obs.data <- data.frame(
    obs = obs.time,
    delta = obs.delta,
    covs = obs.covs
  )
  
  time.point = 0.7046207
  tmp <- create.crossval.data(obs.data, v=5, time.point=time.point)
  
  x_train <- tmp$x_train
  y_train <- tmp$y_train
  y_train.dr <- tmp$y_train.dr
  y_train.bj <- tmp$y_train.bj
  x_test <- tmp$x_test
  y_test <- tmp$y_test
  data.surv <- tmp$data.surv
  ipcw_weights <- tmp$data.ipcw
  time.point <- tmp$time.point
  
  comp <- compare.methods.brier(x_train = x_train, y_train = y_train, x_test = x_test, y_test = y_test, data_surv = data.surv, ipcw_weights = ipcw_weights, time.point = time.point)
  comp_cox_pen[i] = comp$pred.err.cox.pen
  comp_cox[i] = comp$pred.err.cox
  comp_rand[i] = comp$pred.err.rand
  
  comp_bj <- compare.methods.brier(x_train = x_train, y_train = y_train.bj, x_test = x_test, y_test = y_test, data_surv = data.surv, ipcw_weights = ipcw_weights, time.point = time.point)
  comp_bj_cox_pen[i] = comp_bj$pred.err.cox.pen
  comp_bj_cox[i] = comp_bj$pred.err.cox
  comp_bj_rand[i] = comp_bj$pred.err.rand
  
  comp_dr <- compare.methods.brier(x_train = x_train, y_train = y_train.dr, x_test = x_test, y_test = y_test, data_surv = data.surv, ipcw_weights = ipcw_weights, time.point = time.point)
  comp_dr_cox_pen[i] = comp_bj$pred.err.cox.pen
  comp_dr_cox[i] = comp_bj$pred.err.cox
  comp_dr_rand[i] = comp_bj$pred.err.rand
  
  comp_1 <- one.sim.dl.dr(x_train = x_train, y_train.dr = y_train.dr, x_test = x_test, y_test = y_test, data_surv = data.surv, ipcw_weights = ipcw_weights, time.point = time.point)
  comp_dl_dr[i] = comp_1$pred.err.dl
  comp_dl_dr_crossval[i] = comp_1$pred.err.dl.cv
  
  comp_2 <- one.sim.dl.bj(x_train = x_train, y_train.bj = y_train.bj, x_test = x_test, y_test = y_test, data_surv = data.surv, ipcw_weights = ipcw_weights, time.point = time.point)
  comp_dl_bj[i] = comp_2$pred.err.dl
  comp_dl_bj_crossval[i] = comp_2$pred.err.dl.cv
  
  #store results
  df = bind_rows(data.frame(x_train=x_train, y_train.bj = y_train.bj, y_train.dr=y_train.dr),
                 data.frame(x_test=x_test,y_test=y_test), data.frame(data.surv=data.surv),
                 data.frame(ipcw_weights = ipcw_weights),data.frame(time.point=time.point),
                 data.frame(comp = comp, comp_bj=comp_bj, comp_dr=comp_dr, comp_dl_bj=comp_1, comp_dl_dr=comp_2))
  data_name = paste(path,"setting2-iter",i,"data",".csv",sep="")
  write.csv(df, data_name)
}

save.image(file="setting2.RData")

#########Setting3
comp_cox_pen=rep(0,numb)
comp_cox=rep(0,numb)
comp_rand=rep(0,numb)

comp_bj_cox_pen=rep(0,numb)
comp_bj_cox=rep(0,numb)
comp_bj_rand=rep(0,numb)

comp_dr_cox_pen=rep(0,numb)
comp_dr_cox=rep(0,numb)
comp_dr_rand=rep(0,numb)

comp_dl_dr=rep(0,numb)
comp_dl_dr_crossval=rep(0,numb)
comp_dl_bj=rep(0,numb)
comp_dl_bj_crossval=rep(0,numb)

cens_rate = rep(0,numb)


for (i in 1:numb){
  tmp = create_data_setting3(nobs,ncov)
  obs.time = tmp$obs
  obs.delta = tmp$delta
  obs.covs = subset(tmp, select = -c(obs, delta))
  colnames(obs.covs) <- paste0("W", colnames(obs.covs), "_")
  
  cens_rate[i] = length(which(obs.delta==0))/length(obs.time)
  
  obs.data <- data.frame(
    obs = obs.time,
    delta = obs.delta,
    covs = obs.covs
  )
  
  time.point = 2.90718
  tmp <- create.crossval.data(obs.data, v=5, time.point=time.point)
  
  x_train <- tmp$x_train
  y_train <- tmp$y_train
  y_train.dr <- tmp$y_train.dr
  y_train.bj <- tmp$y_train.bj
  x_test <- tmp$x_test
  y_test <- tmp$y_test
  data.surv <- tmp$data.surv
  ipcw_weights <- tmp$data.ipcw
  time.point <- tmp$time.point
  
  comp <- compare.methods.brier(x_train = x_train, y_train = y_train, x_test = x_test, y_test = y_test, data_surv = data.surv, ipcw_weights = ipcw_weights, time.point = time.point)
  comp_cox_pen[i] = comp$pred.err.cox.pen
  comp_cox[i] = comp$pred.err.cox
  comp_rand[i] = comp$pred.err.rand
  
  comp_bj <- compare.methods.brier(x_train = x_train, y_train = y_train.bj, x_test = x_test, y_test = y_test, data_surv = data.surv, ipcw_weights = ipcw_weights, time.point = time.point)
  comp_bj_cox_pen[i] = comp_bj$pred.err.cox.pen
  comp_bj_cox[i] = comp_bj$pred.err.cox
  comp_bj_rand[i] = comp_bj$pred.err.rand
  
  comp_dr <- compare.methods.brier(x_train = x_train, y_train = y_train.dr, x_test = x_test, y_test = y_test, data_surv = data.surv, ipcw_weights = ipcw_weights, time.point = time.point)
  comp_dr_cox_pen[i] = comp_bj$pred.err.cox.pen
  comp_dr_cox[i] = comp_bj$pred.err.cox
  comp_dr_rand[i] = comp_bj$pred.err.rand
  
  comp_1 <- one.sim.dl.dr(x_train = x_train, y_train.dr = y_train.dr, x_test = x_test, y_test = y_test, data_surv = data.surv, ipcw_weights = ipcw_weights, time.point = time.point)
  comp_dl_dr[i] = comp_1$pred.err.dl
  comp_dl_dr_crossval[i] = comp_1$pred.err.dl.cv
  
  comp_2 <- one.sim.dl.bj(x_train = x_train, y_train.bj = y_train.bj, x_test = x_test, y_test = y_test, data_surv = data.surv, ipcw_weights = ipcw_weights, time.point = time.point)
  comp_dl_bj[i] = comp_2$pred.err.dl
  comp_dl_bj_crossval[i] = comp_2$pred.err.dl.cv
  
  #store results
  df = bind_rows(data.frame(x_train=x_train, y_train.bj = y_train.bj, y_train.dr=y_train.dr),
                 data.frame(x_test=x_test,y_test=y_test), data.frame(data.surv=data.surv),
                 data.frame(ipcw_weights = ipcw_weights),data.frame(time.point=time.point),
                 data.frame(comp = comp, comp_bj=comp_bj, comp_dr=comp_dr, comp_dl_bj=comp_1, comp_dl_dr=comp_2))
  data_name = paste(path,"setting3-iter",i,"data",".csv",sep="")
  write.csv(df, data_name)
}

save.image(file="setting3.RData")

#########Setting4
comp_cox_pen=rep(0,numb)
comp_cox=rep(0,numb)
comp_rand=rep(0,numb)

comp_bj_cox_pen=rep(0,numb)
comp_bj_cox=rep(0,numb)
comp_bj_rand=rep(0,numb)

comp_dr_cox_pen=rep(0,numb)
comp_dr_cox=rep(0,numb)
comp_dr_rand=rep(0,numb)

comp_dl_dr=rep(0,numb)
comp_dl_dr_crossval=rep(0,numb)
comp_dl_bj=rep(0,numb)
comp_dl_bj_crossval=rep(0,numb)

cens_rate = rep(0,numb)


for (i in 1:numb){
  tmp = create_data_setting4(nobs,ncov)
  obs.time = tmp$obs
  obs.delta = tmp$delta
  obs.covs = subset(tmp, select = -c(obs, delta))
  colnames(obs.covs) <- paste0("W", colnames(obs.covs), "_")
  
  cens_rate[i] = length(which(obs.delta==0))/length(obs.time)
  
  obs.data <- data.frame(
    obs = obs.time,
    delta = obs.delta,
    covs = obs.covs
  )
  
  time.point = 1.034881
  tmp <- create.crossval.data(obs.data, v=5, time.point=time.point)
  
  x_train <- tmp$x_train
  y_train <- tmp$y_train
  y_train.dr <- tmp$y_train.dr
  y_train.bj <- tmp$y_train.bj
  x_test <- tmp$x_test
  y_test <- tmp$y_test
  data.surv <- tmp$data.surv
  ipcw_weights <- tmp$data.ipcw
  time.point <- tmp$time.point
  
  comp <- compare.methods.brier(x_train = x_train, y_train = y_train, x_test = x_test, y_test = y_test, data_surv = data.surv, ipcw_weights = ipcw_weights, time.point = time.point)
  comp_cox_pen[i] = comp$pred.err.cox.pen
  comp_cox[i] = comp$pred.err.cox
  comp_rand[i] = comp$pred.err.rand
  
  comp_bj <- compare.methods.brier(x_train = x_train, y_train = y_train.bj, x_test = x_test, y_test = y_test, data_surv = data.surv, ipcw_weights = ipcw_weights, time.point = time.point)
  comp_bj_cox_pen[i] = comp_bj$pred.err.cox.pen
  comp_bj_cox[i] = comp_bj$pred.err.cox
  comp_bj_rand[i] = comp_bj$pred.err.rand
  
  comp_dr <- compare.methods.brier(x_train = x_train, y_train = y_train.dr, x_test = x_test, y_test = y_test, data_surv = data.surv, ipcw_weights = ipcw_weights, time.point = time.point)
  comp_dr_cox_pen[i] = comp_bj$pred.err.cox.pen
  comp_dr_cox[i] = comp_bj$pred.err.cox
  comp_dr_rand[i] = comp_bj$pred.err.rand
  
  comp_1 <- one.sim.dl.dr(x_train = x_train, y_train.dr = y_train.dr, x_test = x_test, y_test = y_test, data_surv = data.surv, ipcw_weights = ipcw_weights, time.point = time.point)
  comp_dl_dr[i] = comp_1$pred.err.dl
  comp_dl_dr_crossval[i] = comp_1$pred.err.dl.cv
  
  comp_2 <- one.sim.dl.bj(x_train = x_train, y_train.bj = y_train.bj, x_test = x_test, y_test = y_test, data_surv = data.surv, ipcw_weights = ipcw_weights, time.point = time.point)
  comp_dl_bj[i] = comp_2$pred.err.dl
  comp_dl_bj_crossval[i] = comp_2$pred.err.dl.cv
  
  #store results
  df = bind_rows(data.frame(x_train=x_train, y_train.bj = y_train.bj, y_train.dr=y_train.dr),
                 data.frame(x_test=x_test,y_test=y_test), data.frame(data.surv=data.surv),
                 data.frame(ipcw_weights = ipcw_weights),data.frame(time.point=time.point),
                 data.frame(comp = comp, comp_bj=comp_bj, comp_dr=comp_dr, comp_dl_bj=comp_1, comp_dl_dr=comp_2))
  data_name = paste(path,"setting4-iter",i,"data",".csv",sep="")
  write.csv(df, data_name)
}

save.image(file="setting4.RData")




##################################################################
#run code with existing data
data.set = pbc
#change to event (=death) indicator 
data.set$status[data.set$status == 1] = 0
data.set$status[data.set$status == 2] = 1
#exclude id from set
data.set = subset(data.set, select = -(1))

#delete rows with NA entries
data.set = data.set[complete.cases(data.set),]

#for KM estimators: need distinct survival times
data.set$time = data.set$time + runif(length(data.set$time),min=0, max=0.1)

obs.time = data.set$time
obs.delta = data.set$status
obs.covs = subset(data.set, select = -c(time,status))

obs.data <- data.frame(
  obs = obs.time,
  delta = obs.delta,
  covs = obs.covs
)

time.point = c(1:100)
time.point = time.point * 45
numb=length(time.point)

comp_cox_pen=rep(0,numb)
comp_cox=rep(0,numb)
comp_rand=rep(0,numb)

comp_bj_cox_pen=rep(0,numb)
comp_bj_cox=rep(0,numb)
comp_bj_rand=rep(0,numb)

comp_dr_cox_pen=rep(0,numb)
comp_dr_cox=rep(0,numb)
comp_dr_rand=rep(0,numb)

comp_dl_dr=rep(0,numb)
comp_dl_dr_crossval=rep(0,numb)
comp_dl_bj=rep(0,numb)
comp_dl_bj_crossval=rep(0,numb)

for (t in time.point){
  tmp <- create.crossval.data(obs.data, v=5, time.point=t)
  
  x_train <- tmp$x_train
  y_train <- tmp$y_train
  y_train.dr <- tmp$y_train.dr
  y_train.bj <- tmp$y_train.bj
  x_test <- tmp$x_test
  y_test <- tmp$y_test
  data.surv <- tmp$data.surv
  ipcw_weights <- tmp$data.ipcw
  time.point <- tmp$time.point
  
  comp <- compare.methods.brier(x_train = x_train, y_train = y_train, x_test = x_test, y_test = y_test, data_surv = data.surv, ipcw_weights = ipcw_weights, time.point = time.point)
  comp_cox_pen[i] = comp$pred.err.cox.pen
  comp_cox[i] = comp$pred.err.cox
  comp_rand[i] = comp$pred.err.rand
  
  comp_bj <- compare.methods.brier(x_train = x_train, y_train = y_train.bj, x_test = x_test, y_test = y_test, data_surv = data.surv, ipcw_weights = ipcw_weights, time.point = time.point)
  comp_bj_cox_pen[i] = comp_bj$pred.err.cox.pen
  comp_bj_cox[i] = comp_bj$pred.err.cox
  comp_bj_rand[i] = comp_bj$pred.err.rand
  
  comp_dr <- compare.methods.brier(x_train = x_train, y_train = y_train.dr, x_test = x_test, y_test = y_test, data_surv = data.surv, ipcw_weights = ipcw_weights, time.point = time.point)
  comp_dr_cox_pen[i] = comp_bj$pred.err.cox.pen
  comp_dr_cox[i] = comp_bj$pred.err.cox
  comp_dr_rand[i] = comp_bj$pred.err.rand
  
  comp_1 <- one.sim.dl.dr(x_train = x_train, y_train.dr = y_train.dr, x_test = x_test, y_test = y_test, data_surv = data.surv, ipcw_weights = ipcw_weights, time.point = time.point)
  comp_dl_dr[i] = comp_1$pred.err.dl
  comp_dl_dr_crossval[i] = comp_1$pred.err.dl.cv
  
  comp_2 <- one.sim.dl.bj(x_train = x_train, y_train.bj = y_train.bj, x_test = x_test, y_test = y_test, data_surv = data.surv, ipcw_weights = ipcw_weights, time.point = time.point)
  comp_dl_bj[i] = comp_2$pred.err.dl
  comp_dl_bj_crossval[i] = comp_2$pred.err.dl.cv
  
  #store results
  df = bind_rows(data.frame(x_train=x_train, y_train.bj = y_train.bj, y_train.dr=y_train.dr),
                 data.frame(x_test=x_test,y_test=y_test), data.frame(data.surv=data.surv),
                 data.frame(ipcw_weights = ipcw_weights),data.frame(time.point=time.point),
                 data.frame(comp = comp, comp_bj=comp_bj, comp_dr=comp_dr, comp_dl_bj=comp_1, comp_dl_dr=comp_2))
  data_name = paste(path,"realdata-timepoint",t,".csv",sep="")
  write.csv(df, data_name)
}

save.image(file="realData.RData")
