library(randomForestSRC)

#run code with simulated data


nobs = 300
ncov = 25

obs.data = create_data_setting1(nobs,ncov)
obs = obs.data$obs.time
delta = obs.data$obs.status
covs = subset(obs.data, select = -c(obs, delta))
colnames(covs) <- paste0("W", colnames(covs), "_")

tmp <- create.crossval.data(obs.data, v=5, time.point=6)
