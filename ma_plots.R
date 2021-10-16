library(ggplot2)
library(ggdist)

#first: run ma_runSimulation to define the required parameters

######### Consider bj-trafo
#prepare data, numb is number of iterations performed in ma_runSimulation
tmp_cox <- data.frame(
  mse = comp_bj_cox,
  type = rep("cox", numb)
)

tmp_pencox <- data.frame(
  mse = comp_bj_cox_pen,
  type = rep("pencox", numb)
)

tmp_rand <- data.frame(
  mse = comp_bj_rand,
  type = rep("rand", numb)
)

tmp_bjdl <- data.frame(
  mse = comp_dl_bj,
  type = rep("BJCUDL",numb)
)

tmp_bjcvdl <- data.frame(
  mse = comp_dl_bj_crossval,
  type = rep("BJcvCUDL",numb)
)

comp_total = rbind(tmp_cox,tmp_pencox,tmp_rand, tmp_bjdl, tmp_bjcvdl)
comp_total$type = as.factor(comp_total$type)

#boxplot
p <- ggplot(comp_total, aes(x=type, y=mse, color = type)) + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA) + 
  geom_boxplot(
    width = .25, 
    outlier.shape = NA
  ) +
  geom_point(
    size = 1.3,
    alpha = .3,
    position = position_jitter(
      seed = 1, width = .1
    )
  )
p

######### Consider dr-trafo
#prepare data, numb is number of iterations performed in ma_runSimulation
tmp_cox <- data.frame(
  mse = comp_dr_cox,
  type = rep("cox", numb)
)

tmp_pencox <- data.frame(
  mse = comp_dr_cox_pen,
  type = rep("pencox", numb)
)

tmp_rand <- data.frame(
  mse = comp_dr_rand,
  type = rep("rand", numb)
)

tmp_bjdl <- data.frame(
  mse = comp_dl_dr,
  type = rep("drCUDL", numb)
)

tmp_bjcvdl <- data.frame(
  mse = comp_dl_dr_crossval,
  type = rep("drcvCUDL", numb)
)

comp_total = rbind(tmp_cox,tmp_pencox,tmp_rand, tmp_bjdl, tmp_bjcvdl)
comp_total$type = as.factor(comp_total$type)

#boxplot
p <- ggplot(comp_total, aes(x=type, y=mse, color = type)) + 
  ggdist::stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA) + 
  geom_boxplot(
    width = .25, 
    outlier.shape = NA
  ) +
  geom_point(
    size = 1.3,
    alpha = .3,
    position = position_jitter(
      seed = 1, width = .1
    )
  )
p
