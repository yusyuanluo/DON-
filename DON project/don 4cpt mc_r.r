# load the required package
library(hexbin)
library(mcmcr)
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)    # alternatively, this also loads %>%
library(tidyverse)
library(rstan)
library(bayesplot)
library(corrplot)
library(sensitivity)
library(pksensi)
theme_set(theme_light())
#-------------------MC setpoint check if the parameter distribution can capture the dataset--------------------
# Define the input variable
mName <- "don4cpt.model.r" # the model file put in the model folder
inName <- "don4cpt.mc.in.r" # the input file put in the infile folder

# Create the executable file
makemcsim(mName)

# Run!!
set.seed(1111)
mcsim(mName, inName)
out<-read.delim("sim.out")

vars <- names(out)
index <- which(vars == "Pct_d15g_ex_1.1" | vars == "Pct_d15g_ex_1.24")
X <- apply(out[index[1]:index[2]], 2, quantile,  c(0.5, 0.025, 0.975)) %>% t()
colnames(X) <- c("median", "LCL", "UCL")
df <- as.data.frame(X)
####create a column to assign time points
x <- seq(1, 24, 1)
df$x <- x
df
####obs data
y<-c(0.9448, 8.9734, 15.9371,21.4656, 29.9572, 33.7263, 37.9584, 41.6349, 42.9038, 43.8949, 44.9786, 45.5993, 46.0348, 46.9333, 48.5726, 50.2119, 53.4254, 54.185, 55.8706, 56.5839, 57.2972, 57.5475,58.0062, 58.2565)

obs_data<- data.frame(x,y)
ggplot(df, aes(x = x, y = median)) +
  geom_ribbon(aes(ymin = LCL, ymax = UCL), fill = "grey70", alpha = 0.5) + 
  geom_line() +
  geom_point(data = obs_data, x=x, y=y) +
  labs(x = "Time", y = "%excreted cumulative_d15g")

#--------------------------------------------------------MCMC simulation-------

# Define the input variable
mName <- "don4cpt.model.R" # the model file put in the model folder
inName <- "don4cptv4.mcmc.in.R" # the input file put in the infile folder

# Create the executable file
makemcsim(mName)

# Run!!
set.seed(1111)
out1 <- mcsim(mName, inName)
#chain 2
set.seed(2111)
out2 <- mcsim(mName, inName)
#chain 3
set.seed(3111)
out3 <- mcsim(mName, inName)
#chain 4
set.seed(4111)
out4 <- mcsim(mName, inName)

###define function: mcmc_array###
mcmc_array <- function(data, start_sampling = 0){
  n_chains <- length(data)
  sample_number <- dim(data[[1]])[1] - start_sampling
  dim <- c(sample_number, n_chains, dim(data[[1]])[2])
  n_iter <- dim(data[[1]])[1]
  n_param <- dim(data[[1]])[2]
  x <- array(sample_number:(n_iter * n_chains * n_param), dim = dim)
  for (i in 1:n_chains) {
    x[, i, ] <- as.matrix(data[[i]][(start_sampling + 1):n_iter, ])
  }
  dimnames(x)[[3]] <- names(data[[1]])
  x
}
#########
####organize 4 chains into one list####
sims <- mcmc_array(data = list(out1,out2,out3,out4))
dim(sims)

#####parameter examination of convergence
parms_name <- c("km_d15g.1.", "km_d3g.1.", "ku_d15g.1.", "ku_d3g.1.", "ke_don.1.")
color_scheme_set("mix-blue-red")
mcmc_trace(sims, pars = parms_name, facet_args = list(ncol = 1, strip.position = "left"))

parms_name <- c("kbile_d15g.1.", "kbile_d3g.1.","kgut_liv_d3g.1.", "kgut_liv_d15g.1.", "kgut_elim_d3g.1.", "kgut_elim_d15g.1.")
color_scheme_set("mix-blue-red")
mcmc_trace(sims, pars = parms_name, facet_args = list(ncol = 1, strip.position = "left"))


####density plot 
mcmc_dens_overlay(x = sims[5001:10000,,], pars = parms_name)

mcmc_dens_overlay(x = sims[5001:10000,,], pars = "LnData")

mcmc_pairs(sims[5001:10000,,], pars = parms_name, off_diag_fun = "hex")


#RSTAN PACKAGE parameter distribution output
monitor(sims[,,parms_name], digit=4)

#Take the last 2000 iterations
X <- sims[8001:10000,,] %>% matrix(nrow = 2000*4) 
colnames(X)<-colnames(out1)
write.table(X, file = "setpts.out", row.names = F , sep = "\t")
X_setpts <- mcsim("don4cpt.model.R", "don4cpt_setpts.in.R")
mcmc <- read.delim("simmc.out")

vars<-names(mcmc)
##return position of the variable
index <- which(vars == "Pct_d15g_ex_1.1" | vars == "Pct_d15g_ex_1.24")
##specifiy the column range and 2= "column [1=row]
d15g <- apply(mcmc[index[1]:index[2]], 2, quantile,  c(0.05, 0.5, 0.95)) %>% t()
d15g <- as.data.frame(d15g)
####create a column to assign time points
x <- seq(1, 24, 1)
d15g$x <- x
colnames(d15g)<-c("LCL","mean","HCL", "x")
####obs data

y<-c(0.9448, 8.9734, 15.9371,21.4656, 29.9572, 33.7263, 37.9584, 41.6349, 42.9038, 43.8949, 44.9786, 45.5993, 46.0348, 46.9333, 48.5726, 50.2119, 53.4254, 54.185, 55.8706, 56.5839, 57.2972, 57.5475,58.0062, 58.2565)

####plot pred vs obs _d15g
obs_data<- data.frame(x,y)
ggplot(d15g, aes(x = x, y = mean)) +
  geom_ribbon(aes(ymin = LCL, ymax = HCL), fill = "grey70", alpha = 0.5) + 
  geom_line() +
  geom_point(data = obs_data, x=x, y=y) +
  labs(x = "Time", y = "% cumulative d15g excreted")

##obs_interval _d15g
z<-c(0.9448, 8.0286, 6.9637, 5.5285, 8.4916, 3.7691, 4.2321, 3.6765, 1.2689, 0.9911, 1.0837, 0.6207, 0.4355, 0.8985, 1.6393, 1.6393, 3.2135, 0.7596, 1.6856, 0.7133, 0.7133, 0.2503, 0.4587, 0.2503)
obs_data2<-data.frame(x,z)
write.csv(d15g, "d15g.csv")

b<-read.csv("d15g_pred_4cpt.csv")
b<-as.data.frame(b)
ggplot(b, aes(x = time, y = pred)) +
  geom_line() +
  geom_point(data=obs_data2, x=x, y=z) +
  labs(x = "Time", y = "%excreted")

####--------------------plot pred vs obs don excreted
vars<-names(mcmc)
##return position of the variable
index <- which(vars == "Pct_don_ex_1.1" | vars == "Pct_don_ex_1.24")
##specifiy the column range and 2= "column [1=row]
don_a <- apply(mcmc[index[1]:index[2]], 2, quantile,  c(0.05, 0.5, 0.95)) %>% t()
don_a <- as.data.frame(don_a)
x <- seq(1, 24, 1)
don_a$x <- x
colnames(don_a)<-c("LCL","mean","HCL", "x")

####obs data_don_ex_cumulative
y<-c(3.4232, 7.1248, 10.1299, 13.598, 15.0598, 17.0232, 18.7165, 20.024, 20.5599, 20.7485, 21.2458, 21.473, 21.5459, 21.6959, 21.8459, 22.2274, 22.2275, 22.3389, 22.4118, 22.8319, 23.0591, 23.4483, 23.7527, 23.7528)
obs_data<- data.frame(x,y)
ggplot(don_a, aes(x = x, y = mean)) +
  geom_ribbon(aes(ymin = LCL, ymax = HCL), fill = "grey70", alpha = 0.5) + 
  geom_line() +
  geom_point(data = obs_data, x=x, y=y) +
  labs(x = "Time", y = "% excreted DON (cumulative)")

###obs data_don_ex_individual

z<-c(3.0051, 3.6996,3.0051, 3.4681, 1.4618, 1.9634, 1.6933, 1.3075, 0.5359, 0.1886, 0.4973, 0.2272, 0.0729, 0.15,0.15, 0.3815, 0.0001, 0.1114, 0.0729, 0.4201, 0.2272,0.3892,0.3044,0.0001)
obs_data<- data.frame(x,z)
write.csv(don_a, ' don_a.csv')
indi_pred<-read.csv('don_pred_4cpt.csv')
indi_pred<-as.data.frame(indi_pred)
ggplot(indi_pred, aes(x = time, y = pred)) +
  geom_line() +
  geom_point(data=obs_data, x=x, y=z) +
  labs(x = "Time", y = "% DON excreted")

##-------------Plot pred vs obs d3g_excreted ------------

##return position of the variable
index <- which(vars == "Pct_d3g_ex_1.1" | vars == "Pct_d3g_ex_1.24")
##specifiy the column range and 2= "column [1=row]
d3g_a <- apply(mcmc[index[1]:index[2]], 2, quantile,  c(0.05, 0.5, 0.95)) %>% t()
d3g_a <- as.data.frame(d3g_a)
x <- seq(1, 24, 1)
d3g_a$x <- x
colnames(d3g_a)<-c("LCL","mean","HCL", "x")

####obs data_d3g_ex_cumulative
y<-c(0.0651, 1.7044, 3.7141, 6.1868, 8.9373, 10.6692, 11.8918, 12.8829, 13.411, 13.8697, 14.2126, 14.8642, 14.8685, 15.0571, 15.7859, 16.0131, 17.6679, 18.2423, 19.4726, 19.4727, 19.4728, 19.7, 20.3902, 20.656)
obs_data<- data.frame(x,y)
ggplot(d3g_a, aes(x = x, y = mean)) +
  geom_ribbon(aes(ymin = LCL, ymax = HCL), fill = "grey70", alpha = 0.5) + 
  geom_line() +
  geom_point(data = obs_data, x=x, y=y) +
  labs(x = "Time", y = "% excreted d3g (cumulative)")

#####obsdata_d3g_ex interval
z<-c(0.0651, 1.6393, 2.0097, 2.4727, 2.7505, 1.7319, 1.2226, 0.9911, 0.5281, 0.4587, 0.3429, 0.6516, 0.0043, 0.1886, 0.7288, 0.2272, 1.6548, 0.5744, 1.2303, 0.0001, 0.0001, 0.2272, 0.6902, 0.2658)
obs_data<- data.frame(x,z)
write.csv(d3g_a, ' d3g_a.csv')
indi_pred<-read.csv('d3g_pred_4cpt.csv')
indi_pred<-as.data.frame(indi_pred)
ggplot(indi_pred, aes(x = time, y = pred)) +
  geom_line() +
  geom_point(data=obs_data, x=x, y=z) +
  labs(x = "Time", y = "% d3g excreted")
