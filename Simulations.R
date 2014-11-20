rm(list=c(ls()))

library(RColorBrewer)
library(Rcpp)
library(igraph)
library(beanplot)
library(lattice)
library(emdbook)
source("src/filled_contour.r")
source("src/smooth_pal.R")

colors <- brewer.pal(5,"Spectral")
colors <- smooth_pal(colors, 5)

source("R/SDP_initialcond.R")
sourceCpp("src/SDP_single.cpp")
sourceCpp("src/SDP_single_foreq.cpp")


#Initial Conditions
Rout <- SDP_initialcond(
  max.res.bs <- 1500, 
  num.res <- 15, 
  cons.bs <- 20, 
  tmax <- 500, 
  hab.het.var <- 0, 
  max.enc <- 20, 
  learn.mag <- 0)

res.bs <- Rout[[1]]
xc <- Rout[[2]]
rep.gain <- Rout[[3]]
f.m <- Rout[[4]]
mort <- Rout[[5]]
dec.ls <- Rout[[6]]
rho.vec <- Rout[[7]]
c.learn <- Rout[[8]]
g.forage <- Rout[[9]]
c.forage <- Rout[[10]]
eta <- Rout[[11]]

#Run Backwards simulation
cout <- SDP_single(
  tmax = tmax, 
  res_bs = res.bs, 
  cons_bs = cons.bs,
  xc = xc,
  rep_gain = rep.gain, 
  f_m = f.m, 
  mort = mort, 
  dec_ls = dec.ls, 
  rho_vec = rho.vec, 
  c_learn = c.learn, 
  g_forage = g.forage, 
  c_forage = c.forage
)
W.nr <- cout[[1]]
istar.nr <- cout[[2]]

#Plot the Decision matrix
####### 


col <- RColorBrewer::brewer.pal(6, "Spectral")
col <- rev(smooth_pal(col, 4))

op <- par(mfrow = c(2,2),
          oma = c(6,6,0,0) + 0.1,
          mar = c(0,0,2,2) + 0.1,
          mgp = c(2, 1, 0))

s <- round(seq(1,tmax-1,length.out=4))
for (i in 1:4) {
  time.stamp <- s[i]
  istar.node <- do.call(cbind,lapply(istar.nr,function(x){x[,time.stamp]}))
  #Eliminate the <xc rows
  istar.node <- istar.node[-seq(1,xc-1,1),]
  filled_contour(seq(xc, cons.bs, length.out = nrow(istar.node)), 
                 res.bs, 
                 istar.node,
                 levels = seq(1, max(istar.node)),col = col,
                 lwd = 0.1,xlab="Energetic Reserves",ylab="Resource Size")
}
par(op)

# Has the simulation reached a steady state?
D12 <- numeric(tmax-1)
for (t in 1:(tmax-2)) {
  D1 <- do.call(cbind,lapply(istar.nr,function(x){x[,t]}))
  D2 <- do.call(cbind,lapply(istar.nr,function(x){x[,t+1]}))
  D12[t] <- sqrt(sum((D1 - D2)^2))
}
plot(D12,type="l",ylab="Matrix Diff (t[y] - t[y-1])",xlab="t")





#####
#Starting values
N <- 20 #starting Num. of individuals
tsim <- 1000 #Simulation time
r.alpha <- 2 #Recruit density-independent growth rate
r.beta <- 1/80 #Recruit density-dependent growth rate
r.comp <- 1 #degree of compensation

cout.foreq <- SDP_single_foreq(
  tmax = tmax, 
  res_bs = res.bs, 
  cons_bs = cons.bs, 
  xc = xc, 
  rep_gain = rep.gain, 
  f_m = f.m, 
  mort = mort, 
  dec_ls = dec.ls, 
  rho_vec = rho.vec, 
  c_learn = c.learn, 
  g_forage = g.forage,  
  c_forage = c.forage, 
  W_nr = W.nr, 
  istar_nr = istar.nr, 
  N_init = N, 
  tsim = tsim, 
  alpha = r.alpha, 
  beta = r.beta, 
  comp = r.comp
)

pop.traj <- cout.foreq[[1]]
x.traj <- cout.foreq[[2]]; x.mean <- unlist(lapply(x.traj,median))
t.traj <- cout.foreq[[3]]; t.mean <- unlist(lapply(t.traj,median))
d.traj <- cout.foreq[[4]]; d.mean <- unlist(lapply(d.traj,median))
#Trophic ID varies between 0 and num.res-1 bc they are Cpp indices
trophic.traj <- cout.foreq[[5]]; trophic.mean <- unlist(lapply(trophic.traj,median))
fitness <- cout.foreq[[6]]; fit.mean <- unlist(lapply(fitness,median))

#Note: there is some funky stuff going on with d.traj... Nov 18,2014
#This is because it is recording individuals at time tmax 
#(there are no decisions to make there)
#Set death to occur at tmax-1 for now, and it fixes issue

burnin <- 200
plot(pop.traj,type="l")

#proportional trophic interactions
trophic.v <- unlist(trophic.traj[burnin:tsim])
trophic.prop <- numeric(num.res)
for (i in 1:num.res) {
  trophic.prop[i] <- length(which(trophic.v == i-1))/length(trophic.v)
}
#plot(res.bs/cons.bs,trophic.prop,xlab="Resource/Consumer body size",ylab="Probability")
plot(log(res.bs,10),trophic.prop,pch=16,
     xlab="Log10 resource body size",ylab="Probability",xlim=c(1,4))

#Plot fitness
plot(fit.mean,type="l")
plot(fit.mean[burnin:tsim],pop.traj[burnin:tsim],pch='.',cex=3,col="darkgray")
#This is negative bc higher population sizes are due to accumulated older individuals
#and older individuals have lower fitness values ::
plot(t.mean[burnin:tsim],fit.mean[burnin:tsim],pch='.',cex=3,col="darkgray")







