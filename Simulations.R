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
  max.res.bs <- 1000, 
  num.res <- 20, 
  cons.bs <- 20, 
  tmax <- 80, 
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
cout <- SDP_single(tmax = tmax, 
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


time.stamp <- 1
istar.node <- do.call(cbind,lapply(istar.nr,function(x){x[,time.stamp]}))
#Eliminate the <xc rows
istar.node <- istar.node[-seq(1,xc-1,1),]

col <- RColorBrewer::brewer.pal(6, "Spectral")
col <- rev(smooth_pal(col, 4))

op <- par(mfrow = c(1,1),
          oma = c(6,6,0,0) + 0.1,
          mar = c(0,0,1,1) + 0.1,
          mgp = c(2, 1, 0))

filled_contour(seq(xc, cons.bs, length.out = nrow(istar.node)), 
               res.bs, 
               istar.node,
               levels = seq(1, max(istar.node)),col = col,
               lwd = 0.1,xlab="Energetic Reserves",ylab="Resource Size")
par(op)



#Starting values
N <- 20 #starting Num. of individuals
tsim <- 1000 #Simulation time
r.alpha <- 2 #Recruit density-independent growth rate
r.beta <- 1/80 #Recruit density-dependent growth rate
r.comp <- 1 #degree of compensation

cout.foreq <- SDP_single_foreq(tmax = tmax, 
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
x.traj <- cout.foreq[[2]]; x.mean <- unlist(lapply(x.traj,mean))
t.traj <- cout.foreq[[3]]; t.mean <- unlist(lapply(t.traj,mean))
d.traj <- cout.foreq[[4]]
trophic.traj <- cout.foreq[[5]]
fitness <- cout.foreq[[6]]; fit.mean <- unlist(lapply(fitness,mean))

#Note: there is some funky stuff going on with d.traj... Nov 18,2014


plot(pop.traj,type="l")

#proportional trophic interactions
trophic.v <- unlist(trophic.traj[1:tsim])
trophic.prop <- numeric(20)
for (i in 1:20) {
  trophic.prop[i] <- length(which(trophic.v == i))/length(trophic.v)
}
plot(res.bs,trophic.prop)

plot(rep(1,length(trophic.traj[[1]])),trophic.traj[[1]],xlim=c(0,tsim),ylim=c(0,num.res),pch='.')
for (i in 2:tsim) {
  points(rep(i,length(trophic.traj[[i]])),trophic.traj[[i]],pch='.')
}

#Plot fitness
plot(fit.mean,type="l")


