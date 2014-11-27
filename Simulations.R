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

W.nr.list <- list()
istar.nr.list <- list()
cout.foreq <- list()
num.its <- 10
land.hetero <- seq(0,1,length.out=num.its)
for (j in 4:num.its) {
  print(paste("j=",j,sep=""))
  #Initial Conditions
  Rout <- SDP_initialcond(
    max.res.bs <- 2000, 
    num.res <- 15, 
    cons.bs <- 20, 
    tmax <- 100, 
    hab.het.var <- round(land.hetero[j],2), 
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
  
  W.nr.list[[j]] <- W.nr
  istar.nr.list[[j]] <- istar.nr
  
  #Starting values
  N <- 20 #starting Num. of individuals
  tsim <- 5000 #Simulation time
  r.alpha <- 2 #Recruit density-independent growth rate
  r.beta <- 1/100 #Recruit density-dependent growth rate
  r.comp <- 1 #degree of compensation
  
  cout.foreq[[j]] <- SDP_single_foreq(
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
  
} 



save.image(file="SDP.RData")
  #   istar.node <- do.call(cbind,lapply(istar.nr,function(x){x[,time.stamp]}))
  
  #Plot the Decision matrix
  ####### 
#   col <- RColorBrewer::brewer.pal(6, "Spectral")
#   col <- rev(smooth_pal(col, 4))
#   
#   op <- par(mfrow = c(2,2),
#             oma = c(6,6,0,0) + 0.1,
#             mar = c(0,0,2,2) + 0.1,
#             mgp = c(2, 1, 0))
#   
#   s <- round(seq(1,tmax-1,length.out=4))
#   for (i in 1:4) {
#     time.stamp <- s[i]
#     istar.node <- do.call(cbind,lapply(istar.nr,function(x){x[,time.stamp]}))
#     #Eliminate the <xc rows
#     istar.node <- istar.node[-seq(1,xc-1,1),]
#     filled_contour(seq(xc, cons.bs, length.out = nrow(istar.node)), 
#                    res.bs, 
#                    istar.node,
#                    levels = seq(1, max(istar.node)),col = col,
#                    lwd = 0.1,xlab="Energetic Reserves",ylab="Resource Size")
#   }
#   par(op)
#   
#   # Has the simulation reached a steady state?
#   D12 <- numeric(tmax-1)
#   for (t in 1:(tmax-2)) {
#     D1 <- do.call(cbind,lapply(istar.nr,function(x){x[,t]}))
#     D2 <- do.call(cbind,lapply(istar.nr,function(x){x[,t+1]}))
#     D12[t] <- sqrt(sum((D1 - D2)^2))
#   }
#   plot(D12,type="l",ylab="Matrix Diff (t[y] - t[y-1])",xlab="t")
  #####


load("SDP.Rdata")
#Analysis of Decision matrices
#Ontogenetic analysis
dec.dist.list <- list()
for (j in 1:num.its) {
  dec.dist <- list()
  for (i in 1:(tmax-1)) {
    istar.node <- do.call(cbind,lapply(istar.nr.list[[j]],function(x){x[,i]}))
    dec.numeric <- as.numeric(istar.node)
    dec.dist[[i]] <- dec.numeric[which(dec.numeric != 0)]
  }
  dec.dist.list[[j]] <- dec.dist
}
col.disc <- brewer.pal(10,"RdYlGn")
par(mfrow=c(1,2))
plot(unlist(lapply(dec.dist.list[[1]],mean)),type="l",lwd=3,col=col.disc[1],
     xlab="Time",ylab="Prey-switching mean")
for (j in 2:num.its) {
  lines(unlist(lapply(dec.dist.list[[j]],mean)),type="l",lwd=3,col=col.disc[j])
} 
plot(unlist(lapply(dec.dist.list[[1]],sd)),type="l",lwd=3,col=col.disc[1],
     xlab="Time",ylab="Prey-switching SD")
for (j in 2:num.its) {
  lines(unlist(lapply(dec.dist.list[[j]],sd)),type="l",lwd=3,col=col.disc[j])
} 


# 
# 
# cout.foreq <- list()
# for (j in 1:num.its) {
#   #####
#   #Starting values
#   N <- 20 #starting Num. of individuals
#   tsim <- 5000 #Simulation time
#   r.alpha <- 2 #Recruit density-independent growth rate
#   r.beta <- 1/100 #Recruit density-dependent growth rate
#   r.comp <- 1 #degree of compensation
#   W.nr.h <- W.nr[[j]]
#   istar.nr.h <- istar.nr[[j]]
#   
#   cout.foreq[[j]] <- SDP_single_foreq(
#     tmax = tmax, 
#     res_bs = res.bs, 
#     cons_bs = cons.bs, 
#     xc = xc, 
#     rep_gain = rep.gain, 
#     f_m = f.m, 
#     mort = mort, 
#     dec_ls = dec.ls, 
#     rho_vec = rho.vec, 
#     c_learn = c.learn, 
#     g_forage = g.forage,  
#     c_forage = c.forage, 
#     W_nr = W.nr.h, 
#     istar_nr = istar.nr.h, 
#     N_init = N, 
#     tsim = tsim, 
#     alpha = r.alpha, 
#     beta = r.beta, 
#     comp = r.comp
#   )
#   
# }



#Note: there is some funky stuff going on with d.traj... Nov 18,2014
#This is because it is recording individuals at time tmax 
#(there are no decisions to make there)
#Set death to occur at tmax-1 for now, and it fixes issue
trophic.prop.list <- list()
pop.traj <- list()
pop.mean <- numeric(num.its)
pop.sd <- numeric(num.its)
x.traj <- list(); x.mean <- list()
t.traj <- list(); t.mean <- list()
d.traj <- list(); d.mean <- list()
trophic.traj <- list(); trophic.mean <- list()
fitness <- list(); fit.mean <- list()
for (j in 1:num.its) {
  burnin <- 1000
  pop.traj[[j]] <- cout.foreq[[j]][[1]]
  pop.mean[[j]] <- median(pop.traj[[j]][burnin:tsim])
  pop.sd[[j]] <- sd(pop.traj[[j]][burnin:tsim])
  x.traj[[j]] <- cout.foreq[[j]][[2]]; x.mean[[j]] <- unlist(lapply(x.traj[[j]],median))
  t.traj[[j]] <- cout.foreq[[j]][[3]]; t.mean[[j]] <- unlist(lapply(t.traj[[j]],median))
  d.traj[[j]] <- cout.foreq[[j]][[4]]; d.mean[[j]] <- unlist(lapply(d.traj[[j]],median))
  #Trophic ID varies between 0 and num.res-1 bc they are Cpp indices
  trophic.traj[[j]] <- cout.foreq[[j]][[5]]; trophic.mean[[j]] <- unlist(lapply(trophic.traj[[j]],median))
  fitness[[j]] <- cout.foreq[[j]][[6]]; fit.mean[[j]] <- unlist(lapply(fitness[[j]],median))

  
  #proportional trophic interactions
  trophic.v <- unlist(trophic.traj[[j]][burnin:tsim])
  trophic.prop <- numeric(num.res)
  for (i in 1:num.res) {
    trophic.prop[i] <- length(which(trophic.v == i-1))/length(trophic.v)
  }
  trophic.prop.list[[j]] <- trophic.prop
}

#plot(res.bs/cons.bs,trophic.prop,xlab="Resource/Consumer body size",ylab="Probability")
par(mfrow=c(1,1))
plot(log(res.bs,10),trophic.prop.list[[1]],pch=16,type="b",col="black",
     xlab="Log10 resource body size",ylab="Probability",xlim=c(1,4))
for (j in 2:5) {
  points(log(res.bs,10),trophic.prop.list[[j]],pch=16,type="b",col=col.disc[j])
}

# burnin <- 200
# plot(pop.traj,type="l")
# 
# #proportional trophic interactions
# trophic.v <- unlist(trophic.traj[burnin:tsim])
# trophic.prop <- numeric(num.res)
# for (i in 1:num.res) {
#   trophic.prop[i] <- length(which(trophic.v == i-1))/length(trophic.v)
# }
# #plot(res.bs/cons.bs,trophic.prop,xlab="Resource/Consumer body size",ylab="Probability")
# plot(log(res.bs,10),trophic.prop,pch=16,
#      xlab="Log10 resource body size",ylab="Probability",xlim=c(1,4))

#Plot fitness
plot(fit.mean[[1]],type="l")
plot(fit.mean[[1]][burnin:tsim],pop.traj[[1]][burnin:tsim],pch='.',cex=5,col="darkgray")
#This is slightly negative bc higher population sizes are due to accumulated older individuals
#and older individuals have lower fitness values ::
plot(t.mean[[1]][burnin:tsim],fit.mean[[1]][burnin:tsim],pch='.',cex=5,col=paste(col.disc[1],"50",sep=""),
     xlim=c(18,30),ylim=c(5.5,6.5),xlab="Age mean",ylab="Fitness mean")
for (j in 2:num.its) {
  points(t.mean[[j]][burnin:tsim],fit.mean[[j]][burnin:tsim],pch='.',cex=5,col=paste(col.disc[j],"50",sep=""))
}


#Mean fitness vs. mean & sd population size
par(mfrow=c(1,2))
plot(land.hetero,pop.mean,col=col.disc,pch=16,xlab="Heterogeneity",ylab="Median population size")
plot(unlist(lapply(fit.mean,median)),pop.mean,xlab="Median fitness",ylab="Median population size",pch=16,col=col.disc)

plot(land.hetero,unlist(lapply(t.mean,sd)),xlab="Heterogeneity",ylab="Median age",pch=16,col=col.disc)


plot(unlist(lapply(t.mean,median)),pop.mean,xlab="Median age",ylab="Median population")
plot(unlist(lapply(fit.mean,sd)),pop.sd,xlab = "Fitness SD", ylab = "Population SD")



