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

#Maximum resource body mass
max.res.bs <- 1000
#Number of resources
num.res <- 20

#Define resource body mass vector
res.bs <- round(seq(1,max.res.bs,length.out=num.res),0)

#Logged body mass sequence version
#res.bs <- round(lseq(1,max.res.bs,length.out=num.res),0); res.bs <- unique(res.bs); num.res <- length(res.bs)

#Define consumer body mass
cons.bs <- 20

#Define state matrices for consumer...
tmax <- 70
state.matrix <-  matrix(0,cons.bs,tmax)

#Define critical threshold for consumer
xc <- round((1/4)*cons.bs,0)
xmax <- cons.bs

rep.gain <- numeric(xmax)
for (i in xc:xmax) {
  rep.gain[i] <- i/(10*cons.bs)
}

#Define terminal fitness function
#Currently just a increasing function with respect to x... NEED REFINED
for (i in xc:cons.bs) {
    state.matrix[i,tmax] <- i/cons.bs
}


#Set Habitat Heterogeneity
#seq(0,1):: 0 is an even landscape, 1 is a patchy landscape

#Single value patchiness
hab.het <- rep(0,num.res)

#hab.het <- rep(1,num.res)

#Body size dependent patchiness
#Small animals are uniform; Large animals are patchy
#hab.het <- seq(0,1,length.out=num.res)

#Small animals are patchy; Large animals are uniform
#hab.het <- seq(1,0,length.out=num.res)



max.nu <- 10
min.nu <- 1

#Mean encounter rates of prey in region i should be scaled to body size of resource
#Dispersion of prey in region i should be scaled to the variability in nearest neighbors to nearest neighbors in i
#This is essentially the "Local Heterogeneity"

xi <- seq(10,1,length.out=num.res)
#nu[i,j] <- (1/res.bs[j])*10*(nn2.sd[i]+2)
nu <- max.nu - (hab.het*(max.nu-min.nu))

max.enc <- 20
f.res.patch <- matrix(0,(max.enc+1),num.res)
for (j in 1:num.res) {
  for (k in 1:(max.enc+1)) {
    f.res.patch[k,j] <-  dnbinom(k,mu = xi[j], size = nu[j])
  }
  f.res.patch[,j] <- f.res.patch[,j] / sum(f.res.patch[,j])
}
f.m <- f.res.patch

#plot negative binomial distributions
plot(seq(0,max.enc,1),f.res.patch[,1],type="l",col=colors[1],lwd=3,ylim=c(0,max(f.res.patch)))
for (i in 2:num.res) {
 lines(seq(0,max.enc,1),f.res.patch[,i],col=colors[i],lwd=3)
}

#Establish the probability of successfully capturing a single resource
#The 4 in the denominator moves the half-saturation point of the curve to the right.
encounters <- seq(0,max.enc,1)
rho.vec <- 1 - exp(-encounters^2/(5*max(encounters)))

#Consumer-resource mortality rates
Ratio.RC <- res.bs/cons.bs
mp1 <- 0.003
max.mort <- 0.01
mort <- max.mort - max.mort*(1 - 2*mp1)^(Ratio.RC^2)


#Foraging Gains and Costs (allometric and stoichiometric)
c.rep <- 0.05*cons.bs
eta <- numeric(num.res)
g.forage <- numeric(num.res)
for (i in 1:num.res) {
  #Foraging costs/gains conditional on consumer AND resource
  eta[i] <- 0.1
  #Gains are simply a proportion of the prey's body size (generally 10%)
  #This may become more complicated if stoichiometry is introduced
  g.forage[i] <- eta[i]*res.bs[i]
}
#Costs should scale allometrically :: Shouldn't costs be a function of the resource??? ~equiv. to speed differential?
#Latent trait probability
#The a1, a2, and a3 parameters are generally site-specific
# a1 <- 1.41
# a2 <- 3.73
# a3 <- -1.87
# p <- exp(a1 + a2*log(res.bs/cons.bs) + a3*log(res.bs/cons.bs)^2)
# pr.link <- p/(1+p)
#The probability that the interaction should NOT occur
#We will use this to weight the c.forage term
# pr.nolink <- 1 - pr.link
const <- 3.98/1000 *5 * 12 #3.98 mL O2 /1000 *5kcal * 12 hours = 0.2388 Joules
#Relatively high costs

#Prey-dependent costs
#c.forage <- pr.nolink*(const*cons.bs^(0.66))

#Single cost
c.forage <- rep(const*cons.bs^(0.66),num.res)



#Gradiant of decision possibilities based on resource similarities
num.dec <- 20
a.beta <- seq(1,5,length.out = num.dec)
b.beta <- seq(5,1,length.out = num.dec)
mean.beta <- a.beta/(a.beta + b.beta)
param.beta <- cbind(a.beta,b.beta)


# DOES NOT WORK ~ ASSYMMETRIC
# #Resource similarity ~ cosine similarity index
# #0 = entirely dissimilar
# #1 = entirely similar
# res.sim <- matrix(0,num.res,num.res)
# for (i in 1:num.res) {
#   for (j in 1:num.res) {
#     resi.att <- c(1,res.bs[i]) #in multidimensional trait space, this would be a vector of attributes
#     resj.att <- c(1,res.bs[j])
#     res.sim[i,j] <- (resi.att %*% resj.att)/ sqrt((resi.att %*% resi.att) * (resj.att %*% resj.att))
#   }
# }


#Resource similarity ~ Euclidean Distance ~ This will only work with numeric measures! Right now, just body size
#0 = entirely dissimilar
#1 = entirely similar
res.sim <- matrix(0,num.res,num.res)
max.diff <- sqrt((min(res.bs) - max(res.bs))^2)
for (i in 1:num.res) {
  for (j in 1:num.res) {
    resi.att <- res.bs[i] #in multidimensional trait space, this would be a vector of attributes
    resj.att <- res.bs[j]
    res.sim[i,j] <- 1 - (sqrt((resi.att - resj.att)^2) / max.diff)
  }
}

#Build decision probabilities
#For every 'current' resource, what is the preference probability based on each decision beta dist?

dec.ls <- list()
#Across the focal resources
for (j in 1:num.res) {
  dec.beta <- matrix(0,num.res,num.dec)
  #Across the different decision beta distributions...
  for (i in 1:num.dec) {
    dec.beta[,i] <- dbeta(res.sim[,j],shape1=param.beta[i,1],shape2=param.beta[i,2])
  }
  #Normalize to sum to 1
  #Each column is a prob. set... all choices are possible, but only one is chosen. Must sum to 1
  dec.beta.n <- apply(dec.beta,2,function(x){x/sum(x)})
  dec.ls[[j]] <- dec.beta.n
}
#plot(dec.ls[[1]][,1],type="l",col=colors[1], ylim=c(0,0.25),lwd=3)
#for (i in 2:10) {lines(dec.ls[[1]][,i],type="l",col=colors[i],lwd=3)}


###########################
# Additional Learning costs
###########################
learn.mag <- 0
learn.skew <- 0.5
c.learn <- learn.mag*(1-res.sim)^learn.skew


#Fitness matrices
#Build core consumer fitness matrix
W.xt <- matrix(0,cons.bs,tmax)
istar.xt <-  matrix(0,cons.bs,(tmax-1))
#Build nested fitness lists with the core being W.xt
W.nr <- list()
istar.nr <- list()
for (i in 1:num.res) {
  W.nr[[i]] <- list()
  istar.nr[[i]] <- list()
}
for (i in 1:num.res) {
  W.nr[[i]] <- W.xt
  istar.nr[[i]] <- istar.xt
}
#Build terminal fitness function
for (i in 1:num.res) {
  for (x in xc:xmax) {
    
    #NOTE: I change this linear x/cons.bs to rep.gain[x]... Might be more correct 10/14/14
    term.fit <- rep.gain[x] #x/cons.bs
    
    #Terminal fitness is independent of focal resource
    W.nr[[i]][x,tmax] <- term.fit
    
    #Or Terminal fitness is also dependent on mortality risk of focal resource
    #W.nr[[i]][x,tmax] <- term.fit * (1 - mort[i])
    
  }
}



#Compute fitness matrices
#time
#current resource
#energetic state

#Loop over 'focal resource'
for (r in 1:num.res) {
  
  #Define decision matrix for forcal resource r
  dec.m <- dec.ls[[r]]
  
  #Loop over time
  for (t in seq(tmax-1,1,-1)) {
    
    #Loop over energetic state
    for (x in seq(xc,xmax,1)) {
      
      value.max <- -10
      value <- numeric(num.dec)
      
      #Loop over decisions
      for (i in 1:num.dec) {
        
        #Define vector of preference probabilities across resources corresp. to given decision possibility
        #Each element is the preference probability for food resource i, given the current focal resource r
        pref.vec <- dec.m[,i]
        
        #Loop across resources
        #xp <- numeric(num.res)
        xp <- matrix(0,num.res,max.enc+1)
        for (rr in 1:num.res) {
          #determine the change in x for each amt k resources (which defines the rho.vec)
          delta.x <- x + rho.vec*(g.forage[rr] - (c.forage[rr] + c.learn[r,rr])) - (1-rho.vec)*(c.forage[rr] + c.learn[r,rr])
          #Multiply the prob(k)*delta x and sum
          #xp[rr] <- f.m[,rr] %*% delta.x
          xp[rr,] <- delta.x  
        }
        
        #We must establish boundary conditions
        xp[which(xp<xc)] <- xc
        xp[which(xp>xmax)] <- xmax

        W <- t(
          apply(xp,c(1,2),function(j){
          if ((j != xmax) && (j != xc)) {
            xp.low <- as.integer(j)
            xp.high <- xp.low + 1
            q <- xp.high - j
            q*W.nr[[r]][xp.low,t+1] + (1-q)*W.nr[[r]][xp.high,t+1]
          } else {
            if (j == xc) {
              xp.low <- xc
              xp.high <- xc+1
              q <- 1
              q*W.nr[[r]][xp.low,t+1] + (1-q)*W.nr[[r]][xp.high,t+1]
            } else {
              xp.low <- xmax-1
              xp.high <- xmax
              q <- 0
              q*W.nr[[r]][xp.low,t+1] + (1-q)*W.nr[[r]][xp.high,t+1]
            }
          }
        })
        )
        
        #Apply probabilities of finding k resources
        Wk <- numeric(num.res)
        for (rr in 1:num.res) {
          Wk[rr] <- W[,rr] %*% f.m[,rr]
        }

        Fx <- as.numeric(pref.vec %*% (rep.gain[x] + (1-mort)*Wk))
        
        value[i] <- Fx
        
      }#end decision loop
      
      #What is the fitness-maximizing decision?
      max.value <- which(value == max(value))
      
      #Is there more than one? I hope not!
      if (length(max.value) > 1) {
        print(paste(">1 at r=",r,"x=",x,"t=",t,sep=""))
        istar.nr[[r]][x,t] <- 0
        
        #I don't know if this is legit... 
        #it would assume that the organism would randomly choose a maximization strategy
        #and I don't know what istar would then be
        W.nr[[r]][x,t] <- mean(value[max.value])
        
      } else {
        #Record the fitness-maximizing decision
        istar.nr[[r]][x,t] <- max.value
        #Record the maximum fitness in the fitness matrix
        W.nr[[r]][x,t] <- value[max.value]
      }
      
    }#end energetic state loop
    
    
  }#end time loop
  
  
}#end current resource loop


####
# Alternative Cpp version

sourceCpp("src/SDP_single.cpp")

cout <- SDP_single(tmax, res.bs, cons.bs, rep.gain, 
           f.m, mort, dec.ls, rho.vec, c.learn, g.forage, c.forage)
W.nr <- cout[[1]]
istar.nr <- cout[[2]]


#Time-invariant analysis

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

filled_contour(seq(xc, xmax, length.out = nrow(istar.node)), 
               res.bs, 
               istar.node,
               levels = seq(1, max(istar.node)),col = col,
               lwd = 0.1,xlab="Energetic Reserves",ylab="Resource Size")
par(op)


#Exporting functions
#################


#Uniform Habitat
istarhab0 <- istar.node
save.image("Hab0.RData")
write.table(istarhab0,"istarhab0.csv",col.names=FALSE,row.names=FALSE,sep=",")

istarhab0t1 <- istar.node
write.table(istarhab0t1,"istarhab0t1.csv",col.names=FALSE,row.names=FALSE,sep=",")



#Patchy Habitat
istarhab1 <- istar.node
save.image("Hab1.RData")
write.table(istarhab1,"istarhab1.csv",col.names=FALSE,row.names=FALSE,sep=",")

istarhab1t1 <- istar.node
write.table(istarhab1t1,"istarhab1t1.csv",col.names=FALSE,row.names=FALSE,sep=",")




#Patchines increasing with body size
istarhab_inc <- istar.node
save.image("Hab_inc.RData")
write.table(istarhab_inc,"istarhab_inc.csv",col.names=FALSE,row.names=FALSE,sep=",")

#Patchiness decreasing with body size
istarhab_dec <- istar.node
save.image("Hab_dec.RData")
write.table(istarhab_dec,"istarhab_dec.csv",col.names=FALSE,row.names=FALSE,sep=",")



#Some analytics on the solutions

# Has the simulation reached a steady state?
D12 <- numeric(tmax-1)
for (t in 1:(tmax-2)) {
  D1 <- do.call(cbind,lapply(istar.nr,function(x){x[,t]}))
  D2 <- do.call(cbind,lapply(istar.nr,function(x){x[,t+1]}))
  D12[t] <- sqrt(sum((D1 - D2)^2))
}
pdf("/Users/justinyeakel/Dropbox/PostDoc/2014_Foodweb_SDP/Manuscript/Figure_stationarity.pdf", width=7, height=4)
plot(D12,type="l",ylab="Matrix Diff (t[y] - t[y-1])",xlab="t")
dev.off()

#Plotting Fitness values
r <- 5
plot(W.nr[[r]][1,],type="l",ylim=c(0,3))
for (i in 2:20) {lines(W.nr[[r]][i,])}


# Forward equation


#####################
# Forward equation

#Starting values
N <- 20 #starting Num. of individuals
tsim <- 1000 #Simulation time
r.alpha <- 3 #Recruit density-independent growth rate
r.beta <- 1/80 #Recruit density-dependent growth rate
r.comp <- 1 #degree of compensation

sourceCpp("src/SDP_single_foreq.cpp")

cout.foreq <- SDP_single_foreq(tmax, res.bs, cons.bs, xc, rep.gain, 
                 f.m, mort, dec.ls, rho.vec, c.learn, g.forage, c.forage, 
                 W.nr, istar.nr, N, tsim, eta,
                 r.alpha, r.beta, r.comp)

pop.traj <- cout.foreq[[1]]
x.traj <- cout.foreq[[2]]; x.mean <- unlist(lapply(x.traj,mean))
t.traj <- cout.foreq[[3]]; t.mean <- unlist(lapply(t.traj,mean))
d.traj <- cout.foreq[[4]]
trophic.traj <- cout.foreq[[5]]
fitness <- cout.foreq[[6]]; fit.mean <- unlist(lapply(fitness,mean))

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








