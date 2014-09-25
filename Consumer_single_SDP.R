rm(list=c(ls()))

library(RColorBrewer)
library(Rcpp)
library(igraph)
library(beanplot)
library(lattice)
source("src/filled_contour.r")
colors <- brewer.pal(10,"Spectral")

num.res <- 10
res.bs <- round(rgamma(num.res,shape=5,scale=2),0)
res.bs <- seq(1,num.res,length.out=num.res)
cons.bs <- 10

#Define state matrices for consumer
tmax <- 10
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
#seq(0,1):: 0 is an even landscape, 1 is an open landscape
hab.het <- 0

max.nu <- 5
min.nu <- 1

#Mean encounter rates of prey in region i should be scaled to body size of resource
#Dispersion of prey in region i should be scaled to the variability in nearest neighbors to nearest neighbors in i
#This is essentially the "Local Heterogeneity"
xi <- numeric(num.res)
nu <- numeric(num.res)
for (j in 1:num.res) {
  xi[j] <- (1/res.bs[j])*10
  #nu[i,j] <- (1/res.bs[j])*10*(nn2.sd[i]+2)
  nu[j] <- max.nu - (hab.het*(max.nu-min.nu))
}
#nu <- nu + 1 #Ensure that there are non-zero nu's


max.enc <- 20
f.res.patch <- matrix(0,(max.enc+1),num.res)
for (j in 1:num.res) {
  for (k in 1:(max.enc+1)) {
    f.res.patch[k,j] <-  dnbinom(k,mu = xi[j], size = nu[j])
  }
  f.res.patch[,j] <- f.res.patch[,j] / sum(f.res.patch[,j])
}

#plot negative binomial distributions
#plot(seq(0,max.enc,1),f.res.patch[,1],type="l",col=colors[1],lwd=3,ylim=c(0,0.5))
#for (i in 2:num.res) {
#  lines(seq(0,max.enc,1),f.res.patch[,i],col=colors[i],lwd=3)
#}

#Establish the probability of successfully capturing a single resource
encounters <- seq(0,max.enc,1)
rho.vec <- 1 - exp(-encounters^2/max(encounters))

#Consumer-resource mortality rates
Ratio.RC <- res.bs/cons.bs
mp1 <- 0.1
mort <- 0.5 - 0.5*(1 - 2*mp1)^(Ratio.RC^2)

#Foraging Gains and Costs (allometric and stoichiometric)
c.rep <- 0.05*cons.bs
eta <- numeric(num.res)
g.forage <- numeric(num.res)
for (i in 1:num.res) {
  #Foraging costs/gains conditional on consumer AND resource
  eta[i] <- 0.50
  #Gains are simply a proportion of the prey's body size (generally 10%)
  #This may become more complicated if stoichiometry is introduced
  g.forage[i] <- eta[i]*res.bs[i]
}
#Costs should scale allometrically :: Shouldn't costs be a function of the resource??? ~equiv. to speed differential?
#Latent trait probability
#The a1, a2, and a3 parameters are generally site-specific
a1 <- 1.41
a2 <- 3.73
a3 <- -1.87
p <- exp(a1 + a2*log(res.bs/cons.bs) + a3*log(res.bs/cons.bs)^2)
pr.link <- p/(1+p)
#The probability that the interaction should NOT occur
#We will use this to weight the c.forage term
pr.nolink <- 1 - pr.link
#c.forage <- 0.5*cons.bs^(1/4)
#c.forage <- 0.8*min(g.forage)
const <- 0.2388 #3.98 mL O2 /1000 *5kcal * 12 hours = 0.2388 Joules
#Relatively high costs
c.forage <- pr.nolink*(const*cons.bs^(0.66))
#Costs lower than the lowest gain
#c.forage <- 0.8*min(g.forage)

#Gradiant of decision possibilities based on resource similarities
num.dec <- 20
a.beta <- seq(1,5,length.out = num.dec)
b.beta <- seq(5,1,length.out = num.dec)
mean.beta <- a.beta/(a.beta + b.beta)
param.beta <- cbind(a.beta,b.beta)

#Resource similarity ~ cosine similarity index
#0 = entirely dissimilar
#1 = entirely similar
res.sim <- matrix(0,num.res,num.res)
for (i in 1:num.res) {
  for (j in 1:num.res) {
    resi.att <- c(1,res.bs[i]) #in multidimensional trait space, this would be a vector of attributes
    resj.att <- c(1,res.bs[j])
    res.sim[i,j] <- (resi.att %*% resj.att)/ sqrt((resi.att %*% resi.att) * (resj.att %*% resj.att))
  }
}
#Right now, having trouble with res.sim = 1. It is generally extremely low for all beta dist. other than
#beta (5,1) ~ the last. Rather, the prob(1) should increase as mean beta increases. 
#res.sim <- res.sim - 0.1

# #stretch between 0 and 1 for each similarity column
# for (i in 1:num.res) {
#   res.sim.low <- res.sim[,i] - (min(res.sim[,i]))
#   res.sim[,i] <- res.sim.low/max(res.sim.low)
# }


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
    term.fit <- x/cons.bs
    W.nr[[i]][x,tmax] <- term.fit
  }
}



#Compute fitness matrices
#time
#current resource
#energetic state

#Define the encounter prob. matrix (k by res)
f.m <- f.res.patch

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
        pref.vec <- dec.m[,i]
        
        #Loop across resources
        xp <- numeric(num.res)
        for (rr in 1:num.res) {
          delta.x <- x + rho.vec*(g.forage[rr] - c.forage[rr]) - (1-rho.vec)*c.forage[rr]
          xp[rr] <- f.m[,rr] %*% delta.x
          #We must establish boundary conditions
          # xc <= xp <= xmax
          if (xp[rr] < xc) {xp[rr] <- xc}
          if (xp[rr] > xmax) {xp[rr] <- xmax}
        }
        
        #Interpolation function
        xp.low <- numeric(num.res); xp.high <- numeric(num.res); q <- numeric(num.res)
        for (j in 1:num.res) {
          if ((xp[j] != xmax) && (xp[j] != xc)) {
            xp.low[j] <- as.integer(xp[j])
            xp.high[j] <- xp.low[j] + 1
            q[j] <- xp.high[j] - xp[j]
          }
          if (xp[j] == xc) {
            xp.low[j] <- xc
            xp.high[j] <- xc+1
            q[j] <- 1
          }
          if (xp[j] == xmax) {
            xp.low[j] <- xmax-1
            xp.high[j] <- xmax
            q[j] <- 0
          }
        }
        #q*xp.low + (1-q)*xp.high ;; if q is 1, all weight on xp.low; vice versa
        
        #W is now grabbed from the ascribed fitness value at time t+1 & interpolated
        W <- q*W.nr[[r]][xp.low,t+1] + (1-q)*W.nr[[r]][xp.high,t+1]
        
        Fx <- as.numeric(pref.vec %*% (rep.gain[x] + (1-mort)*W))
        
        value[i] <- Fx
        
      }#end decision loop
      
      #Record the fitness-maximizing decision
      istar.nr[[r]][x,t] <- which(value == max(value))
      
      if (length(which(value==max(value))) > 1) {
        print(paste(">1 at r=",r,"x=",x,"t=",t,sep=""))
      }
      
      #Record the maximum fitness in the fitness matrix
      W.nr[[r]][x,t] <- max(value)
      
    }#end energetic state loop
    
    
  }#end time loop
  
  
}#end current resource loop

#Time-invariant analysis

istar.node <- do.call(cbind,lapply(istar.nr,function(x){x[,1]}))
#Eliminate the <xc rows
istar.node <- istar.node[-seq(1,xc-1,1),]

col <- RColorBrewer::brewer.pal(9, "Blues")
source("src/smooth_pal.R")
col <- smooth_pal(col, 5)

op <- par(mfrow = c(1,1),
          oma = c(5,4,0,0) + 0.1,
          mar = c(0,3,1,1) + 0.1,
          mgp = c(2, 1, 0))

filled_contour(seq(xc, xmax, length.out = nrow(istar.node)), 
               seq(1, num.res,length.out = ncol(istar.node)), 
               istar.node,
               levels = seq(1, max(istar.node)),col = col,
               lwd = 0.1,xlab="Energetic Reserves",ylab="Resource Size")
par(op)


write.table(istarhab0,"istarhab0.csv",col.names=FALSE,row.names=FALSE,sep=",")
write.table(istarhab1,"istarhab1.csv",col.names=FALSE,row.names=FALSE,sep=",")






