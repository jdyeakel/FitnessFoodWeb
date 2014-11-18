SDP_initialcond <- function(max.res.bs, num.res, cons.bs, tmax, hab.het.var, 
                            max.enc, learn.mag) {
  
  #library(RColorBrewer)
  #library(Rcpp)
  #library(igraph)
  #library(beanplot)
  #library(lattice)
  #library(emdbook)
  #source("src/filled_contour.r")
  #source("src/smooth_pal.R")
  
  #colors <- brewer.pal(5,"Spectral")
  #colors <- smooth_pal(colors, 5)
  
  #Maximum resource body mass
  ## max.res.bs <- 1000
  #Number of resources
  ## num.res <- 20
  
  #Define resource body mass vector
  res.bs <- round(seq(1,max.res.bs,length.out=num.res),0)
  
  #Logged body mass sequence version
  #res.bs <- round(lseq(1,max.res.bs,length.out=num.res),0); res.bs <- unique(res.bs); num.res <- length(res.bs)
  
  #Define consumer body mass
  ## cons.bs <- 20
  
  #Define state matrices for consumer...
  ## tmax <- 70
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
  ##hab.het.var <- 0
  hab.het <- rep(hab.het.var,num.res)
  
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
  
  ##max.enc <- 20
  f.res.patch <- matrix(0,(max.enc+1),num.res)
  for (j in 1:num.res) {
    for (k in 1:(max.enc+1)) {
      f.res.patch[k,j] <-  dnbinom(k,mu = xi[j], size = nu[j])
    }
    f.res.patch[,j] <- f.res.patch[,j] / sum(f.res.patch[,j])
  }
  f.m <- f.res.patch
  
  #plot negative binomial distributions
  #   plot(seq(0,max.enc,1),f.res.patch[,1],type="l",col=colors[1],lwd=3,ylim=c(0,max(f.res.patch)))
  #   for (i in 2:num.res) {
  #     lines(seq(0,max.enc,1),f.res.patch[,i],col=colors[i],lwd=3)
  #   }
  
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
  
  
  #   #Fitness matrices
  #   #Build core consumer fitness matrix
  #   W.xt <- matrix(0,cons.bs,tmax)
  #   istar.xt <-  matrix(0,cons.bs,(tmax-1))
  #   #Build nested fitness lists with the core being W.xt
  #   W.nr <- list()
  #   istar.nr <- list()
  #   for (i in 1:num.res) {
  #     W.nr[[i]] <- list()
  #     istar.nr[[i]] <- list()
  #   }
  #   for (i in 1:num.res) {
  #     W.nr[[i]] <- W.xt
  #     istar.nr[[i]] <- istar.xt
  #   }
  #   #Build terminal fitness function
  #   for (i in 1:num.res) {
  #     for (x in xc:xmax) {
  #       
  #       #NOTE: I change this linear x/cons.bs to rep.gain[x]... Might be more correct 10/14/14
  #       term.fit <- rep.gain[x] #x/cons.bs
  #       
  #       #Terminal fitness is independent of focal resource
  #       W.nr[[i]][x,tmax] <- term.fit
  #       
  #       #Or Terminal fitness is also dependent on mortality risk of focal resource
  #       #W.nr[[i]][x,tmax] <- term.fit * (1 - mort[i])
  #       
  #     }
  #   }
  
#   int tmax, NumericVector res_bs, int cons_bs, int xc, NumericVector rep_gain, 
#   NumericMatrix f_m, NumericVector mort, List dec_ls, NumericVector rho_vec, 
#   NumericMatrix c_learn, NumericVector g_forage, NumericVector c_forage
  
  output <- list()
  
  output[[1]] <- res.bs
  output[[2]] <- xc
  output[[3]] <- rep.gain
  output[[4]] <- f.m
  output[[5]] <- mort
  output[[6]] <- dec.ls
  output[[7]] <- rho.vec
  output[[8]] <- c.learn
  output[[9]] <- g.forage
  output[[10]] <- c.forage
  output[[11]] <- eta

  
  return(output)
  
} #End function
