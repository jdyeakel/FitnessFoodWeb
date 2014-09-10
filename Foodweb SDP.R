
library(RColorBrewer)

###################
# Community Metrics
###################

#Establish Species richness
sp_rich <- 10

#Establish species body sizes from a body size distribution
bs_mu <- 50 # Body size mean
bs_sd <- 10 # Body size standard deviation
sp_bs <- round(rnorm(sp_rich,bs_mu,bs_sd),0)

#Define interaction probability for each species
int_matrix <- matrix(0,sp_rich,sp_rich)
a1 <- 1.41
a2 <- 3.73
a3 <- -1.87
for (i in 1:sp_rich) {
  for (j in 1:sp_rich) {
    p <- exp(a1 + a2*log(sp_bs[i]/sp_bs[j]) + a3*log(sp_bs[i]/sp_bs[j])^2)
    int_matrix[i,j] <- p / (1 + p)
  }
}

#Establish potential foods for each species
pot_foods <- list()
cutoff <- 0.7 #probability cutoff for determining potential foods
for (i in 1:sp_rich) {
  pot_foods[[i]] <- which(int_matrix[i,] > cutoff)
}

#Define state matrices for each species
tmax_sp <- rep(10,sp_rich)
state_matrix <- list()
for (i in 1:sp_rich) {
  state_matrix[[i]] <- matrix(0,sp_bs[i],tmax_sp[i])
}
#Define critical threshold for each species
xc <- round((5/8)*sp_bs,0)

#Define terminal fitness function
#Currently just a increasing function with respect to x... NEED REFINED
for (i in 1:sp_rich) {
  for (j in xc[i]:sp_bs[i]) {
    state_matrix[[i]][j,tmax_sp[i]] <- j/sp_bs[i]
  }
}

####################################
# Define Properties of the resources
####################################

#Mean encounter rate
xi <- 1/sp_bs * 100 #The bigger the species, the more rare it is
#Dispersion parameter... Larger = less patchy. Smaller = more patchy.
nu <- 1/sp_bs * 100 #The bigger the species, the more patchy it is... currently unbounded
max_enc <- 20
f_vector <- matrix(0,max_enc+1,sp_rich)
for (i in 1:sp_rich) {
  for (j in 0:max_enc) {
    f_vector[j,i] <- dnbinom(j,mu = xi[i], size = nu[i])
  }
  f_vector[,i] <- f_vector[,i] / sum(f_vector[,i])
}
# Plot probability that a consumer finds k items  of resource i
colors <- brewer.pal(10,"Spectral")
plot(seq(0,max_enc,1),f_vector[,1],type="l",ylim=c(0,max(f_vector)),col=colors[1],lwd=3)
for(i in 2:sp_rich) {
  lines(seq(0,max_enc,1),f_vector[,i],col=colors[i],lwd=3)
} 

#Consumer Resource mortalities (per-encounter mortalities)
mort <- matrix(0,sp_rich,sp_rich)
#What does probability of mortality equal when R/C ratio is 1?
mp1 <- 0.1
for (i in 1:sp_rich) {
  for (j in 1:sp_rich) {
    #Consumer i foraging on Resource j
    #Ratio of Resource to Consumer body size
    Ratio_RC <- sp_bs[j] / sp_bs[i]
    mort[j,i] <- 0.5 - 0.5*(1 - 2*mp1)^(Ratio_RC^2)
  }
}

#Foraging Gains and Costs (allometric and stoichiometric)
eta <- matrix(0,sp_rich,sp_rich)
g_forage <- matrix(0,sp_rich,sp_rich)
c_forage <- matrix(0,sp_rich,sp_rich)
c_rep <- numeric(sp_rich)
for (i in 1:sp_rich) {
  for (j in 1:sp_rich) {
    #Foraging costs/gains conditional on consumer AND resource
    eta[i,j] <- 0.5
    #Gains/costs scale with prob of interaction
    p <- exp(a1 + a2*log(sp_bs[i]/sp_bs[j]) + a3*log(sp_bs[i]/sp_bs[j])^2)
    g_forage[i,j] <- (p/(1+p)) * sp_bs[j]/20
    c_forage[i,j] <- (p/(1+p)) * sp_bs[j]/20
  }
  #Reproductive costs conditional only on consumer
  c_rep[i] <- 0.05*sp_bs[i]
}

##########################
# Build Fitness Matrices
##########################
W_x <- list()
istar <- list()
for (sp in 1:sp_rich) {
  W_x[[sp]] <- matrix(0,sp_bs[sp],tmax_sp[sp])
  istar[[sp]] <-  matrix(0,sp_bs[sp],(tmax_sp[sp]-1))
}




for (sp in 1:sp_rich) {
  #Maximum energetic reserves
  xmax <- sp_bs[sp]
  #Minimum energetic reserves
  xcrit <- xc[sp]
  #Stomach capacity (max amt consumed in one 'temporal window')
  xs <- 100
  tmax <- tmax_sp[sp]
  
  res <- pot_foods[[sp]]
  num_foods <- length(res)
  num_decisions <- num_foods*6 #number of foods * 2 (r/f) * 3 (low,med,high competition)
  
  decisions <- matrix(0,num_decisions,4)
  colnames(decisions) <- c("Resource","Comp","Rep","FValue")
  decisions[,1] <- rep(unlist(lapply(res,function(x){rep(x,3)})),2)
  decisions[,2] <- rep(rep(c(1,2,3),num_foods),2)
  decisions[,3] <- c(rep(0,num_decisions/2),rep(1,num_decisions/2))
  
  comp_vector <- rep(c(1,2,3),num_foods)
  rep_vector <- rep(c(0,1),num_foods)
  
  for (t in seq((tmax-1),1,-1)) {
    
    for (x in xcrit:xmax) {
      
      value_max <- -10
      
      # Iterate across decisions
      for (i in 1:num_decisions) {
        d_i <- decisions[i,]
        res_i <- as.numeric(d_i[1])
        chi <- as.numeric(d_i[2])
        rho <- as.numeric(d_i[3])
        
        xp <- numeric(max_enc+1)
        tic <- 0
        for (k in 0:max_enc) {
          tic <- tic + 1
          
          #ENERGETIC DYNAMICS
          xp[tic] <- x + min((eta[sp,res_i]*g_forage[sp,res_i]*k)/chi, xs) - 
            c_forage[sp,res_i]*chi - rho*c_rep[sp]
          
          #Determine Boundary Conditions
          
          #Define upper and lower for interpolation
          
          
          #Define k Fitness Value
          Wxp[k] <- 
          
        } # end k
        
        # Define Cumulative Fitness Value
        
      } #end i
      
      #Save decision list and fitness values (per x, per t)
      
    } #end x
  } #end t
} #end sp












