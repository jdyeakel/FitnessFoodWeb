library(RColorBrewer)
library(Rcpp)
library(igraph)
library(beanplot)
colors <- brewer.pal(10,"Spectral")

num.res <- 5
res.bs <- round(rgamma(num.res,shape=5,scale=2),0)
cons.bs <- max(res.bs + 5)

#Define state matrices for consumer
tmax <- 10
state.matrix <-  matrix(0,cons.bs,tmax)

#Define critical threshold for consumer
xc <- round((5/8)*cons.bs,0)
xmax <- cons.bs

#Define terminal fitness function
#Currently just a increasing function with respect to x... NEED REFINED
for (i in xc:cons.bs) {
    state.matrix[i,tmax] <- i/cons.bs
}

#Define Spatial landscape & attributes
#Nunmber of spatial patches
n <- 25
lattice.net <- graph.lattice(dimvector=c(sqrt(n),sqrt(n)),directed=FALSE,mutual=TRUE)
adj.m <- get.adjacency(lattice.net)
plot(lattice.net,vertex.size=degree(lattice.net),vertex.color=colors[1],edge.color="lightblue",edge.arrow.size=0.6)

#Find the nearest neighbors of each node
nn <- list()
for (j in 1:n) {
  arow <- adj.m[j,]
  nn[[j]] <- which(arow==1)
}
num.nn <- unlist(lapply(nn,length))
num.nn2 <- num.nn

#How many nearest neighbors does each nearest neighbor have? (minus the node of interest)
nn2 <- list()
for (i in 1:n) {
  nn2[[i]] <- unlist(lapply(nn[c(nn[[i]])],function(x){length(x)}))
}

#Variability in the local region
nn2.sd <- unlist(lapply(nn2,sd))

#Plot Graph where the vertex size is scaled to the 'local landscape variability"
plot(lattice.net,vertex.size=(nn2.sd+1)*3,vertex.label=NA,vertex.color=colors[1],edge.color="lightblue",edge.arrow.size=0.6)


#Mean encounter rates of prey in region i should be scaled to the number of nearest neighbors to i
#Dispersion of prey in region i should be scaled to the variability in nearest neighbors to nearest neighbors in i
#This is essentially the "Local Heterogeneity"
xi <- matrix(0,n,num.res)
nu <- matrix(0,n,num.res)
for (i in 1:n) {
  for (j in 1:num.res) {
    xi[i,j] <- (1/res.bs[j])*10*num.nn[i]
    #nu[i,j] <- (1/res.bs[j])*10*(nn2.sd[i]+2)
    nu[i,j] <- nn2.sd[i]*2
  }
}
nu <- nu + 1 #Ensure that there are non-zero nu's

max.enc <- 20
f.patch <- vector("list",n)
f.res.patch <- lapply(f.patch,function(x){matrix(0,(max.enc+1),num.res)})
for (i in 1:n) {
  for (j in 1:num.res) {
    for (k in 1:(max.enc+1)) {
      f.res.patch[[i]][k,j] <-  dnbinom(k,mu = xi[i,j], size = nu[i,j])
    }
    f.res.patch[[i]][,j] <- f.res.patch[[i]][,j] / sum(f.res.patch[[i]][,j])
  }
}

#Consumer-resource mortality rates
Ratio.RC <- res.bs/cons.bs
mp1 <- 0.2
mort <- 0.5 - 0.5*(1 - 2*mp1)^(Ratio.RC^2)

#Foraging Gains and Costs (allometric and stoichiometric)
c.forage <- cons.bs^(1/4)
c.rep <- 0.05*cons.bs
eta <- numeric(num.res)
g.forage <- numeric(num.res)
for (i in 1:num.res) {
  #Foraging costs/gains conditional on consumer AND resource
  eta[i] <- 0.1
  #Gains/costs scale allometrically
  g.forage[i] <- eta[i]*res.bs[i]
}

#Build consumer fitness matrix
W.xt <- matrix(0,cons.bs,tmax)
istar <-  matrix(0,cons.bs,(tmax-1))

#Gradiant of decision possibilities based on resource similarities
num.dec <- 10
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
#stretch between 0 and 1 for each similarity column
for (i in 1:num.res) {
  res.sim.low <- res.sim[,i] - (min(res.sim[,i]))
  res.sim[,i] <- res.sim.low/max(res.sim.low)
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
plot(dec.ls[[1]][,1],type="l",col=colors[1], ylim=c(0,1))
for (i in 2:10) {lines(dec.ls[[1]][,i],type="l",col=colors[i])}

#Compute fitness matrices
#node
#time
#current resource
#energetic state

for (node in 1:n) {
  
  #Loop over time
  for (t in seq(tmax-1,1,-1)) {
    
    #Loop over 'current resource'
    for (r in 1:num.res) {
      
      #Loop over energetic state
      for (x in seq(xc,xmax,1)) {
        
        value_max <- -10
        
        
        #Loop over decisions
        for (i in 1:num.dec) {
          
          xp <- numeric(max.enc+1)
          
          
          
          
          
          
          
          
          
        }#end decision loop
        
        
        
      }#end energetic state loop
      
      
      
    }#end current resource loop
    
    
    
    
    
  }#end time loop
  
  
}#end node loop









