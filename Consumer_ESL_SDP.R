library(RColorBrewer)
library(Rcpp)
library(igraph)
library(beanplot)
colors <- brewer.pal(3,"Set1")

num.res <- 5
res.bs <- round(rgamma(num.res,shape=5,scale=2),0)
cons.bs <- max(res.bs + 5)

#Define state matrices for consumer
tmax <- 10
state.matrix <-  matrix(0,cons.bs,tmax)

#Define critical threshold for consumer
xc <- round((5/8)*cons.bs,0)

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
xi <- matrix(0,n,num.res)
nu <- matrix(0,n,num.res)
for (i in 1:n) {
  for (j in 1:num.res) {
    xi[i,j] <- (1/res.bs[j])*10*num.nn[i]
    #nu[i,j] <- (1/res.bs[j])*10*(nn2.sd[i]+2)
    nu[i,j] <- nn2.sd[i]*2
  }
}

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










