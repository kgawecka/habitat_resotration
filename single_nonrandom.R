# stochastoc cellular automata
# single species spatially explicit model
# non-random habitat resotration - patches adjacent to occupied patches resotred first

# patch values:
# 1 - occupied
# 0 - empty
# -1 - destroyed

rm(list=ls())


# HABITAT DESTRUCTION ----

# input - habitat destruction

pc <- 0.2         # probablity of colonisation
pe <- 0.2         # probability of extinction
p0 <- 0.1         # initial fraction of occupied sites
n <- 100          # grid size
tmax <- 50        # maximum number of timesteps
dD <- 0.01        # fraction of patches destoryed at a time


# dynamics with habitat destruction

dynamicsD <- function(pc, pe, p0, n, tmax, dD) {
  
  t <- 0:tmax
  pn0 <- round(p0*n*n, digits=0)   # initial number of occupied sites
  if(pn0<1) { return("ERROR: initial number of occupied sites = 0") }
  
  if(dD*n*n<1) { dn <- 1 } 
  else { dn <- round(dD*n*n, digits=0) }     # number of patches destroyed at a time
  d <- seq(0, (n*n), dn)    
  D <- d/(n*n)                               # fraction of habitat destoryed
  
  # intial grid occupancy
  xt <- array(0, dim=c(n, n, length(t)))
  for(i in 1:pn0){
    # select random patch
    i_o <- sample(1:n,1)
    j_o <- sample(1:n,1)
    while(xt[i_o,j_o,1]!=0){
      i_o <- sample(1:n,1)
      j_o <- sample(1:n,1)
    }
    xt[i_o,j_o,1] <- 1
  }
  xt0 <- xt[ , ,1]
  
  pdt <- matrix(0, nrow=length(d), ncol=length(t))    # fraction of occupied sites as function of time and D
  pdt[1,1] <- sum(xt0)/(n*n)
  Pdt <- matrix(0, nrow=length(d), ncol=length(t))    # fraction of occupied sites out of non-destroyed sites as function of time and D
  Pdt[1,1] <- sum(xt0)/(n*n)
  xdeq <- array(0, dim=c(n, n, length(d)))
  
  for(k in 1:length(d)){
    
    # habitat destruction only
    if(k>1){
      
      # initialise grid
      xdeq[ , ,k] <- xdeq[ , ,k-1]
      
      # select random patch
      for(h in 1:dn) {
        i_d <- sample(1:n,1)
        j_d <- sample(1:n,1)
        while(xdeq[i_d,j_d,k]==-1){
          i_d <- sample(1:n,1)
          j_d <- sample(1:n,1)
        }
        xdeq[i_d,j_d,k] <- -1
      }
      
      # grid at t=0 for d[k]
      xt[ , ,1] <- xdeq[ , ,k]
      pdt[k,1] <- (sum(xt[ , , 1])+(k-1)*dn)/(n*n)
      Pdt[k,1] <- (sum(xt[ , , 1])+(k-1)*dn)/(n*n-d[k])
    }
    
    # iterate until steady state
    for(g in 2:length(t)){
      for(i in 1:n){
        for(j in 1:n){
          
          # destoryed sites remain destroyed
          if(xt[i,j,g-1]==-1) { xt[i,j,g] <- -1 }
          
          # occupied site has a probablity of extinction = pe
          else if(xt[i,j,g-1]==1) { xt[i,j,g] <- sample(c(0,1),1, prob=c(pe, 1-pe)) }
          
          # empty site has probability of colonisation = pc, if a neighbour site is occupied
          else if(xt[i,j,g-1]==0) {
            
            # von neumann neighborhood
            neigh <- matrix(c(-1, 0, 1, 0, 0, 1, 0, -1), nrow=4, ncol=2)
            for(h in 1:4){
              i_n <- i+neigh[h,1]
              j_n <- j+neigh[h,2]
              
              if(i_n < 1 || i_n > n || j_n < 1 || j_n > n) { next } # neighbour outside of grid
              
              if      (xt[i_n,j_n,g-1]==1)  { dx <- sample(c(0,1),1, prob=c(1-pc, pc)) }
              else if (xt[i_n,j_n,g-1]==0)  { dx <- 0 }
              else if (xt[i_n,j_n,g-1]==-1) { dx <- 0 }
              
              xt[i,j,g] <- xt[i,j,g-1] + dx
              if(xt[i,j,g]==1) { break }
            }
          }
        }
      }
      pdt[k,g] <- (sum(xt[ , , g])+(k-1)*dn)/(n*n)
      Pdt[k,g] <- (sum(xt[ , , g])+(k-1)*dn)/(n*n-d[k])
    }
    xdeq[ , ,k] <- xt[ , ,length(t)]
  }
  output <- list("D"=D, "pdt"=pdt, "Pdt"=Pdt, "xdeq"=xdeq)
  return(output)
}

outD <- dynamicsD(pc, pe, p0, n, tmax, dD)  # habitat destruction output

D <- outD$D           # fraction of habitat destoryed [0,1]
pdt <- outD$pdt       # fraction of occupied sites as function of time and D
Pdt <- outD$Pdt       # fraction of occupied sites out of non-destroyed sites as function of time and D
xdeq <- outD$xdeq     # equilibrium grid as function of D

plot(x=D, y=pdt[,tmax+1], type="b", col="black", xlim=c(0,1), ylim=c(0,1), xlab="D", ylab="p", pch = 20)
abline(h=0, v=0)


# HABITAT RESTORATION ----

# input - habitat restoration

pc <- 0.2         # probablity of colonisation
pe <- 0.2         # probability of extinction
tmax <- 50        # maximum number of timesteps
DR <- 0.60        # fraction of patches destroyed at the start of restoration
dR <- 0.01        # fraction of patches restored at a time

Re  <- rep(NA,length(D))     # initialise restoration efficiency vector

# initialise matrices for storing results
R_all <- matrix(NA, nrow=length(D), ncol=length(D))      # fraction of habitat loss during resotration
pr_all <- matrix(NA, nrow=length(D), ncol=length(D))     # fraction of occupied sites (mean of replicates)
Pr_all <- matrix(NA, nrow=length(D), ncol=length(D))     # fraction of occupied sites out of non-destroyed sites (mean of replicates)
sd_all <- matrix(NA, nrow=length(D), ncol=length(D))     # st dev of fraction of occupied sites of replicates


# dynamics with habitat restoration

dynamicsR <- function(pc, pe, tmax, DR, dR) {
  
  t <- 0:tmax
  
  for(i in 1:length(D)){
    if(D[i]==DR) { 
      Di <- i 
      break }
    if(i==length(D)) { return("ERROR: DR does not exist") }
  }
  
  if(dR*n*n<1) { dn <- 1 
  } else { dn <- round(dR*n*n, digits=0) }     # number of patches restored at a time
  r <- seq(DR*(n*n), 0, -dn)    
  R <- r/(n*n)                                 # fraction of habitat loss during resotration
  
  xreq <- array(0, dim=c(n, n, length(r)))
  xreq[ , , 1] <- xdeq[ , , Di]
  xt <- array(0, dim=c(n, n, length(t)))
  
  prt <- matrix(0, nrow=length(r), ncol=length(t))    # fraction of occupied sites as function of time and R
  prt[1, ] <- pdt[Di, ]
  Prt <- matrix(0, nrow=length(r), ncol=length(t))    # fraction of occupied sites out of non-destroyed sites as function of time and R
  Prt[1, ] <- Pdt[Di, ]
  
  for(k in 2:length(r)){
    
    # initialise grid
    xreq[ , ,k] <- xreq[ , ,k-1]
    
    n_occ <- 0  # number of destroyed patches with occupied neighbours
    n_noc <- 0  # number of destroyed patches without occupied neighbours
    d_occ <- matrix(NA,nrow=sum(xreq[ , ,k]==-1),ncol=2) # matrix with coordinates of destroyed patches with occupied neighbours
    d_noc <- matrix(NA,nrow=sum(xreq[ , ,k]==-1),ncol=2) # matrix with coordinates of destroyed patches without occupied neighbours
    
    for(i in 1:n){
      for(j in 1:n){
        
        if(xreq[i,j,k]!=-1) { next }
        
        else{
          # check if any adjecent site is occupied
          # von neumann neighborhood
          n1 <- 0
          neigh <- matrix(c(-1, 0, 1, 0, 0, 1, 0, -1), nrow=4, ncol=2)
          for(h in 1:4){
            i_n <- i+neigh[h,1]
            j_n <- j+neigh[h,2]
            
            if(i_n < 1 || i_n > n || j_n < 1 || j_n > n) { next } # neighbour outside of grid
            if(xreq[i_n,j_n,k]==1) {
              n_occ <- n_occ+1
              d_occ[n_occ,1] <- i
              d_occ[n_occ,2] <- j
              n1 <- 1
              break }               # if neighbour occupied, move to next patch
          }
          if(n1==0){
            n_noc <- n_noc+1
            d_noc[n_noc,1] <- i
            d_noc[n_noc,2] <- j
          }
        }
      }
    }
    d_occ <- na.omit(d_occ)
    d_noc <- na.omit(d_noc)
    
    # restore patches adjecent to occupied patches only
    if(n_occ>=dn) {
      rp <- sample(1:n_occ,dn)
      for(h in 1:dn) {
        i_d <- d_occ[rp[h],1]
        j_d <- d_occ[rp[h],2]
        xreq[i_d,j_d,k] <- 0
      }
    }
    
    # restore patches adjecent to occupied patches first
    if(n_occ<dn) {
      if(n_occ>0){
        for(h in 1:n_occ) {
          i_d <- d_occ[h,1]
          j_d <- d_occ[h,2]
          xreq[i_d,j_d,k] <- 0
        }
      }
      rp <- sample(1:n_noc,(dn-n_occ))
      for(h in 1:(dn-n_occ)) {
        i_d <- d_noc[rp[h],1]
        j_d <- d_noc[rp[h],2]
        xreq[i_d,j_d,k] <- 0
      }
    }

    # grid at t=0 for d[k]
    xt[ , ,1] <- xreq[ , ,k]
    pdt[k,1] <- sum(xreq[ , , 1]==1)/(n*n)
    Pdt[k,1] <- sum(xreq[ , , 1]==1)/(n*n-r[k])
    
    # iterate until steady state
    for(g in 2:length(t)){
      for(i in 1:n){
        for(j in 1:n){
          
          # destoryed sites remain destroyed
          if(xt[i,j,g-1]==-1) { xt[i,j,g] <- -1 }
          
          # occupied site has a probablity of extinction = pe
          else if(xt[i,j,g-1]==1) { xt[i,j,g] <- sample(c(0,1),1, prob=c(pe, 1-pe)) }
          
          # empty site has probability of colonisation = pc, if a neighbour site is occupied
          else if(xt[i,j,g-1]==0) {
            
            # von neumann neighborhood
            neigh <- matrix(c(-1, 0, 1, 0, 0, 1, 0, -1), nrow=4, ncol=2)
            for(h in 1:4){
              i_n <- i+neigh[h,1]
              j_n <- j+neigh[h,2]
              
              if(i_n < 1 || i_n > n || j_n < 1 || j_n > n) { next } # neighbour outside of grid
              
              if      (xt[i_n,j_n,g-1]==1)  { dx <- sample(c(0,1),1, prob=c(1-pc, pc)) }
              else if (xt[i_n,j_n,g-1]==0)  { dx <- 0 }
              else if (xt[i_n,j_n,g-1]==-1) { dx <- 0 }
              
              xt[i,j,g] <- xt[i,j,g-1] + dx
              if(xt[i,j,g]==1) { break }
            }
          }
        }
      }
      prt[k,g] <- sum(xt[ , , g]==1)/(n*n)
      Prt[k,g] <- sum(xt[ , , g]==1)/(n*n-r[k])
    }
    xreq[ , ,k] <- xt[ , ,length(t)]
  }
  output <- list("R"=R, "prt"=prt, "Prt"=Prt, "xreq"=xreq)
  return(output)
}

outR <- dynamicsR(pc, pe, tmax, DR, dR)   # habitat restoration output

R <- outR$R                # fraction of habitat loss [DR,0]
prt <- outR$prt            # fraction of occupied sites as function of time and D
Prt <- outR$Prt            # fraction of occupied sites out of non-destroyed sites as function of time and D
xreq <- outR$xreq          # equilibrium grid as function of D

for(k in 1:length(R)){
  image(x=1:n, y=1:n, z=t(apply(xreq[ , ,k], 2, rev)), main=paste("D =", R[k]), xlab="", ylab="", zlim=c(-1,1),
        col=gray.colors(n=3, start = 0, end = 1, rev = FALSE))
  legend(grconvertX(1, "device"), grconvertY(1, "device"), 
         c("occupied", "empty", "destroyed"), fill = gray.colors(n=3, start = 0, end = 1, rev = TRUE), xpd = NA)
}

pr <- prt[, tmax+1]   # equilibrium values only
Pr <- Prt[, tmax+1]
pd <- pdt[, tmax+1]
Pd <- Pdt[, tmax+1]

plot(x=D, y=pd, type="b", col="black", xlim=c(0,1), ylim=c(0,1), xlab="D", ylab="p", pch = 20)
lines(x=R, y=pr, type="b", col="green", pch = 20)
abline(h=0, v=0)


# replicates
nr <- 5        # no of replicates
DR <- 0.10     # fraction of patches destroyed at the start of restoration

R <- seq(DR, 0, -dR)
pr <- matrix(0, nrow=nr, ncol=length(R))
Pr <- matrix(0, nrow=nr, ncol=length(R))
for(i in 1:nr){
  outR <- dynamicsR(pc, pe, tmax, DR, dR)
  prt <- outR$prt
  Prt <- outR$Prt
  pr[i, ] <- prt[ ,tmax+1]
  Pr[i, ] <- Prt[ ,tmax+1]
}

# mean and st dev of replicates
pr_mean <- numeric(length(R))
Pr_mean <- numeric(length(R))
pr_sd <- numeric(length(R))
Pr_sd <- numeric(length(R))
for(i in 1:length(R)) {
  pr_mean[i] <- mean(pr[ , i])
  Pr_mean[i] <- mean(Pr[ , i])
  pr_sd[i] <- sd(pr[ , i])
  Pr_sd[i] <- sd(Pr[ , i])
}

plot(x=D, y=pd, type="b", col="grey", xlim=c(0,1), ylim=c(0,1), xlab="D", ylab="p", pch = 20, main=paste("DR =", DR))
lines(x=R, y=pr_mean, type="p", col="green", pch = 20)
abline(h=0, v=0)
arrows(x0=R, y0=pr_mean-pr_sd, x1=R, y1=pr_mean+pr_sd, length=0.05, angle=90, code=3)


# calculate restoration efficiency

# area below curves
for(i in 1:length(D)){
  if(D[i]==DR) { 
    Di <- i 
    break }
}
AD <- 0
AR <- 0
for(i in 2:Di){
  AD <- AD + 0.5*(pd[i-1]+pd[i])*dD
  AR <- AR + 0.5*(pr_mean[i-1]+pr_mean[i])*dR
}

# restoration efficiency
Re[Di] <- -(AD-AR)/AD

plot(x=D, y=Re, type="b", col="black", xlim=c(0,1), xlab="DR", ylab="Re", pch = 19)
abline(h=0, v=0)


# store results for each DR
R_all[Di,1:length(R)] <- R
pr_all[Di,1:length(pr_mean)] <- pr_mean
Pr_all[Di,1:length(Pr_mean)] <- Pr_mean
sd_all[Di,1:length(pr_sd)] <- pr_sd
