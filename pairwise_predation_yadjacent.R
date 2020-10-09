# stochastoc cellular automata
# pairwise predation spatially explicit model
# non-random habitat resotration - y-adjacent method
#                                - sites adjecent to occupied by predator restored first, 
#                                  then sites adjecent to occupied by prey restored, 
#                                  then empty restored

# species:
# x - prey
# y - predator

# patch values:
# 0 - empty
# 1 - occupied by prey (x) only
# 2 - occupied by predator (y) only
# 3 - occupied by prey (x) and predator (y)
# -1 - destroyed

rm(list=ls())


# HABITAT DESTRUCTION ----

# input - habitat destruction
cx <- 0.4         # probablity of colonisation - prey
ex <- 0.2         # probability of extinction - prey
cy <- 0.4         # probablity of colonisation - predator
ey <- 0.2         # probability of extinction - predator
mu <- 0.5         # probablilty of prey extinction due to predator presence
psi <- 0.8        # probablilty of predator extinction due lack of prey
px0 <- 0.5        # initial fraction of occupied sites - prey
py0 <- 0.5        # initial fraction of occupied sites - predator
n <- 100          # grid size
tmax <- 50        # maximum number of timesteps
dD <- 0.01        # fraction of patches destoryed at a time


# dynamics with habitat destruction
dynamicsD <- function(cx, cy, ex, ey, mu, psi, px0, py0, n, tmax, dD) {
  
  t <- 0:tmax
  pnx0 <- round(px0*n*n, digits=0)   # initial number of occupied sites - prey
  if(pnx0<1) { return("ERROR: initial number of prey occupied sites = 0") }
  pny0 <- round(py0*n*n, digits=0)   # initial number of occupied sites - predator
  if(pny0<1) { return("ERROR: initial number of predator occupied sites = 0") }
  
  if(dD*n*n<1) { dn <- 1 } 
  else { dn <- round(dD*n*n, digits=0) }     # number of patches destroyed at a time
  d <- seq(0, (n*n), dn)    
  D <- d/(n*n)                               # fraction of habitat destoryed
  
  # intial grid occupancy
  xt <- array(0, dim=c(n, n, length(t)))
  # prey
  for(i in 1:pnx0){
    # select random patch
    i_o <- sample(1:n,1)
    j_o <- sample(1:n,1)
    while(xt[i_o,j_o,1]!=0){
      i_o <- sample(1:n,1)
      j_o <- sample(1:n,1)
    }
    xt[i_o,j_o,1] <- 1
  }
  # predator
  for(i in 1:pny0){
    # select random patch
    i_o <- sample(1:n,1)
    j_o <- sample(1:n,1)
    while(xt[i_o,j_o,1]==2 || xt[i_o,j_o,1]==3){
      i_o <- sample(1:n,1)
      j_o <- sample(1:n,1)
    }
    if     (xt[i_o,j_o,1]==0) { xt[i_o,j_o,1] <- 2 }
    else if(xt[i_o,j_o,1]==1) { xt[i_o,j_o,1] <- 3 }
  }
  
  pxdt <- matrix(0, nrow=length(d), ncol=length(t))    # fraction of occupied sites as function of time and D
  pxdt[1,1] <- (sum(xt[ , , 1]==1)+sum(xt[ , , 1]==3))/(n*n)
  Pxdt <- matrix(0, nrow=length(d), ncol=length(t))    # fraction of occupied sites out of non-destroyed sites as function of time and D
  Pxdt[1,1] <- (sum(xt[ , , 1]==1)+sum(xt[ , , 1]==3))/(n*n)
  pydt <- matrix(0, nrow=length(d), ncol=length(t))
  pydt[1,1] <- (sum(xt[ , , 1]==2)+sum(xt[ , , 1]==3))/(n*n)
  Pydt <- matrix(0, nrow=length(d), ncol=length(t))
  Pydt[1,1] <- (sum(xt[ , , 1]==2)+sum(xt[ , , 1]==3))/(n*n)
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
      pxdt[k,1] <- (sum(xt[ , , 1]==1)+sum(xt[ , , 1]==3))/(n*n)
      Pxdt[k,1] <- (sum(xt[ , , 1]==1)+sum(xt[ , , 1]==3))/(n*n-d[k])
      pydt[k,1] <- (sum(xt[ , , 1]==2)+sum(xt[ , , 1]==3))/(n*n)
      Pydt[k,1] <- (sum(xt[ , , 1]==2)+sum(xt[ , , 1]==3))/(n*n-d[k])
    }
    
    # iterate until steady state
    for(g in 2:length(t)){
      
      for(i in 1:n){
        for(j in 1:n){
          
          # von neumann neighborhood
          nx <- 0    # number of neighbours with x
          ny <- 0    # number of neighbours with y
          neigh <- matrix(c(-1, 0, 1, 0, 0, 1, 0, -1), nrow=4, ncol=2)
          for(h in 1:4){
            i_n <- i+neigh[h,1]
            j_n <- j+neigh[h,2]
            
            if(i_n < 1 || i_n > n || j_n < 1 || j_n > n) { next } # neighbour outside of grid
            
            if (xt[i_n,j_n,g-1]==1 || xt[i_n,j_n,g-1]==3) { nx <- nx + 1 }
            if (xt[i_n,j_n,g-1]==2 || xt[i_n,j_n,g-1]==3) { ny <- ny + 1 }
          }
          
          # transition probabilities
          p00 <- (1-cx)^nx * (1-cy)^ny
          p01 <- (1-(1-cx)^nx) * (1-cy)^ny
          p02 <- (1-cx)^nx * (1-(1-cy)^ny)
          p03 <- (1-(1-cx)^nx) * (1-(1-cy)^ny)
          
          p10 <- ex * (1-cy)^ny
          p11 <- (1-ex) * (1-cy)^ny
          p12 <- ex * (1-(1-cy)^ny)
          p13 <- (1-ex) * (1-(1-cy)^ny)
          
          p20 <- (1-cx)^nx * (ey+psi)
          p21 <- (1-(1-cx)^nx) * ey
          p22 <- (1-cx)^nx * (1-ey-psi)
          p23 <- (1-(1-cx)^nx) * (1-ey)
          
          p30 <- (ex+mu) * (ey+psi)
          p31 <- (1-ex-mu) * ey
          p32 <- (ex+mu) * (1-ey-psi)
          p33 <- (1-ex-mu) * (1-ey)
          
          # transitions
          if(xt[i,j,g-1]==-1) { xt[i,j,g] <- -1 }    # destroyed sites remain destroyed
          
          else if(xt[i,j,g-1]==0) { xt[i,j,g] <- sample(c(0,1,2,3),1, prob=c(p00, p01, p02, p03)) }
          
          else if(xt[i,j,g-1]==1) { xt[i,j,g] <- sample(c(0,1,2,3),1, prob=c(p10, p11, p12, p13)) }
          
          else if(xt[i,j,g-1]==2) { xt[i,j,g] <- sample(c(0,1,2,3),1, prob=c(p20, p21, p22, p23)) }
          
          else if(xt[i,j,g-1]==3) { xt[i,j,g] <- sample(c(0,1,2,3),1, prob=c(p30, p31, p32, p33)) }
        }
      }
      
      pxdt[k,g] <- (sum(xt[ , , g]==1)+sum(xt[ , , g]==3))/(n*n)
      Pxdt[k,g] <- (sum(xt[ , , g]==1)+sum(xt[ , , g]==3))/(n*n-d[k])
      pydt[k,g] <- (sum(xt[ , , g]==2)+sum(xt[ , , g]==3))/(n*n)
      Pydt[k,g] <- (sum(xt[ , , g]==2)+sum(xt[ , , g]==3))/(n*n-d[k])
    }
    xdeq[ , ,k] <- xt[ , ,length(t)]
  }
  output <- list("D"=D, "pxdt"=pxdt, "Pxdt"=Pxdt, "pydt"=pydt, "Pydt"=Pydt, "xdeq"=xdeq)
  return(output)
}

outD <- dynamicsD(cx, cy, ex, ey, mu, psi, px0, py0, n, tmax, dD)  # habitat destruction output

D <- outD$D           # fraction of habitat destoryed [0,1]
pxdt <- outD$pxdt     # fraction of x-occupied sites as function of time and D
Pxdt <- outD$Pxdt     # fraction of x-occupied sites out of non-destroyed sites as function of time and D
pydt <- outD$pydt     # fraction of y-occupied sites as function of time and D
Pydt <- outD$Pydt     # fraction of y-occupied sites out of non-destroyed sites as function of time and D
xdeq <- outD$xdeq     # equilibrium grid as function of D

plot(x=D, y=pxdt[ ,tmax+1], type="p", col="black", ylim=c(0,1), ylab="p", pch=20)
lines(x=D, y=pydt[ ,tmax+1], type="p", col="blue", pch=20)
legend(0.7, 1, legend=c("x - prey", "y - predator"), col=c("black", "blue"), pch=20, cex=1)
abline(h=0, v=0)


# HABITAT RESTORATION ----

# input - habitat restoration
cx <- 0.4         # probablity of colonisation - prey
ex <- 0.2         # probability of extinction - prey
cy <- 0.4         # probablity of colonisation - predator
ey <- 0.2         # probability of extinction - predator
mu <- 0.5         # probablilty of prey extinction due to predator presence
psi <- 0.8        # probablilty of predator extinction due lack of prey
tmax <- 50        # maximum number of timesteps
DR <- 0.50        # fraction of patches destroyed at the start of restoration
dR <- 0.01        # fraction of patches restored at a time

Rex  <- rep(NA,length(D))     # initialise restoration efficiency vector
Rey  <- rep(NA,length(D))

# initialise matrices for storing results
R_all <- matrix(NA, nrow=length(D), ncol=length(D))      # fraction of habitat loss during resotration
pxr_all <- matrix(NA, nrow=length(D), ncol=length(D))    # fraction of x-occupied sites (mean of replicates)
Pxr_all <- matrix(NA, nrow=length(D), ncol=length(D))    # fraction of x-occupied sites out of non-destroyed sites (mean of replicates)
pyr_all <- matrix(NA, nrow=length(D), ncol=length(D))    # fraction of y-occupied sites (mean of replicates)
Pyr_all <- matrix(NA, nrow=length(D), ncol=length(D))    # fraction of y-occupied sites out of non-destroyed sites (mean of replicates)
sdx_all <- matrix(NA, nrow=length(D), ncol=length(D))     # st dev of fraction of x-occupied sites of replicates
sdy_all <- matrix(NA, nrow=length(D), ncol=length(D))     # st dev of fraction of y-occupied sites of replicates


# dynamics with habitat restoration
dynamicsR <- function(cx, cy, ex, ey, mu, psi, tmax, DR, dR) {
  
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

  pxrt <- matrix(0, nrow=length(r), ncol=length(t))    # fraction of occupied sites as function of time and R
  pxrt[1, ] <- pxdt[Di, ]
  Pxrt <- matrix(0, nrow=length(r), ncol=length(t))    # fraction of occupied sites out of non-destroyed sites as function of time and R
  Pxrt[1, ] <- Pxdt[Di, ]
  pyrt <- matrix(0, nrow=length(r), ncol=length(t))    
  pyrt[1, ] <- pydt[Di, ]
  Pyrt <- matrix(0, nrow=length(r), ncol=length(t))
  Pyrt[1, ] <- Pydt[Di, ]
  
  for(k in 2:length(r)){
    
    # initialise grid
    xreq[ , ,k] <- xreq[ , ,k-1]
    
    n_x <- 0  # number of destroyed patches with x neighbours only
    n_y <- 0  # number of destroyed patches with y neighbours only
    n_e <- 0  # number of destroyed patches with empty neighbours only
    d_x <- matrix(NA,nrow=sum(xreq[ , ,k]==-1),ncol=2) # matrix with coordinates of destroyed patches with x neighbours
    d_y <- matrix(NA,nrow=sum(xreq[ , ,k]==-1),ncol=2) # matrix with coordinates of destroyed patches with y neighbours
    d_e <- matrix(NA,nrow=sum(xreq[ , ,k]==-1),ncol=2) # matrix with coordinates of destroyed patches with empty neighbours only

    for(i in 1:n){
      for(j in 1:n){
        
        if(xreq[i,j,k]!=-1) { next }
        
        else{
          # check if any adjecent site is occupied
          # von neumann neighborhood
          n1 <- 0
          n2 <- 0
          neigh <- matrix(c(-1, 0, 1, 0, 0, 1, 0, -1), nrow=4, ncol=2)
          for(h in 1:4){
            i_n <- i+neigh[h,1]
            j_n <- j+neigh[h,2]
            
            if(i_n < 1 || i_n > n || j_n < 1 || j_n > n) { next } # neighbour outside of grid
            if(xreq[i_n,j_n,k]==2 || xreq[i_n,j_n,k]==3) {
              n_y <- n_y+1
              d_y[n_y,1] <- i
              d_y[n_y,2] <- j
              n1 <- 1
              break }               # if neighbour occupied by y, move to next patch
            else if(xreq[i_n,j_n,k]==1) {
              n2 <- 1 }
          }
          if(n1==0 && n2==1) {
            n_x <- n_x+1
            d_x[n_x,1] <- i
            d_x[n_x,2] <- j
          }
          else if(n1==0 && n2==0){
            n_e <- n_e+1
            d_e[n_e,1] <- i
            d_e[n_e,2] <- j
          }
        }
      }
    }
    d_x <- na.omit(d_x)
    d_y <- na.omit(d_y)
    d_e <- na.omit(d_e)

    # restore patches adjecent to y-occupied patches
    if(n_y>=dn) {
      rp <- sample(1:n_y,dn)
      for(h in 1:dn) {
        i_d <- d_y[rp[h],1]
        j_d <- d_y[rp[h],2]
        xreq[i_d,j_d,k] <- 0
      }
    }
    
    # restore patches adjecent to y-occupied patches first, then x-occupied, then empty
    if(n_y<dn) {
      if(n_y>0) {
        for(h in 1:n_y) {
          i_d <- d_y[h,1]
          j_d <- d_y[h,2]
          xreq[i_d,j_d,k] <- 0
        }
      }
      if(n_x>=(dn-n_y)) {
        rp <- sample(1:n_x,(dn-n_y))
        for(h in 1:(dn-n_y)) {
          i_d <- d_x[rp[h],1]
          j_d <- d_x[rp[h],2]
          xreq[i_d,j_d,k] <- 0
        }
      }
      else if(n_x<(dn-n_y)) {
        if(n_x>0) {
          for(h in 1:n_x) {
            i_d <- d_x[h,1]
            j_d <- d_x[h,2]
            xreq[i_d,j_d,k] <- 0
          }
        }
        rp <- sample(1:n_e,(dn-n_y-n_x))
        for(h in 1:(dn-n_y-n_x)) {
          i_d <- d_e[rp[h],1]
          j_d <- d_e[rp[h],2]
          xreq[i_d,j_d,k] <- 0
        }
      }
    }
    
    # grid at t=0 for d[k]
    xt[ , ,1] <- xreq[ , ,k]
    pxrt[k,1] <- (sum(xt[ , , 1]==1)+sum(xt[ , , 1]==3))/(n*n)
    Pxrt[k,1] <- (sum(xt[ , , 1]==1)+sum(xt[ , , 1]==3))/(n*n-r[k])
    pyrt[k,1] <- (sum(xt[ , , 1]==2)+sum(xt[ , , 1]==3))/(n*n)
    Pyrt[k,1] <- (sum(xt[ , , 1]==2)+sum(xt[ , , 1]==3))/(n*n-r[k])
    
    # iterate until steady state
    for(g in 2:length(t)){
      for(i in 1:n){
        for(j in 1:n){
          
          # von neumann neighborhood
          nx <- 0    # number of neighbours with x
          ny <- 0    # number of neighbours with y
          neigh <- matrix(c(-1, 0, 1, 0, 0, 1, 0, -1), nrow=4, ncol=2)
          for(h in 1:4){
            i_n <- i+neigh[h,1]
            j_n <- j+neigh[h,2]
            
            if(i_n < 1 || i_n > n || j_n < 1 || j_n > n) { next } # neighbour outside of grid
            
            if (xt[i_n,j_n,g-1]==1 || xt[i_n,j_n,g-1]==3) { nx <- nx + 1 }
            if (xt[i_n,j_n,g-1]==2 || xt[i_n,j_n,g-1]==3) { ny <- ny + 1 }
          }
          
          # transition probabilities
          p00 <- (1-cx)^nx * (1-cy)^ny
          p01 <- (1-(1-cx)^nx) * (1-cy)^ny
          p02 <- (1-cx)^nx * (1-(1-cy)^ny)
          p03 <- (1-(1-cx)^nx) * (1-(1-cy)^ny)
          
          p10 <- ex * (1-cy)^ny
          p11 <- (1-ex) * (1-cy)^ny
          p12 <- ex * (1-(1-cy)^ny)
          p13 <- (1-ex) * (1-(1-cy)^ny)
          
          p20 <- (1-cx)^nx * (ey+psi)
          p21 <- (1-(1-cx)^nx) * ey
          p22 <- (1-cx)^nx * (1-ey-psi)
          p23 <- (1-(1-cx)^nx) * (1-ey)
          
          p30 <- (ex+mu) * (ey+psi)
          p31 <- (1-ex-mu) * ey
          p32 <- (ex+mu) * (1-ey-psi)
          p33 <- (1-ex-mu) * (1-ey)
          
          # transitions
          if(xt[i,j,g-1]==-1) { xt[i,j,g] <- -1 }    # destroyed sites remain destroyed
          
          else if(xt[i,j,g-1]==0) { xt[i,j,g] <- sample(c(0,1,2,3),1, prob=c(p00, p01, p02, p03)) }
          
          else if(xt[i,j,g-1]==1) { xt[i,j,g] <- sample(c(0,1,2,3),1, prob=c(p10, p11, p12, p13)) }
          
          else if(xt[i,j,g-1]==2) { xt[i,j,g] <- sample(c(0,1,2,3),1, prob=c(p20, p21, p22, p23)) }
          
          else if(xt[i,j,g-1]==3) { xt[i,j,g] <- sample(c(0,1,2,3),1, prob=c(p30, p31, p32, p33)) }
        }
      }
      pxrt[k,g] <- (sum(xt[ , , g]==1)+sum(xt[ , , g]==3))/(n*n)
      Pxrt[k,g] <- (sum(xt[ , , g]==1)+sum(xt[ , , g]==3))/(n*n-r[k])
      pyrt[k,g] <- (sum(xt[ , , g]==2)+sum(xt[ , , g]==3))/(n*n)
      Pyrt[k,g] <- (sum(xt[ , , g]==2)+sum(xt[ , , g]==3))/(n*n-r[k])
    }
    xreq[ , ,k] <- xt[ , ,length(t)]
  }
  output <- list("R"=R, "pxrt"=pxrt, "Pxrt"=Pxrt, "pyrt"=pyrt, "Pyrt"=Pyrt, "xreq"=xreq)
  return(output)
}

outR <- dynamicsR(cx, cy, ex, ey, mu, psi, tmax, DR, dR)   # habitat restoration output

R <- outR$R                # fraction of habitat loss [DR,0]
pxrt <- outR$pxrt          # fraction of x-occupied sites as function of time and D
Pxrt <- outR$Pxrt          # fraction of x-occupied sites out of non-destroyed sites as function of time and D
pyrt <- outR$pyrt          # fraction of y-occupied sites as function of time and D
Pyrt <- outR$Pyrt          # fraction of y-occupied sites out of non-destroyed sites as function of time and D
xreq <- outR$xreq          # equilibrium grid as function of D

pxr <- pxrt[, tmax+1]      # equilibrium values only
Pxr <- Pxrt[, tmax+1]
pxd <- pxdt[, tmax+1]
Pxd <- Pxdt[, tmax+1]
pyr <- pyrt[, tmax+1]
Pyr <- Pyrt[, tmax+1]
pyd <- pydt[, tmax+1]
Pyd <- Pydt[, tmax+1]

plot(x=D, y=pxd, type="b", col="grey", xlim=c(0,1), ylim=c(0,1), xlab="D", ylab="p", pch = 20)
lines(x=R, y=pxr, type="b", col="black", pch = 20)
lines(x=D, y=pyd, type="b", col="grey")
lines(x=R, y=pyr, type="b", col="blue")
legend(0.7, 1, legend=c("x - prey", "y - predator"), col=c("black", "blue"), pch=c(20,1), cex=1)
abline(h=0, v=0)

for(k in 1:length(R)){
  image(x=1:n, y=1:n, z=t(apply(xreq[ , ,k], 2, rev)), main=paste("D =", R[k]), xlab="", ylab="", zlim=c(-1,3),
        col=gray.colors(n=5, start = 0, end = 1, rev = FALSE))
  legend(grconvertX(1, "device"), grconvertY(1, "device"), 
         c("destroyed", "empty", "occupied-x", "occupied-y", "occupied-x&y"), fill = gray.colors(n=5, start = 0, end = 1, rev = FALSE), xpd = NA)
}


# replicates
nr <- 5        # no of replicates
DR <- 0.10     # fraction of patches destroyed at the start of restoration

R <- seq(DR, 0, -dR)
pxr <- matrix(0, nrow=nr, ncol=length(R))
Pxr <- matrix(0, nrow=nr, ncol=length(R))
pyr <- matrix(0, nrow=nr, ncol=length(R))
Pyr <- matrix(0, nrow=nr, ncol=length(R))
for(i in 1:nr){
  outR <- dynamicsR(cx, cy, ex, ey, mu, psi, tmax, DR, dR)
  pxrt <- outR$pxrt
  Pxrt <- outR$Pxrt
  pxr[i, ] <- pxrt[ ,tmax+1]
  Pxr[i, ] <- Pxrt[ ,tmax+1]
  pyrt <- outR$pyrt
  Pyrt <- outR$Pyrt
  pyr[i, ] <- pyrt[ ,tmax+1]
  Pyr[i, ] <- Pyrt[ ,tmax+1]
}

# mean and st dev of replicates
pxr_mean <- numeric(length(R))
Pxr_mean <- numeric(length(R))
pxr_sd <- numeric(length(R))
Pxr_sd <- numeric(length(R))
pyr_mean <- numeric(length(R))
Pyr_mean <- numeric(length(R))
pyr_sd <- numeric(length(R))
Pyr_sd <- numeric(length(R))
for(i in 1:length(R)) {
  pxr_mean[i] <- mean(pxr[ , i])
  Pxr_mean[i] <- mean(Pxr[ , i])
  pxr_sd[i] <- sd(pxr[ , i])
  Pxr_sd[i] <- sd(Pxr[ , i])
  pyr_mean[i] <- mean(pyr[ , i])
  Pyr_mean[i] <- mean(Pyr[ , i])
  pyr_sd[i] <- sd(pyr[ , i])
  Pyr_sd[i] <- sd(Pyr[ , i])
}

plot(x=D, y=pxd, type="b", col="grey", xlim=c(0,1), ylim=c(0,1), xlab="D", ylab="p", pch = 20, main=paste("DR =", DR))
lines(x=R, y=pxr_mean, type="b", col="black", pch = 20)
lines(x=D, y=pyd, type="b", col="grey", pch = 20)
lines(x=R, y=pyr_mean, type="b", col="blue", pch = 20)
abline(h=0, v=0)
legend(0.7, 1, legend=c("x - prey", "y - predator"), col=c("black", "blue"), pch=20, cex=1)
arrows(x0=R, y0=pxr_mean-pxr_sd, x1=R, y1=pxr_mean+pxr_sd, length=0.05, angle=90, code=3, col="black")
arrows(x0=R, y0=pyr_mean-pyr_sd, x1=R, y1=pyr_mean+pyr_sd, length=0.05, angle=90, code=3, col="blue")


# calculate restoration efficiency

# area below curves
for(i in 1:length(D)){
  if(D[i]==DR) { 
    Di <- i 
    break }
}
ADx <- 0
ARx <- 0
ADy <- 0
ARy <- 0
for(i in 2:Di){
  ADx <- ADx + 0.5*(pxd[i-1]+pxd[i])*dD
  ARx <- ARx + 0.5*(pxr_mean[i-1]+pxr_mean[i])*dR
  ADy <- ADy + 0.5*(pyd[i-1]+pyd[i])*dD
  ARy <- ARy + 0.5*(pyr_mean[i-1]+pyr_mean[i])*dR
}

# restoration efficiency
Rex[Di] <- -(ADx-ARx)/ADx
Rey[Di] <- -(ADy-ARy)/ADy

plot(x=D, y=Rex, type="b", col="black", xlim=c(0,1), xlab="DR", ylab="Re", pch = 19)
lines(x=D, y=Rey, type="b", col="blue", pch = 19)
legend("bottomleft", c("x", "y"), col=c("black", "blue"), pch=19, cex=1)
abline(h=0, v=0)


# store results for each DR
R_all[Di,1:length(R)] <- R
pxr_all[Di,1:length(pxr_mean)] <- pxr_mean
Pxr_all[Di,1:length(Pxr_mean)] <- Pxr_mean
pyr_all[Di,1:length(pyr_mean)] <- pyr_mean
Pyr_all[Di,1:length(Pyr_mean)] <- Pyr_mean
sdx_all[Di,1:length(pxr_mean)] <- pxr_sd
sdy_all[Di,1:length(pyr_mean)] <- pyr_sd
