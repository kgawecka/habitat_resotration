# stochastoc cellular automata
# trophic APC spatially explicit model
# non-random habitat resotration - z-adjacent method
#                                - sites adjecent to occupied by predator restored first, 
#                                  then sites adjecent to occupied by consumer restored,
#                                  then sites adjecent to occupied by resource restored, 
#                                  then empty restored
# species:
# x - resource
# y - consumer
# z - predator

# patch values:
# -1 - destroyed
# 0 - empty
# 1 - x only
# 2 - y1 only
# 3 - y2 only
# 4 - z only
# 5 - x & y1
# 6 - x & y2
# 7 - x & z
# 8 - y1 & y2
# 9 - y1 & z
# 10 - y2 & z
# 11 - x & y1 & y2
# 12 - x & y1 & z
# 13 - x & y2 & z
# 14 - y1 & y2 & z
# 15 - x & y1 & y2 & z

rm(list=ls())
library(RColorBrewer)


# HABITAT DESTRUCTION ----

# input - habitat destruction
cx  <- 0.7         # probablity of colonisation - resource
cy1 <- 0.7         # probablity of colonisation - consumer 1
cy2 <- 0.9         # probablity of colonisation - consumer 2
cz  <- 0.7         # probablity of colonisation - predator
ex  <- 0.1         # probability of extinction - resource
ey1 <- 0.1         # probability of extinction - consumer 1
ey2 <- 0.1         # probability of extinction - consumer 2
ez  <- 0.1         # probability of extinction - predator

mu1 <- 0.45       # probablilty of resource extinction due to consumer 1 presence
mu2 <- 0.45       # probablilty of resource extinction due to consumer 2 presence
mu3 <- 0.45       # probablilty of consumer 1 extinction due to predator presence
mu4 <- 0.45       # probablilty of consumer 2 extinction due to predator presence
psi1 <- 0.4       # probablilty of consumer 1 extinction due lack of resource
psi2 <- 0.4       # probablilty of consumer 2 extinction due lack of resource
psi3 <- 0.4       # probablilty of predator extinction due lack of consumer 1
psi4 <- 0.4       # probablilty of predator extinction due lack of consumer 2

px0 <- 0.5        # initial fraction of occupied sites - resource
py10 <- 0.5       # initial fraction of occupied sites - consumer 1
py20 <- 0.5       # initial fraction of occupied sites - consumer 2
pz0 <- 0.5        # initial fraction of occupied sites - predator

n <- 100          # grid size
tmax <- 50        # maximum number of timesteps
dD <- 0.01        # fraction of patches destoryed at a time


# dynamics with habitat destruction

dynamicsD <- function(cx, cy1, cy2, cz, ex, ey1, ey2, ez, mu1, mu2, mu3, mu4, psi1, psi2, psi3, psi4, 
                      px0, py10, py20, pz0, n, tmax, dD) {
  
  t <- 0:tmax
  pnx0 <- round(px0*n*n, digits=0)   # initial number of occupied sites - resource
  if(pnx0<1) { return("ERROR: initial number of resource occupied sites = 0") }
  pny10 <- round(py10*n*n, digits=0)   # initial number of occupied sites - consumer 1
  if(pny10<1) { return("ERROR: initial number of consumer 1 occupied sites = 0") }
  pny20 <- round(py20*n*n, digits=0)   # initial number of occupied sites - consumer 2
  if(pny20<1) { return("ERROR: initial number of consumer 2 occupied sites = 0") }
  pnz0 <- round(pz0*n*n, digits=0)   # initial number of occupied sites - predator
  if(pnz0<1) { return("ERROR: initial number of predator occupied sites = 0") }
  
  if(dD*n*n<1) { dn <- 1 
  } else { dn <- round(dD*n*n, digits=0) }     # number of patches destroyed at a time
  d <- seq(0, (n*n), dn)    
  D <- d/(n*n)                                 # fraction of habitat destoryed
  
  # intial grid occupancy
  xt <- array(0, dim=c(n, n, length(t)))
  # resource
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
  # consumer 1
  for(i in 1:pny10){
    # select random patch
    i_o <- sample(1:n,1)
    j_o <- sample(1:n,1)
    while(xt[i_o,j_o,1]==2 || xt[i_o,j_o,1]==5){
      i_o <- sample(1:n,1)
      j_o <- sample(1:n,1)
    }
    if     (xt[i_o,j_o,1]==0) { xt[i_o,j_o,1] <- 2 }
    else if(xt[i_o,j_o,1]==1) { xt[i_o,j_o,1] <- 5 }
  }
  # consumer 2
  for(i in 1:pny20){
    # select random patch
    i_o <- sample(1:n,1)
    j_o <- sample(1:n,1)
    while(xt[i_o,j_o,1]==2 || xt[i_o,j_o,1]==3 || xt[i_o,j_o,1]==5 || xt[i_o,j_o,1]==6){
      i_o <- sample(1:n,1)
      j_o <- sample(1:n,1)
    }
    if     (xt[i_o,j_o,1]==0) { xt[i_o,j_o,1] <- 3 }
    else if(xt[i_o,j_o,1]==1) { xt[i_o,j_o,1] <- 6 }
  }
  # predator
  for(i in 1:pnz0){
    # select random patch
    i_o <- sample(1:n,1)
    j_o <- sample(1:n,1)
    while(xt[i_o,j_o,1]==4 || xt[i_o,j_o,1]==7 || xt[i_o,j_o,1]==9 || xt[i_o,j_o,1]==10 || xt[i_o,j_o,1]==12 || xt[i_o,j_o,1]==13){
      i_o <- sample(1:n,1)
      j_o <- sample(1:n,1)
    }
    if     (xt[i_o,j_o,1]==0) { xt[i_o,j_o,1] <- 4 }
    else if(xt[i_o,j_o,1]==1) { xt[i_o,j_o,1] <- 7 }
    else if(xt[i_o,j_o,1]==2) { xt[i_o,j_o,1] <- 9 }
    else if(xt[i_o,j_o,1]==3) { xt[i_o,j_o,1] <- 10 }
    else if(xt[i_o,j_o,1]==5) { xt[i_o,j_o,1] <- 12 }
    else if(xt[i_o,j_o,1]==6) { xt[i_o,j_o,1] <- 13 }
    
  }
  
  pxdt <- matrix(0, nrow=length(d), ncol=length(t))    # fraction of occupied sites as function of time and D
  pxdt[1,1] <- (sum(xt[ , , 1]==1)+sum(xt[ , , 1]==5)+sum(xt[ , , 1]==6)+sum(xt[ , , 1]==7)
                +sum(xt[ , , 1]==11)+sum(xt[ , , 1]==12)+sum(xt[ , , 1]==13)+sum(xt[ , , 1]==15))/(n*n)
  Pxdt <- matrix(0, nrow=length(d), ncol=length(t))    # fraction of occupied sites out of non-destroyed sites as function of time and D
  Pxdt[1,1] <- (sum(xt[ , , 1]==1)+sum(xt[ , , 1]==5)+sum(xt[ , , 1]==6)+sum(xt[ , , 1]==7)
                +sum(xt[ , , 1]==11)+sum(xt[ , , 1]==12)+sum(xt[ , , 1]==13)+sum(xt[ , , 1]==15))/(n*n)
  
  py1dt <- matrix(0, nrow=length(d), ncol=length(t))
  py1dt[1,1] <- (sum(xt[ , , 1]==2)+sum(xt[ , , 1]==5)+sum(xt[ , , 1]==8)+sum(xt[ , , 1]==9)
                 +sum(xt[ , , 1]==11)+sum(xt[ , , 1]==12)+sum(xt[ , , 1]==14)+sum(xt[ , , 1]==15))/(n*n)
  Py1dt <- matrix(0, nrow=length(d), ncol=length(t))
  Py1dt[1,1] <- (sum(xt[ , , 1]==2)+sum(xt[ , , 1]==5)+sum(xt[ , , 1]==8)+sum(xt[ , , 1]==9)
                 +sum(xt[ , , 1]==11)+sum(xt[ , , 1]==12)+sum(xt[ , , 1]==14)+sum(xt[ , , 1]==15))/(n*n)
  
  py2dt <- matrix(0, nrow=length(d), ncol=length(t))
  py2dt[1,1] <- (sum(xt[ , , 1]==3)+sum(xt[ , , 1]==6)+sum(xt[ , , 1]==8)+sum(xt[ , , 1]==10)
                 +sum(xt[ , , 1]==11)+sum(xt[ , , 1]==13)+sum(xt[ , , 1]==14)+sum(xt[ , , 1]==15))/(n*n)
  Py2dt <- matrix(0, nrow=length(d), ncol=length(t))
  Py2dt[1,1] <- (sum(xt[ , , 1]==3)+sum(xt[ , , 1]==6)+sum(xt[ , , 1]==8)+sum(xt[ , , 1]==10)
                 +sum(xt[ , , 1]==11)+sum(xt[ , , 1]==13)+sum(xt[ , , 1]==14)+sum(xt[ , , 1]==15))/(n*n)
  
  pzdt <- matrix(0, nrow=length(d), ncol=length(t))
  pzdt[1,1] <- (sum(xt[ , , 1]==4)+sum(xt[ , , 1]==7)+sum(xt[ , , 1]==9)+sum(xt[ , , 1]==10)
                +sum(xt[ , , 1]==12)+sum(xt[ , , 1]==13)+sum(xt[ , , 1]==14)+sum(xt[ , , 1]==15))/(n*n)
  Pzdt <- matrix(0, nrow=length(d), ncol=length(t))
  Pzdt[1,1] <- (sum(xt[ , , 1]==4)+sum(xt[ , , 1]==7)+sum(xt[ , , 1]==9)+sum(xt[ , , 1]==10)
                +sum(xt[ , , 1]==12)+sum(xt[ , , 1]==13)+sum(xt[ , , 1]==14)+sum(xt[ , , 1]==15))/(n*n)
  
  xdeq <- array(0, dim=c(n, n, length(d)))
  
  # calculate terms for transition probabilities
  EX1 <- c(ex, (1-ex), ex, ex, ex, (1-ex), (1-ex), (1-ex), 0, ex, ex, 0, (1-ex), (1-ex), 0, 0)
  EX2 <- c((ex+mu1), (1-ex-mu1), (ex+mu1), (ex+mu1), (ex+mu1), (1-ex-mu1), (1-ex-mu1), (1-ex-mu1), 
           0, (ex+mu1), (ex+mu1), 0, (1-ex-mu1), (1-ex-mu1), 0, 0)
  EX3 <- c((ex+mu2), (1-ex-mu2), (ex+mu2), (ex+mu2), (ex+mu2), (1-ex-mu2), (1-ex-mu2), (1-ex-mu2), 
           (ex+mu2), (ex+mu2), (ex+mu2), (1-ex-mu2), (1-ex-mu2), (1-ex-mu2), (ex+mu2), (1-ex-mu2))
  EX4 <- c((ex+mu1+mu2), (1-ex-mu1-mu2), (ex+mu1+mu2), (ex+mu1+mu2), (ex+mu1+mu2), (1-ex-mu1-mu2), (1-ex-mu1-mu2), (1-ex-mu1-mu2), 
           (ex+mu1+mu2), (ex+mu1+mu2), (ex+mu1+mu2), (1-ex-mu1-mu2), (1-ex-mu1-mu2), (1-ex-mu1-mu2), (ex+mu1+mu2), (1-ex-mu1-mu2))
  EY11 <- c((ey1+psi1), ey1, (1-ey1-psi1), (ey1+psi1), (ey1+psi1), (1-ey1), ey1, ey1, 
            (1-ey1-psi1), (1-ey1-psi1), (ey1+psi1), (1-ey1), (1-ey1), ey1, (1-ey1-psi1), (1-ey1))
  EY12 <- c((ey1+psi1+mu3), (ey1+mu3), (1-ey1-psi1-mu3), (ey1+psi1+mu3), (ey1+psi1+mu3), (1-ey1-mu3), (ey1+mu3), (ey1+mu3), 
            (1-ey1-psi1-mu3), (1-ey1-psi1-mu3), (ey1+psi1+mu3), (1-ey1-mu3), (1-ey1-mu3), (ey1+mu3), (1-ey1-psi1-mu3), (1-ey1-mu3))
  EY21 <- c((ey2+psi2), ey2, (ey2+psi2), (1-ey2-psi2), (ey2+psi2), ey2, (1-ey2), ey2, 
            (1-ey2-psi2), (ey2+psi2), (1-ey2-psi2), (1-ey2), ey2, (1-ey2), (1-ey2-psi2), (1-ey2))
  EY22 <- c((ey2+psi2+mu4), (ey2+mu4), (ey2+psi2+mu4), (1-ey2-psi2-mu4), (ey2+psi2+mu4), (ey2+mu4), (1-ey2-mu4), (ey2+mu4), 
            (1-ey2-psi2-mu4), (ey2+psi2+mu4), (1-ey2-psi2-mu4), (1-ey2-mu4), (ey2+mu4), (1-ey2-mu4), (1-ey2-psi2-mu4), (1-ey2-mu4))
  EZ  <- c((ez+psi3/2+psi4/2), (ez+psi3/2+psi4/2), (ez+psi4/2), (ez+psi3/2), 
           (1-ez-psi3/2-psi4/2), (ez+psi4/2), (ez+psi3/2), (1-ez-psi3/2-psi4/2), 
           ez, (1-ez-psi4/2), (1-ez-psi3/2), ez, 
           (1-ez-psi4/2), (1-ez-psi3/2), (1-ez), (1-ez))
  
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
      pxdt[k,1] <- (sum(xt[ , , 1]==1)+sum(xt[ , , 1]==5)+sum(xt[ , , 1]==6)+sum(xt[ , , 1]==7)
                    +sum(xt[ , , 1]==11)+sum(xt[ , , 1]==12)+sum(xt[ , , 1]==13)+sum(xt[ , , 1]==15))/(n*n)
      Pxdt[k,1] <- (sum(xt[ , , 1]==1)+sum(xt[ , , 1]==5)+sum(xt[ , , 1]==6)+sum(xt[ , , 1]==7)
                    +sum(xt[ , , 1]==11)+sum(xt[ , , 1]==12)+sum(xt[ , , 1]==13)+sum(xt[ , , 1]==15))/(n*n-d[k])
      py1dt[k,1] <- (sum(xt[ , , 1]==2)+sum(xt[ , , 1]==5)+sum(xt[ , , 1]==8)+sum(xt[ , , 1]==9)
                     +sum(xt[ , , 1]==11)+sum(xt[ , , 1]==12)+sum(xt[ , , 1]==14)+sum(xt[ , , 1]==15))/(n*n)
      Py1dt[k,1] <- (sum(xt[ , , 1]==2)+sum(xt[ , , 1]==5)+sum(xt[ , , 1]==8)+sum(xt[ , , 1]==9)
                     +sum(xt[ , , 1]==11)+sum(xt[ , , 1]==12)+sum(xt[ , , 1]==14)+sum(xt[ , , 1]==15))/(n*n-d[k])
      py2dt[k,1] <- (sum(xt[ , , 1]==3)+sum(xt[ , , 1]==6)+sum(xt[ , , 1]==8)+sum(xt[ , , 1]==10)
                     +sum(xt[ , , 1]==11)+sum(xt[ , , 1]==13)+sum(xt[ , , 1]==14)+sum(xt[ , , 1]==15))/(n*n)
      Py2dt[k,1] <- (sum(xt[ , , 1]==3)+sum(xt[ , , 1]==6)+sum(xt[ , , 1]==8)+sum(xt[ , , 1]==10)
                     +sum(xt[ , , 1]==11)+sum(xt[ , , 1]==13)+sum(xt[ , , 1]==14)+sum(xt[ , , 1]==15))/(n*n-d[k])
      pzdt[k,1] <- (sum(xt[ , , 1]==4)+sum(xt[ , , 1]==7)+sum(xt[ , , 1]==9)+sum(xt[ , , 1]==10)
                    +sum(xt[ , , 1]==12)+sum(xt[ , , 1]==13)+sum(xt[ , , 1]==14)+sum(xt[ , , 1]==15))/(n*n)
      Pzdt[k,1] <- (sum(xt[ , , 1]==4)+sum(xt[ , , 1]==7)+sum(xt[ , , 1]==9)+sum(xt[ , , 1]==10)
                    +sum(xt[ , , 1]==12)+sum(xt[ , , 1]==13)+sum(xt[ , , 1]==14)+sum(xt[ , , 1]==15))/(n*n-d[k])
    }

    # iterate until steady state
    for(g in 2:length(t)){
      
      for(i in 1:n){
        for(j in 1:n){
          
          # von neumann neighborhood
          nx  <- 0    # number of neighbours with x
          ny1 <- 0    # number of neighbours with y1
          ny2 <- 0    # number of neighbours with y2
          nz  <- 0    # number of neighbours with z
          neigh <- matrix(c(-1, 0, 1, 0, 0, 1, 0, -1), nrow=4, ncol=2)
          for(h in 1:4){
            i_n <- i+neigh[h,1]
            j_n <- j+neigh[h,2]
            
            if(i_n < 1 || i_n > n || j_n < 1 || j_n > n) { next } # neighbour outside of grid
            
            if (xt[i_n,j_n,g-1]==1 || xt[i_n,j_n,g-1]==5 || xt[i_n,j_n,g-1]==6 || xt[i_n,j_n,g-1]==7 || 
                xt[i_n,j_n,g-1]==11 || xt[i_n,j_n,g-1]==12 || xt[i_n,j_n,g-1]==13 || xt[i_n,j_n,g-1]==15) { nx <- nx + 1 }
            if (xt[i_n,j_n,g-1]==2 || xt[i_n,j_n,g-1]==5 || xt[i_n,j_n,g-1]==8 || xt[i_n,j_n,g-1]==9 || 
                xt[i_n,j_n,g-1]==11 || xt[i_n,j_n,g-1]==12 || xt[i_n,j_n,g-1]==14 || xt[i_n,j_n,g-1]==15) { ny1 <- ny1 + 1 }
            if (xt[i_n,j_n,g-1]==3 || xt[i_n,j_n,g-1]==6 || xt[i_n,j_n,g-1]==8 || xt[i_n,j_n,g-1]==10 || 
                xt[i_n,j_n,g-1]==11 || xt[i_n,j_n,g-1]==13 || xt[i_n,j_n,g-1]==14 || xt[i_n,j_n,g-1]==15) { ny2 <- ny2 + 1 }
            if (xt[i_n,j_n,g-1]==4 || xt[i_n,j_n,g-1]==7 || xt[i_n,j_n,g-1]==9 || xt[i_n,j_n,g-1]==10 || 
                xt[i_n,j_n,g-1]==12 || xt[i_n,j_n,g-1]==13 || xt[i_n,j_n,g-1]==14 || xt[i_n,j_n,g-1]==15) { nz <- nz + 1 }
          }
          
          # calculate terms for transition probabilities
          CX <- c((1-cx)^nx, (1-(1-cx)^nx), (1-cx)^nx, (1-cx)^nx, (1-cx)^nx, (1-(1-cx)^nx), (1-(1-cx)^nx), (1-(1-cx)^nx), 
                  (1-cx)^nx, (1-cx)^nx, (1-cx)^nx, (1-(1-cx)^nx), (1-(1-cx)^nx), (1-(1-cx)^nx), (1-cx)^nx, (1-(1-cx)^nx))
          CY1 <- c((1-cy1)^ny1, (1-cy1)^ny1, (1-(1-cy1)^ny1), (1-cy1)^ny1, 
                   (1-cy1)^ny1, (1-(1-cy1)^ny1), (1-cy1)^ny1, (1-cy1)^ny1, 
                   (1-(1-cy1)^ny1), (1-(1-cy1)^ny1), (1-cy1)^ny1, (1-(1-cy1)^ny1), 
                   (1-(1-cy1)^ny1), (1-cy1)^ny1, (1-(1-cy1)^ny1), (1-(1-cy1)^ny1))
          CY2 <- c((1-cy2)^ny2, (1-cy2)^ny2, 1, (1-(1-cy2)^ny2), (1-cy2)^ny2, 1, (1-(1-cy2)^ny2), (1-cy2)^ny2, 
                   0, 1, (1-(1-cy2)^ny2), 0, 1, (1-(1-cy2)^ny2), 0, 0)
          CZ <- c((1-cz)^nz, (1-cz)^nz, (1-cz)^nz, (1-cz)^nz, (1-(1-cz)^nz), (1-cz)^nz, (1-cz)^nz, (1-(1-cz)^nz), 
                  (1-cz)^nz, (1-(1-cz)^nz), (1-(1-cz)^nz), (1-cz)^nz, (1-(1-cz)^nz), (1-(1-cz)^nz), (1-(1-cz)^nz), (1-(1-cz)^nz))
          
          pt <- numeric(length=16)   # vector of transition probabilities
          
          # transitions
          if(xt[i,j,g-1]==-1) { xt[i,j,g] <- -1 }    # destroyed sites remain destroyed
          
          else if(xt[i,j,g-1]==0) { 
            
            for(h in 1:16){ pt[h] <- CX[h] * CY1[h] * CY2[h] * CZ[h] }  # transition probabilities from state 0 to state x
            
            xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), 1, 
                                prob=c(pt[1],pt[2],pt[3],pt[4],pt[5],pt[6],pt[7],pt[8],pt[9],pt[10],pt[11],pt[12],pt[13],pt[14],pt[15],pt[16]))
          }
          
          else if(xt[i,j,g-1]==1) { 
            
            for(h in 1:16){ pt[h] <- EX1[h] * CY1[h] * CY2[h] * CZ[h] }  # transition probabilities from state 1 to state x
            
            xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), 1, 
                                prob=c(pt[1],pt[2],pt[3],pt[4],pt[5],pt[6],pt[7],pt[8],pt[9],pt[10],pt[11],pt[12],pt[13],pt[14],pt[15],pt[16]))
          }
          
          else if(xt[i,j,g-1]==2) { 
            
            for(h in 1:16){ pt[h] <- CX[h] * EY11[h] * CY2[h] * CZ[h] }  # transition probabilities from state 2 to state x
            
            xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), 1, 
                                prob=c(pt[1],pt[2],pt[3],pt[4],pt[5],pt[6],pt[7],pt[8],pt[9],pt[10],pt[11],pt[12],pt[13],pt[14],pt[15],pt[16]))
          }
          
          else if(xt[i,j,g-1]==3) { 
            
            for(h in 1:16){ pt[h] <- CX[h] * CY1[h] * EY21[h] * CZ[h] }  # transition probabilities from state 3 to state x
            
            xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), 1, 
                                prob=c(pt[1],pt[2],pt[3],pt[4],pt[5],pt[6],pt[7],pt[8],pt[9],pt[10],pt[11],pt[12],pt[13],pt[14],pt[15],pt[16]))
          }
          
          else if(xt[i,j,g-1]==4) { 
            
            for(h in 1:16){ pt[h] <- CX[h] * CY1[h] * CY2[h] * EZ[h] }  # transition probabilities from state 4 to state x
            
            xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), 1, 
                                prob=c(pt[1],pt[2],pt[3],pt[4],pt[5],pt[6],pt[7],pt[8],pt[9],pt[10],pt[11],pt[12],pt[13],pt[14],pt[15],pt[16]))
          }
          
          else if(xt[i,j,g-1]==5) { 
            
            for(h in 1:16){ pt[h] <- EX2[h] * EY11[h] * CY2[h] * CZ[h] }  # transition probabilities from state 5 to state x
            
            xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), 1, 
                                prob=c(pt[1],pt[2],pt[3],pt[4],pt[5],pt[6],pt[7],pt[8],pt[9],pt[10],pt[11],pt[12],pt[13],pt[14],pt[15],pt[16]))
          }
          
          else if(xt[i,j,g-1]==6) { 
            
            for(h in 1:16){ pt[h] <- EX3[h] * CY1[h] * EY21[h] * CZ[h] }  # transition probabilities from state 6 to state x
            
            xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), 1, 
                                prob=c(pt[1],pt[2],pt[3],pt[4],pt[5],pt[6],pt[7],pt[8],pt[9],pt[10],pt[11],pt[12],pt[13],pt[14],pt[15],pt[16]))
          }
          
          else if(xt[i,j,g-1]==7) { 
            
            for(h in 1:16){ pt[h] <- EX1[h] * CY1[h] * CY2[h] * EZ[h] }  # transition probabilities from state 7 to state x
            
            xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), 1, 
                                prob=c(pt[1],pt[2],pt[3],pt[4],pt[5],pt[6],pt[7],pt[8],pt[9],pt[10],pt[11],pt[12],pt[13],pt[14],pt[15],pt[16]))
          }
          
          else if(xt[i,j,g-1]==8) { 
            
            for(h in 1:16){ pt[h] <- CX[h] * EY11[h] * EY21[h] * CZ[h] }  # transition probabilities from state 8 to state x
            
            xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), 1, 
                                prob=c(pt[1],pt[2],pt[3],pt[4],pt[5],pt[6],pt[7],pt[8],pt[9],pt[10],pt[11],pt[12],pt[13],pt[14],pt[15],pt[16]))
          }
          
          else if(xt[i,j,g-1]==9) { 
            
            for(h in 1:16){ pt[h] <- CX[h] * EY12[h] * CY2[h] * EZ[h] }  # transition probabilities from state 9 to state x
            
            xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), 1, 
                                prob=c(pt[1],pt[2],pt[3],pt[4],pt[5],pt[6],pt[7],pt[8],pt[9],pt[10],pt[11],pt[12],pt[13],pt[14],pt[15],pt[16]))
          }
          
          else if(xt[i,j,g-1]==10) { 
            
            for(h in 1:16){ pt[h] <- CX[h] * CY1[h] * EY22[h] * EZ[h] }  # transition probabilities from state 10 to state x
            
            xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), 1, 
                                prob=c(pt[1],pt[2],pt[3],pt[4],pt[5],pt[6],pt[7],pt[8],pt[9],pt[10],pt[11],pt[12],pt[13],pt[14],pt[15],pt[16]))
          }
          
          else if(xt[i,j,g-1]==11) { 
            
            for(h in 1:16){ pt[h] <- EX4[h] * EY11[h] * EY21[h] * CZ[h] }  # transition probabilities from state 11 to state x
            
            xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), 1, 
                                prob=c(pt[1],pt[2],pt[3],pt[4],pt[5],pt[6],pt[7],pt[8],pt[9],pt[10],pt[11],pt[12],pt[13],pt[14],pt[15],pt[16]))
          }
          
          else if(xt[i,j,g-1]==12) { 
            
            for(h in 1:16){ pt[h] <- EX2[h] * EY12[h] * CY2[h] * EZ[h] }  # transition probabilities from state 12 to state x
            
            xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), 1, 
                                prob=c(pt[1],pt[2],pt[3],pt[4],pt[5],pt[6],pt[7],pt[8],pt[9],pt[10],pt[11],pt[12],pt[13],pt[14],pt[15],pt[16]))
          }
          
          else if(xt[i,j,g-1]==13) { 
            
            for(h in 1:16){ pt[h] <- EX3[h] * CY1[h] * EY22[h] * EZ[h] }  # transition probabilities from state 13 to state x
            
            xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), 1, 
                                prob=c(pt[1],pt[2],pt[3],pt[4],pt[5],pt[6],pt[7],pt[8],pt[9],pt[10],pt[11],pt[12],pt[13],pt[14],pt[15],pt[16]))
          }
          
          else if(xt[i,j,g-1]==14) { 
            
            for(h in 1:16){ pt[h] <- CX[h] * EY12[h] * EY22[h] * EZ[h] }  # transition probabilities from state 14 to state x
            
            xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), 1, 
                                prob=c(pt[1],pt[2],pt[3],pt[4],pt[5],pt[6],pt[7],pt[8],pt[9],pt[10],pt[11],pt[12],pt[13],pt[14],pt[15],pt[16]))
          }
          
          else if(xt[i,j,g-1]==15) { 
            
            for(h in 1:16){ pt[h] <- EX4[h] * EY12[h] * EY22[h] * EZ[h] }  # transition probabilities from state 15 to state x
            
            xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), 1, 
                                prob=c(pt[1],pt[2],pt[3],pt[4],pt[5],pt[6],pt[7],pt[8],pt[9],pt[10],pt[11],pt[12],pt[13],pt[14],pt[15],pt[16]))
          }
          
        }
      }

      pxdt[k,g] <- (sum(xt[ , , g]==1)+sum(xt[ , , g]==5)+sum(xt[ , , g]==6)+sum(xt[ , , g]==7)
                    +sum(xt[ , , g]==11)+sum(xt[ , , g]==12)+sum(xt[ , , g]==13)+sum(xt[ , , g]==15))/(n*n)
      Pxdt[k,g] <- (sum(xt[ , , g]==1)+sum(xt[ , , g]==5)+sum(xt[ , , g]==6)+sum(xt[ , , g]==7)
                    +sum(xt[ , , g]==11)+sum(xt[ , , g]==12)+sum(xt[ , , g]==13)+sum(xt[ , , g]==15))/(n*n-d[k])
      py1dt[k,g] <- (sum(xt[ , , g]==2)+sum(xt[ , , g]==5)+sum(xt[ , , g]==8)+sum(xt[ , , g]==9)
                     +sum(xt[ , , g]==11)+sum(xt[ , , g]==12)+sum(xt[ , , g]==14)+sum(xt[ , , g]==15))/(n*n)
      Py1dt[k,g] <- (sum(xt[ , , g]==2)+sum(xt[ , , g]==5)+sum(xt[ , , g]==8)+sum(xt[ , , g]==9)
                     +sum(xt[ , , g]==11)+sum(xt[ , , g]==12)+sum(xt[ , , g]==14)+sum(xt[ , , g]==15))/(n*n-d[k])
      py2dt[k,g] <- (sum(xt[ , , g]==3)+sum(xt[ , , g]==6)+sum(xt[ , , g]==8)+sum(xt[ , , g]==10)
                     +sum(xt[ , , g]==11)+sum(xt[ , , g]==13)+sum(xt[ , , g]==14)+sum(xt[ , , g]==15))/(n*n)
      Py2dt[k,g] <- (sum(xt[ , , g]==3)+sum(xt[ , , g]==6)+sum(xt[ , , g]==8)+sum(xt[ , , g]==10)
                     +sum(xt[ , , g]==11)+sum(xt[ , , g]==13)+sum(xt[ , , g]==14)+sum(xt[ , , g]==15))/(n*n-d[k])
      pzdt[k,g] <- (sum(xt[ , , g]==4)+sum(xt[ , , g]==7)+sum(xt[ , , g]==9)+sum(xt[ , , g]==10)
                    +sum(xt[ , , g]==12)+sum(xt[ , , g]==13)+sum(xt[ , , g]==14)+sum(xt[ , , g]==15))/(n*n)
      Pzdt[k,g] <- (sum(xt[ , , g]==4)+sum(xt[ , , g]==7)+sum(xt[ , , g]==9)+sum(xt[ , , g]==10)
                    +sum(xt[ , , g]==12)+sum(xt[ , , g]==13)+sum(xt[ , , g]==14)+sum(xt[ , , g]==15))/(n*n-d[k])
    }
    xdeq[ , ,k] <- xt[ , ,length(t)]
  }
  
  output <- list("D"=D, "pxdt"=pxdt, "Pxdt"=Pxdt, "py1dt"=py1dt, "Py1dt"=Py1dt, "py2dt"=py2dt, "Py2dt"=Py2dt, 
                 "pzdt"=pzdt, "Pzdt"=Pzdt, "xdeq"=xdeq)
  return(output)
}

# habitat destruction output
outD <- dynamicsD(cx, cy1, cy2, cz, ex, ey1, ey2, ez, mu1, mu2, mu3, mu4, psi1, psi2, psi3, psi4, 
                  px0, py10, py20, pz0, n, tmax, dD)

D <- outD$D           # fraction of habitat destoryed [0,1]
pxdt <- outD$pxdt     # fraction of x-occupied sites as function of time and D
Pxdt <- outD$Pxdt     # fraction of x-occupied sites out of non-destroyed sites as function of time and D
py1dt <- outD$py1dt   # fraction of y1-occupied sites as function of time and D
Py1dt <- outD$Py1dt   # fraction of y1-occupied sites out of non-destroyed sites as function of time and D
py2dt <- outD$py2dt   # fraction of y2-occupied sites as function of time and D
Py2dt <- outD$Py2dt   # fraction of y2-occupied sites out of non-destroyed sites as function of time and D
pzdt <- outD$pzdt     # fraction of z-occupied sites as function of time and D
Pzdt <- outD$Pzdt     # fraction of z-occupied sites out of non-destroyed sites as function of time and D
xdeq <- outD$xdeq     # equilibrium grid as function of D

plot(x=D, y=pxdt[ ,tmax+1], type="p", col="green", ylim=c(0,1), ylab="p", pch=20)
lines(x=D, y=py1dt[ ,tmax+1], type="p", col="blue", pch=20)
lines(x=D, y=py2dt[ ,tmax+1], type="p", col="black", pch=20)
lines(x=D, y=pzdt[ ,tmax+1], type="p", col="red", pch=20)
legend(0.7, 1, legend=c("x-resource", "y1-consumer1", "y2-consumer2", "z-predator"), 
       col=c("green", "blue", "black", "red"), pch=20, cex=1)
abline(h=0, v=0)


# HABITAT RESTORATION ----

# input - habitat restoration
cx  <- 0.7         # probablity of colonisation - resource
cy1 <- 0.7         # probablity of colonisation - consumer 1
cy2 <- 0.9         # probablity of colonisation - consumer 2
cz  <- 0.7         # probablity of colonisation - predator
ex  <- 0.1         # probability of extinction - resource
ey1 <- 0.1         # probability of extinction - consumer 1
ey2 <- 0.1         # probability of extinction - consumer 2
ez  <- 0.1         # probability of extinction - predator

mu1 <- 0.45       # probablilty of resource extinction due to consumer 1 presence
mu2 <- 0.45       # probablilty of resource extinction due to consumer 2 presence
mu3 <- 0.45       # probablilty of consumer 1 extinction due to predator presence
mu4 <- 0.45       # probablilty of consumer 2 extinction due to predator presence
psi1 <- 0.4       # probablilty of consumer 1 extinction due lack of resource
psi2 <- 0.4       # probablilty of consumer 2 extinction due lack of resource
psi3 <- 0.4       # probablilty of predator extinction due lack of consumer 1
psi4 <- 0.4       # probablilty of predator extinction due lack of consumer 2

tmax <- 50        # maximum number of timesteps
DR <- 0.10        # fraction of patches destroyed at the start of restoration
dR <- 0.01        # fraction of patches restored at a time

Rex  <- rep(NA,length(D))     # initialise restoration efficiency vector
Rey1  <- rep(NA,length(D))
Rey2  <- rep(NA,length(D))
Rez  <- rep(NA,length(D))

# initialise matrices for storing results
R_all <- matrix(NA, nrow=length(D), ncol=length(D))      # fraction of habitat loss during resotration
pxr_all <- matrix(NA, nrow=length(D), ncol=length(D))    # fraction of x-occupied sites (mean of replicates)
Pxr_all <- matrix(NA, nrow=length(D), ncol=length(D))    # fraction of x-occupied sites out of non-destroyed sites (mean of replicates)
py1r_all <- matrix(NA, nrow=length(D), ncol=length(D))   # fraction of y1-occupied sites (mean of replicates
Py1r_all <- matrix(NA, nrow=length(D), ncol=length(D))   # fraction of y1-occupied sites out of non-destroyed sites (mean of replicates)
py2r_all <- matrix(NA, nrow=length(D), ncol=length(D))   # fraction of y2-occupied sites (mean of replicates
Py2r_all <- matrix(NA, nrow=length(D), ncol=length(D))   # fraction of y2-occupied sites out of non-destroyed sites (mean of replicates)
pzr_all <- matrix(NA, nrow=length(D), ncol=length(D))    # fraction of z-occupied sites (mean of replicates
Pzr_all <- matrix(NA, nrow=length(D), ncol=length(D))    # fraction of z-occupied sites out of non-destroyed sites (mean of replicates)
sdx_all <- matrix(NA, nrow=length(D), ncol=length(D))    # st dev of fraction of x-occupied sites of replicates
sdy1_all <- matrix(NA, nrow=length(D), ncol=length(D))   # st dev of fraction of y1-occupied sites of replicates
sdy2_all <- matrix(NA, nrow=length(D), ncol=length(D))   # st dev of fraction of y2-occupied sites of replicates
sdz_all <- matrix(NA, nrow=length(D), ncol=length(D))    # st dev of fraction of z-occupied sites of replicates


# dynamics with habitat restoration

dynamicsR <- function(cx, cy1, cy2, cz, ex, ey1, ey2, ez, mu1, mu2, mu3, mu4, psi1, psi2, psi3, psi4, 
                      tmax, DR, dR) {
  
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
  py1rt <- matrix(0, nrow=length(r), ncol=length(t))    
  py1rt[1, ] <- py1dt[Di, ]
  Py1rt <- matrix(0, nrow=length(r), ncol=length(t))
  Py1rt[1, ] <- Py1dt[Di, ]
  py2rt <- matrix(0, nrow=length(r), ncol=length(t))    
  py2rt[1, ] <- py2dt[Di, ]
  Py2rt <- matrix(0, nrow=length(r), ncol=length(t))
  Py2rt[1, ] <- Py2dt[Di, ]
  pzrt <- matrix(0, nrow=length(r), ncol=length(t))    
  pzrt[1, ] <- pzdt[Di, ]
  Pzrt <- matrix(0, nrow=length(r), ncol=length(t))
  Pzrt[1, ] <- Pzdt[Di, ]
  
  # calculate terms for transition probabilities
  EX1 <- c(ex, (1-ex), ex, ex, ex, (1-ex), (1-ex), (1-ex), 0, ex, ex, 0, (1-ex), (1-ex), 0, 0)
  EX2 <- c((ex+mu1), (1-ex-mu1), (ex+mu1), (ex+mu1), (ex+mu1), (1-ex-mu1), (1-ex-mu1), (1-ex-mu1), 
           0, (ex+mu1), (ex+mu1), 0, (1-ex-mu1), (1-ex-mu1), 0, 0)
  EX3 <- c((ex+mu2), (1-ex-mu2), (ex+mu2), (ex+mu2), (ex+mu2), (1-ex-mu2), (1-ex-mu2), (1-ex-mu2), 
           (ex+mu2), (ex+mu2), (ex+mu2), (1-ex-mu2), (1-ex-mu2), (1-ex-mu2), (ex+mu2), (1-ex-mu2))
  EX4 <- c((ex+mu1+mu2), (1-ex-mu1-mu2), (ex+mu1+mu2), (ex+mu1+mu2), (ex+mu1+mu2), (1-ex-mu1-mu2), (1-ex-mu1-mu2), (1-ex-mu1-mu2), 
           (ex+mu1+mu2), (ex+mu1+mu2), (ex+mu1+mu2), (1-ex-mu1-mu2), (1-ex-mu1-mu2), (1-ex-mu1-mu2), (ex+mu1+mu2), (1-ex-mu1-mu2))
  EY11 <- c((ey1+psi1), ey1, (1-ey1-psi1), (ey1+psi1), (ey1+psi1), (1-ey1), ey1, ey1, 
            (1-ey1-psi1), (1-ey1-psi1), (ey1+psi1), (1-ey1), (1-ey1), ey1, (1-ey1-psi1), (1-ey1))
  EY12 <- c((ey1+psi1+mu3), (ey1+mu3), (1-ey1-psi1-mu3), (ey1+psi1+mu3), (ey1+psi1+mu3), (1-ey1-mu3), (ey1+mu3), (ey1+mu3), 
            (1-ey1-psi1-mu3), (1-ey1-psi1-mu3), (ey1+psi1+mu3), (1-ey1-mu3), (1-ey1-mu3), (ey1+mu3), (1-ey1-psi1-mu3), (1-ey1-mu3))
  EY21 <- c((ey2+psi2), ey2, (ey2+psi2), (1-ey2-psi2), (ey2+psi2), ey2, (1-ey2), ey2, 
            (1-ey2-psi2), (ey2+psi2), (1-ey2-psi2), (1-ey2), ey2, (1-ey2), (1-ey2-psi2), (1-ey2))
  EY22 <- c((ey2+psi2+mu4), (ey2+mu4), (ey2+psi2+mu4), (1-ey2-psi2-mu4), (ey2+psi2+mu4), (ey2+mu4), (1-ey2-mu4), (ey2+mu4), 
            (1-ey2-psi2-mu4), (ey2+psi2+mu4), (1-ey2-psi2-mu4), (1-ey2-mu4), (ey2+mu4), (1-ey2-mu4), (1-ey2-psi2-mu4), (1-ey2-mu4))
  EZ  <- c((ez+psi3/2+psi4/2), (ez+psi3/2+psi4/2), (ez+psi4/2), (ez+psi3/2), 
           (1-ez-psi3/2-psi4/2), (ez+psi4/2), (ez+psi3/2), (1-ez-psi3/2-psi4/2), 
           ez, (1-ez-psi4/2), (1-ez-psi3/2), ez, 
           (1-ez-psi4/2), (1-ez-psi3/2), (1-ez), (1-ez))
  
  for(k in 2:length(r)){
    
    # initialise grid
    xreq[ , ,k] <- xreq[ , ,k-1]
    
    n_x <- 0  # number of destroyed patches with x neighbours only
    n_y1 <- 0  # number of destroyed patches with y1 neighbours only
    n_y2 <- 0  # number of destroyed patches with y2 neighbours only
    n_z <- 0  # number of destroyed patches with z neighbours only
    n_e <- 0  # number of destroyed patches with empty neighbours only
    d_x <- matrix(NA,nrow=sum(xreq[ , ,k]==-1),ncol=2) # matrix with coordinates of destroyed patches with x neighbours
    d_y1 <- matrix(NA,nrow=sum(xreq[ , ,k]==-1),ncol=2) # matrix with coordinates of destroyed patches with y1 neighbours
    d_y2 <- matrix(NA,nrow=sum(xreq[ , ,k]==-1),ncol=2) # matrix with coordinates of destroyed patches with y2 neighbours
    d_z <- matrix(NA,nrow=sum(xreq[ , ,k]==-1),ncol=2) # matrix with coordinates of destroyed patches with z neighbours
    d_e <- matrix(NA,nrow=sum(xreq[ , ,k]==-1),ncol=2) # matrix with coordinates of destroyed patches with empty neighbours only
    
    for(i in 1:n){
      for(j in 1:n){
        
        if(xreq[i,j,k]!=-1) { next }
        
        else{
          # check if any adjecent site is occupied
          # von neumann neighborhood
          n1 <- 0
          n2 <- 0
          n3 <- 0
          n4 <- 0
          neigh <- matrix(c(-1, 0, 1, 0, 0, 1, 0, -1), nrow=4, ncol=2)
          for(h in 1:4){
            i_n <- i+neigh[h,1]
            j_n <- j+neigh[h,2]
            
            if(i_n < 1 || i_n > n || j_n < 1 || j_n > n) { next } # neighbour outside of grid
            if(xreq[i_n,j_n,k]==4 || xreq[i_n,j_n,k]==7 || xreq[i_n,j_n,k]==9 || xreq[i_n,j_n,k]==10 || 
               xreq[i_n,j_n,k]==12 || xreq[i_n,j_n,k]==13 || xreq[i_n,j_n,k]==14 || xreq[i_n,j_n,k]==15) {
              n_z <- n_z+1
              d_z[n_z,1] <- i
              d_z[n_z,2] <- j
              n1 <- 1
              break }               # if neighbour occupied by z, move to next patch
            else if(xreq[i_n,j_n,k]==3 || xreq[i_n,j_n,k]==6 || xreq[i_n,j_n,k]==8 || xreq[i_n,j_n,k]==11) {
              n2 <- 1 }
            else if(xreq[i_n,j_n,k]==2 || xreq[i_n,j_n,k]==5) {
              n3 <- 1 }
            else if(xreq[i_n,j_n,k]==1) {
              n4 <- 1 }
          }
          if(n1==0 && n2==1) {
            n_y2 <- n_y2+1
            d_y2[n_y2,1] <- i
            d_y2[n_y2,2] <- j
          }
          else if(n1==0 && n2==0 && n3==1) {
            n_y1 <- n_y1+1
            d_y1[n_y1,1] <- i
            d_y1[n_y1,2] <- j
          }
          else if(n1==0 && n2==0 && n3==0 && n4==1) {
            n_x <- n_x+1
            d_x[n_x,1] <- i
            d_x[n_x,2] <- j
          }
          else if(n1==0 && n2==0 && n3==0 && n4==0){
            n_e <- n_e+1
            d_e[n_e,1] <- i
            d_e[n_e,2] <- j
          }
        }
      }
    }
    
    d_x <- na.omit(d_x)
    d_y1 <- na.omit(d_y1)
    d_y2 <- na.omit(d_y2)
    d_z <- na.omit(d_z)
    d_e <- na.omit(d_e)
    
    # restore patches adjecent to z-occupied patches only
    if(n_z>=dn) {
      rp <- sample(1:n_z,dn)
      for(h in 1:dn) {
        i_d <- d_z[rp[h],1]
        j_d <- d_z[rp[h],2]
        xreq[i_d,j_d,k] <- 0
      }
    }
    
    # restore patches adjecent to z-occupied patches first, then y2-occupied, then y1-occupied, then x-occupied, then empty
    if(n_z<dn) {
      if(n_z>0) {
        for(h in 1:n_z) {
          i_d <- d_z[h,1]
          j_d <- d_z[h,2]
          xreq[i_d,j_d,k] <- 0
        }
      }
      if(n_y2>=(dn-n_z)) {
        rp <- sample(1:n_y2,(dn-n_z))
        for(h in 1:(dn-n_z)) {
          i_d <- d_y2[rp[h],1]
          j_d <- d_y2[rp[h],2]
          xreq[i_d,j_d,k] <- 0
        }
      }
      else if(n_y2<(dn-n_z)) {
        if(n_y2>0) {
          for(h in 1:n_y2) {
            i_d <- d_y2[h,1]
            j_d <- d_y2[h,2]
            xreq[i_d,j_d,k] <- 0
          }
        }
        if(n_y1>=(dn-n_z-n_y2)) {
          rp <- sample(1:n_y1,(dn-n_z-n_y2))
          for(h in 1:(dn-n_z-n_y2)) {
            i_d <- d_y1[rp[h],1]
            j_d <- d_y1[rp[h],2]
            xreq[i_d,j_d,k] <- 0
          }
        }
        else if(n_y1<(dn-n_z-n_y2)) {
          if(n_y1>0) {
            for(h in 1:n_y1) {
              i_d <- d_y1[h,1]
              j_d <- d_y1[h,2]
              xreq[i_d,j_d,k] <- 0
            }
          }
          if(n_x>=(dn-n_z-n_y2-n_y1)) {
            rp <- sample(1:n_x,(dn-n_z-n_y2-n_y1))
            for(h in 1:(dn-n_z-n_y2-n_y1)) {
              i_d <- d_x[rp[h],1]
              j_d <- d_x[rp[h],2]
              xreq[i_d,j_d,k] <- 0
            }
          }
          else if(n_x<(dn-n_z-n_y2-n_y1)) {
            if(n_x>0) {
              for(h in 1:n_x) {
                i_d <- d_x[h,1]
                j_d <- d_x[h,2]
                xreq[i_d,j_d,k] <- 0
              }
            }
            rp <- sample(1:n_e,(dn-n_z-n_y2-n_y1-n_x))
            for(h in 1:(dn-n_z-n_y2-n_y1-n_x)) {
              i_d <- d_e[rp[h],1]
              j_d <- d_e[rp[h],2]
              xreq[i_d,j_d,k] <- 0
            }
          }
        }
      }
    }
    
    # grid at t=0 for d[k]
    xt[ , ,1] <- xreq[ , ,k]
    pxrt[k,1] <- (sum(xt[ , , 1]==1)+sum(xt[ , , 1]==5)+sum(xt[ , , 1]==6)+sum(xt[ , , 1]==7)
                  +sum(xt[ , , 1]==11)+sum(xt[ , , 1]==12)+sum(xt[ , , 1]==13)+sum(xt[ , , 1]==15))/(n*n)
    Pxrt[k,1] <- (sum(xt[ , , 1]==1)+sum(xt[ , , 1]==5)+sum(xt[ , , 1]==6)+sum(xt[ , , 1]==7)
                  +sum(xt[ , , 1]==11)+sum(xt[ , , 1]==12)+sum(xt[ , , 1]==13)+sum(xt[ , , 1]==15))/(n*n-r[k])
    py1rt[k,1] <- (sum(xt[ , , 1]==2)+sum(xt[ , , 1]==5)+sum(xt[ , , 1]==8)+sum(xt[ , , 1]==9)
                   +sum(xt[ , , 1]==11)+sum(xt[ , , 1]==12)+sum(xt[ , , 1]==14)+sum(xt[ , , 1]==15))/(n*n)
    Py1rt[k,1] <- (sum(xt[ , , 1]==2)+sum(xt[ , , 1]==5)+sum(xt[ , , 1]==8)+sum(xt[ , , 1]==9)
                   +sum(xt[ , , 1]==11)+sum(xt[ , , 1]==12)+sum(xt[ , , 1]==14)+sum(xt[ , , 1]==15))/(n*n-r[k])
    py2rt[k,1] <- (sum(xt[ , , 1]==3)+sum(xt[ , , 1]==6)+sum(xt[ , , 1]==8)+sum(xt[ , , 1]==10)
                   +sum(xt[ , , 1]==11)+sum(xt[ , , 1]==13)+sum(xt[ , , 1]==14)+sum(xt[ , , 1]==15))/(n*n)
    Py2rt[k,1] <- (sum(xt[ , , 1]==3)+sum(xt[ , , 1]==6)+sum(xt[ , , 1]==8)+sum(xt[ , , 1]==10)
                   +sum(xt[ , , 1]==11)+sum(xt[ , , 1]==13)+sum(xt[ , , 1]==14)+sum(xt[ , , 1]==15))/(n*n-r[k])
    pzrt[k,1] <- (sum(xt[ , , 1]==4)+sum(xt[ , , 1]==7)+sum(xt[ , , 1]==9)+sum(xt[ , , 1]==10)
                  +sum(xt[ , , 1]==12)+sum(xt[ , , 1]==13)+sum(xt[ , , 1]==14)+sum(xt[ , , 1]==15))/(n*n)
    Pzrt[k,1] <- (sum(xt[ , , 1]==4)+sum(xt[ , , 1]==7)+sum(xt[ , , 1]==9)+sum(xt[ , , 1]==10)
                  +sum(xt[ , , 1]==12)+sum(xt[ , , 1]==13)+sum(xt[ , , 1]==14)+sum(xt[ , , 1]==15))/(n*n-r[k])
    
    # iterate until steady state
    for(g in 2:length(t)){
      for(i in 1:n){
        for(j in 1:n){
          
          # von neumann neighborhood
          nx  <- 0    # number of neighbours with x
          ny1 <- 0    # number of neighbours with y1
          ny2 <- 0    # number of neighbours with y2
          nz  <- 0    # number of neighbours with z
          neigh <- matrix(c(-1, 0, 1, 0, 0, 1, 0, -1), nrow=4, ncol=2)
          for(h in 1:4){
            i_n <- i+neigh[h,1]
            j_n <- j+neigh[h,2]
            
            if(i_n < 1 || i_n > n || j_n < 1 || j_n > n) { next } # neighbour outside of grid
            
            if (xt[i_n,j_n,g-1]==1 || xt[i_n,j_n,g-1]==5 || xt[i_n,j_n,g-1]==6 || xt[i_n,j_n,g-1]==7 || 
                xt[i_n,j_n,g-1]==11 || xt[i_n,j_n,g-1]==12 || xt[i_n,j_n,g-1]==13 || xt[i_n,j_n,g-1]==15) { nx <- nx + 1 }
            if (xt[i_n,j_n,g-1]==2 || xt[i_n,j_n,g-1]==5 || xt[i_n,j_n,g-1]==8 || xt[i_n,j_n,g-1]==9 || 
                xt[i_n,j_n,g-1]==11 || xt[i_n,j_n,g-1]==12 || xt[i_n,j_n,g-1]==14 || xt[i_n,j_n,g-1]==15) { ny1 <- ny1 + 1 }
            if (xt[i_n,j_n,g-1]==3 || xt[i_n,j_n,g-1]==6 || xt[i_n,j_n,g-1]==8 || xt[i_n,j_n,g-1]==10 || 
                xt[i_n,j_n,g-1]==11 || xt[i_n,j_n,g-1]==13 || xt[i_n,j_n,g-1]==14 || xt[i_n,j_n,g-1]==15) { ny2 <- ny2 + 1 }
            if (xt[i_n,j_n,g-1]==4 || xt[i_n,j_n,g-1]==7 || xt[i_n,j_n,g-1]==9 || xt[i_n,j_n,g-1]==10 || 
                xt[i_n,j_n,g-1]==12 || xt[i_n,j_n,g-1]==13 || xt[i_n,j_n,g-1]==14 || xt[i_n,j_n,g-1]==15) { nz <- nz + 1 }
          }
          
          # calculate terms for transition probabilities
          CX <- c((1-cx)^nx, (1-(1-cx)^nx), (1-cx)^nx, (1-cx)^nx, (1-cx)^nx, (1-(1-cx)^nx), (1-(1-cx)^nx), (1-(1-cx)^nx), 
                  (1-cx)^nx, (1-cx)^nx, (1-cx)^nx, (1-(1-cx)^nx), (1-(1-cx)^nx), (1-(1-cx)^nx), (1-cx)^nx, (1-(1-cx)^nx))
          CY1 <- c((1-cy1)^ny1, (1-cy1)^ny1, (1-(1-cy1)^ny1), (1-cy1)^ny1, 
                   (1-cy1)^ny1, (1-(1-cy1)^ny1), (1-cy1)^ny1, (1-cy1)^ny1, 
                   (1-(1-cy1)^ny1), (1-(1-cy1)^ny1), (1-cy1)^ny1, (1-(1-cy1)^ny1), 
                   (1-(1-cy1)^ny1), (1-cy1)^ny1, (1-(1-cy1)^ny1), (1-(1-cy1)^ny1))
          CY2 <- c((1-cy2)^ny2, (1-cy2)^ny2, 1, (1-(1-cy2)^ny2), (1-cy2)^ny2, 1, (1-(1-cy2)^ny2), (1-cy2)^ny2, 
                   0, 1, (1-(1-cy2)^ny2), 0, 1, (1-(1-cy2)^ny2), 0, 0)
          CZ <- c((1-cz)^nz, (1-cz)^nz, (1-cz)^nz, (1-cz)^nz, (1-(1-cz)^nz), (1-cz)^nz, (1-cz)^nz, (1-(1-cz)^nz), 
                  (1-cz)^nz, (1-(1-cz)^nz), (1-(1-cz)^nz), (1-cz)^nz, (1-(1-cz)^nz), (1-(1-cz)^nz), (1-(1-cz)^nz), (1-(1-cz)^nz))
          
          pt <- numeric(length=16)   # vector of transition probabilities
          
          # transitions
          if(xt[i,j,g-1]==-1) { xt[i,j,g] <- -1 }    # destroyed sites remain destroyed
          
          else if(xt[i,j,g-1]==0) { 
            
            for(h in 1:16){ pt[h] <- CX[h] * CY1[h] * CY2[h] * CZ[h] }  # transition probabilities from state 0 to state x
            
            xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), 1, 
                                prob=c(pt[1],pt[2],pt[3],pt[4],pt[5],pt[6],pt[7],pt[8],pt[9],pt[10],pt[11],pt[12],pt[13],pt[14],pt[15],pt[16]))
          }
          
          else if(xt[i,j,g-1]==1) { 
            
            for(h in 1:16){ pt[h] <- EX1[h] * CY1[h] * CY2[h] * CZ[h] }  # transition probabilities from state 1 to state x
            
            xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), 1, 
                                prob=c(pt[1],pt[2],pt[3],pt[4],pt[5],pt[6],pt[7],pt[8],pt[9],pt[10],pt[11],pt[12],pt[13],pt[14],pt[15],pt[16]))
          }
          
          else if(xt[i,j,g-1]==2) { 
            
            for(h in 1:16){ pt[h] <- CX[h] * EY11[h] * CY2[h] * CZ[h] }  # transition probabilities from state 2 to state x
            
            xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), 1, 
                                prob=c(pt[1],pt[2],pt[3],pt[4],pt[5],pt[6],pt[7],pt[8],pt[9],pt[10],pt[11],pt[12],pt[13],pt[14],pt[15],pt[16]))
          }
          
          else if(xt[i,j,g-1]==3) { 
            
            for(h in 1:16){ pt[h] <- CX[h] * CY1[h] * EY21[h] * CZ[h] }  # transition probabilities from state 3 to state x
            
            xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), 1, 
                                prob=c(pt[1],pt[2],pt[3],pt[4],pt[5],pt[6],pt[7],pt[8],pt[9],pt[10],pt[11],pt[12],pt[13],pt[14],pt[15],pt[16]))
          }
          
          else if(xt[i,j,g-1]==4) { 
            
            for(h in 1:16){ pt[h] <- CX[h] * CY1[h] * CY2[h] * EZ[h] }  # transition probabilities from state 4 to state x
            
            xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), 1, 
                                prob=c(pt[1],pt[2],pt[3],pt[4],pt[5],pt[6],pt[7],pt[8],pt[9],pt[10],pt[11],pt[12],pt[13],pt[14],pt[15],pt[16]))
          }
          
          else if(xt[i,j,g-1]==5) { 
            
            for(h in 1:16){ pt[h] <- EX2[h] * EY11[h] * CY2[h] * CZ[h] }  # transition probabilities from state 5 to state x
            
            xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), 1, 
                                prob=c(pt[1],pt[2],pt[3],pt[4],pt[5],pt[6],pt[7],pt[8],pt[9],pt[10],pt[11],pt[12],pt[13],pt[14],pt[15],pt[16]))
          }
          
          else if(xt[i,j,g-1]==6) { 
            
            for(h in 1:16){ pt[h] <- EX3[h] * CY1[h] * EY21[h] * CZ[h] }  # transition probabilities from state 6 to state x
            
            xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), 1, 
                                prob=c(pt[1],pt[2],pt[3],pt[4],pt[5],pt[6],pt[7],pt[8],pt[9],pt[10],pt[11],pt[12],pt[13],pt[14],pt[15],pt[16]))
          }
          
          else if(xt[i,j,g-1]==7) { 
            
            for(h in 1:16){ pt[h] <- EX1[h] * CY1[h] * CY2[h] * EZ[h] }  # transition probabilities from state 7 to state x
            
            xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), 1, 
                                prob=c(pt[1],pt[2],pt[3],pt[4],pt[5],pt[6],pt[7],pt[8],pt[9],pt[10],pt[11],pt[12],pt[13],pt[14],pt[15],pt[16]))
          }
          
          else if(xt[i,j,g-1]==8) { 
            
            for(h in 1:16){ pt[h] <- CX[h] * EY11[h] * EY21[h] * CZ[h] }  # transition probabilities from state 8 to state x
            
            xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), 1, 
                                prob=c(pt[1],pt[2],pt[3],pt[4],pt[5],pt[6],pt[7],pt[8],pt[9],pt[10],pt[11],pt[12],pt[13],pt[14],pt[15],pt[16]))
          }
          
          else if(xt[i,j,g-1]==9) { 
            
            for(h in 1:16){ pt[h] <- CX[h] * EY12[h] * CY2[h] * EZ[h] }  # transition probabilities from state 9 to state x
            
            xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), 1, 
                                prob=c(pt[1],pt[2],pt[3],pt[4],pt[5],pt[6],pt[7],pt[8],pt[9],pt[10],pt[11],pt[12],pt[13],pt[14],pt[15],pt[16]))
          }
          
          else if(xt[i,j,g-1]==10) { 
            
            for(h in 1:16){ pt[h] <- CX[h] * CY1[h] * EY22[h] * EZ[h] }  # transition probabilities from state 10 to state x
            
            xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), 1, 
                                prob=c(pt[1],pt[2],pt[3],pt[4],pt[5],pt[6],pt[7],pt[8],pt[9],pt[10],pt[11],pt[12],pt[13],pt[14],pt[15],pt[16]))
          }
          
          else if(xt[i,j,g-1]==11) { 
            
            for(h in 1:16){ pt[h] <- EX4[h] * EY11[h] * EY21[h] * CZ[h] }  # transition probabilities from state 11 to state x
            
            xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), 1, 
                                prob=c(pt[1],pt[2],pt[3],pt[4],pt[5],pt[6],pt[7],pt[8],pt[9],pt[10],pt[11],pt[12],pt[13],pt[14],pt[15],pt[16]))
          }
          
          else if(xt[i,j,g-1]==12) { 
            
            for(h in 1:16){ pt[h] <- EX2[h] * EY12[h] * CY2[h] * EZ[h] }  # transition probabilities from state 12 to state x
            
            xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), 1, 
                                prob=c(pt[1],pt[2],pt[3],pt[4],pt[5],pt[6],pt[7],pt[8],pt[9],pt[10],pt[11],pt[12],pt[13],pt[14],pt[15],pt[16]))
          }
          
          else if(xt[i,j,g-1]==13) { 
            
            for(h in 1:16){ pt[h] <- EX3[h] * CY1[h] * EY22[h] * EZ[h] }  # transition probabilities from state 13 to state x
            
            xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), 1, 
                                prob=c(pt[1],pt[2],pt[3],pt[4],pt[5],pt[6],pt[7],pt[8],pt[9],pt[10],pt[11],pt[12],pt[13],pt[14],pt[15],pt[16]))
          }
          
          else if(xt[i,j,g-1]==14) { 
            
            for(h in 1:16){ pt[h] <- CX[h] * EY12[h] * EY22[h] * EZ[h] }  # transition probabilities from state 14 to state x
            
            xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), 1, 
                                prob=c(pt[1],pt[2],pt[3],pt[4],pt[5],pt[6],pt[7],pt[8],pt[9],pt[10],pt[11],pt[12],pt[13],pt[14],pt[15],pt[16]))
          }
          
          else if(xt[i,j,g-1]==15) { 
            
            for(h in 1:16){ pt[h] <- EX4[h] * EY12[h] * EY22[h] * EZ[h] }  # transition probabilities from state 15 to state x
            
            xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), 1, 
                                prob=c(pt[1],pt[2],pt[3],pt[4],pt[5],pt[6],pt[7],pt[8],pt[9],pt[10],pt[11],pt[12],pt[13],pt[14],pt[15],pt[16]))
          }
          
        }
      }
      pxrt[k,g] <- (sum(xt[ , , g]==1)+sum(xt[ , , g]==5)+sum(xt[ , , g]==6)+sum(xt[ , , g]==7)
                    +sum(xt[ , , g]==11)+sum(xt[ , , g]==12)+sum(xt[ , , g]==13)+sum(xt[ , , g]==15))/(n*n)
      Pxrt[k,g] <- (sum(xt[ , , g]==1)+sum(xt[ , , g]==5)+sum(xt[ , , g]==6)+sum(xt[ , , g]==7)
                    +sum(xt[ , , g]==11)+sum(xt[ , , g]==12)+sum(xt[ , , g]==13)+sum(xt[ , , g]==15))/(n*n-r[k])
      py1rt[k,g] <- (sum(xt[ , , g]==2)+sum(xt[ , , g]==5)+sum(xt[ , , g]==8)+sum(xt[ , , g]==9)
                     +sum(xt[ , , g]==11)+sum(xt[ , , g]==12)+sum(xt[ , , g]==14)+sum(xt[ , , g]==15))/(n*n)
      Py1rt[k,g] <- (sum(xt[ , , g]==2)+sum(xt[ , , g]==5)+sum(xt[ , , g]==8)+sum(xt[ , , g]==9)
                     +sum(xt[ , , g]==11)+sum(xt[ , , g]==12)+sum(xt[ , , g]==14)+sum(xt[ , , g]==15))/(n*n-r[k])
      py2rt[k,g] <- (sum(xt[ , , g]==3)+sum(xt[ , , g]==6)+sum(xt[ , , g]==8)+sum(xt[ , , g]==10)
                     +sum(xt[ , , g]==11)+sum(xt[ , , g]==13)+sum(xt[ , , g]==14)+sum(xt[ , , g]==15))/(n*n)
      Py2rt[k,g] <- (sum(xt[ , , g]==3)+sum(xt[ , , g]==6)+sum(xt[ , , g]==8)+sum(xt[ , , g]==10)
                     +sum(xt[ , , g]==11)+sum(xt[ , , g]==13)+sum(xt[ , , g]==14)+sum(xt[ , , g]==15))/(n*n-r[k])
      pzrt[k,g] <- (sum(xt[ , , g]==4)+sum(xt[ , , g]==7)+sum(xt[ , , g]==9)+sum(xt[ , , g]==10)
                    +sum(xt[ , , g]==12)+sum(xt[ , , g]==13)+sum(xt[ , , g]==14)+sum(xt[ , , g]==15))/(n*n)
      Pzrt[k,g] <- (sum(xt[ , , g]==4)+sum(xt[ , , g]==7)+sum(xt[ , , g]==9)+sum(xt[ , , g]==10)
                    +sum(xt[ , , g]==12)+sum(xt[ , , g]==13)+sum(xt[ , , g]==14)+sum(xt[ , , g]==15))/(n*n-r[k])
    }
    xreq[ , ,k] <- xt[ , ,length(t)]
  }
  output <- list("R"=R, "pxrt"=pxrt, "Pxrt"=Pxrt, "py1rt"=py1rt, "Py1rt"=Py1rt, "py2rt"=py2rt, "Py2rt"=Py2rt, 
                 "pzrt"=pzrt, "Pzrt"=Pzrt, "xreq"=xreq)
  return(output)
}

# habitat restoration output
outR <- dynamicsR(cx, cy1, cy2, cz, ex, ey1, ey2, ez, mu1, mu2, mu3, mu4, psi1, psi2, psi3, psi4, 
                  tmax, DR, dR)

R <- outR$R                # fraction of habitat loss [DR,0]
pxrt <- outR$pxrt          # fraction of x-occupied sites as function of time and D
Pxrt <- outR$Pxrt          # fraction of x-occupied sites out of non-destroyed sites as function of time and D
py1rt <- outR$py1rt        # fraction of y1-occupied sites as function of time and D
Py1rt <- outR$Py1rt        # fraction of y1-occupied sites out of non-destroyed sites as function of time and D
py2rt <- outR$py2rt        # fraction of y2-occupied sites as function of time and D
Py2rt <- outR$Py2rt        # fraction of y2-occupied sites out of non-destroyed sites as function of time and D
pzrt <- outR$pzrt          # fraction of z-occupied sites as function of time and D
Pzrt <- outR$Pzrt          # fraction of z-occupied sites out of non-destroyed sites as function of time and D
xreq <- outR$xreq          # equilibrium grid as function of D

pxr <- pxrt[, tmax+1]      # equilibrium values only
Pxr <- Pxrt[, tmax+1]
pxd <- pxdt[, tmax+1]
Pxd <- Pxdt[, tmax+1]
py1r <- py1rt[, tmax+1]
Py1r <- Py1rt[, tmax+1]
py1d <- py1dt[, tmax+1]
Py1d <- Py1dt[, tmax+1]
py2r <- py2rt[, tmax+1]
Py2r <- Py2rt[, tmax+1]
py2d <- py2dt[, tmax+1]
Py2d <- Py2dt[, tmax+1]
pzr <- pzrt[, tmax+1]
Pzr <- Pzrt[, tmax+1]
pzd <- pzdt[, tmax+1]
Pzd <- Pzdt[, tmax+1]

plot(x=D, y=pxd, type="b", col="lightgreen", xlim=c(0,1), ylim=c(0,1), xlab="D", ylab="p", pch = 20)
lines(x=R, y=pxr, type="b", col="green", pch = 20)
lines(x=D, y=py1d, type="b", col="lightblue", pch = 20)
lines(x=R, y=py1r, type="b", col="blue", pch = 20)
lines(x=D, y=py2d, type="b", col="grey", pch = 20)
lines(x=R, y=py2r, type="b", col="black", pch = 20)
lines(x=D, y=pzd, type="b", col="lightpink", pch = 20)
lines(x=R, y=pzr, type="b", col="red", pch = 20)
legend(0.7, 1, legend=c("x-resource", "y1-consumer1", "y2-consumer2", "z-predator"), 
       col=c("green", "blue", "black", "red"), pch = 20, cex=1)
abline(h=0, v=0)

for(k in 1:length(R)){
  image(x=1:n, y=1:n, z=t(apply(xreq[ , ,k], 2, rev)), main=paste("D =", R[k]), xlab="", ylab="", zlim=c(-1,15),
        col=brewer.pal(n=17, name="Set1"))
  legend(grconvertX(1, "device"), grconvertY(1, "device"), 
         c("destroyed","empty","x","y1","y2","z","x&y1","x&y2","x&z","y1&y2","y1&z","y2&z",
           "x&y1&y2","x&y1&z","x&y2&z","y1&y2&z","x&y1&y2&z"), 
         fill=brewer.pal(n=17, name="Set1"), xpd = NA)
}


# replicates
nr <- 5        # no of replicates
DR <- 0.10     # fraction of patches destroyed at the start of restoration

R <- seq(DR, 0, -dR)
pxr <- matrix(0, nrow=nr, ncol=length(R))
Pxr <- matrix(0, nrow=nr, ncol=length(R))
py1r <- matrix(0, nrow=nr, ncol=length(R))
Py1r <- matrix(0, nrow=nr, ncol=length(R))
py2r <- matrix(0, nrow=nr, ncol=length(R))
Py2r <- matrix(0, nrow=nr, ncol=length(R))
pzr <- matrix(0, nrow=nr, ncol=length(R))
Pzr <- matrix(0, nrow=nr, ncol=length(R))
for(i in 1:nr){
  outR <- dynamicsR(cx, cy1, cy2, cz, ex, ey1, ey2, ez, mu1, mu2, mu3, mu4, psi1, psi2, psi3, psi4, 
                    tmax, DR, dR)
  pxrt <- outR$pxrt
  Pxrt <- outR$Pxrt
  pxr[i, ] <- pxrt[ ,tmax+1]
  Pxr[i, ] <- Pxrt[ ,tmax+1]
  py1rt <- outR$py1rt
  Py1rt <- outR$Py1rt
  py1r[i, ] <- py1rt[ ,tmax+1]
  Py1r[i, ] <- Py1rt[ ,tmax+1]
  py2rt <- outR$py2rt
  Py2rt <- outR$Py2rt
  py2r[i, ] <- py2rt[ ,tmax+1]
  Py2r[i, ] <- Py2rt[ ,tmax+1]
  pzrt <- outR$pzrt
  Pzrt <- outR$Pzrt
  pzr[i, ] <- pzrt[ ,tmax+1]
  Pzr[i, ] <- Pzrt[ ,tmax+1]
}

# mean and st dev of replicates
pxr_mean <- numeric(length(R))
Pxr_mean <- numeric(length(R))
pxr_sd <- numeric(length(R))
Pxr_sd <- numeric(length(R))
py1r_mean <- numeric(length(R))
Py1r_mean <- numeric(length(R))
py1r_sd <- numeric(length(R))
Py1r_sd <- numeric(length(R))
py2r_mean <- numeric(length(R))
Py2r_mean <- numeric(length(R))
py2r_sd <- numeric(length(R))
Py2r_sd <- numeric(length(R))
pzr_mean <- numeric(length(R))
Pzr_mean <- numeric(length(R))
pzr_sd <- numeric(length(R))
Pzr_sd <- numeric(length(R))
for(i in 1:length(R)) {
  pxr_mean[i] <- mean(pxr[ , i])
  Pxr_mean[i] <- mean(Pxr[ , i])
  pxr_sd[i] <- sd(pxr[ , i])
  Pxr_sd[i] <- sd(Pxr[ , i])
  py1r_mean[i] <- mean(py1r[ , i])
  Py1r_mean[i] <- mean(Py1r[ , i])
  py1r_sd[i] <- sd(py1r[ , i])
  Py1r_sd[i] <- sd(Py1r[ , i])
  py2r_mean[i] <- mean(py2r[ , i])
  Py2r_mean[i] <- mean(Py2r[ , i])
  py2r_sd[i] <- sd(py2r[ , i])
  Py2r_sd[i] <- sd(Py2r[ , i])
  pzr_mean[i] <- mean(pzr[ , i])
  Pzr_mean[i] <- mean(Pzr[ , i])
  pzr_sd[i] <- sd(pzr[ , i])
  Pzr_sd[i] <- sd(Pzr[ , i])
}

plot(x=D, y=pxd, type="b", col="lightgreen", xlim=c(0,1), ylim=c(0,1), xlab="D", ylab="p", pch = 20, main=paste("DR =", DR))
lines(x=R, y=pxr_mean, type="b", col="green", pch = 20)
lines(x=D, y=py1d, type="b", col="lightblue", pch = 20)
lines(x=R, y=py1r_mean, type="b", col="blue", pch = 20)
lines(x=D, y=py2d, type="b", col="grey", pch = 20)
lines(x=R, y=py2r_mean, type="b", col="black", pch = 20)
lines(x=D, y=pzd, type="b", col="lightpink", pch = 20)
lines(x=R, y=pzr_mean, type="b", col="red", pch = 20)
abline(h=0, v=0)
legend("topright", legend=c("x-resource", "y1-consumer1", "y2-consumer2", "z-predator"), 
       col=c("green", "blue", "black", "red"), pch = 20, cex=1)
arrows(x0=R, y0=pxr_mean-pxr_sd, x1=R, y1=pxr_mean+pxr_sd, length=0.05, angle=90, code=3, col="green")
arrows(x0=R, y0=py1r_mean-py1r_sd, x1=R, y1=py1r_mean+py1r_sd, length=0.05, angle=90, code=3, col="blue")
arrows(x0=R, y0=py2r_mean-py2r_sd, x1=R, y1=py2r_mean+py2r_sd, length=0.05, angle=90, code=3, col="black")
arrows(x0=R, y0=pzr_mean-pzr_sd, x1=R, y1=pzr_mean+pzr_sd, length=0.05, angle=90, code=3, col="red")


# calculate restoration efficiency

# area below curves
for(i in 1:length(D)){
  if(D[i]==DR) { 
    Di <- i 
    break }
}
ADx <- 0
ARx <- 0
ADy1 <- 0
ARy1 <- 0
ADy2 <- 0
ARy2 <- 0
ADz <- 0
ARz <- 0
for(i in 2:Di){
  ADx <- ADx + 0.5*(pxd[i-1]+pxd[i])*dD
  ARx <- ARx + 0.5*(pxr_mean[i-1]+pxr_mean[i])*dR
  ADy1 <- ADy1 + 0.5*(py1d[i-1]+py1d[i])*dD
  ARy1 <- ARy1 + 0.5*(py1r_mean[i-1]+py1r_mean[i])*dR
  ADy2 <- ADy2 + 0.5*(py2d[i-1]+py2d[i])*dD
  ARy2 <- ARy2 + 0.5*(py2r_mean[i-1]+py2r_mean[i])*dR
  ADz <- ADz + 0.5*(pzd[i-1]+pzd[i])*dD
  ARz <- ARz + 0.5*(pzr_mean[i-1]+pzr_mean[i])*dR
}

# restoration efficiency
Rex[Di] <- -(ADx-ARx)/ADx
Rey1[Di] <- -(ADy1-ARy1)/ADy1
Rey2[Di] <- -(ADy2-ARy2)/ADy2
Rez[Di] <- -(ADz-ARz)/ADz

plot(x=D, y=Rex, type="b", col="green", xlim=c(0,1), xlab="DR", ylab="Re", pch = 19)
lines(x=D, y=Rey1, type="b", col="blue", pch = 19)
lines(x=D, y=Rey2, type="b", col="black", pch = 19)
lines(x=D, y=Rez, type="b", col="red", pch = 19)
legend("bottomleft", c("x", "y1", "y2", "z"), col=c("green", "blue", "black", "red"), pch=19, cex=1)
abline(h=0, v=0)


# store results for each DR
R_all[Di,1:length(R)] <- R
pxr_all[Di,1:length(pxr_mean)] <- pxr_mean
Pxr_all[Di,1:length(Pxr_mean)] <- Pxr_mean
py1r_all[Di,1:length(py1r_mean)] <- py1r_mean
Py1r_all[Di,1:length(Py1r_mean)] <- Py1r_mean
py2r_all[Di,1:length(py2r_mean)] <- py2r_mean
Py2r_all[Di,1:length(Py2r_mean)] <- Py2r_mean
pzr_all[Di,1:length(pzr_mean)] <- pzr_mean
Pzr_all[Di,1:length(Pzr_mean)] <- Pzr_mean
sdx_all[Di,1:length(pxr_mean)] <- pxr_sd
sdy1_all[Di,1:length(py1r_mean)] <- py1r_sd
sdy2_all[Di,1:length(py2r_mean)] <- py2r_sd
sdz_all[Di,1:length(pzr_mean)] <- pzr_sd
