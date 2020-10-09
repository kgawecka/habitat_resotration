# stochastoc cellular automata
# trophic OMN spatially explicit model
# random habitat resotration

# species:
# x - resource
# y - consumer
# z - predator

# patch values:
# -1 - destroyed
# 0 - empty
# 1 - x only
# 2 - y only
# 3 - z only
# 4 - x & y
# 5 - x & z
# 6 - y & z
# 7 - x & y & z

rm(list=ls())
library(RColorBrewer)


# HABITAT DESTRUCTION ----

# input - habitat destruction
cx <- 0.7         # probablity of colonisation - resource
cy <- 0.7         # probablity of colonisation - consumer
cz <- 0.7         # probablity of colonisation - predator
ex <- 0.1         # probability of extinction - resource
ey <- 0.1         # probability of extinction - consumer
ez <- 0.1         # probability of extinction - predator

mu1 <- 0.45       # probablilty of resource extinction due to consumer presence
mu2 <- 0.45       # probablilty of resource extinction due to predator presence
mu3 <- 0.45       # probablilty of consumer extinction due to predator presence
psi1 <- 0.4       # probablilty of consumer extinction due lack of resource
psi2 <- 0.4       # probablilty of predator extinction due lack of resource
psi3 <- 0.4       # probablilty of predator extinction due lack of consumer

px0 <- 0.5        # initial fraction of occupied sites - resource
py0 <- 0.5        # initial fraction of occupied sites - consumer
pz0 <- 0.5        # initial fraction of occupied sites - predator

n <- 100          # grid size
tmax <- 50        # maximum number of timesteps
dD <- 0.01        # fraction of patches destoryed at a time


# dynamics with habitat destruction

dynamicsD <- function(cx, cy, cz, ex, ey, ez, mu1, mu2, mu3, psi1, psi2, psi3, px0, py0, pz0, n, tmax, dD) {
  
  t <- 0:tmax
  pnx0 <- round(px0*n*n, digits=0)   # initial number of occupied sites - resource
  if(pnx0<1) { return("ERROR: initial number of resource occupied sites = 0") }
  pny0 <- round(py0*n*n, digits=0)   # initial number of occupied sites - consumer
  if(pny0<1) { return("ERROR: initial number of consumer occupied sites = 0") }
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
  # consumer
  for(i in 1:pny0){
    # select random patch
    i_o <- sample(1:n,1)
    j_o <- sample(1:n,1)
    while(xt[i_o,j_o,1]==2 || xt[i_o,j_o,1]==4){
      i_o <- sample(1:n,1)
      j_o <- sample(1:n,1)
    }
    if     (xt[i_o,j_o,1]==0) { xt[i_o,j_o,1] <- 2 }
    else if(xt[i_o,j_o,1]==1) { xt[i_o,j_o,1] <- 4 }
  }
  # predator
  for(i in 1:pnz0){
    # select random patch
    i_o <- sample(1:n,1)
    j_o <- sample(1:n,1)
    while(xt[i_o,j_o,1]==3 || xt[i_o,j_o,1]==5 || xt[i_o,j_o,1]==6 || xt[i_o,j_o,1]==7){
      i_o <- sample(1:n,1)
      j_o <- sample(1:n,1)
    }
    if     (xt[i_o,j_o,1]==0) { xt[i_o,j_o,1] <- 3 }
    else if(xt[i_o,j_o,1]==1) { xt[i_o,j_o,1] <- 5 }
    else if(xt[i_o,j_o,1]==2) { xt[i_o,j_o,1] <- 6 }
    else if(xt[i_o,j_o,1]==4) { xt[i_o,j_o,1] <- 7 }
  }
  
  pxdt <- matrix(0, nrow=length(d), ncol=length(t))    # fraction of occupied sites as function of time and D
  pxdt[1,1] <- (sum(xt[ , , 1]==1)+sum(xt[ , , 1]==4)+sum(xt[ , , 1]==5)+sum(xt[ , , 1]==7))/(n*n)
  Pxdt <- matrix(0, nrow=length(d), ncol=length(t))    # fraction of occupied sites out of non-destroyed sites as function of time and D
  Pxdt[1,1] <- (sum(xt[ , , 1]==1)+sum(xt[ , , 1]==4)+sum(xt[ , , 1]==5)+sum(xt[ , , 1]==7))/(n*n)
  
  pydt <- matrix(0, nrow=length(d), ncol=length(t))
  pydt[1,1] <- (sum(xt[ , , 1]==2)+sum(xt[ , , 1]==4)+sum(xt[ , , 1]==6)+sum(xt[ , , 1]==7))/(n*n)
  Pydt <- matrix(0, nrow=length(d), ncol=length(t))
  Pydt[1,1] <- (sum(xt[ , , 1]==2)+sum(xt[ , , 1]==4)+sum(xt[ , , 1]==6)+sum(xt[ , , 1]==7))/(n*n)
  
  pzdt <- matrix(0, nrow=length(d), ncol=length(t))
  pzdt[1,1] <- (sum(xt[ , , 1]==3)+sum(xt[ , , 1]==5)+sum(xt[ , , 1]==6)+sum(xt[ , , 1]==7))/(n*n)
  Pzdt <- matrix(0, nrow=length(d), ncol=length(t))
  Pzdt[1,1] <- (sum(xt[ , , 1]==3)+sum(xt[ , , 1]==5)+sum(xt[ , , 1]==6)+sum(xt[ , , 1]==7))/(n*n)
  xdeq <- array(0, dim=c(n, n, length(d)))

  # calculate terms for transition probabilities
  EX1 <- c(ex, (1-ex), ex, ex, (1-ex), (1-ex), ex, (1-ex))
  EX2 <- c((ex+mu1), (1-ex-mu1), (ex+mu1), (ex+mu1), (1-ex-mu1), (1-ex-mu1), (ex+mu1), (1-ex-mu1))
  EX3 <- c((ex+mu2), (1-ex-mu2), (ex+mu2), (ex+mu2), (1-ex-mu2), (1-ex-mu2), (ex+mu2), (1-ex-mu2))
  EX4 <- c((ex+mu1+mu2), (1-ex-mu1-mu2), (ex+mu1+mu2), (ex+mu1+mu2), (1-ex-mu1-mu2), (1-ex-mu1-mu2), (ex+mu1+mu2), (1-ex-mu1-mu2))
  EY1 <- c((ey+psi1), ey, (1-ey-psi1), (ey+psi1), (1-ey), ey, (1-ey-psi1), (1-ey))
  EY2 <- c((ey+psi1+mu3), (ey+mu3), (1-ey-psi1-mu3), (ey+psi1+mu3), (1-ey-mu3), (ey+mu3), (1-ey-psi1-mu3), (1-ey-mu3))
  EZ  <- c((ez+psi2/2+psi3/2), (ez+psi3/2), (ez+psi2/2), (1-ez-psi2/2-psi3/3), ez, (1-ez-psi3/2), (1-ez-psi2/2), (1-ez))
  
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
      pxdt[k,1] <- (sum(xt[ , , 1]==1)+sum(xt[ , , 1]==4)+sum(xt[ , , 1]==5)+sum(xt[ , , 1]==7))/(n*n)
      Pxdt[k,1] <- (sum(xt[ , , 1]==1)+sum(xt[ , , 1]==4)+sum(xt[ , , 1]==5)+sum(xt[ , , 1]==7))/(n*n-d[k])
      pydt[k,1] <- (sum(xt[ , , 1]==2)+sum(xt[ , , 1]==4)+sum(xt[ , , 1]==6)+sum(xt[ , , 1]==7))/(n*n)
      Pydt[k,1] <- (sum(xt[ , , 1]==2)+sum(xt[ , , 1]==4)+sum(xt[ , , 1]==6)+sum(xt[ , , 1]==7))/(n*n-d[k])
      pzdt[k,1] <- (sum(xt[ , , 1]==3)+sum(xt[ , , 1]==5)+sum(xt[ , , 1]==6)+sum(xt[ , , 1]==7))/(n*n)
      Pzdt[k,1] <- (sum(xt[ , , 1]==3)+sum(xt[ , , 1]==5)+sum(xt[ , , 1]==6)+sum(xt[ , , 1]==7))/(n*n-d[k])
    }
    
    # iterate until steady state
    for(g in 2:length(t)){
      
      for(i in 1:n){
        for(j in 1:n){
          
          # von neumann neighborhood
          nx <- 0    # number of neighbours with x
          ny <- 0    # number of neighbours with y
          nz <- 0    # number of neighbours with z
          neigh <- matrix(c(-1, 0, 1, 0, 0, 1, 0, -1), nrow=4, ncol=2)
          for(h in 1:4){
            i_n <- i+neigh[h,1]
            j_n <- j+neigh[h,2]
            
            if(i_n < 1 || i_n > n || j_n < 1 || j_n > n) { next } # neighbour outside of grid
            
            if (xt[i_n,j_n,g-1]==1 || xt[i_n,j_n,g-1]==4 || xt[i_n,j_n,g-1]==5 || xt[i_n,j_n,g-1]==7) { nx <- nx + 1 }
            if (xt[i_n,j_n,g-1]==2 || xt[i_n,j_n,g-1]==4 || xt[i_n,j_n,g-1]==6 || xt[i_n,j_n,g-1]==7) { ny <- ny + 1 }
            if (xt[i_n,j_n,g-1]==3 || xt[i_n,j_n,g-1]==5 || xt[i_n,j_n,g-1]==6 || xt[i_n,j_n,g-1]==7) { nz <- nz + 1 }
          }
          
          # calculate terms for transition probabilities
          CX <- c((1-cx)^nx, (1-(1-cx)^nx), (1-cx)^nx, (1-cx)^nx, (1-(1-cx)^nx), (1-(1-cx)^nx), (1-cx)^nx, (1-(1-cx)^nx))
          CY <- c((1-cy)^ny, (1-cy)^ny, (1-(1-cy)^ny), (1-cy)^ny, (1-(1-cy)^ny), (1-cy)^ny, (1-(1-cy)^ny), (1-(1-cy)^ny))
          CZ <- c((1-cz)^nz, (1-cz)^nz, (1-cz)^nz, (1-(1-cz)^nz), (1-cz)^nz, (1-(1-cz)^nz), (1-(1-cz)^nz), (1-(1-cz)^nz))
          
          # transition probabilities
          p00 <- CX[1] * CY[1] * CZ[1]
          p01 <- CX[2] * CY[2] * CZ[2]
          p02 <- CX[3] * CY[3] * CZ[3]
          p03 <- CX[4] * CY[4] * CZ[4]
          p04 <- CX[5] * CY[5] * CZ[5]
          p05 <- CX[6] * CY[6] * CZ[6]
          p06 <- CX[7] * CY[7] * CZ[7]
          p07 <- CX[8] * CY[8] * CZ[8]
          
          p10 <- EX1[1] * CY[1] * CZ[1]
          p11 <- EX1[2] * CY[2] * CZ[2]
          p12 <- EX1[3] * CY[3] * CZ[3]
          p13 <- EX1[4] * CY[4] * CZ[4]
          p14 <- EX1[5] * CY[5] * CZ[5]
          p15 <- EX1[6] * CY[6] * CZ[6]
          p16 <- EX1[7] * CY[7] * CZ[7]
          p17 <- EX1[8] * CY[8] * CZ[8]
          
          p20 <- CX[1] * EY1[1] * CZ[1]
          p21 <- CX[2] * EY1[2] * CZ[2]
          p22 <- CX[3] * EY1[3] * CZ[3]
          p23 <- CX[4] * EY1[4] * CZ[4]
          p24 <- CX[5] * EY1[5] * CZ[5]
          p25 <- CX[6] * EY1[6] * CZ[6]
          p26 <- CX[7] * EY1[7] * CZ[7]
          p27 <- CX[8] * EY1[8] * CZ[8]
          
          p30 <- CX[1] * CY[1] * EZ[1]
          p31 <- CX[2] * CY[2] * EZ[2]
          p32 <- CX[3] * CY[3] * EZ[3]
          p33 <- CX[4] * CY[4] * EZ[4]
          p34 <- CX[5] * CY[5] * EZ[5]
          p35 <- CX[6] * CY[6] * EZ[6]
          p36 <- CX[7] * CY[7] * EZ[7]
          p37 <- CX[8] * CY[8] * EZ[8]
          
          p40 <- EX2[1] * EY1[1] * CZ[1]
          p41 <- EX2[2] * EY1[2] * CZ[2]
          p42 <- EX2[3] * EY1[3] * CZ[3]
          p43 <- EX2[4] * EY1[4] * CZ[4]
          p44 <- EX2[5] * EY1[5] * CZ[5]
          p45 <- EX2[6] * EY1[6] * CZ[6]
          p46 <- EX2[7] * EY1[7] * CZ[7]
          p47 <- EX2[8] * EY1[8] * CZ[8]
          
          p50 <- EX3[1] * CY[1] * EZ[1]
          p51 <- EX3[2] * CY[2] * EZ[2]
          p52 <- EX3[3] * CY[3] * EZ[3]
          p53 <- EX3[4] * CY[4] * EZ[4]
          p54 <- EX3[5] * CY[5] * EZ[5]
          p55 <- EX3[6] * CY[6] * EZ[6]
          p56 <- EX3[7] * CY[7] * EZ[7]
          p57 <- EX3[8] * CY[8] * EZ[8]
          
          p60 <- CX[1] * EY2[1] * EZ[1]
          p61 <- CX[2] * EY2[2] * EZ[2]
          p62 <- CX[3] * EY2[3] * EZ[3]
          p63 <- CX[4] * EY2[4] * EZ[4]
          p64 <- CX[5] * EY2[5] * EZ[5]
          p65 <- CX[6] * EY2[6] * EZ[6]
          p66 <- CX[7] * EY2[7] * EZ[7]
          p67 <- CX[8] * EY2[8] * EZ[8]
          
          p70 <- EX4[1] * EY2[1] * EZ[1]
          p71 <- EX4[2] * EY2[2] * EZ[2]
          p72 <- EX4[3] * EY2[3] * EZ[3]
          p73 <- EX4[4] * EY2[4] * EZ[4]
          p74 <- EX4[5] * EY2[5] * EZ[5]
          p75 <- EX4[6] * EY2[6] * EZ[6]
          p76 <- EX4[7] * EY2[7] * EZ[7]
          p77 <- EX4[8] * EY2[8] * EZ[8]
          
          # transitions
          if(xt[i,j,g-1]==-1) { xt[i,j,g] <- -1 }    # destroyed sites remain destroyed
          
          else if(xt[i,j,g-1]==0) { xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7),1, prob=c(p00,p01,p02,p03,p04,p05,p06,p07)) }
          
          else if(xt[i,j,g-1]==1) { xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7),1, prob=c(p10,p11,p12,p13,p14,p15,p16,p17)) }
          
          else if(xt[i,j,g-1]==2) { xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7),1, prob=c(p20,p21,p22,p23,p24,p25,p26,p27)) }
          
          else if(xt[i,j,g-1]==3) { xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7),1, prob=c(p30,p31,p32,p33,p34,p35,p36,p37)) }
          
          else if(xt[i,j,g-1]==4) { xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7),1, prob=c(p40,p41,p42,p43,p44,p45,p46,p47)) }
          
          else if(xt[i,j,g-1]==5) { xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7),1, prob=c(p50,p51,p52,p53,p54,p55,p56,p57)) }
          
          else if(xt[i,j,g-1]==6) { xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7),1, prob=c(p60,p61,p62,p63,p64,p65,p66,p67)) }
          
          else if(xt[i,j,g-1]==7) { xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7),1, prob=c(p70,p71,p72,p73,p74,p75,p76,p77)) }
        }
      }
      
      pxdt[k,g] <- (sum(xt[ , , g]==1)+sum(xt[ , , g]==4)+sum(xt[ , , g]==5)+sum(xt[ , , g]==7))/(n*n)
      Pxdt[k,g] <- (sum(xt[ , , g]==1)+sum(xt[ , , g]==4)+sum(xt[ , , g]==5)+sum(xt[ , , g]==7))/(n*n-d[k])
      
      pydt[k,g] <- (sum(xt[ , , g]==2)+sum(xt[ , , g]==4)+sum(xt[ , , g]==6)+sum(xt[ , , g]==7))/(n*n)
      Pydt[k,g] <- (sum(xt[ , , g]==2)+sum(xt[ , , g]==4)+sum(xt[ , , g]==6)+sum(xt[ , , g]==7))/(n*n-d[k])
      
      pzdt[k,g] <- (sum(xt[ , , g]==3)+sum(xt[ , , g]==5)+sum(xt[ , , g]==6)+sum(xt[ , , g]==7))/(n*n)
      Pzdt[k,g] <- (sum(xt[ , , g]==3)+sum(xt[ , , g]==5)+sum(xt[ , , g]==6)+sum(xt[ , , g]==7))/(n*n-d[k])
    }
    xdeq[ , ,k] <- xt[ , ,length(t)]
  }
  
  output <- list("D"=D, "pxdt"=pxdt, "Pxdt"=Pxdt, "pydt"=pydt, "Pydt"=Pydt, "pzdt"=pzdt, "Pzdt"=Pzdt, "xdeq"=xdeq)
  return(output)
}

# habitat destruction output
outD <- dynamicsD(cx, cy, cz, ex, ey, ez, mu1, mu2, mu3, psi1, psi2, psi3, px0, py0, pz0, n, tmax, dD)

D <- outD$D           # fraction of habitat destoryed [0,1]
pxdt <- outD$pxdt     # fraction of x-occupied sites as function of time and D
Pxdt <- outD$Pxdt     # fraction of x-occupied sites out of non-destroyed sites as function of time and D
pydt <- outD$pydt     # fraction of y-occupied sites as function of time and D
Pydt <- outD$Pydt     # fraction of y-occupied sites out of non-destroyed sites as function of time and D
pzdt <- outD$pzdt     # fraction of z-occupied sites as function of time and D
Pzdt <- outD$Pzdt     # fraction of z-occupied sites out of non-destroyed sites as function of time and D
xdeq <- outD$xdeq     # equilibrium grid as function of D

plot(x=D, y=pxdt[ ,tmax+1], type="p", col="green", ylim=c(0,1), ylab="p", pch=20)
lines(x=D, y=pydt[ ,tmax+1], type="p", col="blue", pch=20)
lines(x=D, y=pzdt[ ,tmax+1], type="p", col="red", pch=20)
legend(0.7, 1, legend=c("x-resource", "y-consumer", "z-predator"), col=c("green", "blue", "red"), pch=20, cex=1)
abline(h=0, v=0)


# HABITAT RESTORATION ----

# input - habitat restoration
cx <- 0.7         # probablity of colonisation - resource
cy <- 0.7         # probablity of colonisation - consumer
cz <- 0.7         # probablity of colonisation - predator
ex <- 0.1         # probability of extinction - resource
ey <- 0.1         # probability of extinction - consumer
ez <- 0.1         # probability of extinction - predator

mu1 <- 0.45       # probablilty of resource extinction due to consumer presence
mu2 <- 0.45       # probablilty of resource extinction due to predator presence
mu3 <- 0.45       # probablilty of consumer extinction due to predator presence
psi1 <- 0.4       # probablilty of consumer extinction due lack of resource
psi2 <- 0.4       # probablilty of predator extinction due lack of resource
psi3 <- 0.4       # probablilty of predator extinction due lack of consumer

tmax <- 50        # maximum number of timesteps
DR <- 0.10        # fraction of patches destroyed at the start of restoration
dR <- 0.01        # fraction of patches restored at a time

Rex  <- rep(NA,length(D))     # initialise restoration efficiency vector
Rey  <- rep(NA,length(D))
Rez  <- rep(NA,length(D))

# initialise matrices for storing results
R_all <- matrix(NA, nrow=length(D), ncol=length(D))      # fraction of habitat loss during resotration
pxr_all <- matrix(NA, nrow=length(D), ncol=length(D))    # fraction of x-occupied sites (mean of replicates)
Pxr_all <- matrix(NA, nrow=length(D), ncol=length(D))    # fraction of x-occupied sites out of non-destroyed sites (mean of replicates)
pyr_all <- matrix(NA, nrow=length(D), ncol=length(D))    # fraction of y-occupied sites (mean of replicates)
Pyr_all <- matrix(NA, nrow=length(D), ncol=length(D))    # fraction of y-occupied sites out of non-destroyed sites (mean of replicates)
pzr_all <- matrix(NA, nrow=length(D), ncol=length(D))    # fraction of z-occupied sites (mean of replicates)
Pzr_all <- matrix(NA, nrow=length(D), ncol=length(D))    # fraction of z-occupied sites out of non-destroyed sites (mean of replicates)
sdx_all <- matrix(NA, nrow=length(D), ncol=length(D))    # st dev of fraction of x-occupied sites of replicates
sdy_all <- matrix(NA, nrow=length(D), ncol=length(D))    # st dev of fraction of y-occupied sites of replicates
sdz_all <- matrix(NA, nrow=length(D), ncol=length(D))    # st dev of fraction of z-occupied sites of replicates


# dynamics with habitat restoration

dynamicsR <- function(cx, cy, cz, ex, ey, ez, mu1, mu2, mu3, psi1, psi2, psi3, tmax, DR, dR) {
  
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
  pzrt <- matrix(0, nrow=length(r), ncol=length(t))    
  pzrt[1, ] <- pzdt[Di, ]
  Pzrt <- matrix(0, nrow=length(r), ncol=length(t))
  Pzrt[1, ] <- Pzdt[Di, ]
  
  # calculate terms for transition probabilities
  EX1 <- c(ex, (1-ex), ex, ex, (1-ex), (1-ex), ex, (1-ex))
  EX2 <- c((ex+mu1), (1-ex-mu1), (ex+mu1), (ex+mu1), (1-ex-mu1), (1-ex-mu1), (ex+mu1), (1-ex-mu1))
  EX3 <- c((ex+mu2), (1-ex-mu2), (ex+mu2), (ex+mu2), (1-ex-mu2), (1-ex-mu2), (ex+mu2), (1-ex-mu2))
  EX4 <- c((ex+mu1+mu2), (1-ex-mu1-mu2), (ex+mu1+mu2), (ex+mu1+mu2), (1-ex-mu1-mu2), (1-ex-mu1-mu2), (ex+mu1+mu2), (1-ex-mu1-mu2))
  EY1 <- c((ey+psi1), ey, (1-ey-psi1), (ey+psi1), (1-ey), ey, (1-ey-psi1), (1-ey))
  EY2 <- c((ey+psi1+mu3), (ey+mu3), (1-ey-psi1-mu3), (ey+psi1+mu3), (1-ey-mu3), (ey+mu3), (1-ey-psi1-mu3), (1-ey-mu3))
  EZ  <- c((ez+psi2/2+psi3/2), (ez+psi3/2), (ez+psi2/2), (1-ez-psi2/2-psi3/3), ez, (1-ez-psi3/2), (1-ez-psi2/2), (1-ez))
  
  for(k in 2:length(r)){
    
    # initialise grid
    xreq[ , ,k] <- xreq[ , ,k-1]
    
    # select random patch from destroyed patches
    for(h in 1:dn) {
      i_d <- sample(1:n,1)
      j_d <- sample(1:n,1)
      while(xreq[i_d,j_d,k]!=-1){
        i_d <- sample(1:n,1)
        j_d <- sample(1:n,1)
      }
      xreq[i_d,j_d,k] <- 0
    }
    
    # grid at t=0 for d[k]
    xt[ , ,1] <- xreq[ , ,k]
    pxrt[k,1] <- (sum(xt[ , , 1]==1)+sum(xt[ , , 1]==4)+sum(xt[ , , 1]==5)+sum(xt[ , , 1]==7))/(n*n)
    Pxrt[k,1] <- (sum(xt[ , , 1]==1)+sum(xt[ , , 1]==4)+sum(xt[ , , 1]==5)+sum(xt[ , , 1]==7))/(n*n-r[k])
    pyrt[k,1] <- (sum(xt[ , , 1]==2)+sum(xt[ , , 1]==4)+sum(xt[ , , 1]==6)+sum(xt[ , , 1]==7))/(n*n)
    Pyrt[k,1] <- (sum(xt[ , , 1]==2)+sum(xt[ , , 1]==4)+sum(xt[ , , 1]==6)+sum(xt[ , , 1]==7))/(n*n-r[k])
    pzrt[k,1] <- (sum(xt[ , , 1]==3)+sum(xt[ , , 1]==5)+sum(xt[ , , 1]==6)+sum(xt[ , , 1]==7))/(n*n)
    Pzrt[k,1] <- (sum(xt[ , , 1]==3)+sum(xt[ , , 1]==5)+sum(xt[ , , 1]==6)+sum(xt[ , , 1]==7))/(n*n-r[k])
    
    # iterate until steady state
    for(g in 2:length(t)){
      for(i in 1:n){
        for(j in 1:n){
          
          # von neumann neighborhood
          nx <- 0    # number of neighbours with x
          ny <- 0    # number of neighbours with y
          nz <- 0    # number of neighbours with z
          neigh <- matrix(c(-1, 0, 1, 0, 0, 1, 0, -1), nrow=4, ncol=2)
          for(h in 1:4){
            i_n <- i+neigh[h,1]
            j_n <- j+neigh[h,2]
            
            if(i_n < 1 || i_n > n || j_n < 1 || j_n > n) { next } # neighbour outside of grid
            
            if (xt[i_n,j_n,g-1]==1 || xt[i_n,j_n,g-1]==4 || xt[i_n,j_n,g-1]==5 || xt[i_n,j_n,g-1]==7) { nx <- nx + 1 }
            if (xt[i_n,j_n,g-1]==2 || xt[i_n,j_n,g-1]==4 || xt[i_n,j_n,g-1]==6 || xt[i_n,j_n,g-1]==7) { ny <- ny + 1 }
            if (xt[i_n,j_n,g-1]==3 || xt[i_n,j_n,g-1]==5 || xt[i_n,j_n,g-1]==6 || xt[i_n,j_n,g-1]==7) { nz <- nz + 1 }
          }
          
          # calculate terms for transition probabilities
          CX <- c((1-cx)^nx, (1-(1-cx)^nx), (1-cx)^nx, (1-cx)^nx, (1-(1-cx)^nx), (1-(1-cx)^nx), (1-cx)^nx, (1-(1-cx)^nx))
          CY <- c((1-cy)^ny, (1-cy)^ny, (1-(1-cy)^ny), (1-cy)^ny, (1-(1-cy)^ny), (1-cy)^ny, (1-(1-cy)^ny), (1-(1-cy)^ny))
          CZ <- c((1-cz)^nz, (1-cz)^nz, (1-cz)^nz, (1-(1-cz)^nz), (1-cz)^nz, (1-(1-cz)^nz), (1-(1-cz)^nz), (1-(1-cz)^nz))
          
          # transition probabilities
          p00 <- CX[1] * CY[1] * CZ[1]
          p01 <- CX[2] * CY[2] * CZ[2]
          p02 <- CX[3] * CY[3] * CZ[3]
          p03 <- CX[4] * CY[4] * CZ[4]
          p04 <- CX[5] * CY[5] * CZ[5]
          p05 <- CX[6] * CY[6] * CZ[6]
          p06 <- CX[7] * CY[7] * CZ[7]
          p07 <- CX[8] * CY[8] * CZ[8]
          
          p10 <- EX1[1] * CY[1] * CZ[1]
          p11 <- EX1[2] * CY[2] * CZ[2]
          p12 <- EX1[3] * CY[3] * CZ[3]
          p13 <- EX1[4] * CY[4] * CZ[4]
          p14 <- EX1[5] * CY[5] * CZ[5]
          p15 <- EX1[6] * CY[6] * CZ[6]
          p16 <- EX1[7] * CY[7] * CZ[7]
          p17 <- EX1[8] * CY[8] * CZ[8]
          
          p20 <- CX[1] * EY1[1] * CZ[1]
          p21 <- CX[2] * EY1[2] * CZ[2]
          p22 <- CX[3] * EY1[3] * CZ[3]
          p23 <- CX[4] * EY1[4] * CZ[4]
          p24 <- CX[5] * EY1[5] * CZ[5]
          p25 <- CX[6] * EY1[6] * CZ[6]
          p26 <- CX[7] * EY1[7] * CZ[7]
          p27 <- CX[8] * EY1[8] * CZ[8]
          
          p30 <- CX[1] * CY[1] * EZ[1]
          p31 <- CX[2] * CY[2] * EZ[2]
          p32 <- CX[3] * CY[3] * EZ[3]
          p33 <- CX[4] * CY[4] * EZ[4]
          p34 <- CX[5] * CY[5] * EZ[5]
          p35 <- CX[6] * CY[6] * EZ[6]
          p36 <- CX[7] * CY[7] * EZ[7]
          p37 <- CX[8] * CY[8] * EZ[8]
          
          p40 <- EX2[1] * EY1[1] * CZ[1]
          p41 <- EX2[2] * EY1[2] * CZ[2]
          p42 <- EX2[3] * EY1[3] * CZ[3]
          p43 <- EX2[4] * EY1[4] * CZ[4]
          p44 <- EX2[5] * EY1[5] * CZ[5]
          p45 <- EX2[6] * EY1[6] * CZ[6]
          p46 <- EX2[7] * EY1[7] * CZ[7]
          p47 <- EX2[8] * EY1[8] * CZ[8]
          
          p50 <- EX3[1] * CY[1] * EZ[1]
          p51 <- EX3[2] * CY[2] * EZ[2]
          p52 <- EX3[3] * CY[3] * EZ[3]
          p53 <- EX3[4] * CY[4] * EZ[4]
          p54 <- EX3[5] * CY[5] * EZ[5]
          p55 <- EX3[6] * CY[6] * EZ[6]
          p56 <- EX3[7] * CY[7] * EZ[7]
          p57 <- EX3[8] * CY[8] * EZ[8]
          
          p60 <- CX[1] * EY2[1] * EZ[1]
          p61 <- CX[2] * EY2[2] * EZ[2]
          p62 <- CX[3] * EY2[3] * EZ[3]
          p63 <- CX[4] * EY2[4] * EZ[4]
          p64 <- CX[5] * EY2[5] * EZ[5]
          p65 <- CX[6] * EY2[6] * EZ[6]
          p66 <- CX[7] * EY2[7] * EZ[7]
          p67 <- CX[8] * EY2[8] * EZ[8]
          
          p70 <- EX4[1] * EY2[1] * EZ[1]
          p71 <- EX4[2] * EY2[2] * EZ[2]
          p72 <- EX4[3] * EY2[3] * EZ[3]
          p73 <- EX4[4] * EY2[4] * EZ[4]
          p74 <- EX4[5] * EY2[5] * EZ[5]
          p75 <- EX4[6] * EY2[6] * EZ[6]
          p76 <- EX4[7] * EY2[7] * EZ[7]
          p77 <- EX4[8] * EY2[8] * EZ[8]
          
          # transitions
          if(xt[i,j,g-1]==-1) { xt[i,j,g] <- -1 }    # destroyed sites remain destroyed
          
          else if(xt[i,j,g-1]==0) { xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7),1, prob=c(p00,p01,p02,p03,p04,p05,p06,p07)) }
          
          else if(xt[i,j,g-1]==1) { xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7),1, prob=c(p10,p11,p12,p13,p14,p15,p16,p17)) }
          
          else if(xt[i,j,g-1]==2) { xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7),1, prob=c(p20,p21,p22,p23,p24,p25,p26,p27)) }
          
          else if(xt[i,j,g-1]==3) { xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7),1, prob=c(p30,p31,p32,p33,p34,p35,p36,p37)) }
          
          else if(xt[i,j,g-1]==4) { xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7),1, prob=c(p40,p41,p42,p43,p44,p45,p46,p47)) }
          
          else if(xt[i,j,g-1]==5) { xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7),1, prob=c(p50,p51,p52,p53,p54,p55,p56,p57)) }
          
          else if(xt[i,j,g-1]==6) { xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7),1, prob=c(p60,p61,p62,p63,p64,p65,p66,p67)) }
          
          else if(xt[i,j,g-1]==7) { xt[i,j,g] <- sample(c(0,1,2,3,4,5,6,7),1, prob=c(p70,p71,p72,p73,p74,p75,p76,p77)) }
        }
      }
      pxrt[k,g] <- (sum(xt[ , , g]==1)+sum(xt[ , , g]==4)+sum(xt[ , , g]==5)+sum(xt[ , , g]==7))/(n*n)
      Pxrt[k,g] <- (sum(xt[ , , g]==1)+sum(xt[ , , g]==4)+sum(xt[ , , g]==5)+sum(xt[ , , g]==7))/(n*n-r[k])
      
      pyrt[k,g] <- (sum(xt[ , , g]==2)+sum(xt[ , , g]==4)+sum(xt[ , , g]==6)+sum(xt[ , , g]==7))/(n*n)
      Pyrt[k,g] <- (sum(xt[ , , g]==2)+sum(xt[ , , g]==4)+sum(xt[ , , g]==6)+sum(xt[ , , g]==7))/(n*n-r[k])
      
      pzrt[k,g] <- (sum(xt[ , , g]==3)+sum(xt[ , , g]==5)+sum(xt[ , , g]==6)+sum(xt[ , , g]==7))/(n*n)
      Pzrt[k,g] <- (sum(xt[ , , g]==3)+sum(xt[ , , g]==5)+sum(xt[ , , g]==6)+sum(xt[ , , g]==7))/(n*n-r[k])
    }
    xreq[ , ,k] <- xt[ , ,length(t)]
  }
  output <- list("R"=R, "pxrt"=pxrt, "Pxrt"=Pxrt, "pyrt"=pyrt, "Pyrt"=Pyrt, "pzrt"=pzrt, "Pzrt"=Pzrt, "xreq"=xreq)
  return(output)
}

# habitat restoration output
outR <- dynamicsR(cx, cy, cz, ex, ey, ez, mu1, mu2, mu3, psi1, psi2, psi3, tmax, DR, dR)

R <- outR$R                # fraction of habitat loss [DR,0]
pxrt <- outR$pxrt          # fraction of x-occupied sites as function of time and D
Pxrt <- outR$Pxrt          # fraction of x-occupied sites out of non-destroyed sites as function of time and D
pyrt <- outR$pyrt          # fraction of y-occupied sites as function of time and D
Pyrt <- outR$Pyrt          # fraction of y-occupied sites out of non-destroyed sites as function of time and D
pzrt <- outR$pzrt          # fraction of z-occupied sites as function of time and D
Pzrt <- outR$Pzrt          # fraction of z-occupied sites out of non-destroyed sites as function of time and D
xreq <- outR$xreq          # equilibrium grid as function of D

pxr <- pxrt[, tmax+1]      # equilibrium values only
Pxr <- Pxrt[, tmax+1]
pxd <- pxdt[, tmax+1]
Pxd <- Pxdt[, tmax+1]
pyr <- pyrt[, tmax+1]
Pyr <- Pyrt[, tmax+1]
pyd <- pydt[, tmax+1]
Pyd <- Pydt[, tmax+1]
pzr <- pzrt[, tmax+1]
Pzr <- Pzrt[, tmax+1]
pzd <- pzdt[, tmax+1]
Pzd <- Pzdt[, tmax+1]

plot(x=D, y=pxd, type="b", col="lightgreen", xlim=c(0,1), ylim=c(0,1), xlab="D", ylab="p", pch = 20)
lines(x=R, y=pxr, type="b", col="green", pch = 20)
lines(x=D, y=pyd, type="b", col="lightblue", pch = 20)
lines(x=R, y=pyr, type="b", col="blue", pch = 20)
lines(x=D, y=pzd, type="b", col="lightpink", pch = 20)
lines(x=R, y=pzr, type="b", col="red", pch = 20)
legend(0.7, 1, legend=c("x-resource", "y-consumer", "z-predator"), col=c("green", "blue", "red"), pch = 20, cex=1)
abline(h=0, v=0)

for(k in 1:length(R)){
  image(x=1:n, y=1:n, z=t(apply(xreq[ , ,k], 2, rev)), main=paste("D =", R[k]), xlab="", ylab="", zlim=c(-1,7),
        col=brewer.pal(n=9, name="Set1"))
  legend(grconvertX(1, "device"), grconvertY(1, "device"), 
         c("destroyed","empty","x","y","z","x&y","x&z","y&z","x&y&z"), 
         fill=brewer.pal(n=9, name="Set1"), xpd = NA)
}


# replicates
nr <- 5        # no of replicates
DR <- 0.10     # fraction of patches destroyed at the start of restoration

R <- seq(DR, 0, -dR)
pxr <- matrix(0, nrow=nr, ncol=length(R))
Pxr <- matrix(0, nrow=nr, ncol=length(R))
pyr <- matrix(0, nrow=nr, ncol=length(R))
Pyr <- matrix(0, nrow=nr, ncol=length(R))
pzr <- matrix(0, nrow=nr, ncol=length(R))
Pzr <- matrix(0, nrow=nr, ncol=length(R))
for(i in 1:nr){
  outR <- dynamicsR(cx, cy, cz, ex, ey, ez, mu1, mu2, mu3, psi1, psi2, psi3, tmax, DR, dR)
  pxrt <- outR$pxrt
  Pxrt <- outR$Pxrt
  pxr[i, ] <- pxrt[ ,tmax+1]
  Pxr[i, ] <- Pxrt[ ,tmax+1]
  pyrt <- outR$pyrt
  Pyrt <- outR$Pyrt
  pyr[i, ] <- pyrt[ ,tmax+1]
  Pyr[i, ] <- Pyrt[ ,tmax+1]
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
pyr_mean <- numeric(length(R))
Pyr_mean <- numeric(length(R))
pyr_sd <- numeric(length(R))
Pyr_sd <- numeric(length(R))
pzr_mean <- numeric(length(R))
Pzr_mean <- numeric(length(R))
pzr_sd <- numeric(length(R))
Pzr_sd <- numeric(length(R))
for(i in 1:length(R)) {
  pxr_mean[i] <- mean(pxr[ , i])
  Pxr_mean[i] <- mean(Pxr[ , i])
  pxr_sd[i] <- sd(pxr[ , i])
  Pxr_sd[i] <- sd(Pxr[ , i])
  pyr_mean[i] <- mean(pyr[ , i])
  Pyr_mean[i] <- mean(Pyr[ , i])
  pyr_sd[i] <- sd(pyr[ , i])
  Pyr_sd[i] <- sd(Pyr[ , i])
  pzr_mean[i] <- mean(pzr[ , i])
  Pzr_mean[i] <- mean(Pzr[ , i])
  pzr_sd[i] <- sd(pzr[ , i])
  Pzr_sd[i] <- sd(Pzr[ , i])
}

plot(x=D, y=pxd, type="b", col="lightgreen", xlim=c(0,1), ylim=c(0,1), xlab="D", ylab="p", pch = 20, main=paste("DR =", DR))
lines(x=R, y=pxr_mean, type="b", col="green", pch = 20)
lines(x=D, y=pyd, type="b", col="lightblue", pch = 20)
lines(x=R, y=pyr_mean, type="b", col="blue", pch = 20)
lines(x=D, y=pzd, type="b", col="lightpink", pch = 20)
lines(x=R, y=pzr_mean, type="b", col="red", pch = 20)
abline(h=0, v=0)
legend(0.7, 1, legend=c("x-resource", "y-consumer", "z-predator"), col=c("green", "blue", "red"), pch = 20, cex=1)
arrows(x0=R, y0=pxr_mean-pxr_sd, x1=R, y1=pxr_mean+pxr_sd, length=0.05, angle=90, code=3, col="green")
arrows(x0=R, y0=pyr_mean-pyr_sd, x1=R, y1=pyr_mean+pyr_sd, length=0.05, angle=90, code=3, col="blue")
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
ADy <- 0
ARy <- 0
ADz <- 0
ARz <- 0
for(i in 2:Di){
  ADx <- ADx + 0.5*(pxd[i-1]+pxd[i])*dD
  ARx <- ARx + 0.5*(pxr_mean[i-1]+pxr_mean[i])*dR
  ADy <- ADy + 0.5*(pyd[i-1]+pyd[i])*dD
  ARy <- ARy + 0.5*(pyr_mean[i-1]+pyr_mean[i])*dR
  ADz <- ADz + 0.5*(pzd[i-1]+pzd[i])*dD
  ARz <- ARz + 0.5*(pzr_mean[i-1]+pzr_mean[i])*dR
}

# restoration efficiency
Rex[Di] <- -(ADx-ARx)/ADx
Rey[Di] <- -(ADy-ARy)/ADy
Rez[Di] <- -(ADz-ARz)/ADz

plot(x=D, y=Rex, type="b", col="green", xlim=c(0,1), xlab="DR", ylab="Re", pch = 19)
lines(x=D, y=Rey, type="b", col="blue", pch = 19)
lines(x=D, y=Rez, type="b", col="red", pch = 19)
legend("bottomleft", c("x", "y", "z"), col=c("green", "blue", "red"), pch=19, cex=1)
abline(h=0, v=0)


# store results for each DR
R_all[Di,1:length(R)] <- R
pxr_all[Di,1:length(pxr_mean)] <- pxr_mean
Pxr_all[Di,1:length(Pxr_mean)] <- Pxr_mean
pyr_all[Di,1:length(pyr_mean)] <- pyr_mean
Pyr_all[Di,1:length(Pyr_mean)] <- Pyr_mean
pzr_all[Di,1:length(pzr_mean)] <- pzr_mean
Pzr_all[Di,1:length(Pzr_mean)] <- Pzr_mean
sdx_all[Di,1:length(pxr_mean)] <- pxr_sd
sdy_all[Di,1:length(pyr_mean)] <- pyr_sd
sdz_all[Di,1:length(pzr_mean)] <- pzr_sd
