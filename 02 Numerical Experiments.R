##### Setting working directory
setwd("working/directory")

##### Loading needed packages
library(nvc)
library(latex2exp)

##### User-defined functions for simulation
### For general random vector setting
pcs_gen <- function(d0 = 2,ktau.max=0.90,family_set = c("indep","gaussian","clayton","gumbel","frank","joe")){
  
  pcs <- vector("list",d0-1)
  for(i in 1:length(pcs)){
    pcs[[i]] <- vector("list",d0-i)
  }
  
  for(j in 1:(d0-1)){
    for(i in 1:(d0-j)){
      fam_i <- sample(1:length(family_set),1)
      if(family_set[fam_i] %in% c("indep")){
        pcs[[i]][[j]] <- bicop_dist(family_set[fam_i])
      } else if(family_set[fam_i] %in% c("gaussian","frank")){
        ktau <- runif(1,min=-ktau.max,max=ktau.max)
        cop.par <- ktau_to_par(family = family_set[fam_i],tau = ktau)
        pcs[[i]][[j]] <- bicop_dist(family_set[fam_i], 0, as.numeric(cop.par))
      } else {
        ktau <- runif(1,min=0,max=ktau.max)
        cop.par <- ktau_to_par(family = family_set[fam_i],tau = ktau)
        pcs[[i]][[j]] <- bicop_dist(family_set[fam_i], 0, as.numeric(cop.par))
      }
      
    }
  }
  return(pcs)
}

### For nonlinear spectral vector dependence setting
sim.nl.signal <- function(nsec = 200, samprate = 100, freq = 30,
                          q = 2, noise.sd = 0.5) {
  
  Tlength <- nsec*samprate
  
  X <- matrix(NA, nrow = Tlength, ncol = q)
  Y <- matrix(NA, nrow = Tlength, ncol = q)
  
  ## Latent block-specific spectral powers for X
  PX <- runif(nsec)
  
  ## U-shaped, non-monotone dependence
  PY <- (PX - 0.5)^2
  
  for (b in 1:nsec) {
    
    b.ind <- ((b - 1) * samprate + 1):(b * samprate)
    t.ind <- (1:samprate)/samprate
    
    ## Independent phases between regions
    phi.x <- runif(q, min = 0, max = 2 * pi)
    phi.y <- runif(q, min = 0, max = 2 * pi)
    
    for (j in 1:q) {
      
      X[b.ind, j] <- sqrt(PX[b]) * cos(2 * pi * freq * t.ind + phi.x[j]) +
        rnorm(samprate, sd = noise.sd)
      
      Y[b.ind, j] <- sqrt(PY[b]) * cos(2 * pi * freq * t.ind + phi.y[j]) +
        rnorm(samprate, sd = noise.sd)
    }
    
  }
  
  return(list(X = X,Y = Y))
}



#####################################################
##### Simulation for linear spectral dependence #####
#####################################################

nrep <- 1000
nboot <- 500
ncores <- 8
samprate <- 100
nsec_all <- c(50,100,200)
fband <- c("delta","theta","alpha","beta2","gamma2")


##########
# Case 1 #
##########

c.ind <- 1
rhythms <- c("delta","theta2","alpha","beta2","gamma2")
set.seed(c.ind) # seed for reproducibility

est_comp <- array(NA,dim = c(samprate/2,nrep))
coh_comp <- array(NA,dim = c(samprate/2,nrep))
cancoh_comp <- array(NA,dim = c(samprate/2,nrep))
bandcoh_comp <- array(NA,dim = c(5,nrep))
kencoh_comp <- array(NA,dim = c(5,nrep))
pval_comp <- array(NA,dim = c(samprate/2,nrep))
pval.adj_comp <- array(NA,dim = c(samprate/2,nrep))

for(n.ind in 1:length(nsec_all)){
  nsec <- nsec_all[n.ind]
  Tlength <- nsec*samprate
  
  for(r in 1:nrep){
    
    Z_alpha <- sim.fband(nsec = nsec,samprate = samprate,fband = rep(fband[3],4))
    
    x1 <- 0.5*Z_alpha[,1] + 0.5*rnorm(Tlength)
    x2 <- 0.5*Z_alpha[,2] + 0.5*rnorm(Tlength)
    y1 <- 0.5*Z_alpha[,1] + 0.5*rnorm(Tlength)
    y2 <- 0.5*Z_alpha[,2] + 0.5*rnorm(Tlength)
    x <- cbind(x1,x2)
    y <- cbind(y1,y2)
    
    samp <- nvc(x = x, y = y,blength = samprate,perm.inv = TRUE,
                b.scale = FALSE,nboot = nboot,ncores = ncores)
    est_comp[,r] <- samp$nvc
    pval_comp[,r] <- samp$pval
    pval.adj_comp[,r] <- samp$pval.adj
    
    coh.all <- array(NA,dim = c(samprate/2,ncol(x)*ncol(y)))
    bandcoh.all <- array(NA,dim = c(5,ncol(x)*ncol(y)))
    
    for(ch1 in 1:ncol(x)){
      for(ch2 in 1:ncol(y)){
        coh.all[,ncol(y)*(ch1-1)+ch2] <- coh(x = x[,ch1],y = y[,ch2],blength = samprate)$coh
        bandcoh.all[,ncol(y)*(ch1-1)+ch2] <- band.coh(x = x[,ch1],y = y[,ch2],samprate = samprate,
                                                      fband = rhythms)
      }
    }
    
    coh_comp[,r] <- apply(coh.all,MARGIN = 1,FUN = mean,na.rm = TRUE)
    bandcoh_comp[,r] <- apply(bandcoh.all,MARGIN = 1,FUN = mean,na.rm = TRUE)
    
    cancoh_comp[,r] <- cancoh(x = x,y = y,samprate = samprate,
                              blength = samprate)$cancoh[-1]
    
    kencoh_comp[,r] <- kencoh(x,y,samprate = samprate,fband = rhythms)
    
    cat(paste0("+ case = ",c.ind,", n = ",nsec_all[n.ind],", r = ",r),"\r")
    
  }
  
  saveRDS(est_comp,paste0("est_c",c.ind,"_n",n.ind,".rds"))
  saveRDS(pval_comp,paste0("pval_c",c.ind,"_n",n.ind,".rds"))
  saveRDS(pval.adj_comp,paste0("pvaladj_c",c.ind,"_n",n.ind,".rds"))
  saveRDS(coh_comp,paste0("coh_c",c.ind,"_n",n.ind,".rds"))
  saveRDS(cancoh_comp,paste0("cancoh_c",c.ind,"_n",n.ind,".rds"))
  saveRDS(bandcoh_comp,paste0("bandcoh_c",c.ind,"_n",n.ind,".rds"))
  saveRDS(kencoh_comp,paste0("kencoh_c",c.ind,"_n",n.ind,".rds"))
  
}



##########
# Case 2 #
##########

c.ind <- 2
rhythms <- c("delta","theta2","alpha","beta2","gamma2")
set.seed(c.ind) # seed for reproducibility

est_comp <- array(NA,dim = c(samprate/2,nrep))
coh_comp <- array(NA,dim = c(samprate/2,nrep))
cancoh_comp <- array(NA,dim = c(samprate/2,nrep))
bandcoh_comp <- array(NA,dim = c(5,nrep))
kencoh_comp <- array(NA,dim = c(5,nrep))
pval_comp <- array(NA,dim = c(samprate/2,nrep))
pval.adj_comp <- array(NA,dim = c(samprate/2,nrep))

for(n.ind in 1:length(nsec_all)){
  nsec <- nsec_all[n.ind]
  Tlength <- nsec*samprate
  
  for(r in 1:nrep){
    
    Z_alpha <- sim.fband(nsec = nsec,samprate = samprate,fband = rep(fband[3],4))
    
    x1 <- 0.5*Z_alpha[,1] + 0.5*rnorm(Tlength)
    x2 <- 0.5*Z_alpha[,2] + 0.5*rnorm(Tlength)
    y1 <- 0.5*Z_alpha[,1] + 0.5*rnorm(Tlength)
    y2 <- 0.5*Z_alpha[,3] + 0.5*rnorm(Tlength)
    x <- cbind(x1,x2)
    y <- cbind(y1,y2)
    
    samp <- nvc(x = x, y = y,blength = samprate,perm.inv = TRUE,
                b.scale = FALSE,nboot = nboot,ncores = ncores)
    est_comp[,r] <- samp$nvc
    pval_comp[,r] <- samp$pval
    pval.adj_comp[,r] <- samp$pval.adj
    
    coh.all <- array(NA,dim = c(samprate/2,ncol(x)*ncol(y)))
    bandcoh.all <- array(NA,dim = c(5,ncol(x)*ncol(y)))
    
    for(ch1 in 1:ncol(x)){
      for(ch2 in 1:ncol(y)){
        coh.all[,ncol(y)*(ch1-1)+ch2] <- coh(x = x[,ch1],y = y[,ch2],blength = samprate)$coh
        bandcoh.all[,ncol(y)*(ch1-1)+ch2] <- band.coh(x = x[,ch1],y = y[,ch2],samprate = samprate,
                                                      fband = rhythms)
      }
    }
    
    coh_comp[,r] <- apply(coh.all,MARGIN = 1,FUN = mean,na.rm = TRUE)
    bandcoh_comp[,r] <- apply(bandcoh.all,MARGIN = 1,FUN = mean,na.rm = TRUE)
    
    cancoh_comp[,r] <- cancoh(x = x,y = y,samprate = samprate,
                              blength = samprate)$cancoh[-1]
    
    kencoh_comp[,r] <- kencoh(x,y,samprate = samprate,fband = rhythms)
    
    cat(paste0("+ case = ",c.ind,", n = ",nsec_all[n.ind],", r = ",r),"\r")
    
  }
  
  saveRDS(est_comp,paste0("est_c",c.ind,"_n",n.ind,".rds"))
  saveRDS(pval_comp,paste0("pval_c",c.ind,"_n",n.ind,".rds"))
  saveRDS(pval.adj_comp,paste0("pvaladj_c",c.ind,"_n",n.ind,".rds"))
  saveRDS(coh_comp,paste0("coh_c",c.ind,"_n",n.ind,".rds"))
  saveRDS(cancoh_comp,paste0("cancoh_c",c.ind,"_n",n.ind,".rds"))
  saveRDS(bandcoh_comp,paste0("bandcoh_c",c.ind,"_n",n.ind,".rds"))
  saveRDS(kencoh_comp,paste0("kencoh_c",c.ind,"_n",n.ind,".rds"))
  
}



##########
# Case 3 #
##########

c.ind <- 3
rhythms <- c("delta","theta2","alpha","beta2","gamma2")
set.seed(c.ind) # seed for reproducibility

est_comp <- array(NA,dim = c(samprate/2,nrep))
coh_comp <- array(NA,dim = c(samprate/2,nrep))
cancoh_comp <- array(NA,dim = c(samprate/2,nrep))
bandcoh_comp <- array(NA,dim = c(5,nrep))
kencoh_comp <- array(NA,dim = c(5,nrep))
pval_comp <- array(NA,dim = c(samprate/2,nrep))
pval.adj_comp <- array(NA,dim = c(samprate/2,nrep))

for(n.ind in 1:length(nsec_all)){
  nsec <- nsec_all[n.ind]
  Tlength <- nsec*samprate
  
  for(r in 1:nrep){
    
    Z_alpha <- sim.fband(nsec = nsec,samprate = samprate,fband = rep(fband[3],4))
    
    x1 <- 0.5*Z_alpha[,1] + 0.5*rnorm(Tlength)
    x2 <- 0.5*Z_alpha[,2] + 0.5*rnorm(Tlength)
    y1 <- 0.5*Z_alpha[,3] + 0.5*rnorm(Tlength)
    y2 <- 0.5*Z_alpha[,4] + 0.5*rnorm(Tlength)
    x <- cbind(x1,x2)
    y <- cbind(y1,y2)
    
    samp <- nvc(x = x, y = y,blength = samprate,perm.inv = TRUE,
                b.scale = FALSE,nboot = nboot,ncores = ncores)
    est_comp[,r] <- samp$nvc
    pval_comp[,r] <- samp$pval
    pval.adj_comp[,r] <- samp$pval.adj
    
    coh.all <- array(NA,dim = c(samprate/2,ncol(x)*ncol(y)))
    bandcoh.all <- array(NA,dim = c(5,ncol(x)*ncol(y)))
    
    for(ch1 in 1:ncol(x)){
      for(ch2 in 1:ncol(y)){
        coh.all[,ncol(y)*(ch1-1)+ch2] <- coh(x = x[,ch1],y = y[,ch2],blength = samprate)$coh
        bandcoh.all[,ncol(y)*(ch1-1)+ch2] <- band.coh(x = x[,ch1],y = y[,ch2],samprate = samprate,
                                                      fband = rhythms)
      }
    }
    
    coh_comp[,r] <- apply(coh.all,MARGIN = 1,FUN = mean,na.rm = TRUE)
    bandcoh_comp[,r] <- apply(bandcoh.all,MARGIN = 1,FUN = mean,na.rm = TRUE)
    
    cancoh_comp[,r] <- cancoh(x = x,y = y,samprate = samprate,
                              blength = samprate)$cancoh[-1]
    
    kencoh_comp[,r] <- kencoh(x,y,samprate = samprate,fband = rhythms)
    
    cat(paste0("+ case = ",c.ind,", n = ",nsec_all[n.ind],", r = ",r),"\r")
    
  }
  
  saveRDS(est_comp,paste0("est_c",c.ind,"_n",n.ind,".rds"))
  saveRDS(pval_comp,paste0("pval_c",c.ind,"_n",n.ind,".rds"))
  saveRDS(pval.adj_comp,paste0("pvaladj_c",c.ind,"_n",n.ind,".rds"))
  saveRDS(coh_comp,paste0("coh_c",c.ind,"_n",n.ind,".rds"))
  saveRDS(cancoh_comp,paste0("cancoh_c",c.ind,"_n",n.ind,".rds"))
  saveRDS(bandcoh_comp,paste0("bandcoh_c",c.ind,"_n",n.ind,".rds"))
  saveRDS(kencoh_comp,paste0("kencoh_c",c.ind,"_n",n.ind,".rds"))
  
}

##########
# Case 4 #
##########

c.ind <- 4
rhythms <- c("delta2","theta","alpha2","beta2","gamma2")
set.seed(c.ind) # seed for reproducibility

est_comp <- array(NA,dim = c(samprate/2,nrep))
coh_comp <- array(NA,dim = c(samprate/2,nrep))
cancoh_comp <- array(NA,dim = c(samprate/2,nrep))
bandcoh_comp <- array(NA,dim = c(5,nrep))
kencoh_comp <- array(NA,dim = c(5,nrep))
pval_comp <- array(NA,dim = c(samprate/2,nrep))
pval.adj_comp <- array(NA,dim = c(samprate/2,nrep))

for(n.ind in 1:length(nsec_all)){
  nsec <- nsec_all[n.ind]
  Tlength <- nsec*samprate
  
  for(r in 1:nrep){
    
    Z_theta <- sim.fband(nsec = nsec,samprate = samprate,fband = rep(fband[2],6))
    Z_gamma <- sim.fband(nsec = nsec,samprate = samprate,fband = rep(fband[5],6))
    
    x1 <- 0.35*Z_theta[,1] + 0.35*Z_gamma[,1] + 0.3*rnorm(Tlength)
    x2 <- 0.35*Z_theta[,2] + 0.35*Z_gamma[,2] + 0.3*rnorm(Tlength)
    x3 <- 0.35*Z_theta[,3] + 0.35*Z_gamma[,3] + 0.3*rnorm(Tlength)
    y1 <- 0.35*Z_theta[,1] + 0.35*Z_gamma[,1] + 0.3*rnorm(Tlength)
    y2 <- 0.35*Z_theta[,2] + 0.35*Z_gamma[,2] + 0.3*rnorm(Tlength)
    y3 <- 0.35*Z_theta[,3] + 0.35*Z_gamma[,4] + 0.3*rnorm(Tlength)
    x <- cbind(x1,x2,x3)
    y <- cbind(y1,y2,y3)
    
    samp <- nvc(x = x, y = y,blength = samprate,perm.inv = TRUE,
                b.scale = FALSE,nboot = nboot,ncores = ncores)
    est_comp[,r] <- samp$nvc
    pval_comp[,r] <- samp$pval
    pval.adj_comp[,r] <- samp$pval.adj
    
    coh.all <- array(NA,dim = c(samprate/2,ncol(x)*ncol(y)))
    bandcoh.all <- array(NA,dim = c(5,ncol(x)*ncol(y)))
    
    for(ch1 in 1:ncol(x)){
      for(ch2 in 1:ncol(y)){
        coh.all[,ncol(y)*(ch1-1)+ch2] <- coh(x = x[,ch1],y = y[,ch2],blength = samprate)$coh
        bandcoh.all[,ncol(y)*(ch1-1)+ch2] <- band.coh(x = x[,ch1],y = y[,ch2],samprate = samprate,
                                                      fband = rhythms)
      }
    }
    
    coh_comp[,r] <- apply(coh.all,MARGIN = 1,FUN = mean,na.rm = TRUE)
    bandcoh_comp[,r] <- apply(bandcoh.all,MARGIN = 1,FUN = mean,na.rm = TRUE)
    
    cancoh_comp[,r] <- cancoh(x = x,y = y,samprate = samprate,
                              blength = samprate)$cancoh[-1]
    
    kencoh_comp[,r] <- kencoh(x,y,samprate = samprate,fband = rhythms)
    
    cat(paste0("+ case = ",c.ind,", n = ",nsec_all[n.ind],", r = ",r),"\r")
    
  }
  
  saveRDS(est_comp,paste0("est_c",c.ind,"_n",n.ind,".rds"))
  saveRDS(pval_comp,paste0("pval_c",c.ind,"_n",n.ind,".rds"))
  saveRDS(pval.adj_comp,paste0("pvaladj_c",c.ind,"_n",n.ind,".rds"))
  saveRDS(coh_comp,paste0("coh_c",c.ind,"_n",n.ind,".rds"))
  saveRDS(cancoh_comp,paste0("cancoh_c",c.ind,"_n",n.ind,".rds"))
  saveRDS(bandcoh_comp,paste0("bandcoh_c",c.ind,"_n",n.ind,".rds"))
  saveRDS(kencoh_comp,paste0("kencoh_c",c.ind,"_n",n.ind,".rds"))
  
}

##########
# Case 5 #
##########

c.ind <- 5
rhythms <- c("delta2","theta","alpha2","beta2","gamma2")
set.seed(c.ind) # seed for reproducibility

est_comp <- array(NA,dim = c(samprate/2,nrep))
coh_comp <- array(NA,dim = c(samprate/2,nrep))
cancoh_comp <- array(NA,dim = c(samprate/2,nrep))
bandcoh_comp <- array(NA,dim = c(5,nrep))
kencoh_comp <- array(NA,dim = c(5,nrep))
pval_comp <- array(NA,dim = c(samprate/2,nrep))
pval.adj_comp <- array(NA,dim = c(samprate/2,nrep))

for(n.ind in 1:length(nsec_all)){
  nsec <- nsec_all[n.ind]
  Tlength <- nsec*samprate
  
  for(r in 1:nrep){
    
    Z_theta <- sim.fband(nsec = nsec,samprate = samprate,fband = rep(fband[2],6))
    Z_gamma <- sim.fband(nsec = nsec,samprate = samprate,fband = rep(fband[5],6))
    
    x1 <- 0.35*Z_theta[,1] + 0.35*Z_gamma[,1] + 0.3*rnorm(Tlength)
    x2 <- 0.35*Z_theta[,2] + 0.35*Z_gamma[,2] + 0.3*rnorm(Tlength)
    x3 <- 0.35*Z_theta[,3] + 0.35*Z_gamma[,3] + 0.3*rnorm(Tlength)
    y1 <- 0.35*Z_theta[,1] + 0.35*Z_gamma[,1] + 0.3*rnorm(Tlength)
    y2 <- 0.35*Z_theta[,2] + 0.35*Z_gamma[,2] + 0.3*rnorm(Tlength)
    y3 <- 0.35*Z_theta[,4] + 0.35*Z_gamma[,3] + 0.3*rnorm(Tlength)
    x <- cbind(x1,x2,x3)
    y <- cbind(y1,y2,y3)
    
    samp <- nvc(x = x, y = y,blength = samprate,perm.inv = TRUE,
                b.scale = FALSE,nboot = nboot,ncores = ncores)
    est_comp[,r] <- samp$nvc
    pval_comp[,r] <- samp$pval
    pval.adj_comp[,r] <- samp$pval.adj
    
    coh.all <- array(NA,dim = c(samprate/2,ncol(x)*ncol(y)))
    bandcoh.all <- array(NA,dim = c(5,ncol(x)*ncol(y)))
    
    for(ch1 in 1:ncol(x)){
      for(ch2 in 1:ncol(y)){
        coh.all[,ncol(y)*(ch1-1)+ch2] <- coh(x = x[,ch1],y = y[,ch2],blength = samprate)$coh
        bandcoh.all[,ncol(y)*(ch1-1)+ch2] <- band.coh(x = x[,ch1],y = y[,ch2],samprate = samprate,
                                                      fband = rhythms)
      }
    }
    
    coh_comp[,r] <- apply(coh.all,MARGIN = 1,FUN = mean,na.rm = TRUE)
    bandcoh_comp[,r] <- apply(bandcoh.all,MARGIN = 1,FUN = mean,na.rm = TRUE)
    
    cancoh_comp[,r] <- cancoh(x = x,y = y,samprate = samprate,
                              blength = samprate)$cancoh[-1]
    
    kencoh_comp[,r] <- kencoh(x,y,samprate = samprate,fband = rhythms)
    
    cat(paste0("+ case = ",c.ind,", n = ",nsec_all[n.ind],", r = ",r),"\r")
    
  }
  
  saveRDS(est_comp,paste0("est_c",c.ind,"_n",n.ind,".rds"))
  saveRDS(pval_comp,paste0("pval_c",c.ind,"_n",n.ind,".rds"))
  saveRDS(pval.adj_comp,paste0("pvaladj_c",c.ind,"_n",n.ind,".rds"))
  saveRDS(coh_comp,paste0("coh_c",c.ind,"_n",n.ind,".rds"))
  saveRDS(cancoh_comp,paste0("cancoh_c",c.ind,"_n",n.ind,".rds"))
  saveRDS(bandcoh_comp,paste0("bandcoh_c",c.ind,"_n",n.ind,".rds"))
  saveRDS(kencoh_comp,paste0("kencoh_c",c.ind,"_n",n.ind,".rds"))
  
}


### Importing results
samprate <- 100
freq <- 1:50

est.mean <- array(NA,dim = c(length(freq),5,3))
est.lo <- array(NA,dim = c(length(freq),5,3))
est.hi <- array(NA,dim = c(length(freq),5,3))
pval.mean <- array(NA,dim = c(length(freq),5,3))
cancoh.mean <- array(NA,dim = c(length(freq),5,3))
cancoh.lo <- array(NA,dim = c(length(freq),5,3))
cancoh.hi <- array(NA,dim = c(length(freq),5,3))
est.comp <- array(NA,dim = c(length(freq),1000,5,3))
bandcoh.comp <- array(NA,dim = c(5,1000,5,3))
kencoh.comp <- array(NA,dim = c(5,1000,5,3))


for(c in 1:5){
  for(n in 1:3){
    est <- readRDS(paste("01 Simulations/NVC JASA/Frequency/est_c",c,"_n",n,".rds",sep = ""))
    est.comp[,,c,n] <- est
    
    pval <- readRDS(paste("01 Simulations/NVC JASA/Frequency/pval_c",c,"_n",n,".rds",sep = ""))
    
    cancoh <- readRDS(paste("01 Simulations/NVC JASA/Frequency/cancoh_c",c,"_n",n,".rds",sep = ""))
    
    bandcoh <- readRDS(paste("01 Simulations/NVC JASA/Frequency/bandcoh_c",c,"_n",n,".rds",sep = ""))
    bandcoh.comp[,,c,n] <- bandcoh
    
    kencoh <- readRDS(paste("01 Simulations/NVC JASA/Frequency/kencoh_c",c,"_n",n,".rds",sep = ""))
    kencoh.comp[,,c,n] <- kencoh
    
    est.mean[,c,n] <- apply(est,MARGIN = 1,FUN = mean,na.rm = TRUE)
    est.hi[,c,n] <- apply(est,MARGIN = 1,FUN = quantile,probs = 0.975,na.rm = TRUE)
    est.lo[,c,n] <- apply(est,MARGIN = 1,FUN = quantile,probs = 0.025,na.rm = TRUE)
    
    pval.mean[,c,n] <- apply(pval<0.05,MARGIN = 1,FUN = mean,na.rm = TRUE)
    
    cancoh.mean[,c,n] <- apply(cancoh,MARGIN = 1,FUN = mean,na.rm = TRUE)
    cancoh.hi[,c,n] <- apply(cancoh,MARGIN = 1,FUN = quantile,probs = 0.975,na.rm = TRUE)
    cancoh.lo[,c,n] <- apply(cancoh,MARGIN = 1,FUN = quantile,probs = 0.025,na.rm = TRUE)
    
  }
}


### Sampling rate and frequency bands
samprate <- 100
freq <- 1:50
ind_delta <- freq[freq > 0.5 & freq <= 4]
ind_theta <- freq[freq > 4 & freq < 8]
ind_alpha <- freq[freq >= 8 & freq <= 12]
ind_delta2 <- freq[freq > 0.5 & freq < 4]
ind_theta2 <- freq[freq >= 4 & freq <= 8]
ind_alpha2 <- freq[freq > 8 & freq <= 12]
ind_beta <- freq[freq >= 21 & freq <= 25]
ind_gamma <- freq[freq >= 36 & freq <= 40]


### Organizing results
est.sd <- apply(est.comp,MARGIN = c(1,3,4),FUN = sd,na.rm = TRUE)

c.ind <- 1:3
power1 <- apply(pval.mean[ind_alpha,c.ind,],MARGIN = c(2,3),FUN = mean,na.rm = TRUE)
size1 <- apply(pval.mean[setdiff(freq,ind_alpha),c.ind,],MARGIN = c(2,3),FUN = mean,na.rm = TRUE)
std.p1 <- apply(est.sd[ind_alpha,c.ind,],MARGIN = c(2,3),FUN = mean,na.rm = TRUE)
std.s1 <- apply(est.sd[setdiff(freq,ind_alpha),c.ind,],MARGIN = c(2,3),FUN = mean,na.rm = TRUE)

c.ind <- 4:5
power2 <- apply(pval.mean[ind_theta2,c.ind,],MARGIN = c(2,3),FUN = mean,na.rm = TRUE)
power3 <- apply(pval.mean[ind_gamma,c.ind,],MARGIN = c(2,3),FUN = mean,na.rm = TRUE)
size2 <- apply(pval.mean[setdiff(freq,c(ind_theta2,ind_gamma)),c.ind,],MARGIN = c(2,3),FUN = mean,na.rm = TRUE)
std.p2 <- apply(est.sd[ind_theta2,c.ind,],MARGIN = c(2,3),FUN = mean,na.rm = TRUE)
std.p3 <- apply(est.sd[ind_gamma,c.ind,],MARGIN = c(2,3),FUN = mean,na.rm = TRUE)
std.s2 <- apply(est.sd[setdiff(freq,c(ind_theta2,ind_gamma)),c.ind,],MARGIN = c(2,3),FUN = mean,na.rm = TRUE)

power_and_size <- rbind(power1[1,],size1[1,],power1[2,],size1[2,],power1[3,],size1[3,],
                        power2[1,],power3[1,],size2[1,],power2[2,],power3[2,],size2[2,])
stdev <- rbind(std.p1[1,],std.s1[1,],std.p1[2,],std.s1[2,],std.p1[3,],std.s1[3,],
               std.p2[1,],std.p3[1,],std.s2[1,],std.p2[2,],std.p3[2,],std.s2[2,])
cases <- c(rep(paste0("Case ",1:3),each = 2),rep(paste0("Case ",4:5),each = 3))
fbands <- c(rep(c("alpha","not alpha"),3),rep(c("theta","gamma","not theta-gamma"),2))

### Table of size and power per frequency band
cbind(cases,fbands,round(stdev,4),round(power_and_size,4))


### Visualizing results
png("sim_linear_c123.png",width = 1500, height = 1500)
par(mfcol = c(3,2),mar = c(8,5,2,2),oma = c(10,8,8,0),xpd = FALSE)

for(n in 1:3){
  plot(NA,NA,type = "l",ylim = c(-0.2,0.6),xlim = c(0,50),xlab = "",ylab = "",cex.axis = 2.75)
  box(lwd = 3)
  mtext(TeX("Frequency $\\omega$ (in Hz)"),side = 1,cex = 2.5,line = 5)
  abline(h = 0,col = "black",lty = 2,lwd = 2)
  
  polygon(x = c(freq, rev(freq)),
          y = c(est.hi[,1,n], rev(est.lo[,1,n])),
          col = adjustcolor("darkviolet",alpha.f = 0.15), border = NA)
  lines(freq,est.mean[,1,n],col = "darkviolet",lwd = 3)
  lines(freq,est.lo[,1,n],lty = 3,col = "darkviolet",lwd = 3)
  lines(freq,est.hi[,1,n],lty = 3,col = "darkviolet",lwd = 3)
  
  polygon(x = c(freq, rev(freq)),
          y = c(est.hi[,2,n], rev(est.lo[,2,n])),
          col = adjustcolor("darkorange",alpha.f = 0.15), border = NA)
  lines(freq,est.mean[,2,n],col = "darkorange",lwd = 3)
  lines(freq,est.lo[,2,n],lty = 3,col = "darkorange",lwd = 3)
  lines(freq,est.hi[,2,n],lty = 3,col = "darkorange",lwd = 3)
  
  polygon(x = c(freq, rev(freq)),
          y = c(est.hi[,3,n], rev(est.lo[,3,n])),
          col = adjustcolor("limegreen",alpha.f = 0.15), border = NA)
  lines(freq,est.mean[,3,n],col = "limegreen",lwd = 3)
  lines(freq,est.lo[,3,n],lty = 3,col = "limegreen",lwd = 3)
  lines(freq,est.hi[,3,n],lty = 3,col = "limegreen",lwd = 3)
  
}

for(n in 1:3){
  plot(NA,NA,type = "l",ylim = c(0,1),xlim = c(0,50),xlab = "",ylab = "",cex.axis = 2.75)
  box(lwd = 3)
  mtext(TeX("Frequency $\\omega$ (in Hz)"),side = 1,cex = 2.5,line = 5)
  abline(h = 0.05,col = "black",lty = 2,lwd = 2)
  
  lines(freq,pval.mean[,1,n],col = "darkviolet",lwd = 5)
  
  lines(freq,pval.mean[,2,n],col = "darkorange",lwd = 5)
  
  lines(freq,pval.mean[,3,n],col = "limegreen",lwd = 5)
  
  
}

mtext(TeX("Estimates",bold = TRUE,italic = TRUE),outer = TRUE,
      side = 3,cex = 3,line = 2,adj = 0.225)
mtext(TeX("Prob(Rej. H0)",bold = TRUE,italic = TRUE),outer = TRUE,
      side = 3,cex = 3,line = 2,adj = 0.79)

mtext(TeX("n = 50",bold = TRUE,italic = TRUE),outer = TRUE,
      side = 2,cex = 3,line = 2,adj = 0.875)
mtext(TeX("n = 100",bold = TRUE,italic = TRUE),outer = TRUE,
      side = 2,cex = 3,line = 2,adj = 0.5)
mtext(TeX("n = 200",bold = TRUE,italic = TRUE),outer = TRUE,
      side = 2,cex = 3,line = 2,adj = 0.15)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend(x = -0.4,y = -0.95,legend = c("Case 1  ","Case 2  ","Case 3"),
       col = c("darkviolet","darkorange", "limegreen"),
       lwd = 5, xpd = TRUE, horiz = TRUE, cex = 4,
       pt.bg = c(adjustcolor("darkviolet",alpha.f = 0.15),
                 adjustcolor("darkorange",alpha.f = 0.15),
                 adjustcolor("limegreen",alpha.f = 0.15)),
       pch = 22,pt.cex = 10,bty = "n")

dev.off()



png("sim_linear_c45.png",width = 1500, height = 1500)
par(mfcol = c(3,2),mar = c(8,5,2,2),oma = c(10,8,8,0),xpd = FALSE)

for(n in 1:3){
  plot(NA,NA,type = "l",ylim = c(-0.2,0.6),xlim = c(0,50),xlab = "",ylab = "",cex.axis = 2.75)
  box(lwd = 3)
  mtext(TeX("Frequency $\\omega$ (in Hz)"),side = 1,cex = 2.5,line = 5)
  abline(h = 0,col = "black",lty = 2,lwd = 2)
  
  polygon(x = c(freq, rev(freq)),
          y = c(est.hi[,4,n], rev(est.lo[,4,n])),
          col = adjustcolor("darkblue",alpha.f = 0.15), border = NA)
  lines(freq,est.mean[,4,n],col = "darkblue",lwd = 3)
  lines(freq,est.lo[,4,n],lty = 3,col = "darkblue",lwd = 3)
  lines(freq,est.hi[,4,n],lty = 3,col = "darkblue",lwd = 3)
  
  polygon(x = c(freq, rev(freq)),
          y = c(est.hi[,5,n], rev(est.lo[,5,n])),
          col = adjustcolor("red",alpha.f = 0.15), border = NA)
  lines(freq,est.mean[,5,n],col = "red",lwd = 3)
  lines(freq,est.lo[,5,n],lty = 3,col = "red",lwd = 3)
  lines(freq,est.hi[,5,n],lty = 3,col = "red",lwd = 3)
  
}

for(n in 1:3){
  plot(NA,NA,type = "l",ylim = c(0,1),xlim = c(0,50),xlab = "",ylab = "",cex.axis = 2.75)
  box(lwd = 3)
  mtext(TeX("Frequency $\\omega$ (in Hz)"),side = 1,cex = 2.5,line = 5)
  abline(h = 0.05,col = "black",lty = 2,lwd = 2)
  
  lines(freq,pval.mean[,4,n],col = "darkblue",lwd = 5)
  
  lines(freq,pval.mean[,5,n],col = "red",lwd = 5)
  
}

mtext(TeX("Estimates",bold = TRUE,italic = TRUE),outer = TRUE,
      side = 3,cex = 3,line = 2,adj = 0.225)
mtext(TeX("Prob(Rej. H0)",bold = TRUE,italic = TRUE),outer = TRUE,
      side = 3,cex = 3,line = 2,adj = 0.79)

mtext(TeX("n = 50",bold = TRUE,italic = TRUE),outer = TRUE,
      side = 2,cex = 3,line = 2,adj = 0.875)
mtext(TeX("n = 100",bold = TRUE,italic = TRUE),outer = TRUE,
      side = 2,cex = 3,line = 2,adj = 0.5)
mtext(TeX("n = 200",bold = TRUE,italic = TRUE),outer = TRUE,
      side = 2,cex = 3,line = 2,adj = 0.15)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend(x = -0.2,y = -0.95,legend = c("Case 4  ","Case 5"),
       col = c("darkblue", "red"),
       lwd = 5, xpd = TRUE, horiz = TRUE, cex = 4,
       pt.bg = c(adjustcolor("darkblue",alpha.f = 0.15),
                 adjustcolor("red",alpha.f = 0.15)),
       pch = 22,pt.cex = 10,bty = "n")

dev.off()


### Organizing results
# Case 1-3
c.ind <- 1:3
est.fband <- array(NA,dim = c(5,1000,5,3))
est.fband[1,,c.ind,] <- apply(est.comp[ind_delta,,c.ind,],MARGIN = c(2,3,4),FUN = mean,na.rm = TRUE)
est.fband[2,,c.ind,] <- apply(est.comp[ind_theta,,c.ind,],MARGIN = c(2,3,4),FUN = mean,na.rm = TRUE)
est.fband[3,,c.ind,] <- apply(est.comp[ind_alpha,,c.ind,],MARGIN = c(2,3,4),FUN = mean,na.rm = TRUE)
est.fband[4,,c.ind,] <- apply(est.comp[ind_beta,,c.ind,],MARGIN = c(2,3,4),FUN = mean,na.rm = TRUE)
est.fband[5,,c.ind,] <- apply(est.comp[ind_gamma,,c.ind,],MARGIN = c(2,3,4),FUN = mean,na.rm = TRUE)

# Case 4-5
c.ind <- 4:5
est.fband[1,,c.ind,] <- apply(est.comp[ind_delta2,,c.ind,],MARGIN = c(2,3,4),FUN = mean,na.rm = TRUE)
est.fband[2,,c.ind,] <- apply(est.comp[ind_theta2,,c.ind,],MARGIN = c(2,3,4),FUN = mean,na.rm = TRUE)
est.fband[3,,c.ind,] <- apply(est.comp[ind_alpha2,,c.ind,],MARGIN = c(2,3,4),FUN = mean,na.rm = TRUE)
est.fband[4,,c.ind,] <- apply(est.comp[ind_beta,,c.ind,],MARGIN = c(2,3,4),FUN = mean,na.rm = TRUE)
est.fband[5,,c.ind,] <- apply(est.comp[ind_gamma,,c.ind,],MARGIN = c(2,3,4),FUN = mean,na.rm = TRUE)


### Visualizing results
png("sim_linear_c123_with_cancoh.png",width = 2250, height = 1500)
par(mfcol = c(3,3),mar = c(8,5,2,2),oma = c(10,8,8,0),xpd = FALSE)

color_all <- c("darkviolet","darkorange","limegreen","darkblue","red")

for(c in 1:3){
  for(n in 1:3){
    plot(NA,NA,type = "l",ylim = c(-0.2,1),xlim = c(0,50),xlab = "",ylab = "",cex.axis = 2.75)
    box(lwd = 4)
    mtext(TeX("Frequency $\\omega$ (in Hz)"),side = 1,cex = 2.5,line = 5)
    abline(h = 0,col = "black",lty = 2,lwd = 2)
    
    polygon(x = c(freq, rev(freq)),
            y = c(est.hi[,c,n], rev(est.lo[,c,n])),
            col = adjustcolor(color_all[c],alpha.f = 0.15), border = NA)
    lines(freq,est.mean[,c,n],col = color_all[c],lwd = 4)
    lines(freq,est.lo[,c,n],lty = 3,col = color_all[c],lwd = 4)
    lines(freq,est.hi[,c,n],lty = 3,col = color_all[c],lwd = 4)
    
    polygon(x = c(freq, rev(freq)),
            y = c(cancoh.hi[,c,n], rev(cancoh.lo[,c,n])),
            col = adjustcolor("black",alpha.f = 0.15), border = NA)
    lines(freq,cancoh.mean[,c,n],col = "black",lwd = 4)
    lines(freq,cancoh.lo[,c,n],lty = 3,col = "black",lwd = 4)
    lines(freq,cancoh.hi[,c,n],lty = 3,col = "black",lwd = 4)
    
    
  }
}


mtext(TeX("Case 1",bold = TRUE,italic = TRUE),outer = TRUE,
      side = 3,cex = 3,line = 1,adj = 0.15)
mtext(TeX("Case 2",bold = TRUE,italic = TRUE),outer = TRUE,
      side = 3,cex = 3,line = 1,adj = 0.5)
mtext(TeX("Case 3",bold = TRUE,italic = TRUE),outer = TRUE,
      side = 3,cex = 3,line = 1,adj = 0.85)


mtext(TeX("n = 50",bold = TRUE,italic = TRUE),outer = TRUE,
      side = 2,cex = 3,line = 2,adj = 0.875)
mtext(TeX("n = 100",bold = TRUE,italic = TRUE),outer = TRUE,
      side = 2,cex = 3,line = 2,adj = 0.5)
mtext(TeX("n = 200",bold = TRUE,italic = TRUE),outer = TRUE,
      side = 2,cex = 3,line = 2,adj = 0.15)


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend(x = -0.855,y = -0.95,legend = c("NVC  ", "CanCoh"),
       col = c("darkviolet","black"),
       lwd = 4, xpd = TRUE, horiz = TRUE, cex = 4,
       pt.bg = c(adjustcolor("darkviolet",alpha.f = 0.15),
                 adjustcolor("black",alpha.f = 0.15)),
       pch = 22,pt.cex = 10,bty = "n")
legend(x = -0.155,y = -0.95,legend = c("NVC  ","CanCoh"),
       col = c("darkorange", "black"),
       lwd = 4, xpd = TRUE, horiz = TRUE, cex = 4,
       pt.bg = c(adjustcolor("darkorange",alpha.f = 0.15),
                 adjustcolor("black",alpha.f = 0.15)),
       pch = 22,pt.cex = 10,bty = "n")
legend(x = 0.54,y = -0.95,legend = c("NVC  ","CanCoh"),
       col = c("limegreen", "black"),
       lwd = 4, xpd = TRUE, horiz = TRUE, cex = 4,
       pt.bg = c(adjustcolor("limegreen",alpha.f = 0.15),
                 adjustcolor("black",alpha.f = 0.15)),
       pch = 22,pt.cex = 10,bty = "n")

dev.off()




png("sim_linear_c45_with_cancoh.png",width = 1500, height = 1500)
par(mfcol = c(3,2),mar = c(8,5,2,2),oma = c(10,8,8,0),xpd = FALSE)

color_all <- c("darkviolet","darkorange","limegreen","darkblue","red")

for(c in 4:5){
  for(n in 1:3){
    plot(NA,NA,type = "l",ylim = c(-0.2,1),xlim = c(0,50),xlab = "",ylab = "",cex.axis = 2.75)
    box(lwd = 4)
    mtext(TeX("Frequency $\\omega$ (in Hz)"),side = 1,cex = 2.5,line = 5)
    abline(h = 0,col = "black",lty = 2,lwd = 2)
    
    polygon(x = c(freq, rev(freq)),
            y = c(est.hi[,c,n], rev(est.lo[,c,n])),
            col = adjustcolor(color_all[c],alpha.f = 0.15), border = NA)
    lines(freq,est.mean[,c,n],col = color_all[c],lwd = 4)
    lines(freq,est.lo[,c,n],lty = 3,col = color_all[c],lwd = 4)
    lines(freq,est.hi[,c,n],lty = 3,col = color_all[c],lwd = 4)
    
    polygon(x = c(freq, rev(freq)),
            y = c(cancoh.hi[,c,n], rev(cancoh.lo[,c,n])),
            col = adjustcolor("black",alpha.f = 0.15), border = NA)
    lines(freq,cancoh.mean[,c,n],col = "black",lwd = 4)
    lines(freq,cancoh.lo[,c,n],lty = 3,col = "black",lwd = 4)
    lines(freq,cancoh.hi[,c,n],lty = 3,col = "black",lwd = 4)
    
    
  }
}


mtext(TeX("Case 4",bold = TRUE,italic = TRUE),outer = TRUE,
      side = 3,cex = 3,line = 1,adj = 0.225)
mtext(TeX("Case 5",bold = TRUE,italic = TRUE),outer = TRUE,
      side = 3,cex = 3,line = 1,adj = 0.775)


mtext(TeX("n = 50",bold = TRUE,italic = TRUE),outer = TRUE,
      side = 2,cex = 3,line = 2,adj = 0.875)
mtext(TeX("n = 100",bold = TRUE,italic = TRUE),outer = TRUE,
      side = 2,cex = 3,line = 2,adj = 0.5)
mtext(TeX("n = 200",bold = TRUE,italic = TRUE),outer = TRUE,
      side = 2,cex = 3,line = 2,adj = 0.15)


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend(x = -0.75,y = -0.95,legend = c("NVC  ", "CanCoh"),
       col = c("darkblue","black"),
       lwd = 4, xpd = TRUE, horiz = TRUE, cex = 4,
       pt.bg = c(adjustcolor("darkblue",alpha.f = 0.15),
                 adjustcolor("black",alpha.f = 0.15)),
       pch = 22,pt.cex = 10,bty = "n")
legend(x = 0.275,y = -0.95,legend = c("NVC  ","CanCoh"),
       col = c("red", "black"),
       lwd = 4, xpd = TRUE, horiz = TRUE, cex = 4,
       pt.bg = c(adjustcolor("red",alpha.f = 0.15),
                 adjustcolor("black",alpha.f = 0.15)),
       pch = 22,pt.cex = 10,bty = "n")

dev.off()



### Sampling rate and frequency bands
samprate <- 100
freq <- 1:50
ind_delta <- freq[freq > 0.5 & freq <= 4]
ind_theta <- freq[freq > 4 & freq < 8]
ind_alpha <- freq[freq >= 8 & freq <= 12]
ind_delta2 <- freq[freq > 0.5 & freq < 4]
ind_theta2 <- freq[freq >= 4 & freq <= 8]
ind_alpha2 <- freq[freq > 8 & freq <= 12]
ind_beta <- freq[freq >= 21 & freq <= 25]
ind_gamma <- freq[freq >= 36 & freq <= 40]

### Organizing results
# Case 1-3
c.ind <- 1:3
est.fband <- array(NA,dim = c(5,1000,5,3))
est.fband[1,,c.ind,] <- apply(est.comp[ind_delta,,c.ind,],MARGIN = c(2,3,4),FUN = mean,na.rm = TRUE)
est.fband[2,,c.ind,] <- apply(est.comp[ind_theta,,c.ind,],MARGIN = c(2,3,4),FUN = mean,na.rm = TRUE)
est.fband[3,,c.ind,] <- apply(est.comp[ind_alpha,,c.ind,],MARGIN = c(2,3,4),FUN = mean,na.rm = TRUE)
est.fband[4,,c.ind,] <- apply(est.comp[ind_beta,,c.ind,],MARGIN = c(2,3,4),FUN = mean,na.rm = TRUE)
est.fband[5,,c.ind,] <- apply(est.comp[ind_gamma,,c.ind,],MARGIN = c(2,3,4),FUN = mean,na.rm = TRUE)

# Case 4-5
c.ind <- 4:5
est.fband[1,,c.ind,] <- apply(est.comp[ind_delta2,,c.ind,],MARGIN = c(2,3,4),FUN = mean,na.rm = TRUE)
est.fband[2,,c.ind,] <- apply(est.comp[ind_theta2,,c.ind,],MARGIN = c(2,3,4),FUN = mean,na.rm = TRUE)
est.fband[3,,c.ind,] <- apply(est.comp[ind_alpha2,,c.ind,],MARGIN = c(2,3,4),FUN = mean,na.rm = TRUE)
est.fband[4,,c.ind,] <- apply(est.comp[ind_beta,,c.ind,],MARGIN = c(2,3,4),FUN = mean,na.rm = TRUE)
est.fband[5,,c.ind,] <- apply(est.comp[ind_gamma,,c.ind,],MARGIN = c(2,3,4),FUN = mean,na.rm = TRUE)


### Visualizing results
png("sim_linear_c123_with_pbc_kencoh.png",width = 1200,height = 1500)
par(mfrow = c(3,1),mar = c(10,7,1,1),oma = c(1,5,1,1))
boxplot(cbind(est.fband[3,,1,],NA,est.fband[3,,2,],NA,est.fband[3,,3,]),ylim = c(-0.1,1),
        xaxt = "n",yaxt = "n",cex = 1,lwd = 3.5)
axis(side = 2, at = 0.2*(-1:5),labels = 0.2*(-1:5),lwd = 5,cex.axis = 3.5,font = 2)
axis(side = 1, at = c(1:3,5:7,9:11),labels = rep(c(50,100,200),3),lwd = 5,cex.axis = 3.5,padj = 0.5,font = 2)
box(lwd = 5)
abline(h = 0,lty = 2,lwd = 3)

axis(side = 3, at = c(1:3),labels = c(NA,NA),lwd = 5,line = -45.5)
axis(side = 3, at = 2,labels = "Case 1",lwd = 5,cex.axis = 4,font = 2,line = -50,tick = FALSE)
axis(side = 3, at = c(5:7),labels = c(NA,NA),lwd = 5,line = -45.5)
axis(side = 3, at = 6,labels = "Case 2",lwd = 5,cex.axis = 4,font = 2,line = -50,tick = FALSE)
axis(side = 3, at = c(9:11),labels = c(NA,NA),lwd = 5,line = -45.5)
axis(side = 3, at = 10,labels = "Case 3",lwd = 5,cex.axis = 4,font = 2,line = -50,tick = FALSE)

mtext("n = ",side = 1,cex = 3,font = 2,line = 2,adj = -0.01)
mtext("NVC",font = 2,side = 2,cex = 5,line = 6)

boxplot(cbind(bandcoh.comp[3,,1,],NA,bandcoh.comp[3,,2,],NA,bandcoh.comp[3,,3,]),ylim = c(-0.1,1),
        xaxt = "n",yaxt = "n",cex = 1,lwd = 3.5)
axis(side = 2, at = 0.2*(-1:5),labels = 0.2*(-1:5),lwd = 5,cex.axis = 3.5,font = 2)
axis(side = 1, at = c(1:3,5:7,9:11),labels = rep(c(50,100,200),3),lwd = 5,cex.axis = 3.5,padj = 0.5,font = 2)
box(lwd = 5)
abline(h = 0,lty = 2,lwd = 3)

axis(side = 3, at = c(1:3),labels = c(NA,NA),lwd = 5,line = -45.5)
axis(side = 3, at = 2,labels = "Case 1",lwd = 5,cex.axis = 4,font = 2,line = -50,tick = FALSE)
axis(side = 3, at = c(5:7),labels = c(NA,NA),lwd = 5,line = -45.5)
axis(side = 3, at = 6,labels = "Case 2",lwd = 5,cex.axis = 4,font = 2,line = -50,tick = FALSE)
axis(side = 3, at = c(9:11),labels = c(NA,NA),lwd = 5,line = -45.5)
axis(side = 3, at = 10,labels = "Case 3",lwd = 5,cex.axis = 4,font = 2,line = -50,tick = FALSE)

mtext("n = ",side = 1,cex = 3,font = 2,line = 2,adj = -0.01)
mtext("PBC",font = 2,side = 2,cex = 5,line = 6)

boxplot(cbind(kencoh.comp[3,,1,],NA,kencoh.comp[3,,2,],NA,kencoh.comp[3,,3,]),ylim = c(-0.1,1),
        xaxt = "n",yaxt = "n",cex = 1,lwd = 3.5)
axis(side = 2, at = 0.2*(-1:5),labels = 0.2*(-1:5),lwd = 5,cex.axis = 3.5,font = 2)
axis(side = 1, at = c(1:3,5:7,9:11),labels = rep(c(50,100,200),3),lwd = 5,cex.axis = 3.5,padj = 0.5,font = 2)
box(lwd = 5)
abline(h = 0,lty = 2,lwd = 3)

axis(side = 3, at = c(1:3),labels = c(NA,NA),lwd = 5,line = -45.5)
axis(side = 3, at = 2,labels = "Case 1",lwd = 5,cex.axis = 4,font = 2,line = -50,tick = FALSE)
axis(side = 3, at = c(5:7),labels = c(NA,NA),lwd = 5,line = -45.5)
axis(side = 3, at = 6,labels = "Case 2",lwd = 5,cex.axis = 4,font = 2,line = -50,tick = FALSE)
axis(side = 3, at = c(9:11),labels = c(NA,NA),lwd = 5,line = -45.5)
axis(side = 3, at = 10,labels = "Case 3",lwd = 5,cex.axis = 4,font = 2,line = -50,tick = FALSE)

mtext("n = ",side = 1,cex = 3,font = 2,line = 2,adj = -0.01)
mtext("KenCoh",font = 2,side = 2,cex = 5,line = 6)

dev.off()

png("sim_linear_c45_with_pbc_kencoh.png",width = 1200,height = 1500)
par(mfrow = c(3,1),mar = c(12,7,1,1),oma = c(1,5,1,1))
boxplot(cbind(est.fband[2,,4,],est.fband[5,,4,],NA,est.fband[2,,5,],est.fband[5,,5,]),ylim = c(-0.1,1),
        xaxt = "n",yaxt = "n",cex = 1,lwd = 3.5)
axis(side = 2, at = 0.2*(-1:5),labels = 0.2*(-1:5),lwd = 5,cex.axis = 3,font = 2)
axis(side = 1, at = c(1:6,8:13),labels = rep(c(50,100,200,50,100,200),2),lwd = 5,cex.axis = 3,padj = 0.5,font = 2)
box(lwd = 5)
abline(h = 0,lty = 2,lwd = 3)

axis(side = 3, at = c(1:3),labels = c(NA,NA),lwd = 5,line = -42.5)
axis(side = 3, at = c(4:6),labels = c(NA,NA),lwd = 5,line = -42.5)
axis(side = 3, at = 2,labels = "Theta",lwd = 5,cex.axis = 3.5,font = 2,line = -46.5,tick = FALSE)
axis(side = 3, at = 5,labels = "Gamma",lwd = 5,cex.axis = 3.5,font = 2,line = -46.5,tick = FALSE)
axis(side = 3, at = c(2,5),labels = c(NA,NA),lwd = 5,line = -46.5)
axis(side = 3, at = 3.5,labels = "Case 4",lwd = 5,cex.axis = 4,font = 2,line = -50.5,tick = FALSE)

axis(side = 3, at = c(8:10),labels = c(NA,NA),lwd = 5,line = -42.5)
axis(side = 3, at = c(11:13),labels = c(NA,NA),lwd = 5,line = -42.5)
axis(side = 3, at = 9,labels = "Theta",lwd = 5,cex.axis = 3.5,font = 2,line = -46.5,tick = FALSE)
axis(side = 3, at = 12,labels = "Gamma",lwd = 5,cex.axis = 3.5,font = 2,line = -46.5,tick = FALSE)
axis(side = 3, at = c(9,12),labels = c(NA,NA),lwd = 5,line = -46.5)
axis(side = 3, at = 10.5,labels = "Case 5",lwd = 5,cex.axis = 4,font = 2,line = -50.5,tick = FALSE)

mtext("n = ",side = 1,cex = 2.5,font = 2,line = 2,adj = 0)
mtext("NVC",font = 2,side = 2,cex = 5,line = 6)



boxplot(cbind(bandcoh.comp[2,,4,],bandcoh.comp[5,,4,],NA,bandcoh.comp[2,,5,],bandcoh.comp[5,,5,]),ylim = c(-0.1,1),
        xaxt = "n",yaxt = "n",cex = 1,lwd = 3.5)
axis(side = 2, at = 0.2*(-1:5),labels = 0.2*(-1:5),lwd = 5,cex.axis = 3,font = 2)
axis(side = 1, at = c(1:6,8:13),labels = rep(c(50,100,200,50,100,200),2),lwd = 5,cex.axis = 3,padj = 0.5,font = 2)
box(lwd = 5)
abline(h = 0,lty = 2,lwd = 3)

axis(side = 3, at = c(1:3),labels = c(NA,NA),lwd = 5,line = -42.5)
axis(side = 3, at = c(4:6),labels = c(NA,NA),lwd = 5,line = -42.5)
axis(side = 3, at = 2,labels = "Theta",lwd = 5,cex.axis = 3.5,font = 2,line = -46.5,tick = FALSE)
axis(side = 3, at = 5,labels = "Gamma",lwd = 5,cex.axis = 3.5,font = 2,line = -46.5,tick = FALSE)
axis(side = 3, at = c(2,5),labels = c(NA,NA),lwd = 5,line = -46.5)
axis(side = 3, at = 3.5,labels = "Case 4",lwd = 5,cex.axis = 4,font = 2,line = -50.5,tick = FALSE)

axis(side = 3, at = c(8:10),labels = c(NA,NA),lwd = 5,line = -42.5)
axis(side = 3, at = c(11:13),labels = c(NA,NA),lwd = 5,line = -42.5)
axis(side = 3, at = 9,labels = "Theta",lwd = 5,cex.axis = 3.5,font = 2,line = -46.5,tick = FALSE)
axis(side = 3, at = 12,labels = "Gamma",lwd = 5,cex.axis = 3.5,font = 2,line = -46.5,tick = FALSE)
axis(side = 3, at = c(9,12),labels = c(NA,NA),lwd = 5,line = -46.5)
axis(side = 3, at = 10.5,labels = "Case 5",lwd = 5,cex.axis = 4,font = 2,line = -50.5,tick = FALSE)

mtext("n = ",side = 1,cex = 2.5,font = 2,line = 2,adj = 0)
mtext("PBC",font = 2,side = 2,cex = 5,line = 6)



boxplot(cbind(kencoh.comp[2,,4,],kencoh.comp[5,,4,],NA,kencoh.comp[2,,5,],kencoh.comp[5,,5,]),ylim = c(-0.1,1),
        xaxt = "n",yaxt = "n",cex = 1,lwd = 3.5)
axis(side = 2, at = 0.2*(-1:5),labels = 0.2*(-1:5),lwd = 5,cex.axis = 3,font = 2)
axis(side = 1, at = c(1:6,8:13),labels = rep(c(50,100,200,50,100,200),2),lwd = 5,cex.axis = 3,padj = 0.5,font = 2)
box(lwd = 5)
abline(h = 0,lty = 2,lwd = 3)

axis(side = 3, at = c(1:3),labels = c(NA,NA),lwd = 5,line = -42.5)
axis(side = 3, at = c(4:6),labels = c(NA,NA),lwd = 5,line = -42.5)
axis(side = 3, at = 2,labels = "Theta",lwd = 5,cex.axis = 3.5,font = 2,line = -46.5,tick = FALSE)
axis(side = 3, at = 5,labels = "Gamma",lwd = 5,cex.axis = 3.5,font = 2,line = -46.5,tick = FALSE)
axis(side = 3, at = c(2,5),labels = c(NA,NA),lwd = 5,line = -46.5)
axis(side = 3, at = 3.5,labels = "Case 4",lwd = 5,cex.axis = 4,font = 2,line = -50.5,tick = FALSE)

axis(side = 3, at = c(8:10),labels = c(NA,NA),lwd = 5,line = -42.5)
axis(side = 3, at = c(11:13),labels = c(NA,NA),lwd = 5,line = -42.5)
axis(side = 3, at = 9,labels = "Theta",lwd = 5,cex.axis = 3.5,font = 2,line = -46.5,tick = FALSE)
axis(side = 3, at = 12,labels = "Gamma",lwd = 5,cex.axis = 3.5,font = 2,line = -46.5,tick = FALSE)
axis(side = 3, at = c(9,12),labels = c(NA,NA),lwd = 5,line = -46.5)
axis(side = 3, at = 10.5,labels = "Case 5",lwd = 5,cex.axis = 4,font = 2,line = -50.5,tick = FALSE)

mtext("n = ",side = 1,cex = 2.5,font = 2,line = 2,adj = 0)
mtext("KenCoh",font = 2,side = 2,cex = 5,line = 6)

dev.off()



#################################################
##### Simulation for general random vectors #####
#################################################

nrep <- 1000
nboot <- 500
ncores <- 8

####################
# Independent Case #
####################

set.seed(1) # seed for reproducibility
family_set_all <- list("indep","gaussian","gumbel","clayton","frank","joe",
                       c("indep","gaussian","gumbel","clayton","frank","joe"))
margin_all <- c("unif","norm")
n_all <- c(50,100,200)
q_all <- c(2,3,4)

pval_comp <- array(NA,dim = c(nrep,length(n_all),length(q_all),
                              length(family_set_all),length(margin_all)))
est_comp <- array(NA,dim = c(nrep,length(n_all),length(q_all),
                             length(family_set_all),length(margin_all)))

for(m.ind in 1:length(margin_all)){
  margin <- margin_all[m.ind]
  
  for(n.ind in 1:length(n_all)){
    n <- n_all[n.ind]
    
    for(q.ind in 1:length(q_all)){
      q <- q_all[q.ind]
      
      for(f.ind in 1:length(family_set_all)){
        family_set <- family_set_all[[f.ind]]
        
        for(b in 1:nrep){
          pcs <- pcs_gen(d0 = q,family_set = family_set)
          if(margin == "unif"){
            cop <- vinecop_dist(pair_copulas = pcs,structure = cvine_structure(q:1))
            x <- rvinecop(n=n,cop)
            y <- rvinecop(n=n,cop)
          } else if(margin == "norm"){
            cop <- vine_dist(margins = list(list(distr = "norm")), pair_copulas = pcs,
                             structure = cvine_structure(q:1))
            x <- rvine(n=n,cop)
            y <- rvine(n=n,cop)
          }
          
          res <- xi.vec(x,y,perm.inv = TRUE,nboot = nboot,ncores = ncores)
          est_comp[b,n.ind,q.ind,f.ind,m.ind] <- res$xi.vec
          pval_comp[b,n.ind,q.ind,f.ind,m.ind] <- res$pval
          cat(b,"\r")
        }
        
        saveRDS(est_comp,"est_comp_indep.rds")
        saveRDS(pval_comp,"pval_comp_indep.rds")
        cat(paste0("\n+ margin = ",margin,", n = ",n,", q = ",q,
                   ", family = ",ifelse(f.ind != 7,family_set,"mixed"),
                   ", size = ", mean(pval_comp[,n.ind,q.ind,f.ind,m.ind]<0.05)),"\n")
        
      } 
    }
  }
}


##########################
# Linear Dependence Case #
##########################

set.seed(2) # seed for reproducibility
family_set_all <- list("indep","gaussian","gumbel","clayton","frank","joe",
                       c("indep","gaussian","gumbel","clayton","frank","joe"))
margin_all <- c("unif","norm")
n_all <- c(50,100,200)
q_all <- c(2,3,4)

pval_comp <- array(NA,dim = c(nrep,length(n_all),length(q_all),
                              length(family_set_all),length(margin_all)))
est_comp <- array(NA,dim = c(nrep,length(n_all),length(q_all),
                             length(family_set_all),length(margin_all)))

for(m.ind in 1:length(margin_all)){
  margin <- margin_all[m.ind]
  
  for(n.ind in 1:length(n_all)){
    n <- n_all[n.ind]
    
    for(q.ind in 1:length(q_all)){
      q <- q_all[q.ind]
      
      for(f.ind in 1:length(family_set_all)){
        family_set <- family_set_all[[f.ind]]
        
        for(b in 1:nrep){
          pcs <- pcs_gen(d0 = q,family_set = family_set)
          if(margin == "unif"){
            cop <- vinecop_dist(pair_copulas = pcs,structure = cvine_structure(q:1))
            x <- rvinecop(n=n,cop)
            y <- 0.7*x + 0.3*rvinecop(n=n,cop)
          } else if(margin == "norm"){
            cop <- vine_dist(margins = list(list(distr = "norm")), pair_copulas = pcs,
                             structure = cvine_structure(q:1))
            x <- rvine(n=n,cop)
            y <- 0.7*x + 0.3*rvine(n=n,cop)
          }
          
          res <- xi.vec(x,y,perm.inv = TRUE,nboot = nboot,ncores = ncores)
          est_comp[b,n.ind,q.ind,f.ind,m.ind] <- res$xi.vec
          pval_comp[b,n.ind,q.ind,f.ind,m.ind] <- res$pval
          cat(b,"\r")
        }
        
        saveRDS(est_comp,"est_comp_linear.rds")
        saveRDS(pval_comp,"pval_comp_linear.rds")
        cat(paste0("\n+ margin = ",margin,", n = ",n,", q = ",q,
                   ", family = ",ifelse(f.ind != 7,family_set,"mixed"),
                   ", size = ", mean(pval_comp[,n.ind,q.ind,f.ind,m.ind]<0.05)),"\n")
        
      } 
    }
  }
}

#############################
# Quadratic Dependence Case #
#############################

set.seed(3) # seed for reproducibility
family_set_all <- list("indep","gaussian","gumbel","clayton","frank","joe",
                       c("indep","gaussian","gumbel","clayton","frank","joe"))
margin_all <- c("unif","norm")
n_all <- c(50,100,200)
q_all <- c(2,3,4)

pval_comp <- array(NA,dim = c(nrep,length(n_all),length(q_all),
                              length(family_set_all),length(margin_all)))
est_comp <- array(NA,dim = c(nrep,length(n_all),length(q_all),
                             length(family_set_all),length(margin_all)))

for(m.ind in 1:length(margin_all)){
  margin <- margin_all[m.ind]
  
  for(n.ind in 1:length(n_all)){
    n <- n_all[n.ind]
    
    for(q.ind in 1:length(q_all)){
      q <- q_all[q.ind]
      
      for(f.ind in 1:length(family_set_all)){
        family_set <- family_set_all[[f.ind]]
        
        for(b in 1:nrep){
          pcs <- pcs_gen(d0 = q,family_set = family_set)
          if(margin == "unif"){
            cop <- vinecop_dist(pair_copulas = pcs,structure = cvine_structure(q:1))
            x <- rvinecop(n=n,cop)
            y <- 0.7*x^2 + 0.3*rvinecop(n=n,cop)
          } else if(margin == "norm"){
            cop <- vine_dist(margins = list(list(distr = "norm")), pair_copulas = pcs,
                             structure = cvine_structure(q:1))
            x <- rvine(n=n,cop)
            y <- 0.7*x^2 + 0.3*rvine(n=n,cop)
          }
          
          res <- xi.vec(x,y,perm.inv = TRUE,nboot = nboot,ncores = ncores)
          est_comp[b,n.ind,q.ind,f.ind,m.ind] <- res$xi.vec
          pval_comp[b,n.ind,q.ind,f.ind,m.ind] <- res$pval
          cat(b,"\r")
        }
        
        saveRDS(est_comp,"est_comp_quad.rds")
        saveRDS(pval_comp,"pval_comp_quad.rds")
        cat(paste0("\n+ margin = ",margin,", n = ",n,", q = ",q,
                   ", family = ",ifelse(f.ind != 7,family_set,"mixed"),
                   ", size = ", mean(pval_comp[,n.ind,q.ind,f.ind,m.ind]<0.05)),"\n")
        
      } 
    }
  }
}

############################
# Sinusoid Dependence Case #
############################

set.seed(4) # seed for reproducibility
family_set_all <- list("indep","gaussian","gumbel","clayton","frank","joe",
                       c("indep","gaussian","gumbel","clayton","frank","joe"))
margin_all <- c("unif","norm")
n_all <- c(50,100,200)
q_all <- c(2,3,4)

pval_comp <- array(NA,dim = c(nrep,length(n_all),length(q_all),
                              length(family_set_all),length(margin_all)))
est_comp <- array(NA,dim = c(nrep,length(n_all),length(q_all),
                             length(family_set_all),length(margin_all)))

for(m.ind in 1:length(margin_all)){
  margin <- margin_all[m.ind]
  
  for(n.ind in 1:length(n_all)){
    n <- n_all[n.ind]
    
    for(q.ind in 1:length(q_all)){
      q <- q_all[q.ind]
      
      for(f.ind in 1:length(family_set_all)){
        family_set <- family_set_all[[f.ind]]
        
        for(b in 1:nrep){
          pcs <- pcs_gen(d0 = q,family_set = family_set)
          if(margin == "unif"){
            cop <- vinecop_dist(pair_copulas = pcs,structure = cvine_structure(q:1))
            x <- rvinecop(n=n,cop)
            y <- 0.7*sin(2*pi*x) + 0.3*rvinecop(n=n,cop)
          } else if(margin == "norm"){
            cop <- vine_dist(margins = list(list(distr = "norm")), pair_copulas = pcs,
                             structure = cvine_structure(q:1))
            x <- rvine(n=n,cop)
            y <- 0.7*sin(2*pi*x/3) + 0.3*rvine(n=n,cop)
          }
          
          res <- xi.vec(x,y,perm.inv = TRUE,nboot = nboot,ncores = ncores)
          est_comp[b,n.ind,q.ind,f.ind,m.ind] <- res$xi.vec
          pval_comp[b,n.ind,q.ind,f.ind,m.ind] <- res$pval
          cat(b,"\r")
        }
        
        saveRDS(est_comp,"est_comp_sin.rds")
        saveRDS(pval_comp,"pval_comp_sin.rds")
        cat(paste0("\n+ margin = ",margin,", n = ",n,", q = ",q,
                   ", family = ",ifelse(f.ind != 7,family_set,"mixed"),
                   ", size = ", mean(pval_comp[,n.ind,q.ind,f.ind,m.ind]<0.05)),"\n")
        
      } 
    }
  }
}


### Importing results
est_indep <- readRDS("01 Simulations/NVC JASA/Size and Power/est_comp_indep.rds")
pval_indep <- readRDS("01 Simulations/NVC JASA/Size and Power/pval_comp_indep.rds")
est_linear <- readRDS("01 Simulations/NVC JASA/Size and Power/est_comp_linear.rds")
pval_linear <- readRDS("01 Simulations/NVC JASA/Size and Power/pval_comp_linear.rds")
est_quad <- readRDS("01 Simulations/NVC JASA/Size and Power/est_comp_quad.rds")
pval_quad <- readRDS("01 Simulations/NVC JASA/Size and Power/pval_comp_quad.rds")
est_sin <- readRDS("01 Simulations/NVC JASA/Size and Power/est_comp_sin.rds")
pval_sin <- readRDS("01 Simulations/NVC JASA/Size and Power/pval_comp_sin.rds")

### Organizing results
pval.mean_indep <- apply(pval_indep<0.05,MARGIN = c(2,3,4,5),FUN = mean,na.rm = TRUE)
pval.mean_linear <- apply(pval_linear<0.05,MARGIN = c(2,3,4,5),FUN = mean,na.rm = TRUE)
pval.mean_quad <- apply(pval_quad<0.05,MARGIN = c(2,3,4,5),FUN = mean,na.rm = TRUE)
pval.mean_sin <- apply(pval_sin<0.05,MARGIN = c(2,3,4,5),FUN = mean,na.rm = TRUE)

# Uniform margins
ps_indep1 <- NULL
ps_linear1 <- NULL
ps_quad1 <- NULL
ps_sin1 <- NULL
for(t.ind in 1:7){
  ps_indep1 <- rbind(ps_indep1,as.vector(pval.mean_indep[,,t.ind,1]))
  ps_linear1 <- rbind(ps_linear1,as.vector(pval.mean_linear[,,t.ind,1]))
  ps_quad1 <- rbind(ps_quad1,as.vector(pval.mean_quad[,,t.ind,1]))
  ps_sin1 <- rbind(ps_sin1,as.vector(pval.mean_sin[,,t.ind,1]))
}

# Gaussian margins
ps_indep2 <- NULL
ps_linear2 <- NULL
ps_quad2 <- NULL
ps_sin2 <- NULL
for(t.ind in 1:7){
  ps_indep2 <- rbind(ps_indep2,as.vector(pval.mean_indep[,,t.ind,2]))
  ps_linear2 <- rbind(ps_linear2,as.vector(pval.mean_linear[,,t.ind,2]))
  ps_quad2 <- rbind(ps_quad2,as.vector(pval.mean_quad[,,t.ind,2]))
  ps_sin2 <- rbind(ps_sin2,as.vector(pval.mean_sin[,,t.ind,2]))
}

# Table of size and power
types <- rep(c("Independence","Gaussian","Gumbel","Clayton","Frank","Joe","Mixed"),4)
func <- rep(c("Independent","Linear","Quadratic","Sinusoid"),each = 7)
cbind(func,types,round(rbind(ps_indep1,ps_linear1,ps_quad1,ps_sin1),3)) # uniform margins
cbind(func,types,round(rbind(ps_indep2,ps_linear2,ps_quad2,ps_sin2),3)) # gaussian margins


################################################################
##### Illustration on nonlinear spectral vector dependence #####
################################################################

set.seed(1) # seed for reproducibility
samprate <- 100
nsec_all <- c(50,100,200)
rhythms <- c("delta","theta","alpha","beta2","gamma2")

est_comp <- array(NA,dim = c(samprate/2,1000,3))
coh_comp <- array(NA,dim = c(samprate/2,1000,3))
cancoh_comp <- array(NA,dim = c(samprate/2,1000,3))
bandcoh_comp <- array(NA,dim = c(5,1000,3))
kencoh_comp <- array(NA,dim = c(5,1000,3))

for(n in 1:3){
  nsec <- nsec_all[n]
  for(r in 1:1000){
    
    dat <- sim.nl.signal(nsec = nsec,samprate = 100,freq = 21,
                         q = 2,noise.sd = 0.1)
    dat2 <- sim.nl.signal(nsec = nsec,samprate = 100,freq = 22,
                          q = 2,noise.sd = 0.1)
    dat3 <- sim.nl.signal(nsec = nsec,samprate = 100,freq = 23,
                          q = 2,noise.sd = 0.1)
    dat4 <- sim.nl.signal(nsec = nsec,samprate = 100,freq = 24,
                          q = 2,noise.sd = 0.1)
    dat5 <- sim.nl.signal(nsec = nsec,samprate = 100,freq = 25,
                          q = 2,noise.sd = 0.1)
    
    x <-0.9*(scale(dat$X)+scale(dat2$X)+scale(dat3$X)+scale(dat4$X)+scale(dat5$X)) + 0.1*matrix(rnorm(nrow(dat$X)*ncol(dat$X)),ncol = ncol(dat$X))
    y <-0.9*(scale(dat$Y)+scale(dat2$Y)+scale(dat3$Y)+scale(dat4$Y)+scale(dat5$Y)) + 0.1*matrix(rnorm(nrow(dat$Y)*ncol(dat$Y)),ncol = ncol(dat$Y))
    
    samp <- nvc(x = x, y = y,blength = samprate,perm.inv = TRUE,
                b.scale = FALSE,nboot = 0,ncores = 8)
    est_comp[,r,n] <- samp$nvc
    
    coh.all <- array(NA,dim = c(samprate/2,ncol(x)*ncol(y)))
    bandcoh.all <- array(NA,dim = c(5,ncol(x)*ncol(y)))
    
    for(ch1 in 1:ncol(x)){
      for(ch2 in 1:ncol(y)){
        coh.all[,ncol(y)*(ch1-1)+ch2] <- coh(x = x[,ch1],y = y[,ch2],blength = samprate)$coh
        bandcoh.all[,ncol(y)*(ch1-1)+ch2] <- band.coh(x = x[,ch1],y = y[,ch2],samprate = samprate,
                                                      fband = rhythms)
      }
    }
    
    coh_comp[,r,n] <- apply(coh.all,MARGIN = 1,FUN = mean,na.rm = TRUE)
    bandcoh_comp[,r,n] <- apply(bandcoh.all,MARGIN = 1,FUN = mean,na.rm = TRUE)
    
    cancoh_comp[,r,n] <- cancoh(x = x,y = y,samprate = samprate,
                                blength = samprate)$cancoh[-1]
    
    kencoh_comp[,r,n] <- kencoh(x,y,samprate = samprate,fband = rhythms)
    
    cat(r,"\r")
    
  }
}
saveRDS(est_comp,"est_nl.rds")
saveRDS(coh_comp,"coh_nl.rds")
saveRDS(cancoh_comp,"cancoh_nl.rds")
saveRDS(bandcoh_comp,"bandcoh_nl.rds")
saveRDS(kencoh_comp,"kencoh_nl.rds")

### Importing results
est_nl <- readRDS("01 Simulations/NVC JASA/Nonlinear Example/est_nl.rds")
coh_nl <- readRDS("01 Simulations/NVC JASA/Nonlinear Example/coh_nl.rds")
cancoh_nl <- readRDS("01 Simulations/NVC JASA/Nonlinear Example/cancoh_nl.rds")
bandcoh_nl <- readRDS("01 Simulations/NVC JASA/Nonlinear Example/bandcoh_nl.rds")
kencoh_nl <- readRDS("01 Simulations/NVC JASA/Nonlinear Example/kencoh_nl.rds")

### Organizing results
samprate <- 100
freq <- 1:50

ind_delta <- freq[freq > 0.5 & freq <= 4]
ind_theta <- freq[freq > 4 & freq < 8]
ind_alpha <- freq[freq >= 8 & freq <= 12]
ind_beta <- freq[freq >= 21 & freq <= 25]
ind_gamma <- freq[freq >= 36 & freq <= 40]

est.mean <- apply(est_nl,MARGIN = c(1,3),FUN = mean,na.rm = TRUE)
est.hi <- apply(est_nl,MARGIN = c(1,3),FUN = quantile,probs = 0.975,na.rm = TRUE)
est.lo <- apply(est_nl,MARGIN = c(1,3),FUN = quantile,probs = 0.025,na.rm = TRUE)

cancoh.mean <- apply(cancoh_nl,MARGIN = c(1,3),FUN = mean,na.rm = TRUE)
cancoh.hi <- apply(cancoh_nl,MARGIN = c(1,3),FUN = quantile,probs = 0.975,na.rm = TRUE)
cancoh.lo <- apply(cancoh_nl,MARGIN = c(1,3),FUN = quantile,probs = 0.025,na.rm = TRUE)

est.fband_nl <- array(NA,dim = c(5,1000,3))
est.fband_nl[1,,] <- apply(est_nl[ind_delta,,],MARGIN = c(2,3),FUN = mean,na.rm = TRUE)
est.fband_nl[2,,] <- apply(est_nl[ind_theta,,],MARGIN = c(2,3),FUN = mean,na.rm = TRUE)
est.fband_nl[3,,] <- apply(est_nl[ind_alpha,,],MARGIN = c(2,3),FUN = mean,na.rm = TRUE)
est.fband_nl[4,,] <- apply(est_nl[ind_beta,,],MARGIN = c(2,3),FUN = mean,na.rm = TRUE)
est.fband_nl[5,,] <- apply(est_nl[ind_gamma,,],MARGIN = c(2,3),FUN = mean,na.rm = TRUE)

### Visualizing results
png("nonlinear_with_cancoh.png",width = 2250, height = 600)
par(mfcol = c(1,3),mar = c(10,7,5,2),oma = c(4,1,1,1),xpd = FALSE)

nsec_all <- c(50,100,200)

for(n in 1:3){
  plot(NA,NA,type = "l",ylim = c(-0.2,1),xlim = c(0,50),xlab = "",ylab = "",cex.axis = 2.75)
  box(lwd = 4)
  mtext(paste0("n = ",nsec_all[n]),side = 3,cex = 3.5,font = 2,line = 1.5)
  mtext(TeX("Frequency $\\omega$ (in Hz)"),side = 1,cex = 2.5,line = 5)
  mtext("Estimates",side = 2,cex = 3,line = 3.5)
  abline(h = 0,col = "black",lty = 2,lwd = 2)
  
  polygon(x = c(freq, rev(freq)),
          y = c(est.hi[,n], rev(est.lo[,n])),
          col = adjustcolor("darkviolet",alpha.f = 0.15), border = NA)
  lines(freq,est.mean[,n],col = "darkviolet",lwd = 4)
  lines(freq,est.lo[,n],lty = 3,col = "darkviolet",lwd = 4)
  lines(freq,est.hi[,n],lty = 3,col = "darkviolet",lwd = 4)
  
  polygon(x = c(freq, rev(freq)),
          y = c(cancoh.hi[,n], rev(cancoh.lo[,n])),
          col = adjustcolor("black",alpha.f = 0.15), border = NA)
  lines(freq,cancoh.mean[,n],col = "black",lwd = 4)
  lines(freq,cancoh.lo[,n],lty = 3,col = "black",lwd = 4)
  lines(freq,cancoh.hi[,n],lty = 3,col = "black",lwd = 4)
  
  
}

mtext(TeX("n = 50",bold = TRUE,italic = TRUE),outer = TRUE,
      side = 2,cex = 3,line = 2,adj = 0.875)
mtext(TeX("n = 100",bold = TRUE,italic = TRUE),outer = TRUE,
      side = 2,cex = 3,line = 2,adj = 0.5)
mtext(TeX("n = 200",bold = TRUE,italic = TRUE),outer = TRUE,
      side = 2,cex = 3,line = 2,adj = 0.15)


par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom",legend = c("NVC  ", "CanCoh"),
       col = c("darkviolet","black"),
       lwd = 4, xpd = TRUE, horiz = TRUE, cex = 4,
       pt.bg = c(adjustcolor("darkviolet",alpha.f = 0.15),
                 adjustcolor("black",alpha.f = 0.15)),
       pch = 22,pt.cex = 10,bty = "n")

dev.off()


png("nonlinear_with_pbc_kencoh.png",width = 1200,height = 600)
par(mfrow = c(1,1),mar = c(7,6,1,1))
boxplot(cbind(est.fband_nl[4,,],NA,bandcoh_nl[4,,],NA,kencoh_nl[4,,]),ylim = c(-0.1,1),
        xaxt = "n",yaxt = "n",cex = 1,lwd = 2.5)
axis(side = 2, at = 0.2*(-1:5),labels = 0.2*(-1:5),lwd = 3,cex.axis = 2,font = 2)
axis(side = 1, at = c(1:3,5:7,9:11),labels = rep(c(50,100,200),3),lwd = 3,cex.axis = 2,padj = 0.5,font = 2)
box(lwd = 3)
abline(h = 0,lty = 2,lwd = 3)

axis(side = 3, at = c(1:3),labels = c(NA,NA),lwd = 3,line = -37)
axis(side = 3, at = 2,labels = "NVC",lwd = 3,cex.axis = 3,font = 2,line = -40.5,tick = FALSE)
axis(side = 3, at = c(5:7),labels = c(NA,NA),lwd = 3,line = -37)
axis(side = 3, at = 6,labels = "PBC",lwd = 5,cex.axis = 3,font = 2,line = -40.5,tick = FALSE)
axis(side = 3, at = c(9:11),labels = c(NA,NA),lwd = 3,line = -37)
axis(side = 3, at = 10,labels = "KenCoh",lwd = 5,cex.axis = 3,font = 2,line = -40.5,tick = FALSE)

mtext("n = ",side = 1,cex = 2,font = 2,line = 1.5,adj = 0)
mtext("Estimates",font = 2,side = 2,cex = 2.5,line = 4)


dev.off()




#############################################################################
##### Empirical unbiasedness for choice of the number of permutations Q #####
#############################################################################

nsec_all <- c(100,200)
samprate <- 100

freq <- c(0.02,0.22,0.40)
est_comp <- array(NA,dim = c(length(freq),1000,4,2))

set.seed(1) # seed for reproducibility

for(n.ind in 1:length(nsec_all)){
  nsec <- nsec_all[n.ind]
  
  for(r in 1:1000){
    
    Z_beta <- sim.fband(nsec = nsec,samprate = samprate,fband = rep("beta2",6))
    Tlength <- nsec*samprate
    
    x1 <- 0.95*Z_beta[,1] + 0.05*rnorm(Tlength)
    x2 <- 0.95*Z_beta[,2] + 0.05*rnorm(Tlength)
    x3 <- 0.90*Z_beta[,3] + 0.10*rnorm(Tlength)
    x4 <- 0.90*Z_beta[,4] + 0.10*rnorm(Tlength)
    x5 <- 0.85*Z_beta[,5] + 0.15*rnorm(Tlength)
    x6 <- 0.85*Z_beta[,6] + 0.15*rnorm(Tlength)
    
    y1 <- 0.95*Z_beta[,1] + 0.05*rnorm(Tlength)
    y2 <- 0.95*Z_beta[,2] + 0.05*rnorm(Tlength)
    y3 <- 0.90*Z_beta[,3] + 0.10*rnorm(Tlength)
    y4 <- 0.90*Z_beta[,4] + 0.10*rnorm(Tlength)
    y5 <- 0.85*Z_beta[,5] + 0.15*rnorm(Tlength)
    y6 <- 0.85*Z_beta[,6] + 0.15*rnorm(Tlength)
    
    x <- cbind(x1,x2,x3,x4,x5,x6)
    y <- cbind(y1,y2,y3,y4,y5,y6)
    
    samp <- nvc(x = x, y = y,blength = samprate,perm.inv = TRUE,nperm = 30,freq = freq)
    est_comp[,r,1,n.ind] <- samp$nvc
    
    samp <- nvc(x = x, y = y,blength = samprate,perm.inv = TRUE,nperm = 60,freq = freq)
    est_comp[,r,2,n.ind] <- samp$nvc
    
    samp <- nvc(x = x, y = y,blength = samprate,perm.inv = TRUE,nperm = 120,freq = freq)
    est_comp[,r,3,n.ind] <- samp$nvc
    
    samp <- nvc(x = x, y = y,blength = samprate,perm.inv = TRUE,freq = freq)
    est_comp[,r,4,n.ind] <- samp$nvc
    
    saveRDS(est_comp,"simest_permsens.rds")
    
    cat(paste0("+ nsec = ",nsec,", r = ",r),"\r")
    
  }
  
}

### Importing and visualizing results
est_comp <- readRDS("simest_permsens.rds")[,,,1] # n = 100 case

png("permsens_n100.png",width = 1600,height = 700)
par(mar = c(7,7,4,1))
boxplot(cbind(est_comp[1,,],NA,est_comp[2,,],NA,est_comp[3,,]),ylim = c(-0.1,0.7),
        xaxt = "n",yaxt = "n",cex = 1,lwd = 3.5)
title("(a) n = 100",line = 1,adj = 0,lwd = 4,cex.main = 4)
axis(side = 2, at = 0.1*(-1:5),labels = 0.1*(-1:5),lwd = 5,cex.axis = 2,font = 2)
axis(side = 1, at = c(1:4,6:9,11:14),labels = rep(c(30,60,120,720),3),lwd = 5,cex.axis = 2,padj = 0.5,font = 2)
box(lwd = 5)
mtext(TeX("$\\bar{T}_{\\omega_{k},n}$  estimates",bold = TRUE),side = 2,cex = 3,line = 3)
mtext(TeX("Number of random permutations $Q$",bold = TRUE),side = 1,cex = 3,line = 5)

axis(side = 1, at = c(1,4),labels = c(NA,NA),lwd = 5,line = -17)
axis(side = 1, at = 2.5,labels = TeX("$\\omega_k = 0.02$  $\\approx$  at 2 Hz",bold = TRUE),
     lwd = 5,cex.axis = 2,line = -20,tick = FALSE)

axis(side = 1, at = c(6,9),labels = c(NA,NA),lwd = 5,line = -31)
axis(side = 1, at = 7.5,labels = TeX("$\\omega_k = 0.22$  $\\approx$  at 22 Hz",bold = TRUE),
     lwd = 5,cex.axis = 2,line = -34,tick = FALSE)

axis(side = 1, at = c(11,14),labels = c(NA,NA),lwd = 5,line = -17)
axis(side = 1, at = 12.5,labels = TeX("$\\omega_k = 0.40$  $\\approx$  at 40 Hz",bold = TRUE),
     lwd = 5,cex.axis = 2,line = -20,tick = FALSE)

abline(h = 0,lty = 2,lwd = 3)
dev.off()

est_comp <- readRDS("simest_permsens.rds")[,,,2] # n = 200 case

png("permsens_n200.png",width = 1600,height = 700)
par(mar = c(7,7,4,1))
boxplot(cbind(est_comp[1,,],NA,est_comp[2,,],NA,est_comp[3,,]),ylim = c(-0.1,0.7),
        xaxt = "n",yaxt = "n",cex = 1,lwd = 3.5)
title("(b) n = 200",line = 1,adj = 0,lwd = 4,cex.main = 4)
axis(side = 2, at = 0.1*(-1:5),labels = 0.1*(-1:5),lwd = 5,cex.axis = 2,font = 2)
axis(side = 1, at = c(1:4,6:9,11:14),labels = rep(c(30,60,120,720),3),lwd = 5,cex.axis = 2,padj = 0.5,font = 2)
box(lwd = 5)
mtext(TeX("$\\bar{T}_{\\omega_{k},n}$  estimates",bold = TRUE),side = 2,cex = 3,line = 3)
mtext(TeX("Number of random permutations $Q$",bold = TRUE),side = 1,cex = 3,line = 5)

axis(side = 1, at = c(1,4),labels = c(NA,NA),lwd = 5,line = -17)
axis(side = 1, at = 2.5,labels = TeX("$\\omega_k = 0.02$  $\\approx$  at 2 Hz",bold = TRUE),
     lwd = 5,cex.axis = 2,line = -20,tick = FALSE)

axis(side = 1, at = c(6,9),labels = c(NA,NA),lwd = 5,line = -31)
axis(side = 1, at = 7.5,labels = TeX("$\\omega_k = 0.22$  $\\approx$  at 22 Hz",bold = TRUE),
     lwd = 5,cex.axis = 2,line = -34,tick = FALSE)

axis(side = 1, at = c(11,14),labels = c(NA,NA),lwd = 5,line = -17)
axis(side = 1, at = 12.5,labels = TeX("$\\omega_k = 0.40$  $\\approx$  at 40 Hz",bold = TRUE),
     lwd = 5,cex.axis = 2,line = -20,tick = FALSE)

abline(h = 0,lty = 2,lwd = 3)
dev.off()

