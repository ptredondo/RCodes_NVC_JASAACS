##### EEG analysis: NVC, PBC, KenCoh, nonspectral T measure, RBP and sensitivity analyses
# Download OpenNeuro dataset ds004504, version 1.0.9, from:
# https://openneuro.org/datasets/ds004504/versions/1.0.9
# Save its contents in Dataset/ folder (participants.tsv and derivatives/ must be
# directly inside Dataset/).
# install.packages(c("R.matlab", "permuco", "LSTS", "remotes"))
# devtools::install_github("ptredondo/nvc")

##### Setting working directory
setwd("working/directory")

##### Loading needed packages
library(nvc)
library(R.matlab)
library(igraph)
library(caret)
library(e1071)
library(class)
library(randomForest)
library(dplyr)

##### Subject information
subj_info <- read.table("Dataset/participants.tsv",header = TRUE,sep = "\t")
subj_all <- subj_info$participant_id

nsec_all <- rep(NA,length(subj_all))
for(subj in 1:length(subj_all)){
  data <- readMat(paste0("Dataset/derivatives/",subj_all[subj],"/eeg/",subj_all[subj],"_task-eyesclosed_eeg.set"))
  nsec_all[subj] <- ncol(data$data)/500
  cat(paste0("+ Subject = ",subj),"\r")
}
subj_info$nsec <- nsec_all

##### Random start for downsampling
set.seed(1)
start_subj <- sample(1:5,length(subj_all),replace = TRUE)
saveRDS(start_subj,"start_subj_s100.rds")

start_subj <- sample(1:2,length(subj_all),replace = TRUE)
saveRDS(start_subj,"start_subj_s250.rds")


##### Channel names
subj <- 1
data <- readMat(paste0("Dataset/derivatives/",subj_all[subj],"/eeg/",subj_all[subj],"_task-eyesclosed_eeg.set"))
chnames <- rep(NA,19)
for(ch in 1:19){
  chnames[ch] <- data$chanlocs[,,ch]$labels[1,1]
}

##### Regions of interest (ROIs)
ROI <- list(LF = c("Fp1", "F3", "F7"),
            RF = c("Fp2", "F4", "F8"),
            LT = c("T3", "T5"),
            RT = c("T4", "T6"),
            C = c("C3", "C4", "Cz"),
            P = c("P3", "P4", "Pz"),
            O = c("O1", "O2"))

##### ROI and band names
group_names <- c("RF", "RT", "C", "O", "P", "LT", "LF")
band_names <- c("delta", "theta", "alpha", "beta", "gamma")

#########################################
##### Generate results based on NVC #####
#########################################

generate.results.NVC <- function(samprate,blength = samprate,scale.overall = TRUE,
                                 conn_mat.name,pval.name,pval_adj.name){
  
  freq_all <- seq(0,0.5,length = blength/2+1)[-1]
  freq_all <- freq_all[freq_all*samprate<=50]
  
  conn_mat <- array(NA,dim = c(length(ROI),length(ROI),length(freq_all),length(subj_all)))
  
  for(subj in 1:length(subj_all)){
    data <- readMat(paste0("Dataset/derivatives/",subj_all[subj],"/eeg/",subj_all[subj],"_task-eyesclosed_eeg.set"))
    eeg_mat <- t(data$data)
    if(scale.overall){
      eeg_mat <- scale(eeg_mat[5*500 + 1:((floor(nsec_all[subj])-10)*500),])
    } else {
      eeg_mat <- eeg_mat[5*500 + 1:((floor(nsec_all[subj])-10)*500),]
    }
    
    if(samprate == 100){
      start_subj <- readRDS("start_subj_s100.rds")
      eeg_mat <- eeg_mat[start_subj[subj] + (0:((floor(nsec_all[subj])-10)*100-1))*5,]
    }
    
    if(samprate == 250){
      start_subj <- readRDS("start_subj_s250.rds")
      eeg_mat <- eeg_mat[start_subj[subj] + (0:((floor(nsec_all[subj])-10)*250-1))*2,]
    }
    
    for(r1 in 1:length(ROI)){
      for(r2 in 1:length(ROI)){
        if(r1 == r2){next}
        
        elec_r1 <- (1:19)[chnames %in% ROI[[r1]]]
        elec_r2 <- (1:19)[chnames %in% ROI[[r2]]]
        
        y_mat <- eeg_mat[,elec_r1]
        x_mat <- eeg_mat[,elec_r2]
        
        est <- nvc(x = x_mat, y = y_mat,b.scale = FALSE,perm.inv = TRUE,nperm = NULL,
                   blength = blength,freq = freq_all,nboot = 0)
        
        conn_mat[r1,r2,,subj] <- est$nvc
        
        saveRDS(conn_mat,conn_mat.name)
        
        cat(paste("+ DONE: Subject = ",subj,", R1 = ",r1,", R2 = ",r2,sep = ""),"\r")
        
      }  
    }
  }
  
  symm_mat <- array(NA,dim = dim(conn_mat))
  for(subj in 1:length(subj_all)){
    for(r1 in 1:length(ROI)){
      for(r2 in 1:length(ROI)){
        if(r1 == r2){next}
        symm_mat[r1,r2,,subj] <- ifelse(conn_mat[r1,r2,,subj] > conn_mat[r2,r1,,subj],conn_mat[r1,r2,,subj],conn_mat[r2,r1,,subj])
      }
    }
  }
  
  ind_delta <- (1:length(freq_all))[freq_all*samprate > 0.5 & freq_all*samprate <= 4]
  ind_theta <- (1:length(freq_all))[freq_all*samprate > 4 & freq_all*samprate <= 8]
  ind_alpha <- (1:length(freq_all))[freq_all*samprate > 8 & freq_all*samprate <= 12]
  ind_beta <- (1:length(freq_all))[freq_all*samprate > 12 & freq_all*samprate <= 30]
  ind_gamma <- (1:length(freq_all))[freq_all*samprate > 30 & freq_all*samprate <= 45]
  
  symm_delta <- apply(symm_mat[,,ind_delta,],MARGIN = c(1,2,4),FUN = mean,na.rm = TRUE)
  symm_theta <- apply(symm_mat[,,ind_theta,],MARGIN = c(1,2,4),FUN = mean,na.rm = TRUE)
  symm_alpha <- apply(symm_mat[,,ind_alpha,],MARGIN = c(1,2,4),FUN = mean,na.rm = TRUE)
  symm_beta <- apply(symm_mat[,,ind_beta,],MARGIN = c(1,2,4),FUN = mean,na.rm = TRUE)
  symm_gamma <- apply(symm_mat[,,ind_gamma,],MARGIN = c(1,2,4),FUN = mean,na.rm = TRUE)
  
  symm_fband <- array(NA,dim = c(length(ROI),length(ROI),5,length(subj_all)))
  symm_fband[,,1,] <- symm_delta
  symm_fband[,,2,] <- symm_theta
  symm_fband[,,3,] <- symm_alpha
  symm_fband[,,4,] <- symm_beta
  symm_fband[,,5,] <- symm_gamma
  
  ### Testing by pair
  nboot <- 10^5
  pval_pair <- array(NA,dim = c(length(ROI),length(ROI),5,3))
  
  for(r1 in 1:6){
    for(r2 in (r1+1):7){
      for(f in 1:5){
        group <- c(rep("delta_A",36),rep("delta_C",29),rep("delta_F",23))
        data_fband <- data.frame(group = group,x = symm_fband[r1,r2,f,])
        pval_fband_ac <- aovperm(x ~ group, data = data_fband[1:65,], np = nboot)$table[1,5]
        pval_fband_af <- aovperm(x ~ group, data = data_fband[c(1:36,66:88),], np = nboot)$table[1,5]
        pval_fband_cf <- aovperm(x ~ group, data = data_fband[37:88,], np = nboot)$table[1,5]
        pval_pair[r1,r2,f,1] <- pval_fband_ac
        pval_pair[r1,r2,f,2] <- pval_fband_af
        pval_pair[r1,r2,f,3] <- pval_fband_cf
        
      }
      saveRDS(pval_pair,pval.name)
      cat(paste("+ DONE: R1 = ",r1,", R2 = ",r2,sep = ""),"\r")
    }
  }
  
  
  pval_all <- NULL
  for(f in 1:5){
    for(r1 in 1:6){
      for(r2 in (r1+1):7){
        for(p in 1:3){
          pval_all <- rbind(pval_all,data.frame(r1 = r1,r2 = r2,f = f,p = p,pval = pval_pair[r1,r2,f,p]))
        }
      }
    }
  }
  
  pval_all$pval_adj <- p.adjust(pval_all$pval,method = "BH")
  sum(pval_all$pval_adj < 0.1,na.rm = TRUE)
  pval_pair_adj <- array(NA,dim = c(length(ROI),length(ROI),5,3))
  for(freq in 1:5){
    for(roi1 in 1:6){
      for(roi2 in (roi1+1):7){
        for(pair in 1:3){
          pval_pair_adj[roi1,roi2,freq,pair] <- subset(pval_all,r1 == roi1 & r2 == roi2 & f == freq & p == pair)$pval_adj
          pval_pair_adj[roi2,roi1,freq,pair] <- subset(pval_all,r1 == roi1 & r2 == roi2 & f == freq & p == pair)$pval_adj
        }
      }
    }
  }
  saveRDS(pval_pair_adj,pval_adj.name)
}

##### Sensitivitiy analysis for downsampling and standardization effect
set.seed(1) # seed for reproducibility

### 100 Hz, standardized
generate.results.NVC(samprate = 100,scale.overall = TRUE,
                     conn_mat.name = "conn_mat_s100_sct.rds",
                     pval.name = "pval_s100_sct.rds",
                     pval_adj.name = "pval_adj_s100_sct.rds")

### 250 Hz, standardized
generate.results.NVC(samprate = 250,scale.overall = TRUE,
                     conn_mat.name = "conn_mat_s250_sct.rds",
                     pval.name = "pval_s250_sct.rds",
                     pval_adj.name = "pval_adj_s250_sct.rds")

### 500 Hz, standardized
generate.results.NVC(samprate = 500,scale.overall = TRUE,
                     conn_mat.name = "conn_mat_s500_sct.rds",
                     pval.name = "pval_s500_sct.rds",
                     pval_adj.name = "pval_adj_s500_sct.rds")

### 100 Hz, unstandardized
generate.results.NVC(samprate = 100,scale.overall = FALSE,
                     conn_mat.name = "conn_mat_s100_scf.rds",
                     pval.name = "pval_s100_scf.rds",
                     pval_adj.name = "pval_adj_s100_scf.rds")

### 250 Hz, unstandardized
generate.results.NVC(samprate = 250,scale.overall = FALSE,
                     conn_mat.name = "conn_mat_s250_scf.rds",
                     pval.name = "pval_s250_scf.rds",
                     pval_adj.name = "pval_adj_s250_scf.rds")

### 500 Hz, unstandardized
generate.results.NVC(samprate = 500,scale.overall = FALSE,
                     conn_mat.name = "conn_mat_s500_scf.rds",
                     pval.name = "pval_s500_scf.rds",
                     pval_adj.name = "pval_adj_s500_scf.rds")

##### Sensitivitiy analysis for choice of block length effect

### B = 750 (1.5 second blocks)
generate.results.NVC(samprate = 500,blength = 750,scale.overall = FALSE,
                     conn_mat.name = "conn_mat_s500_b750_scf.rds",
                     pval.name = "pval_s500_b750_scf.rds",
                     pval_adj.name = "pval_adj_s500_b750_scf.rds")

### B = 1000 (2 second blocks)
generate.results.NVC(samprate = 500,blength = 1000,scale.overall = FALSE,
                     conn_mat.name = "conn_mat_s500_b1000_scf.rds",
                     pval.name = "pval_s500_b1000_scf.rds",
                     pval_adj.name = "pval_adj_s500_b1000_scf.rds")

#########################################
##### Generate results based on PBC #####
#########################################

generate.results.PBC <- function(samprate,scale.overall = TRUE,
                                 conn_mat.name,pval.name,pval_adj.name,
                                 fband = c("delta","theta","alpha","beta","gamma")){
  
  conn_mat <- array(NA,dim = c(length(ROI),length(ROI),5,length(subj_all)))
  
  for(subj in 1:length(subj_all)){
    data <- readMat(paste0("Dataset/derivatives/",subj_all[subj],"/eeg/",subj_all[subj],"_task-eyesclosed_eeg.set"))
    eeg_mat <- t(data$data)
    if(scale.overall){
      eeg_mat <- scale(eeg_mat[5*500 + 1:((floor(nsec_all[subj])-10)*500),])
    } else {
      eeg_mat <- eeg_mat[5*500 + 1:((floor(nsec_all[subj])-10)*500),]
    }
    
    if(samprate == 100){
      start_subj <- readRDS("start_subj_s100.rds")
      eeg_mat <- eeg_mat[start_subj[subj] + (0:((floor(nsec_all[subj])-10)*100-1))*5,]
    }
    
    if(samprate == 250){
      start_subj <- readRDS("start_subj_s250.rds")
      eeg_mat <- eeg_mat[start_subj[subj] + (0:((floor(nsec_all[subj])-10)*250-1))*2,]
    }
    
    for(r1 in 1:length(ROI)){
      for(r2 in 1:length(ROI)){
        if(r1 >= r2){next}
        
        elec_r1 <- (1:19)[chnames %in% ROI[[r1]]]
        elec_r2 <- (1:19)[chnames %in% ROI[[r2]]]
        
        y_mat <- eeg_mat[,elec_r1]
        x_mat <- eeg_mat[,elec_r2]
        
        bandcoh.all <- array(NA,dim = c(5,ncol(x_mat)*ncol(y_mat)))
        
        for(ch1 in 1:ncol(x_mat)){
          for(ch2 in 1:ncol(y_mat)){
            bandcoh.all[,ncol(y_mat)*(ch1-1)+ch2] <- band.coh(x = x_mat[,ch1],y = y_mat[,ch2],samprate = samprate,
                                                              fband = fband)
          }
        }
        
        conn_mat[r1,r2,,subj] <- conn_mat[r2,r1,,subj] <- apply(bandcoh.all,MARGIN = 1,FUN = mean,na.rm = TRUE)
        
        saveRDS(conn_mat,conn_mat.name)
        
        cat(paste("+ DONE: Subject = ",subj,", R1 = ",r1,", R2 = ",r2,sep = ""),"\r")
        
      }  
    }
  }
  
  symm_fband <- conn_mat
  
  ### Testing by pair
  nboot <- 10^5
  pval_pair <- array(NA,dim = c(length(ROI),length(ROI),5,3))
  
  for(r1 in 1:6){
    for(r2 in (r1+1):7){
      for(f in 1:5){
        group <- c(rep("delta_A",36),rep("delta_C",29),rep("delta_F",23))
        data_fband <- data.frame(group = group,x = symm_fband[r1,r2,f,])
        pval_fband_ac <- aovperm(x ~ group, data = data_fband[1:65,], np = nboot)$table[1,5]
        pval_fband_af <- aovperm(x ~ group, data = data_fband[c(1:36,66:88),], np = nboot)$table[1,5]
        pval_fband_cf <- aovperm(x ~ group, data = data_fband[37:88,], np = nboot)$table[1,5]
        pval_pair[r1,r2,f,1] <- pval_fband_ac
        pval_pair[r1,r2,f,2] <- pval_fband_af
        pval_pair[r1,r2,f,3] <- pval_fband_cf
        
      }
      saveRDS(pval_pair,pval.name)
      cat(paste("+ DONE: R1 = ",r1,", R2 = ",r2,sep = ""),"\r")
    }
  }
  
  
  pval_all <- NULL
  for(f in 1:5){
    for(r1 in 1:6){
      for(r2 in (r1+1):7){
        for(p in 1:3){
          pval_all <- rbind(pval_all,data.frame(r1 = r1,r2 = r2,f = f,p = p,pval = pval_pair[r1,r2,f,p]))
        }
      }
    }
  }
  
  pval_all$pval_adj <- p.adjust(pval_all$pval,method = "BH")
  sum(pval_all$pval_adj < 0.1,na.rm = TRUE)
  pval_pair_adj <- array(NA,dim = c(length(ROI),length(ROI),5,3))
  for(freq in 1:5){
    for(roi1 in 1:6){
      for(roi2 in (roi1+1):7){
        for(pair in 1:3){
          pval_pair_adj[roi1,roi2,freq,pair] <- subset(pval_all,r1 == roi1 & r2 == roi2 & f == freq & p == pair)$pval_adj
          pval_pair_adj[roi2,roi1,freq,pair] <- subset(pval_all,r1 == roi1 & r2 == roi2 & f == freq & p == pair)$pval_adj
        }
      }
    }
  }
  saveRDS(pval_pair_adj,pval_adj.name)
}


set.seed(2) # seed for reproducibility
generate.results.PBC(samprate = 500,scale.overall = FALSE,
                     conn_mat.name = "pbc_conn_mat_s500_scf.rds",
                     pval.name = "pbc_pval_s500_scf.rds",
                     pval_adj.name = "pbc_pval_adj_s500_scf.rds")




############################################
##### Generate results based on KenCoh #####
############################################

generate.results.KenCoh <- function(samprate,scale.overall = TRUE,
                                    conn_mat.name,pval.name,pval_adj.name,
                                    fband = c("delta","theta","alpha","beta","gamma")){
  
  conn_mat <- array(NA,dim = c(length(ROI),length(ROI),5,length(subj_all)))
  
  for(subj in 1:length(subj_all)){
    data <- readMat(paste0("Dataset/derivatives/",subj_all[subj],"/eeg/",subj_all[subj],"_task-eyesclosed_eeg.set"))
    eeg_mat <- t(data$data)
    if(scale.overall){
      eeg_mat <- scale(eeg_mat[5*500 + 1:((floor(nsec_all[subj])-10)*500),])
    } else {
      eeg_mat <- eeg_mat[5*500 + 1:((floor(nsec_all[subj])-10)*500),]
    }
    
    if(samprate == 100){
      start_subj <- readRDS("start_subj_s100.rds")
      eeg_mat <- eeg_mat[start_subj[subj] + (0:((floor(nsec_all[subj])-10)*100-1))*5,]
    }
    
    if(samprate == 250){
      start_subj <- readRDS("start_subj_s250.rds")
      eeg_mat <- eeg_mat[start_subj[subj] + (0:((floor(nsec_all[subj])-10)*250-1))*2,]
    }
    
    for(r1 in 1:length(ROI)){
      for(r2 in 1:length(ROI)){
        if(r1 >= r2){next}
        
        elec_r1 <- (1:19)[chnames %in% ROI[[r1]]]
        elec_r2 <- (1:19)[chnames %in% ROI[[r2]]]
        
        y_mat <- eeg_mat[,elec_r1]
        x_mat <- eeg_mat[,elec_r2]
        
        conn_mat[r1,r2,,subj] <- conn_mat[r2,r1,,subj] <- kencoh(x_mat,y_mat,samprate = samprate,fband = fband)
        
        saveRDS(conn_mat,conn_mat.name)
        
        cat(paste("+ DONE: Subject = ",subj,", R1 = ",r1,", R2 = ",r2,sep = ""),"\r")
        
      }  
    }
  }
  
  symm_fband <- conn_mat
  
  ### Testing by pair
  nboot <- 10^5
  pval_pair <- array(NA,dim = c(length(ROI),length(ROI),5,3))
  
  for(r1 in 1:6){
    for(r2 in (r1+1):7){
      for(f in 1:5){
        group <- c(rep("delta_A",36),rep("delta_C",29),rep("delta_F",23))
        data_fband <- data.frame(group = group,x = symm_fband[r1,r2,f,])
        pval_fband_ac <- aovperm(x ~ group, data = data_fband[1:65,], np = nboot)$table[1,5]
        pval_fband_af <- aovperm(x ~ group, data = data_fband[c(1:36,66:88),], np = nboot)$table[1,5]
        pval_fband_cf <- aovperm(x ~ group, data = data_fband[37:88,], np = nboot)$table[1,5]
        pval_pair[r1,r2,f,1] <- pval_fband_ac
        pval_pair[r1,r2,f,2] <- pval_fband_af
        pval_pair[r1,r2,f,3] <- pval_fband_cf
        
      }
      saveRDS(pval_pair,pval.name)
      cat(paste("+ DONE: R1 = ",r1,", R2 = ",r2,sep = ""),"\r")
    }
  }
  
  
  pval_all <- NULL
  for(f in 1:5){
    for(r1 in 1:6){
      for(r2 in (r1+1):7){
        for(p in 1:3){
          pval_all <- rbind(pval_all,data.frame(r1 = r1,r2 = r2,f = f,p = p,pval = pval_pair[r1,r2,f,p]))
        }
      }
    }
  }
  
  pval_all$pval_adj <- p.adjust(pval_all$pval,method = "BH")
  sum(pval_all$pval_adj < 0.1,na.rm = TRUE)
  pval_pair_adj <- array(NA,dim = c(length(ROI),length(ROI),5,3))
  for(freq in 1:5){
    for(roi1 in 1:6){
      for(roi2 in (roi1+1):7){
        for(pair in 1:3){
          pval_pair_adj[roi1,roi2,freq,pair] <- subset(pval_all,r1 == roi1 & r2 == roi2 & f == freq & p == pair)$pval_adj
          pval_pair_adj[roi2,roi1,freq,pair] <- subset(pval_all,r1 == roi1 & r2 == roi2 & f == freq & p == pair)$pval_adj
        }
      }
    }
  }
  saveRDS(pval_pair_adj,pval_adj.name)
}

set.seed(3) # seed for reproducibility
generate.results.KenCoh(samprate = 500,scale.overall = FALSE,
                        conn_mat.name = "kencoh_conn_mat_s500_scf.rds",
                        pval.name = "kencoh_pval_s500_scf.rds",
                        pval_adj.name = "kencoh_pval_adj_s500_scf.rds")


###########################################################
##### Generate results based on nonspectral T measure #####
###########################################################

generate.results.T <- function(samprate,scale.overall = TRUE,
                               conn_mat.name,pval.name,pval_adj.name){
  
  conn_mat <- array(NA,dim = c(length(ROI),length(ROI),length(subj_all)))
  for(subj in 1:length(subj_all)){
    data <- readMat(paste0("Dataset/derivatives/",subj_all[subj],"/eeg/",subj_all[subj],"_task-eyesclosed_eeg.set"))
    eeg_mat <- t(data$data)
    
    if(scale.overall){
      eeg_mat <- scale(eeg_mat[5*500 + 1:((floor(nsec_all[subj])-10)*500),])
    } else {
      eeg_mat <- eeg_mat[5*500 + 1:((floor(nsec_all[subj])-10)*500),]
    }
    
    if(samprate == 100){
      start_subj <- readRDS("start_subj_s100.rds")
      eeg_mat <- eeg_mat[start_subj[subj] + (0:((floor(nsec_all[subj])-10)*100-1))*5,]
    }
    
    if(samprate == 250){
      start_subj <- readRDS("start_subj_s250.rds")
      eeg_mat <- eeg_mat[start_subj[subj] + (0:((floor(nsec_all[subj])-10)*250-1))*2,]
    }
    
    for(r1 in 1:length(ROI)){
      for(r2 in 1:length(ROI)){
        if(r1 == r2){next}
        
        elec_r1 <- (1:19)[chnames %in% ROI[[r1]]]
        elec_r2 <- (1:19)[chnames %in% ROI[[r2]]]
        
        y_mat <- eeg_mat[,elec_r1]
        x_mat <- eeg_mat[,elec_r2]
        
        est <- xi.vec(x = x_mat, y = y_mat,perm.inv = TRUE,nperm = NULL,nboot = 0)
        
        conn_mat[r1,r2,subj] <- est
        
        saveRDS(conn_mat,conn_mat.name)
        
        cat(paste("+ DONE: Subject = ",subj,", R1 = ",r1,", R2 = ",r2,sep = ""),"\r")
        
      }  
    }
  }
  
  symm_mat <- array(NA,dim = dim(conn_mat))
  for(subj in 1:length(subj_all)){
    for(r1 in 1:length(ROI)){
      for(r2 in 1:length(ROI)){
        if(r1 == r2){next}
        symm_mat[r1,r2,subj] <- ifelse(conn_mat[r1,r2,subj] > conn_mat[r2,r1,subj],conn_mat[r1,r2,subj],conn_mat[r2,r1,subj])
      }
    }
  }
  
  ### Testing by pair
  nboot <- 10^5
  pval_pair <- array(NA,dim = c(length(ROI),length(ROI),3))
  
  for(r1 in 1:6){
    for(r2 in (r1+1):7){
      group <- c(rep("delta_A",36),rep("delta_C",29),rep("delta_F",23))
      data_mat <- data.frame(group = group,x = symm_mat[r1,r2,])
      pval_ac <- aovperm(x ~ group, data = data_mat[1:65,], np = nboot)$table[1,5]
      pval_af <- aovperm(x ~ group, data = data_mat[c(1:36,66:88),], np = nboot)$table[1,5]
      pval_cf <- aovperm(x ~ group, data = data_mat[37:88,], np = nboot)$table[1,5]
      pval_pair[r1,r2,1] <- pval_ac
      pval_pair[r1,r2,2] <- pval_af
      pval_pair[r1,r2,3] <- pval_cf
      saveRDS(pval_pair,pval.name)
      cat(paste("+ DONE: R1 = ",r1,", R2 = ",r2,sep = ""),"\r")
    }
  }
  
  
  pval_all <- NULL
  for(r1 in 1:6){
    for(r2 in (r1+1):7){
      for(p in 1:3){
        pval_all <- rbind(pval_all,data.frame(r1 = r1,r2 = r2,p = p,pval = pval_pair[r1,r2,p]))
      }
    }
  }
  
  pval_all$pval_adj <- p.adjust(pval_all$pval,method = "BH")
  
  pval_pair_adj <- array(NA,dim = c(length(ROI),length(ROI),3))
  for(roi1 in 1:6){
    for(roi2 in (roi1+1):7){
      for(pair in 1:3){
        pval_pair_adj[roi1,roi2,pair] <- subset(pval_all,r1 == roi1 & r2 == roi2 & p == pair)$pval_adj
        pval_pair_adj[roi2,roi1,pair] <- subset(pval_all,r1 == roi1 & r2 == roi2 & p == pair)$pval_adj
      }
    }
  }
  saveRDS(pval_pair_adj,pval_adj.name)
  
}

set.seed(4) # seed for reproducibility
generate.results.T(samprate = 500,scale.overall = FALSE,
                   conn_mat.name = "conn_mat_nf.rds",
                   pval.name = "pval_nf.rds",
                   pval_adj.name = "pval_adj_nf.rds")


##########################################################
##### Generate results for relative band power (RBP) #####
##########################################################

rel.power <- function(y,blength = 500,samprate = 500){
  
  Tlength <- nrow(y)
  bnum <- Tlength/blength
  
  freq_all <- seq(0,0.5,length = samprate/2+1)[-1]
  freq_all <- freq_all[freq_all*samprate<=50]
  
  ind_delta <- (1:length(freq_all))[freq_all*samprate > 0.5 & freq_all*samprate <= 4]
  ind_theta <- (1:length(freq_all))[freq_all*samprate > 4 & freq_all*samprate <= 8]
  ind_alpha <- (1:length(freq_all))[freq_all*samprate > 8 & freq_all*samprate <= 12]
  ind_beta <- (1:length(freq_all))[freq_all*samprate > 12 & freq_all*samprate <= 30]
  ind_gamma <- (1:length(freq_all))[freq_all*samprate > 30 & freq_all*samprate <= 45]
  
  y.power <- array(NA,dim = c(bnum,ncol(y),5))
  for(ch in 1:ncol(y)){
    for(b in 1:bnum){
      perio <- LSTS::periodogram(y[blength*(b-1) + 1:blength,ch],plot = FALSE)[[1]]
      y.power[b,ch,1] <- sum(perio[ind_delta])/sum(perio[1:length(freq_all)])
      y.power[b,ch,2] <- sum(perio[ind_theta])/sum(perio[1:length(freq_all)])
      y.power[b,ch,3] <- sum(perio[ind_alpha])/sum(perio[1:length(freq_all)])
      y.power[b,ch,4] <- sum(perio[ind_beta])/sum(perio[1:length(freq_all)])
      y.power[b,ch,5] <- sum(perio[ind_gamma])/sum(perio[1:length(freq_all)])
    }
  }
  return(apply(y.power,MARGIN = 3,FUN = mean,na.rm = TRUE))
}



generate.results.relp <- function(samprate,blength,scale.overall = TRUE,
                                  relp_mat.name,pval_relp.name,pval_relp_adj.name){
  
  freq_all <- seq(0,0.5,length = blength/2+1)[-1]
  freq_all <- freq_all[freq_all*samprate<=50]
  
  relp_mat <- array(NA,dim = c(length(ROI),5,length(subj_all)))
  
  for(subj in 1:length(subj_all)){
    data <- readMat(paste0("Dataset/derivatives/",subj_all[subj],"/eeg/",subj_all[subj],"_task-eyesclosed_eeg.set"))
    eeg_mat <- t(data$data)
    if(scale.overall){
      eeg_mat <- scale(eeg_mat[5*500 + 1:((floor(nsec_all[subj])-10)*500),])
    } else {
      eeg_mat <- eeg_mat[5*500 + 1:((floor(nsec_all[subj])-10)*500),]
    }
    
    if(samprate == 100){
      start_subj <- readRDS("start_subj_s100.rds")
      eeg_mat <- eeg_mat[start_subj[subj] + (0:((floor(nsec_all[subj])-10)*100-1))*5,]
    }
    
    if(samprate == 250){
      start_subj <- readRDS("start_subj_s250.rds")
      eeg_mat <- eeg_mat[start_subj[subj] + (0:((floor(nsec_all[subj])-10)*250-1))*2,]
    }
    
    for(r1 in 1:length(ROI)){
      
      elec_r1 <- (1:19)[chnames %in% ROI[[r1]]]
      
      y_mat <- eeg_mat[,elec_r1]
      
      relp_mat[r1,,subj] <- rel.power(y_mat,blength = blength,samprate = samprate)
      
      saveRDS(relp_mat,relp_mat.name)
      
      cat(paste("+ DONE: Subject = ",subj,", R1 = ",r1,sep = ""),"\r")
      
    }
    
  }
  
  ### Testing by pair
  nboot <- 10^5
  pval_relp <- array(NA,dim = c(length(ROI),5,3))
  
  for(r1 in 1:length(ROI)){
    for(f in 1:5){
      group <- c(rep("delta_A",36),rep("delta_C",29),rep("delta_F",23))
      data_fband <- data.frame(group = group,x = relp_mat[r1,f,])
      pval_fband_ac <- aovperm(x ~ group, data = data_fband[1:65,], np = nboot)$table[1,5]
      pval_fband_af <- aovperm(x ~ group, data = data_fband[c(1:36,66:88),], np = nboot)$table[1,5]
      pval_fband_cf <- aovperm(x ~ group, data = data_fband[37:88,], np = nboot)$table[1,5]
      pval_relp[r1,f,1] <- pval_fband_ac
      pval_relp[r1,f,2] <- pval_fband_af
      pval_relp[r1,f,3] <- pval_fband_cf
      
    }
    saveRDS(pval_relp,pval_relp.name)
    cat(paste("+ DONE: R1 = ",r1,sep = ""),"\r")
  }
  
  
  pval_all <- NULL
  for(pr in 1:3){
    for(f in 1:5){
      for(r1 in 1:7){
        pval_all <- rbind(pval_all,data.frame(r1 = r1,f = f,pr = pr,pval = pval_relp[r1,f,pr]))
      }
    }
  }
  pval_all$pval_adj <- p.adjust(pval_all$pval,method = "BH")
  
  pval_relp_adj <- array(NA,dim = c(length(ROI),5,3))
  for(pr2 in 1:3){
    for(freq in 1:5){
      for(roi1 in 1:7){
        pval_relp_adj[roi1,freq,pr2] <- subset(pval_all,r1 == roi1 & f == freq & pr == pr2)$pval_adj
      }
    }
  }
  saveRDS(pval_relp_adj,pval_relp_adj.name)
  
}

set.seed(5) # seed for reproducibility
generate.results.relp(samprate = 500,blength = 500,scale.overall = FALSE,
                      relp_mat.name = "relp_mat_s500_scf.rds",
                      pval_relp.name = "pval_relp_s500_scf.rds",
                      pval_relp_adj.name = "pval_relp_adj_s500_scf.rds")


### Summarize results
relp_mat <- readRDS("relp_mat_s500_scf.rds")
pval_relp <- readRDS("pval_relp_s500_scf.rds")
pval_relp_adj <- readRDS("pval_relp_adj_s500_scf.rds")

pval_all <- NULL
for(f in 1:5){
  for(r1 in 1:length(ROI)){
    
    pval_all <- rbind(pval_all,
                      data.frame(r1 = r1,f = f,p = 1,
                                 pval = pval_relp[r1,f,1],
                                 pval_adj = pval_relp_adj[r1,f,1]))
    
    pval_all <- rbind(pval_all,
                      data.frame(r1 = r1,f = f,p = 2,
                                 pval = pval_relp[r1,f,2],
                                 pval_adj = pval_relp_adj[r1,f,2]))
    
    pval_all <- rbind(pval_all,
                      data.frame(r1 = r1,f = f,p = 3,
                                 pval = pval_relp[r1,f,3],
                                 pval_adj = pval_relp_adj[r1,f,3]))
  }
}

band_names <- c("delta", "theta", "alpha", "beta", "gamma")
pval_all$Comparison <- ifelse(pval_all$p == 1,"CN vs. AD",
                              ifelse(pval_all$p == 2,"AD vs. FTD","CN vs. FTD"))
pval_all$ROI <- names(ROI)[pval_all$r1]
pval_all$'Freq. Band' <- band_names[pval_all$f]

comp <- pval_all[,c(6:8,4:5)] 
comp# full table of significance of average RBP comparisons

comp1 <- subset(comp,Comparison == "CN vs. AD")
comp1[order(comp1$pval_adj),] # significant RBP for CN vs. AD.

comp2 <- subset(comp,Comparison == "CN vs. FTD")
comp2[order(comp2$pval_adj),] # significant RBP for CN vs. FTD.

#################################################
##### Visualization of connectivity results #####
#################################################

band_names <- c("delta", "theta", "alpha", "beta", "gamma")
generate.plots <- function(type = "nvc",samprate,sig.level = 0.05,
                           conn_mat.name,pval_adj.name,
                           plot.name_ad,plot.title_ad,
                           plot.name_ftd,plot.title_ftd){
  
  conn_mat <- readRDS(conn_mat.name)
  
  if(type == "nvc"){
    symm_mat <- array(NA,dim = dim(conn_mat))
    for(subj in 1:length(subj_all)){
      for(r1 in 1:length(ROI)){
        for(r2 in 1:length(ROI)){
          if(r1 == r2){next}
          symm_mat[r1,r2,,subj] <- ifelse(conn_mat[r1,r2,,subj] > conn_mat[r2,r1,,subj],conn_mat[r1,r2,,subj],conn_mat[r2,r1,,subj])
        }
      }
    }
    
    freq_all <- seq(0,0.5,length = samprate/2+1)[-1]
    freq_all <- freq_all[freq_all*samprate<=50]
    
    symm_delta <- apply(symm_mat[,,1:4,],MARGIN = c(1,2,4),FUN = mean,na.rm = TRUE)
    symm_theta <- apply(symm_mat[,,5:8,],MARGIN = c(1,2,4),FUN = mean,na.rm = TRUE)
    symm_alpha <- apply(symm_mat[,,9:12,],MARGIN = c(1,2,4),FUN = mean,na.rm = TRUE)
    symm_beta <- apply(symm_mat[,,13:30,],MARGIN = c(1,2,4),FUN = mean,na.rm = TRUE)
    symm_gamma <- apply(symm_mat[,,31:45,],MARGIN = c(1,2,4),FUN = mean,na.rm = TRUE)
    
    symm_conn <- array(NA,dim = c(7,7,5,88))
    symm_conn[,,1,] <- symm_delta
    symm_conn[,,2,] <- symm_theta
    symm_conn[,,3,] <- symm_alpha
    symm_conn[,,4,] <- symm_beta
    symm_conn[,,5,] <- symm_gamma
  } else if(type %in% c("pbc","kencoh")){
    symm_conn <- conn_mat
  }
  
  symm_agg <- array(NA,dim = c(7,7,5,3))
  symm_agg[,,,1] <- apply(symm_conn[,,,1:36],MARGIN = c(1,2,3),FUN = median,na.rm = TRUE)
  symm_agg[,,,2] <- apply(symm_conn[,,,37:65],MARGIN = c(1,2,3),FUN = median,na.rm = TRUE)
  symm_agg[,,,3] <- apply(symm_conn[,,,66:88],MARGIN = c(1,2,3),FUN = median,na.rm = TRUE)
  symm_agg <- symm_agg[c(2,4,5,7,6,3,1),c(2,4,5,7,6,3,1),,]
  
  pval_pair_adj <- readRDS(pval_adj.name)
  pval_ac <- pval_pair_adj[c(2,4,5,7,6,3,1),c(2,4,5,7,6,3,1),,1] < sig.level
  pval_cf <- pval_pair_adj[c(2,4,5,7,6,3,1),c(2,4,5,7,6,3,1),,3] < sig.level
  
  
  edges_ac <- NULL
  for(ch1 in 1:6){
    for(ch2 in (ch1+1):7){
      for(f in 1:5){
        if(pval_ac[ch1,ch2,f]){
          entry <- data.frame(from = paste0("N",(ch1-1)*5+f),
                              to = paste0("N",(ch2-1)*5+f),
                              band = band_names[f],
                              weight = 7,
                              line_type = ifelse(symm_agg[ch1,ch2,f,1] > symm_agg[ch1,ch2,f,2],1,3))
          edges_ac <- rbind(edges_ac,entry)
        }
      }
    }
  }
  
  edges_ac$band <- factor(edges_ac$band,levels = band_names,ordered = TRUE)
  edges_ac <- edges_ac[order(edges_ac$band),]
  edges_ac$band <- as.character(edges_ac$band)
  # edges_ac <- edges_ac[1:3,]
  
  
  
  edges_cf <- NULL
  for(ch1 in 1:6){
    for(ch2 in (ch1+1):7){
      for(f in 1:5){
        if(pval_cf[ch1,ch2,f]){
          entry <- data.frame(from = paste0("N",(ch1-1)*5+f),
                              to = paste0("N",(ch2-1)*5+f),
                              band = band_names[f],
                              weight = 7,
                              line_type = ifelse(symm_agg[ch1,ch2,f,3] > symm_agg[ch1,ch2,f,2],1,3))
          edges_cf <- rbind(edges_cf,entry)
        }
      }
    }
  }
  
  
  edges_cf$band <- factor(edges_cf$band,levels = band_names,ordered = TRUE)
  edges_cf <- edges_cf[order(edges_cf$band),]
  edges_cf$band <- as.character(edges_cf$band)
  
  generate.connect(edges_input = edges_ac,
                   plot.name = plot.name_ad,
                   plot.title = plot.title_ad)
  
  generate.connect(edges_input = edges_cf,
                   plot.name = plot.name_ftd,
                   plot.title = plot.title_ftd)
  
}

##### Main results
generate.plots(type = "nvc",samprate = 500,
               conn_mat.name = "conn_mat_s500_scf.rds",
               pval_adj.name = "pval_adj_s500_scf.rds",
               plot.name_ad = "NVC_AD.png",
               plot.title_ad = "NVC (CN vs. AD): 500 Hz, Unstandardized",
               plot.name_ftd = "NVC_FTD.png",
               plot.title_ftd = "NVC (CN vs. FTD): 500 Hz, Unstandardized")


generate.plots(type = "pbc",samprate = 500,
               conn_mat.name = "pbc_conn_mat_s500_scf.rds",
               pval_adj.name = "pbc_pval_adj_s500_scf.rds",
               plot.name_ad = "PBC_AD.png",
               plot.title_ad = "PBC (CN vs. AD): 500 Hz, Unstandardized",
               plot.name_ftd = "PBC_FTD.png",
               plot.title_ftd = "PBC (CN vs. FTD): 500 Hz, Unstandardized")


generate.plots(type = "kencoh",samprate = 500,
               conn_mat.name = "kencoh_conn_mat_s500_scf.rds",
               pval_adj.name = "kencoh_pval_adj_s500_scf.rds",
               plot.name_ad = "KenCoh_AD.png",
               plot.title_ad = "KenCoh (CN vs. AD): 500 Hz, Unstandardized",
               plot.name_ftd = "KenCoh_FTD.png",
               plot.title_ftd = "KenCoh (CN vs. FTD): 500 Hz, Unstandardized")

##### Sensitivity analysis results
generate.plots(type = "nvc",samprate = 100,
               conn_mat.name = "conn_mat_s100_scf.rds",
               pval_adj.name = "pval_adj_s100_scf.rds",
               plot.name_ad = "NVC_AD_s100_scf.png",
               plot.title_ad = "NVC (CN vs. AD): 100 Hz, Unstandardized",
               plot.name_ftd = "NVC_FTD_s100_scf.png",
               plot.title_ftd = "NVC (CN vs. FTD): 100 Hz, Unstandardized")

generate.plots(type = "nvc",samprate = 100,
               conn_mat.name = "conn_mat_s100_sct.rds",
               pval_adj.name = "pval_adj_s100_sct.rds",
               plot.name_ad = "NVC_AD_s100_sct.png",
               plot.title_ad = "NVC (CN vs. AD): 100 Hz, Standardized",
               plot.name_ftd = "NVC_FTD_s100_sct.png",
               plot.title_ftd = "NVC (CN vs. FTD): 100 Hz, Standardized")


generate.plots(type = "nvc",samprate = 250,
               conn_mat.name = "conn_mat_s250_scf.rds",
               pval_adj.name = "pval_adj_s250_scf.rds",
               plot.name_ad = "NVC_AD_s250_scf.png",
               plot.title_ad = "NVC (CN vs. AD): 250 Hz, Unstandardized",
               plot.name_ftd = "NVC_FTD_s250_scf.png",
               plot.title_ftd = "NVC (CN vs. FTD): 250 Hz, Unstandardized")

generate.plots(type = "nvc",samprate = 250,
               conn_mat.name = "conn_mat_s250_sct.rds",
               pval_adj.name = "pval_adj_s250_sct.rds",
               plot.name_ad = "NVC_AD_s250_sct.png",
               plot.title_ad = "NVC (CN vs. AD): 250 Hz, Standardized",
               plot.name_ftd = "NVC_FTD_s250_sct.png",
               plot.title_ftd = "NVC (CN vs. FTD): 250 Hz, Standardized")


generate.plots(type = "nvc",samprate = 500,
               conn_mat.name = "conn_mat_s500_sct.rds",
               pval_adj.name = "pval_adj_s500_sct.rds",
               plot.name_ad = "NVC_AD_s500_sct.png",
               plot.title_ad = "NVC (CN vs. AD): 500 Hz, Standardized",
               plot.name_ftd = "NVC_FTD_s500_sct.png",
               plot.title_ftd = "NVC (CN vs. FTD): 500 Hz, Standardized")



generate.plots(type = "nvc",samprate = 500,
               conn_mat.name = "conn_mat_s500_b750_scf.rds",
               pval_adj.name = "pval_adj_s500_b750_scf.rds",
               plot.name_ad = "NVC_AD_s500_b750_scf.png",
               plot.title_ad = "NVC (CN vs. AD): 500 Hz, B = 750, Unstandardized",
               plot.name_ftd = "NVC_FTD_s500_b750_scf.png",
               plot.title_ftd = "NVC (CN vs. FTD): 500 Hz, B = 750, Unstandardized")


generate.plots(type = "nvc",samprate = 500,
               conn_mat.name = "conn_mat_s500_b1000_scf.rds",
               pval_adj.name = "pval_adj_s500_b1000_scf.rds",
               plot.name_ad = "NVC_AD_s500_b1000_scf.png",
               plot.title_ad = "NVC (CN vs. AD): 500 Hz, B = 1000, Unstandardized",
               plot.name_ftd = "NVC_FTD_s500_b1000_scf.png",
               plot.title_ftd = "NVC (CN vs. FTD): 500 Hz, B = 1000, Unstandardized")

############################################
##### Computing classification metrics #####
############################################

compute.metrics <- function(pred_df, model_name,gcomp){
  
  cm <- caret::confusionMatrix(
    factor(pred_df$predicted, levels = gcomp),
    factor(pred_df$actual, levels = gcomp),
    positive = setdiff(gcomp, "CN")
  )
  
  acc  <- cm$overall["Accuracy"]
  sens <- cm$byClass["Sensitivity"]   # TP / (TP + FN)
  spec <- cm$byClass["Specificity"]   # TN / (TN + FP)
  f1   <- cm$byClass["F1"]            # Harmonic mean of precision and recall
  
  metrics <- data.frame(Model = model_name,
                        Accuracy = round(acc, 3),
                        Sensitivity = round(sens, 3),
                        Specificity = round(spec, 3),
                        F1_Score = round(f1, 3))
  return(metrics)
}


generate.metrics <- function(type,samprate = 500,blength = 500,scale.true = "f",
                             conn.num = 10,relp.num = 5,pr = 1){
  
  if(samprate == blength){
    setting <- paste0("s",samprate,"_sc",scale.true)
  } else {
    setting <- paste0("s",samprate,"_b",blength,"_sc",scale.true)
  }
  
  conn_mat <- readRDS(paste0("",type,"conn_mat_",setting,".rds"))
  pval_pair_adj <- readRDS(paste0("",type,"pval_adj_",setting,".rds"))
  relp_mat <- readRDS(paste0("relp_mat_",setting,".rds"))
  pval_relp_adj <- readRDS(paste0("pval_relp_adj_",setting,".rds"))
  
  if(type %in% c("pbc_","kencoh_")){
    
    results.name <- paste0("metrics_",type,"c",conn.num,"_r",relp.num,"_pr",pr,".rds")
    symm_conn <- conn_mat
    
  } else {
    
    results.name <- paste0("metrics_nvc_c",conn.num,"_r",relp.num,"_pr",pr,".rds")
    subj_all <- 1:88
    symm_mat <- array(NA,dim = dim(conn_mat))
    for(subj in 1:length(subj_all)){
      for(r1 in 1:length(ROI)){
        for(r2 in 1:length(ROI)){
          if(r1 == r2){next}
          symm_mat[r1,r2,,subj] <- ifelse(conn_mat[r1,r2,,subj] > conn_mat[r2,r1,,subj],conn_mat[r1,r2,,subj],conn_mat[r2,r1,,subj])
        }
      }
    }
    
    freq_all <- seq(0,0.5,length = blength/2+1)[-1]
    freq_all <- freq_all[freq_all*samprate<=50]
    
    ind_delta <- (1:length(freq_all))[freq_all*samprate > 0.5 & freq_all*samprate <= 4]
    ind_theta <- (1:length(freq_all))[freq_all*samprate > 4 & freq_all*samprate <= 8]
    ind_alpha <- (1:length(freq_all))[freq_all*samprate > 8 & freq_all*samprate <= 12]
    ind_beta <- (1:length(freq_all))[freq_all*samprate > 12 & freq_all*samprate <= 30]
    ind_gamma <- (1:length(freq_all))[freq_all*samprate > 30 & freq_all*samprate <= 45]
    
    symm_conn <- array(NA,dim = c(length(ROI),length(ROI),5,88))
    symm_conn[,,1,] <- apply(symm_mat[,,ind_delta,],MARGIN = c(1,2,4),FUN = mean,na.rm = TRUE)
    symm_conn[,,2,] <- apply(symm_mat[,,ind_theta,],MARGIN = c(1,2,4),FUN = mean,na.rm = TRUE)
    symm_conn[,,3,] <- apply(symm_mat[,,ind_alpha,],MARGIN = c(1,2,4),FUN = mean,na.rm = TRUE)
    symm_conn[,,4,] <- apply(symm_mat[,,ind_beta,],MARGIN = c(1,2,4),FUN = mean,na.rm = TRUE)
    symm_conn[,,5,] <- apply(symm_mat[,,ind_gamma,],MARGIN = c(1,2,4),FUN = mean,na.rm = TRUE)
    
  }
  
  pval_all <- NULL
  for(f in 1:5){
    for(r1 in 1:(length(ROI)-1)){
      for(r2 in (r1+1):length(ROI)){
        
        pval_all <- rbind(pval_all,
                          data.frame(r1 = r1,r2 = r2,f = f,p = 1,
                                     pval = pval_pair_adj[r1,r2,f,1]))
        pval_all <- rbind(pval_all,
                          data.frame(r1 = r1,r2 = r2,f = f,p = 2,
                                     pval = pval_pair_adj[r1,r2,f,2]))
        
        pval_all <- rbind(pval_all,
                          data.frame(r1 = r1,r2 = r2,f = f,p = 3,
                                     pval = pval_pair_adj[r1,r2,f,3]))
        
      }
    }
  }
  if(type == "kencoh_" & pr == 3){
    pval_sub <- subset(pval_all, p == pr)
  } else {
    pval_sub <- subset(pval_all, p == pr & pval < 0.05)
  }
  conn_top <- pval_sub[order(pval_sub$pval),][1:conn.num,]
  
  f_data <- data.frame(subject_id = 1:88,label = c(rep("AD",36),rep("CN",29),rep("FTD",23)))
  if(conn.num>0){
    for(r.ind in 1:nrow(conn_top)){
      sub <- conn_top[r.ind,]
      entry <- data.frame(conn = symm_conn[sub$r1,sub$r2,sub$f,])
      names(entry) <- paste0("C",r.ind)
      f_data <- cbind(f_data,entry)
    }
  }
  
  pval_all <- NULL
  for(f in 1:5){
    for(r1 in 1:length(ROI)){
      
      pval_all <- rbind(pval_all,
                        data.frame(r1 = r1,f = f,p = 1,
                                   pval = pval_relp_adj[r1,f,1]))
      
      pval_all <- rbind(pval_all,
                        data.frame(r1 = r1,f = f,p = 2,
                                   pval = pval_relp_adj[r1,f,2]))
      
      pval_all <- rbind(pval_all,
                        data.frame(r1 = r1,f = f,p = 3,
                                   pval = pval_relp_adj[r1,f,3]))
    }
  }
  
  pval_sub <- subset(pval_all, p == pr & pval < 0.05)
  relp_top <- pval_sub[order(pval_sub$pval),][1:relp.num,]
  
  if(relp.num>0){
    for(r.ind in 1:nrow(relp_top)){
      sub <- relp_top[r.ind,]
      entry <- data.frame(relp = relp_mat[sub$r1,sub$f,])
      names(entry) <- paste0("R",r.ind)
      f_data <- cbind(f_data,entry)
    }
  }
  
  if(pr == 1){
    gcomp <- c("AD","CN")
  } else if(pr == 3){
    gcomp <- c("CN","FTD")
  }
  
  f_data <- subset(f_data,label %in% gcomp)
  f_data$label <- factor(f_data$label,levels = gcomp)
  
  # Storage for predictions
  results_svm <- data.frame(actual = character(), predicted = character())
  results_svm2 <- data.frame(actual = character(), predicted = character())
  results_svm3 <- data.frame(actual = character(), predicted = character())
  results_svm4 <- data.frame(actual = character(), predicted = character())
  results_knn <- data.frame(actual = character(), predicted = character())
  results_knn2 <- data.frame(actual = character(), predicted = character())
  results_rf  <- data.frame(actual = character(), predicted = character())
  results_rf2  <- vector("list",100)
  for(i in 1:100){
    results_rf2[[i]] <- data.frame(actual = character(), predicted = character())
  }
  
  
  # LOSO-CV
  for (subject in unique(f_data$subject_id)) {
    train_data <- subset(f_data, subject_id != subject)
    test_data  <- subset(f_data, subject_id == subject)
    
    train_x <- train_data[, -c(1,2)]
    train_y <- train_data$label
    test_x  <- test_data[, -c(1,2)]
    test_y  <- test_data$label
    
    model_svm <- svm(train_x, train_y, kernel = "linear")
    pred_svm  <- predict(model_svm, test_x)
    
    model_svm2 <- svm(train_x, train_y, kernel = "polynomial")
    pred_svm2  <- predict(model_svm2, test_x)
    
    model_svm3 <- svm(train_x, train_y, kernel = "radial")
    pred_svm3  <- predict(model_svm3, test_x)
    
    model_svm4 <- svm(train_x, train_y, kernel = "sigmoid")
    pred_svm4  <- predict(model_svm4, test_x)
    
    pred_knn  <- class::knn(train = train_x, test = test_x, cl = train_y, k = 5)
    pred_knn2  <- class::knn(train = train_x, test = test_x, cl = train_y, k = 10)
    
    model_rf  <- randomForest(train_x, train_y, ntree = 50)
    pred_rf   <- predict(model_rf, test_x)
    
    results_svm <- rbind(results_svm, data.frame(actual = test_y, predicted = pred_svm))
    results_svm2 <- rbind(results_svm2, data.frame(actual = test_y, predicted = pred_svm2))
    results_svm3 <- rbind(results_svm3, data.frame(actual = test_y, predicted = pred_svm3))
    results_svm4 <- rbind(results_svm4, data.frame(actual = test_y, predicted = pred_svm4))
    results_knn <- rbind(results_knn, data.frame(actual = test_y, predicted = pred_knn))
    results_knn2 <- rbind(results_knn2, data.frame(actual = test_y, predicted = pred_knn2))
    results_rf  <- rbind(results_rf,  data.frame(actual = test_y, predicted = pred_rf))
    
    for(i in 1:100){
      model_rf  <- randomForest(train_x, train_y, ntree = 50)
      pred_rf   <- predict(model_rf, test_x)
      results_rf2[[i]]  <- rbind(results_rf2[[i]],  data.frame(actual = test_y, predicted = pred_rf))
      cat(i,"\r")
    }
    
  }
  
  # Calculate metrics for each model
  metrics_svm <- compute.metrics(results_svm, "SVM (Linear)",gcomp)
  metrics_svm2 <- compute.metrics(results_svm2, "SVM (Polynomial)",gcomp)
  metrics_svm3 <- compute.metrics(results_svm3, "SVM (Radial)",gcomp)
  metrics_svm4 <- compute.metrics(results_svm4, "SVM (Sigmoid)",gcomp)
  metrics_knn <- compute.metrics(results_knn, "kNN (k = 5)",gcomp)
  metrics_knn2 <- compute.metrics(results_knn2, "kNN (k = 10)",gcomp)
  metrics_rf  <- compute.metrics(results_rf,  "Random Forest (Single)",gcomp)
  metrics_rf2 <- NULL
  for(i in 1:100){
    metrics_rf2 <- rbind(metrics_rf2,compute.metrics(results_rf2[[i]],  "Random Forest (Multiple)",gcomp))
  }
  metrics_rf2_summ <- aggregate(cbind(Accuracy,Sensitivity,Specificity,F1_Score)~Model,data = metrics_rf2,FUN = mean,na.rm = TRUE)
  metrics_rf2_summ[,-1] <- round(metrics_rf2_summ[-1],3)
  
  final_results <- bind_rows(metrics_svm,metrics_svm2,metrics_svm3,metrics_svm4,
                             metrics_knn,metrics_knn2,metrics_rf,metrics_rf2_summ)
  row.names(final_results) <- NULL
  
  saveRDS(final_results,results.name)
}

##### Generating classification performance metrics
set.seed(1) # seed for reproducibility

### PBC metrics
generate.metrics(type = "pbc_",samprate = 500,blength = 500,scale.true = "f",
                 conn.num = 5,relp.num = 0,pr = 1)
generate.metrics(type = "pbc_",samprate = 500,blength = 500,scale.true = "f",
                 conn.num = 10,relp.num = 0,pr = 1)
generate.metrics(type = "pbc_",samprate = 500,blength = 500,scale.true = "f",
                 conn.num = 5,relp.num = 5,pr = 1)
generate.metrics(type = "pbc_",samprate = 500,blength = 500,scale.true = "f",
                 conn.num = 10,relp.num = 5,pr = 1)
generate.metrics(type = "pbc_",samprate = 500,blength = 500,scale.true = "f",
                 conn.num = 5,relp.num = 0,pr = 3)
generate.metrics(type = "pbc_",samprate = 500,blength = 500,scale.true = "f",
                 conn.num = 5,relp.num = 5,pr = 3)

### KenCoh metrics
generate.metrics(type = "kencoh_",samprate = 500,blength = 500,scale.true = "f",
                 conn.num = 5,relp.num = 0,pr = 1)
generate.metrics(type = "kencoh_",samprate = 500,blength = 500,scale.true = "f",
                 conn.num = 10,relp.num = 0,pr = 1)
generate.metrics(type = "kencoh_",samprate = 500,blength = 500,scale.true = "f",
                 conn.num = 5,relp.num = 5,pr = 1)
generate.metrics(type = "kencoh_",samprate = 500,blength = 500,scale.true = "f",
                 conn.num = 10,relp.num = 5,pr = 1)
generate.metrics(type = "kencoh_",samprate = 500,blength = 500,scale.true = "f",
                 conn.num = 5,relp.num = 0,pr = 3)
generate.metrics(type = "kencoh_",samprate = 500,blength = 500,scale.true = "f",
                 conn.num = 5,relp.num = 5,pr = 3)

### NVC metrics
generate.metrics(type = "",samprate = 500,blength = 500,scale.true = "f",
                 conn.num = 5,relp.num = 0,pr = 1)
generate.metrics(type = "",samprate = 500,blength = 500,scale.true = "f",
                 conn.num = 10,relp.num = 0,pr = 1)
generate.metrics(type = "",samprate = 500,blength = 500,scale.true = "f",
                 conn.num = 5,relp.num = 5,pr = 1)
generate.metrics(type = "",samprate = 500,blength = 500,scale.true = "f",
                 conn.num = 10,relp.num = 5,pr = 1)
generate.metrics(type = "",samprate = 500,blength = 500,scale.true = "f",
                 conn.num = 5,relp.num = 0,pr = 3)
generate.metrics(type = "",samprate = 500,blength = 500,scale.true = "f",
                 conn.num = 5,relp.num = 5,pr = 3)

##### Metrics
metrics_nvc_c5_r0_pr1 <- readRDS("metrics_nvc_c5_r0_pr1.rds")
metrics_nvc_c10_r0_pr1 <- readRDS("metrics_nvc_c10_r0_pr1.rds")
metrics_nvc_c5_r5_pr1 <- readRDS("metrics_nvc_c5_r5_pr1.rds")
metrics_nvc_c10_r5_pr1 <- readRDS("metrics_nvc_c10_r5_pr1.rds")
metrics_nvc_c5_r0_pr3 <- readRDS("metrics_nvc_c5_r0_pr3.rds")
metrics_nvc_c5_r5_pr3 <- readRDS("metrics_nvc_c5_r5_pr3.rds")

metrics_pbc_c5_r0_pr1 <- readRDS("metrics_pbc_c5_r0_pr1.rds")
metrics_pbc_c10_r0_pr1 <- readRDS("metrics_pbc_c10_r0_pr1.rds")
metrics_pbc_c5_r5_pr1 <- readRDS("metrics_pbc_c5_r5_pr1.rds")
metrics_pbc_c10_r5_pr1 <- readRDS("metrics_pbc_c10_r5_pr1.rds")
metrics_pbc_c5_r0_pr3 <- readRDS("metrics_pbc_c5_r0_pr3.rds")
metrics_pbc_c5_r5_pr3 <- readRDS("metrics_pbc_c5_r5_pr3.rds")

metrics_kencoh_c5_r0_pr1 <- readRDS("metrics_kencoh_c5_r0_pr1.rds")
metrics_kencoh_c10_r0_pr1 <- readRDS("metrics_kencoh_c10_r0_pr1.rds")
metrics_kencoh_c5_r5_pr1 <- readRDS("metrics_kencoh_c5_r5_pr1.rds")
metrics_kencoh_c10_r5_pr1 <- readRDS("metrics_kencoh_c10_r5_pr1.rds")
metrics_kencoh_c5_r0_pr3 <- readRDS("metrics_kencoh_c5_r0_pr3.rds")
metrics_kencoh_c5_r5_pr3 <- readRDS("metrics_kencoh_c5_r5_pr3.rds")


##### SVM Results (Best classification across all four kernels)
types <- rep(c("NVC","PBC","KenCoh"),each = 4)
features <- rep(c("C5","C5M5","C10","C10M5"),3)
metrics.na <- data.frame(Accuracy = NA,Sensitivity = NA,Specificity = NA,F1_Score = NA)
metric.summ <- cbind(types,features,
                     rbind(cbind(rbind(metrics_nvc_c5_r0_pr1[1,-1],
                                       metrics_nvc_c5_r5_pr1[1,-1],
                                       metrics_nvc_c10_r0_pr1[3,-1],
                                       metrics_nvc_c10_r5_pr1[1,-1]),
                                 rbind(metrics_nvc_c5_r0_pr3[4,-1],
                                       metrics_nvc_c5_r5_pr3[3,-1],
                                       metrics.na,metrics.na)),
                           cbind(rbind(metrics_pbc_c5_r0_pr1[4,-1],
                                       metrics_pbc_c5_r5_pr1[1,-1],
                                       metrics_pbc_c10_r0_pr1[2,-1],
                                       metrics_pbc_c10_r5_pr1[2,-1]),
                                 rbind(metrics_pbc_c5_r0_pr3[4,-1],
                                       metrics_pbc_c5_r5_pr3[2,-1],
                                       metrics.na,metrics.na)),
                           cbind(rbind(metrics_kencoh_c5_r0_pr1[3,-1],
                                       metrics_kencoh_c5_r5_pr1[1,-1],
                                       metrics_kencoh_c10_r0_pr1[4,-1],
                                       metrics_kencoh_c10_r5_pr1[1,-1]),
                                 rbind(metrics_kencoh_c5_r0_pr3[4,-1],
                                       metrics_kencoh_c5_r5_pr3[2,-1],
                                       metrics.na,metrics.na))))

metric.summ # best classification performance via SVM with different kernels


##### Results per algorithm
### m.ind == 1 (SVM, linear)
### m.ind == 2 (SVM, polynomial)
### m.ind == 3 (SVM, radial)
### m.ind == 4 (SVM, sigmoid)
### m.ind == 5 (k-NN, k=5)
### m.ind == 6 (k-NN, k=10)
### m.ind == 7 (single random forest, N_trees = 50)
### m.ind == 8 (average of 100 random forest, N_trees = 50)
### run from 5 to 8
m.ind <- 5
metric.summ <- cbind(types,features,
                     rbind(cbind(rbind(metrics_nvc_c5_r0_pr1[m.ind,-1],
                                       metrics_nvc_c5_r5_pr1[m.ind,-1],
                                       metrics_nvc_c10_r0_pr1[m.ind,-1],
                                       metrics_nvc_c10_r5_pr1[m.ind,-1]),
                                 rbind(metrics_nvc_c5_r0_pr3[m.ind,-1],
                                       metrics_nvc_c5_r5_pr3[m.ind,-1],
                                       metrics.na,metrics.na)),
                           cbind(rbind(metrics_pbc_c5_r0_pr1[m.ind,-1],
                                       metrics_pbc_c5_r5_pr1[m.ind,-1],
                                       metrics_pbc_c10_r0_pr1[m.ind,-1],
                                       metrics_pbc_c10_r5_pr1[m.ind,-1]),
                                 rbind(metrics_pbc_c5_r0_pr3[m.ind,-1],
                                       metrics_pbc_c5_r5_pr3[m.ind,-1],
                                       metrics.na,metrics.na)),
                           cbind(rbind(metrics_kencoh_c5_r0_pr1[m.ind,-1],
                                       metrics_kencoh_c5_r5_pr1[m.ind,-1],
                                       metrics_kencoh_c10_r0_pr1[m.ind,-1],
                                       metrics_kencoh_c10_r5_pr1[m.ind,-1]),
                                 rbind(metrics_kencoh_c5_r0_pr3[m.ind,-1],
                                       metrics_kencoh_c5_r5_pr3[m.ind,-1],
                                       metrics.na,metrics.na))))

metric.summ # classification performance per method

