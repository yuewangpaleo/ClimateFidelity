####  The code is to calculate the climate fidelity of plant taxa 
####  during the past 18,000 years in North America by pollen record
####  by Yue Wang, March 2022
####  R Version 4.2.0 (Vigorous Calisthenics)

##  packages
library(ecospat) #ecospat v.3.2
library(viridis) #viridis v.0.6.2
library(abind) #abind v.1.4-5
library(MASS) #MASS v.7.3-55
library(raster)

##  working directory
WD <- "~/Dropbox (GaTech)/BioSci-McGuire/Users/Yue/Writing/ClimateFidelity/v3/ClimateFidelityDataPackage"
WD <- "working directory" 
setwd(WD)


##  1. Load data  ====
pollen <- read.csv("pollen.csv")

##  taxa
taxa_tree <- c("Pinus","Quercus","Picea","Betula",
               "Alnus","Tsuga","Cupress","Fagus",
               "Ulmus","Abies","Fraxinus","Salix")
taxa_herb <- c("Poaceae","Cyperaceae","Artemisia","ChenAm")
taxa_all <- c(taxa_tree,taxa_herb)
taxa_dis <- c("ChenAm","Cupress","Tsuga","Poaceae",
              "Cyperaceae","Salix","Artemisia","Picea",
              "Alnus","Betula","Pinus","Abies",
              "Ulmus","Quercus","Fagus","Fraxinus")


##  2. Calculate climate fidelity by niche overlap, similarity, and equivalency ====
##  sample size:
##  N = 20, 50, 100, deglacial samples, all
##  For the main paper, similarity & equivalency are calculated using background
##  climates that sum all samples from t(n) and t(n+1); whereas overlap is calculated
##  using background climates from across all time periods (to ensure that all
##  overlap values are scaled similarly)

##  N = all ----
sim_all <- list()
eqvl_all <- list()
overlap <- matrix(NA,ncol = length(taxa_all),nrow = 5)
colnames(overlap) <- taxa_all
rownames(overlap) <- c(0,seq(2,14,4))
sim_p <- overlap
eqvl_p <- overlap

for (i in 1:length(taxa_all)) {
  sim_i <- list()
  eqvl_i <- list()
  
  for (k in 1:5) {
    if (k==1) {
      pollen_k <- subset(pollen, Age <= 0)[,c(which(colnames(pollen)==taxa_all[i]),6,7)]
      pollen_kplus <- subset(pollen, Age <= 2000 & Age > 0)[,c(which(colnames(pollen)==taxa_all[i]),6,7)]
    }
    if (k==2) {
      pollen_k <- subset(pollen, Age <= 2000 & Age > 0)[,c(which(colnames(pollen)==taxa_all[i]),6,7)]
      pollen_kplus <- subset(pollen, Age <= (k+1-2)*4000+2000 & Age > (k+1-3)*4000+2000)[,c(which(colnames(pollen)==taxa_all[i]),6,7)]
    }
    if (k > 2) {
      pollen_k <- subset(pollen, Age <= (k-2)*4000+2000 & Age > (k-3)*4000+2000)[,c(which(colnames(pollen)==taxa_all[i]),6,7)]
      pollen_kplus <- subset(pollen, Age <= (k+1-2)*4000+2000 & Age > (k+1-3)*4000+2000)[,c(which(colnames(pollen)==taxa_all[i]),6,7)]
    }
    
    pollen_k_i <- pollen_k[pollen_k[,1]>0,]
    pollen_kplus_i <- pollen_kplus[pollen_kplus[,1]>0,]
    
    if (nrow(pollen_k_i) >= 5 & nrow(pollen_kplus_i) >= 5) {
      pollen_i_grid_k <- ecospat.grid.clim.dyn(glob = rbind(pollen_k[,2:3],pollen_kplus[,2:3]),
                                               glob1 = pollen_k[,2:3],
                                               sp = pollen_k_i[,2:3],
                                               R = 100, th.sp = 0.05, th.env = 0.05)
      pollen_i_grid_kplus <- ecospat.grid.clim.dyn(glob = rbind(pollen_k[,2:3],pollen_kplus[,2:3]),
                                                   glob1 = pollen_kplus[,2:3],
                                                   sp = pollen_kplus_i[,2:3],
                                                   R = 100, th.sp = 0.05, th.env = 0.05)
      
      similarity_k <- ecospat.niche.similarity.test(pollen_i_grid_kplus,pollen_i_grid_k,rep = 100)
      equivalency_k <- ecospat.niche.equivalency.test(pollen_i_grid_kplus,pollen_i_grid_k,rep = 100)
      overlap[k,i] <- similarity_k$obs$D
      sim_p[k,i] <- similarity_k$p.D
      eqvl_p[k,i] <- equivalency_k$p.D
      sim_i[[k]] <- similarity_k
      eqvl_i[[k]] <- equivalency_k
    }
  }
  sim_all[[i]] <- sim_i
  eqvl_all[[i]] <- eqvl_i
}
saveRDS(sim_all,"./Similarity/Similarity_all_all.rds")
saveRDS(eqvl_all,"./Equivalency/Equivalency_all_all.rds")
write.csv(overlap,"./Overlap/Overlap_all.csv")
write.csv(sim_p,"./Similarity/Similarity_p_all.csv")
write.csv(eqvl_p,"./Equivalency/Equivalency_p_all.csv")

##  N = 20  ----
pollen_all <- array(NA,dim = c(20*6,length(taxa_all)*8,0))
sim_all <- list()
eqvl_all <- list()
overlap <- array(NA,dim = c(5,length(taxa_all),0))
sim_p <- array(NA,dim = c(5,length(taxa_all),0))
eqvl_p <- array(NA,dim = c(5,length(taxa_all),0))
for (l in 1:100) {
  pollen_l <- matrix(NA,ncol = 0,nrow = 20*6)
  sim_all_l <- list()
  eqvl_all_l <- list()
  overlap_l <- matrix(NA,ncol = length(taxa_all),nrow = 5)
  colnames(overlap_l) <- taxa_all
  rownames(overlap_l) <- c(0,seq(2,14,4))
  sim_p_l <- overlap_l
  eqvl_p_l <- overlap_l
  
  for (i in 1:length(taxa_all)) {
    ##  produce the randomly sampled pollen data
    pollen_i <- matrix(NA,ncol = 8,nrow = 0)
    for (k in 1:6) {
      if (k==1) 
        pollen_k <- subset(pollen, Age <= 0)[,c(which(colnames(pollen)==taxa_all[i]),1:7)]
      if (k==2) 
        pollen_k <- subset(pollen, Age <= 2000 & Age > 0)[,c(which(colnames(pollen)==taxa_all[i]),1:7)]
      if (k > 2) 
        pollen_k <- subset(pollen, Age <= (k-2)*4000+2000 & Age > (k-3)*4000+2000)[,c(which(colnames(pollen)==taxa_all[i]),1:7)]
      
      pollen_k_i <- pollen_k[pollen_k[,1]>0,]
      if (nrow(pollen_k_i) >= 20) {
        pollen_k_i <- pollen_k_i[sample(1:nrow(pollen_k_i),20),]
        pollen_i <- rbind(pollen_i,pollen_k_i)
      }
      if (nrow(pollen_k_i) < 20) {
        pollen_NA <- matrix(NA,ncol = 8,nrow = 20)
        colnames(pollen_NA) <- colnames(pollen_k_i)
        pollen_i <- rbind(pollen_i,pollen_NA)
      }
    }
    
    ##  calculate climate fidelity
    sim_i <- list()
    eqvl_i <- list()
    for (k in 1:5) {
      if (k==1) {
        pollen_k <- subset(pollen_i, Age <= 0)[,c(which(colnames(pollen_i)==taxa_all[i]),7,8)]
        pollen_kplus <- subset(pollen_i, Age <= 2000 & Age > 0)[,c(which(colnames(pollen_i)==taxa_all[i]),7,8)]
      }
      if (k==2) {
        pollen_k <- subset(pollen_i, Age <= 2000 & Age > 0)[,c(which(colnames(pollen_i)==taxa_all[i]),7,8)]
        pollen_kplus <- subset(pollen_i, Age <= (k+1-2)*4000+2000 & Age > (k+1-3)*4000+2000)[,c(which(colnames(pollen_i)==taxa_all[i]),7,8)]
      }
      if (k > 2) {
        pollen_k <- subset(pollen_i, Age <= (k-2)*4000+2000 & Age > (k-3)*4000+2000)[,c(which(colnames(pollen_i)==taxa_all[i]),7,8)]
        pollen_kplus <- subset(pollen_i, Age <= (k+1-2)*4000+2000 & Age > (k+1-3)*4000+2000)[,c(which(colnames(pollen_i)==taxa_all[i]),7,8)]
      }
      
      pollen_k_i <- pollen_k[pollen_k[,1]>0,]
      pollen_kplus_i <- pollen_kplus[pollen_kplus[,1]>0,]
      
      if (nrow(pollen_k_i) >= 5 & nrow(pollen_kplus_i) >= 5) {
        pollen_i_grid_k <- ecospat.grid.clim.dyn(glob = rbind(pollen_k[,2:3],pollen_kplus[,2:3]),
                                                 glob1 = pollen_k[,2:3],
                                                 sp = pollen_k_i[,2:3],
                                                 R = 100, th.sp = 0.05, th.env = 0.05)
        pollen_i_grid_kplus <- ecospat.grid.clim.dyn(glob = rbind(pollen_k[,2:3],pollen_kplus[,2:3]),
                                                     glob1 = pollen_kplus[,2:3],
                                                     sp = pollen_kplus_i[,2:3],
                                                     R = 100, th.sp = 0.05, th.env = 0.05)
        
        similarity_k <- ecospat.niche.similarity.test(pollen_i_grid_kplus,pollen_i_grid_k,rep = 100)
        equivalency_k <- ecospat.niche.equivalency.test(pollen_i_grid_kplus,pollen_i_grid_k,rep = 100)
        overlap_l[k,i] <- similarity_k$obs$D
        sim_p_l[k,i] <- similarity_k$p.D
        eqvl_p_l[k,i] <- equivalency_k$p.D
        sim_i[[k]] <- similarity_k
        eqvl_i[[k]] <- equivalency_k
      }
    }
    sim_all[[i]] <- sim_i
    eqvl_all[[i]] <- eqvl_i
  }
  pollen_all <- abind(pollen_all,pollen_l,along = 3)
  sim_all[[l]] <- sim_all_l
  eqvl_all[[l]] <- eqvl_all_l
  overlap <- abind(overlap,overlap_l,along = 3)
  sim_p <- abind(sim_p,sim_p_l,along = 3)
  eqvl_p <- abind(eqvl_p,eqvl_p_l,along = 3)
}
saveRDS(pollen_all,"./PollenData/pollen_all_N20.rds")
saveRDS(sim_all,"./Similarity/Similarity_all_N20.rds")
saveRDS(eqvl_all,"./Equivalency/Equivalency_all_N20.rds")
saveRDS(overlap,"./Overlap/Overlap_N20.rds")
saveRDS(sim_p,"./Similarity/Similarity_p_N20.rds")
saveRDS(eqvl_p,"./Equivalency/Equivalency_p_N20.rds")

##  N = 50  ----
pollen_all <- array(NA,dim = c(50*6,length(taxa_all)*8,0))
sim_all <- list()
eqvl_all <- list()
overlap <- array(NA,dim = c(5,length(taxa_all),0))
sim_p <- array(NA,dim = c(5,length(taxa_all),0))
eqvl_p <- array(NA,dim = c(5,length(taxa_all),0))
for (l in 1:100) {
  pollen_l <- matrix(NA,ncol = 0,nrow = 50*6)
  sim_all_l <- list()
  eqvl_all_l <- list()
  overlap_l <- matrix(NA,ncol = length(taxa_all),nrow = 5)
  colnames(overlap_l) <- taxa_all
  rownames(overlap_l) <- c(0,seq(2,14,4))
  sim_p_l <- overlap_l
  eqvl_p_l <- overlap_l
  
  for (i in 1:length(taxa_all)) {
    ##  produce the randomly sampled pollen data
    pollen_i <- matrix(NA,ncol = 8,nrow = 0)
    for (k in 1:6) {
      if (k==1) 
        pollen_k <- subset(pollen, Age <= 0)[,c(which(colnames(pollen)==taxa_all[i]),1:7)]
      if (k==2) 
        pollen_k <- subset(pollen, Age <= 2000 & Age > 0)[,c(which(colnames(pollen)==taxa_all[i]),1:7)]
      if (k > 2) 
        pollen_k <- subset(pollen, Age <= (k-2)*4000+2000 & Age > (k-3)*4000+2000)[,c(which(colnames(pollen)==taxa_all[i]),1:7)]
      
      pollen_k_i <- pollen_k[pollen_k[,1]>0,]
      if (nrow(pollen_k_i) >= 50) {
        pollen_k_i <- pollen_k_i[sample(1:nrow(pollen_k_i),50),]
        pollen_i <- rbind(pollen_i,pollen_k_i)
      }
      if (nrow(pollen_k_i) < 50) {
        pollen_NA <- matrix(NA,ncol = 8,nrow = 50)
        colnames(pollen_NA) <- colnames(pollen_k_i)
        pollen_i <- rbind(pollen_i,pollen_NA)
      }
    }
    
    ##  calculate climate fidelity
    sim_i <- list()
    eqvl_i <- list()
    for (k in 1:5) {
      if (k==1) {
        pollen_k <- subset(pollen_i, Age <= 0)[,c(which(colnames(pollen_i)==taxa_all[i]),7,8)]
        pollen_kplus <- subset(pollen_i, Age <= 2000 & Age > 0)[,c(which(colnames(pollen_i)==taxa_all[i]),7,8)]
      }
      if (k==2) {
        pollen_k <- subset(pollen_i, Age <= 2000 & Age > 0)[,c(which(colnames(pollen_i)==taxa_all[i]),7,8)]
        pollen_kplus <- subset(pollen_i, Age <= (k+1-2)*4000+2000 & Age > (k+1-3)*4000+2000)[,c(which(colnames(pollen_i)==taxa_all[i]),7,8)]
      }
      if (k > 2) {
        pollen_k <- subset(pollen_i, Age <= (k-2)*4000+2000 & Age > (k-3)*4000+2000)[,c(which(colnames(pollen_i)==taxa_all[i]),7,8)]
        pollen_kplus <- subset(pollen_i, Age <= (k+1-2)*4000+2000 & Age > (k+1-3)*4000+2000)[,c(which(colnames(pollen_i)==taxa_all[i]),7,8)]
      }
      
      pollen_k_i <- pollen_k[pollen_k[,1]>0,]
      pollen_kplus_i <- pollen_kplus[pollen_kplus[,1]>0,]
      
      if (nrow(pollen_k_i) >= 5 & nrow(pollen_kplus_i) >= 5) {
        pollen_i_grid_k <- ecospat.grid.clim.dyn(glob = rbind(pollen_k[,2:3],pollen_kplus[,2:3]),
                                                 glob1 = pollen_k[,2:3],
                                                 sp = pollen_k_i[,2:3],
                                                 R = 100, th.sp = 0.05, th.env = 0.05)
        pollen_i_grid_kplus <- ecospat.grid.clim.dyn(glob = rbind(pollen_k[,2:3],pollen_kplus[,2:3]),
                                                     glob1 = pollen_kplus[,2:3],
                                                     sp = pollen_kplus_i[,2:3],
                                                     R = 100, th.sp = 0.05, th.env = 0.05)
        
        similarity_k <- ecospat.niche.similarity.test(pollen_i_grid_kplus,pollen_i_grid_k,rep = 100)
        equivalency_k <- ecospat.niche.equivalency.test(pollen_i_grid_kplus,pollen_i_grid_k,rep = 100)
        overlap_l[k,i] <- similarity_k$obs$D
        sim_p_l[k,i] <- similarity_k$p.D
        eqvl_p_l[k,i] <- equivalency_k$p.D
        sim_i[[k]] <- similarity_k
        eqvl_i[[k]] <- equivalency_k
      }
    }
    sim_all[[i]] <- sim_i
    eqvl_all[[i]] <- eqvl_i
  }
  pollen_all <- abind(pollen_all,pollen_l,along = 3)
  sim_all[[l]] <- sim_all_l
  eqvl_all[[l]] <- eqvl_all_l
  overlap <- abind(overlap,overlap_l,along = 3)
  sim_p <- abind(sim_p,sim_p_l,along = 3)
  eqvl_p <- abind(eqvl_p,eqvl_p_l,along = 3)
}
saveRDS(pollen_all,"./PollenData/pollen_all_N50.rds")
saveRDS(sim_all,"./Similarity/Similarity_all_N50.rds")
saveRDS(eqvl_all,"./Equivalency/Equivalency_all_N50.rds")
saveRDS(overlap,"./Overlap/Overlap_N50.rds")
saveRDS(sim_p,"./Similarity/Similarity_p_N50.rds")
saveRDS(eqvl_p,"./Equivalency/Equivalency_p_N50.rds")

##  N = 100  ----
pollen_all <- array(NA,dim = c(100*6,length(taxa_all)*8,0))
sim_all <- list()
eqvl_all <- list()
overlap <- array(NA,dim = c(5,length(taxa_all),0))
sim_p <- array(NA,dim = c(5,length(taxa_all),0))
eqvl_p <- array(NA,dim = c(5,length(taxa_all),0))
for (l in 1:100) {
  pollen_l <- matrix(NA,ncol = 0,nrow = 100*6)
  sim_all_l <- list()
  eqvl_all_l <- list()
  overlap_l <- matrix(NA,ncol = length(taxa_all),nrow = 5)
  colnames(overlap_l) <- taxa_all
  rownames(overlap_l) <- c(0,seq(2,14,4))
  sim_p_l <- overlap_l
  eqvl_p_l <- overlap_l
  
  for (i in 1:length(taxa_all)) {
    ##  produce the randomly sampled pollen data
    pollen_i <- matrix(NA,ncol = 8,nrow = 0)
    for (k in 1:6) {
      if (k==1) 
        pollen_k <- subset(pollen, Age <= 0)[,c(which(colnames(pollen)==taxa_all[i]),1:7)]
      if (k==2) 
        pollen_k <- subset(pollen, Age <= 2000 & Age > 0)[,c(which(colnames(pollen)==taxa_all[i]),1:7)]
      if (k > 2) 
        pollen_k <- subset(pollen, Age <= (k-2)*4000+2000 & Age > (k-3)*4000+2000)[,c(which(colnames(pollen)==taxa_all[i]),1:7)]
      
      pollen_k_i <- pollen_k[pollen_k[,1]>0,]
      if (nrow(pollen_k_i) >= 100) {
        pollen_k_i <- pollen_k_i[sample(1:nrow(pollen_k_i),100),]
        pollen_i <- rbind(pollen_i,pollen_k_i)
      }
      if (nrow(pollen_k_i) < 100) {
        pollen_NA <- matrix(NA,ncol = 8,nrow = 100)
        colnames(pollen_NA) <- colnames(pollen_k_i)
        pollen_i <- rbind(pollen_i,pollen_NA)
      }
    }
    
    ##  calculate climate fidelity
    sim_i <- list()
    eqvl_i <- list()
    for (k in 1:5) {
      if (k==1) {
        pollen_k <- subset(pollen_i, Age <= 0)[,c(which(colnames(pollen_i)==taxa_all[i]),7,8)]
        pollen_kplus <- subset(pollen_i, Age <= 2000 & Age > 0)[,c(which(colnames(pollen_i)==taxa_all[i]),7,8)]
      }
      if (k==2) {
        pollen_k <- subset(pollen_i, Age <= 2000 & Age > 0)[,c(which(colnames(pollen_i)==taxa_all[i]),7,8)]
        pollen_kplus <- subset(pollen_i, Age <= (k+1-2)*4000+2000 & Age > (k+1-3)*4000+2000)[,c(which(colnames(pollen_i)==taxa_all[i]),7,8)]
      }
      if (k > 2) {
        pollen_k <- subset(pollen_i, Age <= (k-2)*4000+2000 & Age > (k-3)*4000+2000)[,c(which(colnames(pollen_i)==taxa_all[i]),7,8)]
        pollen_kplus <- subset(pollen_i, Age <= (k+1-2)*4000+2000 & Age > (k+1-3)*4000+2000)[,c(which(colnames(pollen_i)==taxa_all[i]),7,8)]
      }
      
      pollen_k_i <- pollen_k[pollen_k[,1]>0,]
      pollen_kplus_i <- pollen_kplus[pollen_kplus[,1]>0,]
      
      if (nrow(pollen_k_i) >= 5 & nrow(pollen_kplus_i) >= 5) {
        pollen_i_grid_k <- ecospat.grid.clim.dyn(glob = rbind(pollen_k[,2:3],pollen_kplus[,2:3]),
                                                 glob1 = pollen_k[,2:3],
                                                 sp = pollen_k_i[,2:3],
                                                 R = 100, th.sp = 0.05, th.env = 0.05)
        pollen_i_grid_kplus <- ecospat.grid.clim.dyn(glob = rbind(pollen_k[,2:3],pollen_kplus[,2:3]),
                                                     glob1 = pollen_kplus[,2:3],
                                                     sp = pollen_kplus_i[,2:3],
                                                     R = 100, th.sp = 0.05, th.env = 0.05)
        
        similarity_k <- ecospat.niche.similarity.test(pollen_i_grid_kplus,pollen_i_grid_k,rep = 100)
        equivalency_k <- ecospat.niche.equivalency.test(pollen_i_grid_kplus,pollen_i_grid_k,rep = 100)
        overlap_l[k,i] <- similarity_k$obs$D
        sim_p_l[k,i] <- similarity_k$p.D
        eqvl_p_l[k,i] <- equivalency_k$p.D
        sim_i[[k]] <- similarity_k
        eqvl_i[[k]] <- equivalency_k
      }
    }
    sim_all[[i]] <- sim_i
    eqvl_all[[i]] <- eqvl_i
  }
  pollen_all <- abind(pollen_all,pollen_l,along = 3)
  sim_all[[l]] <- sim_all_l
  eqvl_all[[l]] <- eqvl_all_l
  overlap <- abind(overlap,overlap_l,along = 3)
  sim_p <- abind(sim_p,sim_p_l,along = 3)
  eqvl_p <- abind(eqvl_p,eqvl_p_l,along = 3)
}
saveRDS(pollen_all,"./PollenData/pollen_all_N100.rds")
saveRDS(sim_all,"./Similarity/Similarity_all_N100.rds")
saveRDS(eqvl_all,"./Equivalency/Equivalency_all_N100.rds")
saveRDS(overlap,"./Overlap/Overlap_N100.rds")
saveRDS(sim_p,"./Similarity/Similarity_p_N100.rds")
saveRDS(eqvl_p,"./Equivalency/Equivalency_p_N100.rds")

##  N = DG samples ----
SampleSize <- vector(length = 6)
for (k in 1:6) {
  if (k==1) 
    pollen_k <- subset(pollen, Age <= 0)
  if (k==2) 
    pollen_k <- subset(pollen, Age <= 2000 & Age > 0)
  if (k > 2) 
    pollen_k <- subset(pollen, Age <= (k-2)*4000+2000 & Age > (k-3)*4000+2000)
  
  SampleSize[k] <- nrow(pollen_k)
}

##  calculation
pollen_all <- array(NA,dim = c(SampleSize[6]*6,ncol(pollen),0))
sim_all <- list()
eqvl_all <- list()
overlap <- array(NA,dim = c(5,length(taxa_all),0))
sim_p <- array(NA,dim = c(5,length(taxa_all),0))
eqvl_p <- array(NA,dim = c(5,length(taxa_all),0))
for (l in 1:100) {
  ##  produce the randomly sampled pollen data
  pollen_l <- subset(pollen, Age <= (6-2)*4000+2000 & Age > (6-3)*4000+2000)
  for (k in 5:1) {
    if (k==1) 
      pollen_k <- subset(pollen, Age <= 0)
    if (k==2) 
      pollen_k <- subset(pollen, Age <= 2000 & Age > 0)
    if (k > 2) 
      pollen_k <- subset(pollen, Age <= (k-2)*4000+2000 & Age > (k-3)*4000+2000)
    
    pollen_k <- pollen_k[sample(1:nrow(pollen_k),SampleSize[6]),]
    pollen_l <- rbind(pollen_l,pollen_k)
  }
  pollen_all <- abind(pollen_all,pollen_l,along = 3)
  
  ##  calculate climate fidelity
  for (i in 1:length(taxa_all)) {
    sim_i <- list()
    eqvl_i <- list()
    for (k in 1:5) {
      if (k==1) {
        pollen_k <- subset(pollen_l, Age <= 0)[,c(which(colnames(pollen_l)==taxa_all[i]),6,7)]
        pollen_kplus <- subset(pollen_l, Age <= 2000 & Age > 0)[,c(which(colnames(pollen_l)==taxa_all[i]),6,7)]
      }
      if (k==2) {
        pollen_k <- subset(pollen_l, Age <= 2000 & Age > 0)[,c(which(colnames(pollen_l)==taxa_all[i]),6,7)]
        pollen_kplus <- subset(pollen_l, Age <= (k+1-2)*4000+2000 & Age > (k+1-3)*4000+2000)[,c(which(colnames(pollen_l)==taxa_all[i]),6,7)]
      }
      if (k > 2) {
        pollen_k <- subset(pollen_l, Age <= (k-2)*4000+2000 & Age > (k-3)*4000+2000)[,c(which(colnames(pollen_l)==taxa_all[i]),6,7)]
        pollen_kplus <- subset(pollen_l, Age <= (k+1-2)*4000+2000 & Age > (k+1-3)*4000+2000)[,c(which(colnames(pollen_l)==taxa_all[i]),6,7)]
      }
      
      pollen_k_i <- pollen_k[pollen_k[,1]>0,]
      pollen_kplus_i <- pollen_kplus[pollen_kplus[,1]>0,]
      
      if (nrow(pollen_k_i) >= 5 & nrow(pollen_kplus_i) >= 5) {
        pollen_i_grid_k <- ecospat.grid.clim.dyn(glob = rbind(pollen_k[,2:3],pollen_kplus[,2:3]),
                                                 glob1 = pollen_k[,2:3],
                                                 sp = pollen_k_i[,2:3],
                                                 R = 100, th.sp = 0.05, th.env = 0.05)
        pollen_i_grid_kplus <- ecospat.grid.clim.dyn(glob = rbind(pollen_k[,2:3],pollen_kplus[,2:3]),
                                                     glob1 = pollen_kplus[,2:3],
                                                     sp = pollen_kplus_i[,2:3],
                                                     R = 100, th.sp = 0.05, th.env = 0.05)
        
        similarity_k <- ecospat.niche.similarity.test(pollen_i_grid_kplus,pollen_i_grid_k,rep = 100)
        equivalency_k <- ecospat.niche.equivalency.test(pollen_i_grid_kplus,pollen_i_grid_k,rep = 100)
        overlap_l[k,i] <- similarity_k$obs$D
        sim_p_l[k,i] <- similarity_k$p.D
        eqvl_p_l[k,i] <- equivalency_k$p.D
        sim_i[[k]] <- similarity_k
        eqvl_i[[k]] <- equivalency_k
      }
    }
    sim_all[[i]] <- sim_i
    eqvl_all[[i]] <- eqvl_i
  }
  sim_all[[l]] <- sim_all_l
  eqvl_all[[l]] <- eqvl_all_l
  overlap <- abind(overlap,overlap_l,along = 3)
  sim_p <- abind(sim_p,sim_p_l,along = 3)
  eqvl_p <- abind(eqvl_p,eqvl_p_l,along = 3)
}
saveRDS(pollen_all,"./Pollen/pollen_all_DG.rds")
saveRDS(sim_all,"./Similarity/Similarity_all_DG.rds")
saveRDS(eqvl_all,"./Equivalency/Equivalency_all_DG.rds")
saveRDS(overlap,"./Overlap/Overlap_DG.rds")
saveRDS(sim_p,"./Similarity/Similarity_p_DG.rds")
saveRDS(eqvl_p,"./Equivalency/Equivalency_p_DG.rds")


##  3. Distance: kernel density ====
##  N = all ----
distance_all <- matrix(NA,nrow=5,ncol=length(taxa_all))
colnames(distance_all) <- taxa_all
rownames(distance_all) <- c(0,seq(2,14,4))
for (i in 1:length(taxa_all)) {
  pollen_i <- pollen[pollen[,which(colnames(pollen)==taxa_all[i])]>0,c(i,2:4)]
  
  for (k in 1:5) {
    if (k==1) {
      pollen_i_k <- subset(pollen_i, Age <= 0)
      pollen_i_kplus <- subset(pollen_i, Age <= 2000 & Age > 0)
    }
    if (k==2) {
      pollen_i_k <- subset(pollen_i, Age <= 2000 & Age > 0)
      pollen_i_kplus <- subset(pollen_i, Age <= (k+1-2)*4000+2000 & Age > (k+1-3)*4000+2000)
    }
    if (k > 2) {
      pollen_i_k <- subset(pollen_i, Age <= (k-2)*4000+2000 & Age > (k-3)*4000+2000)
      pollen_i_kplus <- subset(pollen_i, Age <= (k+1-2)*4000+2000 & Age > (k+1-3)*4000+2000)
    }
    
    pollen_i_k <- unique(pollen_i_k[,2:3])
    pollen_i_kplus <- unique(pollen_i_kplus[,2:3])
    
    if (nrow(pollen_i_k) > 1 & nrow(pollen_i_kplus) > 1) {
      pollen_i_k_kde <- kde2d(pollen_i_k$lon,pollen_i_k$lat,n=200,lims=c(range(pollen$lon),range(pollen$lat)))
      pollen_i_kplus_kde <- kde2d(pollen_i_kplus$lon,pollen_i_kplus$lat,n=200,lims=c(range(pollen$lon),range(pollen$lat)))
      
      pollen_i_k_kde_0.95 <- as.numeric(pollen_i_k_kde[[3]])[which(pollen_i_k_kde[[3]] > quantile(pollen_i_k_kde[[3]],0.05))]
      pollen_i_k_kde_0.95_w <- pollen_i_k_kde_0.95/sum(pollen_i_k_kde_0.95)
      pollen_i_k_kde_0.95_x <- rep(pollen_i_k_kde[[1]],each=200)[which(pollen_i_k_kde[[3]] > quantile(pollen_i_k_kde[[3]],0.05))]
      pollen_i_k_kde_0.95_y <- rep(pollen_i_k_kde[[2]],times=200)[which(pollen_i_k_kde[[3]] > quantile(pollen_i_k_kde[[3]],0.05))]
      pollen_i_k_kde_centroid <- c(sum(pollen_i_k_kde_0.95_x*pollen_i_k_kde_0.95_w),sum(pollen_i_k_kde_0.95_y*pollen_i_k_kde_0.95_w))
      
      pollen_i_kplus_kde_0.95 <- as.numeric(pollen_i_kplus_kde[[3]])[which(pollen_i_kplus_kde[[3]] > quantile(pollen_i_kplus_kde[[3]],0.05))]
      pollen_i_kplus_kde_0.95_w <- pollen_i_kplus_kde_0.95/sum(pollen_i_kplus_kde_0.95)
      pollen_i_kplus_kde_0.95_x <- rep(pollen_i_kplus_kde[[1]],each=200)[which(pollen_i_kplus_kde[[3]] > quantile(pollen_i_kplus_kde[[3]],0.05))]
      pollen_i_kplus_kde_0.95_y <- rep(pollen_i_kplus_kde[[2]],times=200)[which(pollen_i_kplus_kde[[3]] > quantile(pollen_i_kplus_kde[[3]],0.05))]
      pollen_i_kplus_kde_centroid <- c(sum(pollen_i_kplus_kde_0.95_x*pollen_i_kplus_kde_0.95_w),sum(pollen_i_kplus_kde_0.95_y*pollen_i_kplus_kde_0.95_w))
      
      distance_all[k,i] <- pointDistance(pollen_i_k_kde_centroid,pollen_i_kplus_kde_centroid,lonlat = T)
    }
  }
}
write.csv(distance_all,"./Distance/Distance_all_kd.csv")

##  N = 20, 50, 100, DG samples ----
pollen_N20 <- readRDS("./Pollen/pollen_all_N20.rds")
pollen_N50 <- readRDS("./Pollen/pollen_all_N50.rds")
pollen_N100 <- readRDS("./Pollen/pollen_all_N100.rds")
pollen_DG <- readRDS("./Pollen/pollen_all_DG.rds")
distance_N20 <- array(NA,dim = c(5,length(taxa_all),0))
distance_N50 <- array(NA,dim = c(5,length(taxa_all),0))
distance_N100 <- array(NA,dim = c(5,length(taxa_all),0))
distance_DG <- array(NA,dim = c(5,length(taxa_all),0))
for (l in 1:100) {
  distance_N20_l <- matrix(NA,nrow=5,ncol=length(taxa_all))
  distance_N50_l <- matrix(NA,nrow=5,ncol=length(taxa_all))
  distance_N100_l <- matrix(NA,nrow=5,ncol=length(taxa_all))
  distance_DG_l <- matrix(NA,nrow=5,ncol=length(taxa_all))
  
  for (i in 1:length(taxa_all)) {
    taxa_i <- which(colnames(pollen)==taxa_all[i])
    
    pollen_N20_i <- as.data.frame(pollen_N20[,,l][,(i-1)*8+c(1,3:5)])
    pollen_N50_i <- as.data.frame(pollen_N50[,,l][,(i-1)*8+c(1,3:5)])
    pollen_N100_i <- as.data.frame(pollen_N100[,,l][,(i-1)*8+c(1,3:5)])
    pollen_DG_i <- as.data.frame(pollen_DG[,,l][,c(taxa_i,2:4)])
    
    pollen_N20_i <- as.data.frame(sapply(pollen_N20_i,as.numeric))
    pollen_N50_i <- as.data.frame(sapply(pollen_N50_i,as.numeric))
    pollen_N100_i <- as.data.frame(sapply(pollen_N100_i,as.numeric))
    pollen_DG_i <- as.data.frame(sapply(pollen_DG_i,as.numeric))
    pollen_DG_i <- pollen_DG_i[pollen_DG_i[,1]>0,]
    
    for (k in 1:5) {
      if (k==1) {
        pollen_N20_i_k <- subset(pollen_N20_i, Age <= 0)
        pollen_N50_i_k <- subset(pollen_N50_i, Age <= 0)
        pollen_N100_i_k <- subset(pollen_N100_i, Age <= 0)
        pollen_DG_i_k <- subset(pollen_DG_i, Age <= 0)
        pollen_N20_i_kplus <- subset(pollen_N20_i, Age <= 2000 & Age > 0)
        pollen_N50_i_kplus <- subset(pollen_N50_i, Age <= 2000 & Age > 0)
        pollen_N100_i_kplus <- subset(pollen_N100_i, Age <= 2000 & Age > 0)
        pollen_DG_i_kplus <- subset(pollen_DG_i, Age <= 2000 & Age > 0)
      }
      if (k==2) {
        pollen_N20_i_k <- subset(pollen_N20_i, Age <= 2000 & Age > 0)
        pollen_N50_i_k <- subset(pollen_N50_i, Age <= 2000 & Age > 0)
        pollen_N100_i_k <- subset(pollen_N100_i, Age <= 2000 & Age > 0)
        pollen_DG_i_k <- subset(pollen_DG_i, Age <= 2000 & Age > 0)
        pollen_N20_i_kplus <- subset(pollen_N20_i, Age <= (k+1-2)*4000+2000 & Age > (k+1-3)*4000+2000)
        pollen_N50_i_kplus <- subset(pollen_N50_i, Age <= (k+1-2)*4000+2000 & Age > (k+1-3)*4000+2000)
        pollen_N100_i_kplus <- subset(pollen_N100_i, Age <= (k+1-2)*4000+2000 & Age > (k+1-3)*4000+2000)
        pollen_DG_i_kplus <- subset(pollen_DG_i, Age <= (k+1-2)*4000+2000 & Age > (k+1-3)*4000+2000)
      }
      if (k > 2) {
        pollen_N20_i_k <- subset(pollen_N20_i, Age <= (k-2)*4000+2000 & Age > (k-3)*4000+2000)
        pollen_N50_i_k <- subset(pollen_N50_i, Age <= (k-2)*4000+2000 & Age > (k-3)*4000+2000)
        pollen_N100_i_k <- subset(pollen_N100_i, Age <= (k-2)*4000+2000 & Age > (k-3)*4000+2000)
        pollen_DG_i_k <- subset(pollen_DG_i, Age <= (k-2)*4000+2000 & Age > (k-3)*4000+2000)
        pollen_N20_i_kplus <- subset(pollen_N20_i, Age <= (k+1-2)*4000+2000 & Age > (k+1-3)*4000+2000)
        pollen_N50_i_kplus <- subset(pollen_N50_i, Age <= (k+1-2)*4000+2000 & Age > (k+1-3)*4000+2000)
        pollen_N100_i_kplus <- subset(pollen_N100_i, Age <= (k+1-2)*4000+2000 & Age > (k+1-3)*4000+2000)
        pollen_DG_i_kplus <- subset(pollen_DG_i, Age <= (k+1-2)*4000+2000 & Age > (k+1-3)*4000+2000)
      }
      
      pollen_N20_i_k <- unique(pollen_N20_i_k[,2:3])
      pollen_N50_i_k <- unique(pollen_N50_i_k[,2:3])
      pollen_N100_i_k <- unique(pollen_N100_i_k[,2:3])
      pollen_DG_i_k <- unique(pollen_DG_i_k[,2:3])
      pollen_N20_i_kplus <- unique(pollen_N20_i_kplus[,2:3])
      pollen_N50_i_kplus <- unique(pollen_N50_i_kplus[,2:3])
      pollen_N100_i_kplus <- unique(pollen_N100_i_kplus[,2:3])
      pollen_DG_i_kplus <- unique(pollen_DG_i_kplus[,2:3])
      
      for (m in 1:4) {
        if (m == 1) {
          pollen_i_k <- pollen_N20_i_k
          pollen_i_kplus <- pollen_N20_i_kplus
        }
        if (m == 2) {
          pollen_i_k <- pollen_N50_i_k
          pollen_i_kplus <- pollen_N50_i_kplus
        }
        if (m == 3) {
          pollen_i_k <- pollen_N100_i_k
          pollen_i_kplus <- pollen_N100_i_kplus
        }
        if (m == 4) {
          pollen_i_k <- pollen_DG_i_k
          pollen_i_kplus <- pollen_DG_i_kplus
        }
        
        if (nrow(pollen_i_k) > 1 & nrow(pollen_i_kplus) > 1) {
          pollen_i_k_kde <- kde2d(pollen_i_k$lon,pollen_i_k$lat,n=200,lims=c(range(pollen$lon),range(pollen$lat)))
          pollen_i_kplus_kde <- kde2d(pollen_i_kplus$lon,pollen_i_kplus$lat,n=200,lims=c(range(pollen$lon),range(pollen$lat)))
          
          pollen_i_k_kde_0.95 <- as.numeric(pollen_i_k_kde[[3]])[which(pollen_i_k_kde[[3]] > quantile(pollen_i_k_kde[[3]],0.05))]
          pollen_i_k_kde_0.95_w <- pollen_i_k_kde_0.95/sum(pollen_i_k_kde_0.95)
          pollen_i_k_kde_0.95_x <- rep(pollen_i_k_kde[[1]],each=200)[which(pollen_i_k_kde[[3]] > quantile(pollen_i_k_kde[[3]],0.05))]
          pollen_i_k_kde_0.95_y <- rep(pollen_i_k_kde[[2]],times=200)[which(pollen_i_k_kde[[3]] > quantile(pollen_i_k_kde[[3]],0.05))]
          pollen_i_k_kde_centroid <- c(sum(pollen_i_k_kde_0.95_x*pollen_i_k_kde_0.95_w),sum(pollen_i_k_kde_0.95_y*pollen_i_k_kde_0.95_w))
          
          pollen_i_kplus_kde_0.95 <- as.numeric(pollen_i_kplus_kde[[3]])[which(pollen_i_kplus_kde[[3]] > quantile(pollen_i_kplus_kde[[3]],0.05))]
          pollen_i_kplus_kde_0.95_w <- pollen_i_kplus_kde_0.95/sum(pollen_i_kplus_kde_0.95)
          pollen_i_kplus_kde_0.95_x <- rep(pollen_i_kplus_kde[[1]],each=200)[which(pollen_i_kplus_kde[[3]] > quantile(pollen_i_kplus_kde[[3]],0.05))]
          pollen_i_kplus_kde_0.95_y <- rep(pollen_i_kplus_kde[[2]],times=200)[which(pollen_i_kplus_kde[[3]] > quantile(pollen_i_kplus_kde[[3]],0.05))]
          pollen_i_kplus_kde_centroid <- c(sum(pollen_i_kplus_kde_0.95_x*pollen_i_kplus_kde_0.95_w),sum(pollen_i_kplus_kde_0.95_y*pollen_i_kplus_kde_0.95_w))
        }
        
        if (m == 1 & nrow(pollen_i_k) > 1 & nrow(pollen_i_kplus) > 1)
          distance_N20_l[k,i] <- pointDistance(pollen_i_k_kde_centroid,pollen_i_kplus_kde_centroid,lonlat = T)
        if (m == 2 & nrow(pollen_i_k) > 1 & nrow(pollen_i_kplus) > 1)
          distance_N50_l[k,i] <- pointDistance(pollen_i_k_kde_centroid,pollen_i_kplus_kde_centroid,lonlat = T)
        if (m == 3 & nrow(pollen_i_k) > 1 & nrow(pollen_i_kplus) > 1)
          distance_N100_l[k,i] <- pointDistance(pollen_i_k_kde_centroid,pollen_i_kplus_kde_centroid,lonlat = T)
        if (m == 4 & nrow(pollen_i_k) > 1 & nrow(pollen_i_kplus) > 1)
          distance_DG_l[k,i] <- pointDistance(pollen_i_k_kde_centroid,pollen_i_kplus_kde_centroid,lonlat = T)
      }
    }
  }
  distance_N20 <- abind(distance_N20,distance_N20_l,along = 3)
  distance_N50 <- abind(distance_N50,distance_N50_l,along = 3)
  distance_N100 <- abind(distance_N100,distance_N100_l,along = 3)
  distance_DG <- abind(distance_DG,distance_DG_l,along = 3)
}
saveRDS(distance_N20,"./Distance/Distance_N20_kd.rds")
saveRDS(distance_N50,"./Distance/Distance_N50_kd.rds")
saveRDS(distance_N100,"./Distance/Distance_N100_kd.rds")
saveRDS(distance_DG,"./Distance/Distance_DG_kd.rds")


##  4. Figures  ====
##  niche overlap  ----
titie_name <- c("Recent","Late Holocene (LH)","Mid-Holocene (MH)","Early Holocene (EH)","Deglaciation (DG)")
for (i in 1:length(taxa_all)) {
  pdf(paste0("./Overlap/",taxa_all[i],".pdf"),height = 5,width = 8,useDingbats = F)
  layout(matrix(seq(1,6,1),nrow = 2,ncol = 3,byrow = T))
  for (k in 5:1) {
    if (k==1) {
      pollen_k <- subset(pollen, Age <= 0)[,c(which(colnames(pollen)==taxa_all[i]),6,7)]
      pollen_kplus <- subset(pollen, Age <= 2000 & Age > 0)[,c(which(colnames(pollen)==taxa_all[i]),6,7)]
    }
    if (k==2) {
      pollen_k <- subset(pollen, Age <= 2000 & Age > 0)[,c(which(colnames(pollen)==taxa_all[i]),6,7)]
      pollen_kplus <- subset(pollen, Age <= (k+1-2)*4000+2000 & Age > (k+1-3)*4000+2000)[,c(which(colnames(pollen)==taxa_all[i]),6,7)]
    }
    if (k > 2) {
      pollen_k <- subset(pollen, Age <= (k-2)*4000+2000 & Age > (k-3)*4000+2000)[,c(which(colnames(pollen)==taxa_all[i]),6,7)]
      pollen_kplus <- subset(pollen, Age <= (k+1-2)*4000+2000 & Age > (k+1-3)*4000+2000)[,c(which(colnames(pollen)==taxa_all[i]),6,7)]
    }
    
    pollen_k_i <- pollen_k[pollen_k[,1]>0,]
    pollen_kplus_i <- pollen_kplus[pollen_kplus[,1]>0,]
    
    if (nrow(pollen_k_i) >= 5 & nrow(pollen_kplus_i) >= 5) {
      pollen_i_grid_k <- ecospat.grid.clim.dyn(glob = rbind(pollen_k[,2:3],pollen_kplus[,2:3]),
                                               glob1 = pollen_k[,2:3],
                                               sp = pollen_k_i[,2:3],
                                               R = 100, th.sp = 0.05, th.env = 0.05)
      pollen_i_grid_kplus <- ecospat.grid.clim.dyn(glob = rbind(pollen_k[,2:3],pollen_kplus[,2:3]),
                                                   glob1 = pollen_kplus[,2:3],
                                                   sp = pollen_kplus_i[,2:3],
                                                   R = 100, th.sp = 0.05, th.env = 0.05)
      ecospat.plot.niche.dyn(pollen_i_grid_kplus,
                             pollen_i_grid_k,
                             quant = F,
                             name.axis1 = "Mean Annual Temperature (\u00B0C)",
                             name.axis2 = "Mean Annual Precipitation (mm)",
                             title = titie_name[k],
                             interest = 1,
                             colz1 = plasma(6,alpha = 0.8)[6-k],
                             colz2 = plasma(6,alpha = 0.8)[6-(k-1)],
                             colinter = adjustcolor("turquoise4",alpha.f = 0.5),
                             colZ1 = plasma(6)[6-k],
                             colZ2 = plasma(6)[6-(k-1)])
    }
    if (nrow(pollen_k_i) < 5 | nrow(pollen_kplus_i) < 5) 
      plot.new()
  }
  # legend
  plot.new()
  legend("topleft",
         legend = c("18-14 ka","14-10 ka",
                    "10-6 ka","6-2 ka",
                    "2 ka - 1950 AD","after 1950 AD",
                    "Overlap"),
         title = "Age bins",
         fill = c(plasma(6)[1:6],"turquoise4"),
         border = c(plasma(6)[1:6],"turquoise4"),
         cex = 1.4,
         horiz = F,
         box.lty = 0,
         xpd = T)
  dev.off()
}

##  p-values  ----
sim_all <- readRDS("./Similarity/Similarity_all_all.rds")
eqvl_all <- readRDS("./Equivalency/Equivalency_all_all.rds")

title_name <- c("Recent","Late Holocene (LH)","Mid-Holocene (MH)","Early Holocene (EH)","Deglaciation (DG)")

for (i in 1:length(taxa_all)) {
  pdf(paste0("./Similarity/",taxa_all[i],".pdf"),height = 5,width = 8,useDingbats = F)
  layout(matrix(seq(1,6,1),nrow = 2,ncol = 3,byrow = T))
  for (k in 5:1) 
    ecospat.plot.overlap.test(sim_all[[i]][[k]],"D",paste0(title_name[k]))
  dev.off()
  
  pdf(paste0("./Equivalency/",taxa_all[i],".pdf"),height = 5,width = 8,useDingbats = F)
  layout(matrix(seq(1,6,1),nrow = 2,ncol = 3,byrow = T))
  for (k in 5:1) 
    ecospat.plot.overlap.test(eqvl_all[[i]][[k]],"D",paste0(title_name[k]))
  dev.off()
}


##  5. Resilience vs climate fidelity ----
##  climate fidelity map  ----
##  climate fidelity score for each taxa
##  score: significant: +1; marginally significant: +0.5; not significant: 0
NicheConservatism_taxa <- data.frame(taxa=taxa_dis,
                                     score=c(4.5,4.5,4,5,5,4.5,4.5,5,5,4,5,4.5,3.5,5,4,3.5))
write.csv(NicheConservatism_taxa,"./Resilience/NCscore_plants.csv",row.names=F)

##  modern plant site climate fidelity score
pollen_modern <- subset(pollen,Age == -99)
pollen_modern$TotalSum <- rowSums(pollen_modern[,8:23])
pollen_modern <- subset(pollen_modern,TotalSum > 0.8)
NC <- pollen_modern[,c(2,3,5)]
NC$score <- 0

for (i in 1:nrow(NC)) {
  for (j in 8:23) {
    if (pollen_modern[i,j] > 0) {
      NC_posi <- which(taxa_dis==colnames(pollen_modern)[j])
      NC$score[i] <- NC$score[i] + NicheConservatism_taxa$score[NC_posi]
    }
  }
}

NC$score_4 <- NA
NC$score_4[NC$score <= quantile(NC$score,0.25)] <- 1
NC$score_4[NC$score <= quantile(NC$score,0.5) & NC$score > quantile(NC$score,0.25)] <- 2
NC$score_4[NC$score <= quantile(NC$score,0.75) & NC$score > quantile(NC$score,0.5)] <- 3
NC$score_4[NC$score > quantile(NC$score,0.75)] <- 4
write.csv(NC,"./Resilience/NicheConservatismScore.csv",row.names = F)

##  maps
##  get your map
map

NC$color <- mako(60)[60-NC$score]
NC$color_4 <- NA
NC$color_4[NC$score_4==1] <- mako(50)[45]
NC$color_4[NC$score_4==2] <- mako(50)[29]
NC$color_4[NC$score_4==3] <- mako(50)[20]
NC$color_4[NC$score_4==4] <- mako(50)[5]

pdf("./Resilience/map_NC.pdf",width = 8,height = 7,useDingbats=F)
plot(map, xlim=c(-160,-40), ylim=c(20,80), xlab = "Longitude", ylab="Latitude")
points(NC$lon,NC$lat,col=NC$color,pch=16)
axis(1,at=seq(-160,-40,by=20))
axis(2,at=seq(20,80,by=10))
dev.off()

pdf("./Resilience/map_NC_4.pdf",width = 8,height = 7,useDingbats=F)
plot(map, xlim=c(-160,-40), ylim=c(20,80), xlab = "Longitude", ylab="Latitude")
points(NC$lon,NC$lat,col=NC$color_4,pch=16)
axis(1,at=seq(-160,-40,by=20))
axis(2,at=seq(20,80,by=10))
dev.off()

##  relationship b/t climate fidelity and resilience  ----
NC <- read.csv("./Resilience/NicheConservatismScore.csv")
Resilience <- read.csv("./Resilience/TNC.csv")
Resilience <- subset(Resilience,Resilience > -9999 & Terr_Resil >= -3500 & Landscape_ > -9999 & Climate_Fl > -9999)

##  t-test
##  resilience
t.test(Resilience$Terr_Resil[Resilience$score>quantile(Resilience$score,0.75)],
       Resilience$Terr_Resil[Resilience$score<quantile(Resilience$score,0.25)])
##  p<0.001
##  climate connectivity
t.test(Resilience$Climate_Fl[Resilience$score>quantile(Resilience$score,0.75)],
       Resilience$Climate_Fl[Resilience$score<quantile(Resilience$score,0.25)])
##  p=0.008
