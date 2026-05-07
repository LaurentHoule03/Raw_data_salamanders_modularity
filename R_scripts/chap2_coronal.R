
# Chapter 2

# Coronal distance

source(paste0(getwd(),"/functions.R"))

# Librairies ####

library(stringr)
library(readxl)
library(tidyverse)
library(performance)
library(ggrepel)
library(viridis)

# Importation ####

import.csv <- function(path, pattern = "_", i.name = 2, i.ID = 3, sep = ";"){
  
  library(stringr)
  library(readxl)
  
  files <- list.files(path, full.names = T)
  nfiles <- list.files(path, full.names = F)
  
  out <- list()
  
  for(i in 1:length(files)){
    
    impxl <- read.csv(files[i], sep = sep)
    out[[i]] <- impxl
    
    nSP <- str_split_i(i = i.name, pattern = pattern, nfiles[i])
    nID <- str_split_i(i = i.ID, pattern = pattern, nfiles[i])

    if(is.na(nID) == F){
      
      if(nID == "skull.csv")nID <- NA
      
      if(is.na(nID) == F){
        
        RID <- str_split_i(i = 1, pattern = ".csv", nID)
        if(RID != "")nSP <- paste(nSP,RID,sep = "_")
        
      }
    }

    names(out)[i] <- nSP
    
    
  }
  return(out)
}

sk.data <- import.csv(path = "F:/MicroCT_Scan/xx_landmark_brain_salamanders_xx/Chapiters/doc-chap2/coronal/Morpho_skull")
sf.data <- import.csv(path = "F:/MicroCT_Scan/xx_landmark_brain_salamanders_xx/Chapiters/doc-chap2/coronal/Morpho_soft")

# Calculation of euclidian distances between skull and the mid point of coronal suture and between skull and optic chiasm ####

eucl.dist <- function(data, from = "skull", to = "coronal mid", names = 3, x = 4,y = 5,z = 6){
  
  out <- data.frame(individual = NULL, distance = NULL)
  
  for(i in 1:length(data)){
    
    dif <- setdiff(c(from,to), data[[i]][[names]])
    if(is_empty(dif) == T){
      
      for(j in 1:length(data[[i]][[1]])){
        
        if(data[[i]][j,names] == from){
          from.pos <- data[[i]][j,c(x,y,z)]
          
        }
        if(data[[i]][j,names] == to){
          to.pos <- data[[i]][j,c(x,y,z)]
          
        }
        
      }

      distance <- sqrt((to.pos[[1]] - from.pos[[1]])^2 + (to.pos[[2]] - from.pos[[2]])^2 + (to.pos[[3]] - from.pos[[3]])^2)
      out <- rbind(out, data.frame(indiviual = names(data)[i], distance = distance))
    }else{
      
      print(names(data)[i])
      
    }
  }
  
  colnames(out)[2] <- paste(from,to,sep = "_to_")
  return(out)
}

dist.coro.sk <- eucl.dist(data = sk.data, from = "skull", to = "coronal mid")
dist.coro.sf <- eucl.dist(data = sf.data,from = "skull", to = "coronal mid")

dist.opchi.sf <- eucl.dist(data = sf.data,from = "skull", to = "forebrain-midbrain")

merge.info <- function(data1, data2, names){
  
  l1 <- length(data1[[1]])
  l2 <- length(data2[[1]])
  
  names.common <- c()
  
  for(i in 1:length(data1[[1]])){
    
    
    for(j in 1:length(data2[[1]])){
      
      if(data1$indiviual[i] == data2$indiviual[j]){
        
        data2 <- data2[-j,]
        break
        
      }
      
    }
    
    
  }
  
  new <- rbind(data1,data2)
  return(new)
}

new.coro <- merge.info(dist.coro.sk,dist.coro.sf, names = dist.opchi.sf$indiviual)

# Final distance data.frame

distances.ind <- cbind(dist.opchi.sf, skull_to_coronal = new.coro[match(dist.opchi.sf$indiviual,new.coro$indiviual),2])

distances_WNecturus.ind <- distances.ind[-which(distances.ind$indiviual == "NEMA_N01" | distances.ind$indiviual == "NEMA_N02" ),]

# Averaging species ####

average.dist <- function(data,names, pattern = "_", i = 1){
  
  library(stringr)
  library(tidyverse)
  
  name.sp <- factor(str_split_i(i = i, pattern = pattern, names))
  data$SP <- name.sp
  tab.sp <- as.data.frame(table(name.sp))

  for(j in 1:length(levels(name.sp))){
    
    if(tab.sp$Freq[j] > 1){
      

      dat.j <- filter(data, SP == levels(name.sp)[j])
      data <- data[-which(data$SP == levels(name.sp)[j]),]
      
      new <- data.frame(indiviual = dat.j$SP[1], var1 = mean(dat.j[[2]]), var2 = mean(dat.j[[3]]), SP = dat.j$SP[1])
      colnames(new)[2] <- colnames(dat.j)[2]
      colnames(new)[3] <- colnames(dat.j)[3]
      data <- rbind(data, new)
      
    }
    
  }
  out <- data[,-length(data)]
  colnames(out)[1] <- "Species"
  return(out)
}

distances.sp <- average.dist(data = distances.ind, names = distances.ind$indiviual)

distances_WNecturus.sp <- distances.sp[-which(distances.sp$Species == "NEMA" | distances.sp$Species == "NEMA" ),]

# Plot ####

GR <- ggplot(distances.sp, aes(x = log10(`skull_to_forebrain-midbrain`), y = log10(skull_to_coronal)))+geom_point(size = 4, shape = 21, fill = cividis(100)[20])+labs(x = "Log Optic chiasm position", y = "Log Coronal suture position")+theme_bw()
GR

ggsave("./output/Linear_morphometrics/coronal_optic_dist.pdf", height = 6, width = 8)

GR_n <- ggplot(distances.sp, aes(x = log10(`skull_to_forebrain-midbrain`), y = log10(skull_to_coronal)))+geom_point(size = 4, shape = 21, fill = "grey")+labs(x = "Log Optic chiasm position", y = "Log Coronal suture position")+geom_text_repel(aes(label = Species))+theme_bw()
GR_n

GR_f <- ggplot(distances.sp, aes(x = log10(`skull_to_forebrain-midbrain`), y = log10(skull_to_coronal)))+geom_point(size = 4, shape = 21, aes(fill = abs(log10(`skull_to_forebrain-midbrain`) - log10(skull_to_coronal))))+labs(fill = "optic chiasm - coronal suture difference", x = "Log Optic chiasm position", y = "Log Coronal suture position")+scale_fill_viridis_c()+theme_bw() + theme(legend.position = "bottom")
GR_f

# Pearson correlation test ####

ct <- cor.test(log10(distances.sp$`skull_to_forebrain-midbrain`), log10(distances.sp$skull_to_coronal))
ct

check_normality(ct)
