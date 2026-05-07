
# Chapter 2

# Shape PCA - Step 5

# Laurent Houle

# Importing functions ####

source(paste0("F:/MicroCT_Scan/xx_landmark_brain_salamanders_xx/Chapiters/doc-chap2","/functions.R"))

# Librairies ####

library(stringr)
library(readxl)
library(geomorph)
library(Morpho)
library(tidyverse)
library(ape)
library(phytools)
library(viridis)
library(ggplot2)
library(vegan)
library(scales)
library(MKinfer)
library(ggpubr)


# Importation ####
load("./output/RDA/tree.rda")
load("./output/RDA/allom_rm_data.rda")

dimnames(shape.data.both$shape.allom.rm.sp)[[3]] <- sub(" ", "_", shape.data.both$infos.sp$Species)

land.d.brain <- read_excel(path = "./Patches/description_land_brain.xlsx",1)
land.d.endo <- read_excel("./Patches/description_land_endo.xlsx",1)


# Intraspecific PCA ####

## Brain ####

### PCA objects ####

shape.data.both$symm.coords.ind <- unique.names.array(shape = shape.data.both$symm.coords.ind, ID = shape.data.both$infos$ID2)

pca.brain.ind <- gm.prcomp(shape.data.both$symm.coords.ind[br.lands,,])

PCA.plots.brain.ind <- morphospace.visualization(pca.obj = pca.brain.ind, 
                                             text = shape.data.both$infos$acr,
                                             nb.pc = 3, 
                                             PCA.hull = T, 
                                             phylomorpho.hull = F,
                                             col.fill = shape.data.both$infos$Species, 
                                             name.fill = "Species")
PCA.plots.brain.ind$PC1_VS_PC2+theme(legend.position = "none")

data.plot <- PCA.plots.brain.ind$PC1_VS_PC2$data
data.hull <- PCA.plots.brain.ind$PC1_VS_PC2$plot_env$hull1

n.species.mult <- unique(as.character(data.hull$FILL))

for(i in 1:length(data.plot[[1]])){
  ind <- "wrong"
  for(j in 1:length(n.species.mult)){
    if(data.plot$FILL[i] == n.species.mult[j]){
      ind <- "good"
    }
  }
  if(ind == "wrong"){
    data.plot$FILL[i] <- NA
  }
}
rownames(data.plot) <- NULL
data.plot$name.plot <- c(rep(NA,11), "B. pinchonii", rep(NA,9),"B. subpalmata",rep(NA,23),"I. bellii",rep(NA,3),"N. maculosus",rep(NA,2),"N. viridescens",rep(NA,10),"P. vehiculum",rep(NA,9),"T. granulosa",rep(NA,9))

br_ind <- ggplot(data = data.plot, aes(x = X, y = Y,fill = FILL))+geom_polygon(data = data.hull, aes(x = PC1, y = PC2, fill = FILL), color = "black")+geom_point(size = 3, shape = 21) + geom_text_repel(aes(label = name.plot), fontface = "italic")+labs(x= PCA.plots.brain.ind$PC1_VS_PC2$labels$x, y = PCA.plots.brain.ind$PC1_VS_PC2$labels$y, title = "Brain shape")+theme_bw()+theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

## Endocast ####

### PCA object ####

pca.endo.ind <- gm.prcomp(shape.data.both$symm.coords.ind[en.lands,,])

### Symmetric shape ####

PCA.plots.endo.ind <- morphospace.visualization(pca.obj = pca.endo.ind, 
                                             text = shape.data.both$infos$acr,
                                             nb.pc = 3, 
                                             PCA.hull = T, 
                                             phylomorpho.hull = F,
                                             col.fill = shape.data.both$infos$Species, 
                                             name.fill = "Species")
PCA.plots.endo.ind$PC1_VS_PC2 +theme(legend.position = "none")

data.plot <- PCA.plots.endo.ind$PC1_VS_PC2$data
data.hull <- PCA.plots.endo.ind$PC1_VS_PC2$plot_env$hull1

n.species.mult <- unique(as.character(data.hull$FILL))

for(i in 1:length(data.plot[[1]])){
  ind <- "wrong"
  for(j in 1:length(n.species.mult)){
    if(data.plot$FILL[i] == n.species.mult[j]){
      ind <- "good"
    }
  }
  if(ind == "wrong"){
    data.plot$FILL[i] <- NA
  }
}
rownames(data.plot) <- NULL
data.plot$name.plot <- c(rep(NA,11), "B. pinchonii", rep(NA,9),"B. subpalmata",rep(NA,23),"I. bellii",rep(NA,3),"N. maculosus",rep(NA,2),"N. viridescens",rep(NA,10),"P. vehiculum",rep(NA,9),"T. granulosa",rep(NA,9))

en_ind <- ggplot(data = data.plot, aes(x = X, y = Y,fill = FILL))+geom_polygon(data = data.hull, aes(x = PC1, y = PC2, fill = FILL), color = "black")+geom_point(size = 3, shape = 21) + geom_text_repel(aes(label = name.plot), fontface = "italic")+labs(x= PCA.plots.endo.ind$PC1_VS_PC2$labels$x, y = PCA.plots.endo.ind$PC1_VS_PC2$labels$y, title = "Endocast shape")+theme_bw()+theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

## Both ####

### PCA object ####

pca.both.ind <- gm.prcomp(shape.data.both$symm.coords.ind)

### Symmetric shape ####

PCA.plots.both.ind <- morphospace.visualization(pca.obj = pca.both.ind, 
                                             text = shape.data.both$infos$acr,
                                             nb.pc = 3, 
                                             PCA.hull = T, 
                                             phylomorpho.hull = F,
                                             col.fill = shape.data.both$infos$Species, 
                                             name.fill = "Species")
PCA.plots.both.ind$PC1_VS_PC2 +theme(legend.position = "none")

data.plot <- PCA.plots.both.ind$PC1_VS_PC2$data
data.hull <- PCA.plots.both.ind$PC1_VS_PC2$plot_env$hull1

n.species.mult <- unique(as.character(data.hull$FILL))

for(i in 1:length(data.plot[[1]])){
  ind <- "wrong"
  for(j in 1:length(n.species.mult)){
    if(data.plot$FILL[i] == n.species.mult[j]){
      ind <- "good"
    }
  }
  if(ind == "wrong"){
    data.plot$FILL[i] <- NA
  }
}
rownames(data.plot) <- NULL
data.plot$name.plot <- c(rep(NA,11), "B. pinchonii", rep(NA,9),"B. subpalmata",rep(NA,23),"I. bellii",rep(NA,3),"N. maculosus",rep(NA,2),"N. viridescens",rep(NA,10),"P. vehiculum",rep(NA,9),"T. granulosa",rep(NA,9))

bo_ind <- ggplot(data = data.plot, aes(x = X, y = Y,fill = FILL))+geom_polygon(data = data.hull, aes(x = PC1, y = PC2, fill = FILL), color = "black")+geom_point(size = 3, shape = 21) + geom_text_repel(aes(label = name.plot), fontface = "italic")+labs(x= PCA.plots.both.ind$PC1_VS_PC2$labels$x, y = PCA.plots.both.ind$PC1_VS_PC2$labels$y, title = "Brain and endocast shape")+theme_bw()+theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

all_ind <- ggarrange(br_ind, en_ind, bo_ind, ncol= 1, nrow = 3)
all_ind

## Procrustes distance histograms ####

dimnames(shape.data.both$symm.coords.ind)[[3]] <- shape.data.both$infos$Species

PDist.data <- geomorph.data.frame(shape = shape.data.both$symm.coords.ind, group = shape.data.both$infos$Species)

shape.comp <- procD.comp(spec.shape = PDist.data$shape, groups = PDist.data$group, category.names = c("Among all specimens", "Within species", "Among species means"), title = "")
SHAPEC <- shape.comp$Plot+scale_fill_viridis_d(option = "C") + theme(legend.position = "bottom")

SHAPEC


## Saving plots ####

ggsave(all_ind, filename = "./output/PCA/intraspec_PCAs.pdf", height = 15, width = 7)

ggsave(SHAPEC, filename = "ProcD_histogram_ind.pdf", path = "./output/PCA", height = 5, width = 7)

# Among species PCA ####

## Brain ####

### PCA objects ####

pca.brain.sp <- gm.prcomp(shape.data.both$symm.coords.sp[br.lands,,])
pca.brain.sp.allo.rm <- gm.prcomp(shape.data.both$shape.allom.rm.sp[br.lands,,])

### Symmetric shape ####

PCA.plots.brain.sp <- morphospace.visualization(pca.obj = pca.brain.sp,
                                             text = shape.data.both$infos.sp$acr,
                                             nb.pc = 3, 
                                             PCA.hull = T,
                                             phylomorpho = T,
                                             tree = MCC,
                                             text_phylo = MCC$tip.label[match(sub(" ", "_", shape.data.both$infos.sp$Species),MCC$tip.label)],
                                             phylomorpho.hull = T, 
                                             col.fill = shape.data.both$infos.sp$Genus, 
                                             name.fill = "Genus")
# Preliminary plot
PCA.plots.brain.sp$PC1_VS_PC2 +theme(legend.position = "bottom")

# Data for plot

PCA.data.br <- data.frame(PC1 = PCA.plots.brain.sp$PC1_VS_PC2$plot_env$pca.data$PC1, PC2 = PCA.plots.brain.sp$PC1_VS_PC2$plot_env$pca.data$PC2, genus = PCA.plots.brain.sp$PC1_VS_PC2$plot_env$pca.data$FILL)

# Convex hull data

hull1 <- PCA.data.br %>%
  slice(chull(PC1, PC2))
hull1 <- hull1[NULL,]
for(j in 1:length(levels(PCA.plots.brain.sp$PC1_VS_PC2$plot_env$pca.data$FILL))){
  new <- filter(PCA.data.br, genus == levels(PCA.plots.brain.sp$PC1_VS_PC2$plot_env$pca.data$FILL)[j])
  
  if(length(new$PC1) >= 2){
    hull.new <- new %>%
      slice(chull(PC1, PC2))
    hull1 <- rbind(hull1,hull.new)
  }
  
}

# Text

text.genus <- c("Ambystoma",rep(NA,12),"Bolitoglossa",rep(NA,16),"Eurycea",rep(NA,5),"Hynobius",rep(NA,16),"Plethodon",rep(NA,5),"Pseudoeurycea",rep(NA,8),"Thorius",rep(NA,7))

# Phylogenetic tree

R.devices::suppressGraphics({
  lines.b <- phylomorphospace(MCC, matrix(c(PCA.data.br$PC1, PCA.data.br$PC2), ncol = 2, dimnames = list(MCC$tip.label[match(sub(" ", "_", shape.data.both$infos.sp$Species),MCC$tip.label)], c("PC1","PC2"))), nsteps = 200, ftype="off")
})

lines.PM <- data.frame(x = NULL,xend = NULL,y = NULL,yend = NULL)

for(h in 1:length(lines.b$edge[,1])){
  
  new <- data.frame(x = lines.b$xx[[lines.b$edge[h,1]]], xend = lines.b$xx[[lines.b$edge[h,2]]], y = lines.b$yy[[lines.b$edge[h,1]]], yend = lines.b$yy[[lines.b$edge[h,2]]])
  lines.PM <- rbind(lines.PM, new)
}

# Plot

neutral <- "grey"

brain_genus <- ggplot(PCA.data.br, aes(x = PC1, y = PC2, fill = genus)) + geom_polygon(data = hull1, aes(x = PC1, y = PC2, fill = genus), alpha = 0.5)+ geom_segment(data = lines.PM, aes(x=x,y=y,yend=yend,xend=xend),inherit.aes = F) + geom_point(size = 3, shape = 21) + scale_fill_manual(values = c("#444444",neutral,neutral,"#A95C2F",neutral,"#990099",neutral,"#109618",neutral,"#DC3912","#3366CC","#FF9900",rep(neutral,28)))+labs(title = "Brain", x = PCA.plots.brain.sp$PC1_VS_PC2$labels$x, y = PCA.plots.brain.sp$PC1_VS_PC2$labels$y) +geom_text_repel(aes(label = text.genus),inherit.aes = T) + theme_bw() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

### Size corrected shape ####

PCA.plots.brain.sp.allo.rm <- morphospace.visualization(pca.obj = pca.brain.sp.allo.rm, 
                                                text = shape.data.both$infos.sp$acr,
                                                nb.pc = 3, 
                                                PCA.hull = T, 
                                                phylomorpho.hull = F,
                                                col.fill = shape.data.both$infos.sp$Genus, 
                                                name.fill = "Genus")
# Preliminary plot
PCA.plots.brain.sp.allo.rm$PC1_VS_PC2 +theme(legend.position = "bottom")

# Data for plot
PCA.data.br.al <- data.frame(PC1 = PCA.plots.brain.sp.allo.rm$PC1_VS_PC2$plot_env$pca.data$PC1, PC2 = PCA.plots.brain.sp.allo.rm$PC1_VS_PC2$plot_env$pca.data$PC2, genus = PCA.plots.brain.sp.allo.rm$PC1_VS_PC2$plot_env$pca.data$FILL)

PCA.data.br.al <- PCA.data.br.al[order(PCA.data.br.al$genus),]

# Convex hull data

hull1 <- PCA.data.br.al %>%
  slice(chull(PC1, PC2))
hull1 <- hull1[NULL,]
for(j in 1:length(levels(PCA.plots.brain.sp.allo.rm$PC1_VS_PC2$plot_env$pca.data$FILL))){
  new <- filter(PCA.data.br.al, genus == levels(PCA.plots.brain.sp.allo.rm$PC1_VS_PC2$plot_env$pca.data$FILL)[j])
  
  if(length(new$PC1) >= 2){
    hull.new <- new %>%
      slice(chull(PC1, PC2))
    hull1 <- rbind(hull1,hull.new)
  }
  
}

# Text

text.genus <- c("Ambystoma",rep(NA,12),"Bolitoglossa",rep(NA,16),"Eurycea",rep(NA,6),"Hynobius",rep(NA,15),"Plethodon",rep(NA,6),"Pseudoeurycea",rep(NA,7),"Thorius",rep(NA,7))
PCA.data.br.al$text.genus <- text.genus

# Phylogenetic tree

R.devices::suppressGraphics({
  lines.b <- phylomorphospace(MCC, matrix(c(PCA.data.br.al$PC1, PCA.data.br.al$PC2), ncol = 2, dimnames = list(MCC$tip.label[match(sub(" ", "_", shape.data.both$infos.sp$Species),MCC$tip.label)], c("PC1","PC2"))), nsteps = 200, ftype="off")
})

lines.PM <- data.frame(x = NULL,xend = NULL,y = NULL,yend = NULL)

for(h in 1:length(lines.b$edge[,1])){
  
  new <- data.frame(x = lines.b$xx[[lines.b$edge[h,1]]], xend = lines.b$xx[[lines.b$edge[h,2]]], y = lines.b$yy[[lines.b$edge[h,1]]], yend = lines.b$yy[[lines.b$edge[h,2]]])
  lines.PM <- rbind(lines.PM, new)
}

# Plot

neutral <- "grey"

brain_genus.al <- ggplot(PCA.data.br.al, aes(x = PC1, y = PC2, fill = genus)) + geom_polygon(data = hull1, aes(x = PC1, y = PC2, fill = genus), alpha = 0.5)+ geom_segment(data = lines.PM, aes(x=x,y=y,yend=yend,xend=xend),inherit.aes = F) + geom_point(size = 3, shape = 21) + scale_fill_manual(values = c("#444444",neutral,neutral,"#A95C2F",neutral,"#990099",neutral,"#109618",neutral,"#DC3912","#3366CC","#FF9900",rep(neutral,28)))+labs(title = "Brain (size correction)", x = PCA.plots.brain.sp.allo.rm$PC1_VS_PC2$labels$x, y = PCA.plots.brain.sp.allo.rm$PC1_VS_PC2$labels$y) +geom_text_repel(aes(label = text.genus),inherit.aes = T) + theme_bw() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
brain_genus.al


## Endocast ####

### PCA objects ####

pca.endo.sp <- gm.prcomp(shape.data.both$symm.coords.sp[en.lands,,])
pca.endo.sp.allo.rm <- gm.prcomp(shape.data.both$shape.allom.rm.sp[en.lands,,])

### Symmetric shape ####

PCA.plots.endo.sp <- morphospace.visualization(pca.obj = pca.endo.sp,
                                            nb.pc = 3, 
                                            PCA.hull = T, 
                                            phylomorpho.hull = T,
                                            phylomorpho = T,
                                            tree = MCC,
                                            text_phylo = MCC$tip.label[match(sub(" ", "_", shape.data.both$infos.sp$Species),MCC$tip.label)], 
                                            col.fill = shape.data.both$infos.sp$Genus, 
                                            name.fill = "Family")

# Data for plot

PCA.data.en <- data.frame(PC1 = PCA.plots.endo.sp$PC1_VS_PC2$plot_env$pca.data$PC1, PC2 = PCA.plots.endo.sp$PC1_VS_PC2$plot_env$pca.data$PC2, genus = PCA.plots.endo.sp$PC1_VS_PC2$plot_env$pca.data$FILL)

# Convex hull data

hull1 <- PCA.data.en %>%
  slice(chull(PC1, PC2))
hull1 <- hull1[NULL,]
for(j in 1:length(levels(PCA.plots.endo.sp$PC1_VS_PC2$plot_env$pca.data$FILL))){
  new <- filter(PCA.data.en, genus == levels(PCA.plots.endo.sp$PC1_VS_PC2$plot_env$pca.data$FILL)[j])
  
  if(length(new$PC1) >= 2){
    hull.new <- new %>%
      slice(chull(PC1, PC2))
    hull1 <- rbind(hull1,hull.new)
  }
  
}

# Text

text.genus <- c("Ambystoma",rep(NA,12),"Bolitoglossa",rep(NA,16),"Eurycea",rep(NA,5),"Hynobius",rep(NA,16),"Plethodon",rep(NA,5),"Pseudoeurycea",rep(NA,8),"Thorius",rep(NA,7))

# Phylogenetic tree

R.devices::suppressGraphics({
  lines.b <- phylomorphospace(MCC, matrix(c(PCA.data.en$PC1, PCA.data.en$PC2), ncol = 2, dimnames = list(MCC$tip.label[match(sub(" ", "_", shape.data.both$infos.sp$Species),MCC$tip.label)], c("PC1","PC2"))), nsteps = 200, ftype="off")
})

lines.PM <- data.frame(x = NULL,xend = NULL,y = NULL,yend = NULL)

for(h in 1:length(lines.b$edge[,1])){
  
  new <- data.frame(x = lines.b$xx[[lines.b$edge[h,1]]], xend = lines.b$xx[[lines.b$edge[h,2]]], y = lines.b$yy[[lines.b$edge[h,1]]], yend = lines.b$yy[[lines.b$edge[h,2]]])
  lines.PM <- rbind(lines.PM, new)
}

# Plot

neutral <- "grey"

endo_genus <- ggplot(PCA.data.en, aes(x = PC1, y = PC2, fill = genus)) + geom_polygon(data = hull1, aes(x = PC1, y = PC2, fill = genus), alpha = 0.5)+ geom_segment(data = lines.PM, aes(x=x,y=y,yend=yend,xend=xend),inherit.aes = F) + geom_point(size = 3, shape = 21) + scale_fill_manual(values = c("#444444",neutral,neutral,"#A95C2F",neutral,"#990099",neutral,"#109618",neutral,"#DC3912","#3366CC","#FF9900",rep(neutral,28)))+labs(title = "Endocast", x = PCA.plots.endo.sp$PC1_VS_PC2$labels$x, y = PCA.plots.endo.sp$PC1_VS_PC2$labels$y) +geom_text_repel(aes(label = text.genus),inherit.aes = T) + theme_bw() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
endo_genus

### Size corrected shape ####

PCA.plots.endo.sp.allo.rm <- morphospace.visualization(pca.obj = pca.endo.sp.allo.rm, 
                                               text = shape.data.both$infos.sp$acr,
                                               nb.pc = 3, 
                                               PCA.hull = T, 
                                               phylomorpho.hull = F,
                                               col.fill = shape.data.both$infos.sp$Genus, 
                                               name.fill = "Genus")
# Preliminary plot
PCA.plots.endo.sp.allo.rm$PC1_VS_PC2 +theme(legend.position = "bottom")

# Data for plot

PCA.data.en.al <- data.frame(PC1 = PCA.plots.endo.sp.allo.rm$PC1_VS_PC2$plot_env$pca.data$PC1, PC2 = PCA.plots.endo.sp.allo.rm$PC1_VS_PC2$plot_env$pca.data$PC2, genus = PCA.plots.endo.sp.allo.rm$PC1_VS_PC2$plot_env$pca.data$FILL)

PCA.data.en.al <- PCA.data.en.al[order(PCA.data.en.al$genus),]

# Convex hull data

hull1 <- PCA.data.en.al %>%
  slice(chull(PC1, PC2))
hull1 <- hull1[NULL,]
for(j in 1:length(levels(PCA.plots.endo.sp.allo.rm$PC1_VS_PC2$plot_env$pca.data$FILL))){
  new <- filter(PCA.data.en.al, genus == levels(PCA.plots.endo.sp.allo.rm$PC1_VS_PC2$plot_env$pca.data$FILL)[j])
  
  if(length(new$PC1) >= 2){
    hull.new <- new %>%
      slice(chull(PC1, PC2))
    hull1 <- rbind(hull1,hull.new)
  }
  
}

# Text

text.genus <- c("Ambystoma",rep(NA,12),"Bolitoglossa",rep(NA,16),"Eurycea",rep(NA,6),"Hynobius",rep(NA,15),"Plethodon",rep(NA,6),"Pseudoeurycea",rep(NA,7),"Thorius",rep(NA,7))
PCA.data.en.al$text.genus <- text.genus

# Phylogenetic tree

R.devices::suppressGraphics({
  lines.b <- phylomorphospace(MCC, matrix(c(PCA.data.en.al$PC1, PCA.data.en.al$PC2), ncol = 2, dimnames = list(MCC$tip.label[match(sub(" ", "_", shape.data.both$infos.sp$Species),MCC$tip.label)], c("PC1","PC2"))), nsteps = 200, ftype="off")
})

lines.PM <- data.frame(x = NULL,xend = NULL,y = NULL,yend = NULL)

for(h in 1:length(lines.b$edge[,1])){
  
  new <- data.frame(x = lines.b$xx[[lines.b$edge[h,1]]], xend = lines.b$xx[[lines.b$edge[h,2]]], y = lines.b$yy[[lines.b$edge[h,1]]], yend = lines.b$yy[[lines.b$edge[h,2]]])
  lines.PM <- rbind(lines.PM, new)
}

# Plot

neutral <- "grey"

endo_genus.al <- ggplot(PCA.data.en.al, aes(x = PC1, y = PC2, fill = genus)) + geom_polygon(data = hull1, aes(x = PC1, y = PC2, fill = genus), alpha = 0.5)+ geom_segment(data = lines.PM, aes(x=x,y=y,yend=yend,xend=xend),inherit.aes = F) + geom_point(size = 3, shape = 21) + scale_fill_manual(values = c("#444444",neutral,neutral,"#A95C2F",neutral,"#990099",neutral,"#109618",neutral,"#DC3912","#3366CC","#FF9900",rep(neutral,28)))+labs(title = "Endocast (size correction)", x = PCA.plots.endo.sp.allo.rm$PC1_VS_PC2$labels$x, y = PCA.plots.endo.sp.allo.rm$PC1_VS_PC2$labels$y) +geom_text_repel(aes(label = text.genus),inherit.aes = T) + theme_bw() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
endo_genus.al

both_genus <- ggarrange(brain_genus,endo_genus,brain_genus.al,endo_genus.al, ncol = 2, nrow = 2)
both_genus

## Both ####

### PCA objects ####

pca.both.sp <- gm.prcomp(shape.data.both$symm.coords.sp)
pca.both.sp.allo.rm <- gm.prcomp(shape.data.both$shape.allom.rm.sp)



MCC$tip.label[match(sub(" ", "_", shape.data.both$infos.sp$Species),MCC$tip.label)]
match(sub(" ", "_", shape.data.both$infos.sp$Species),MCC$tip.label)

### Symmetric shape ####

PCA.plots.both.sp <- morphospace.visualization(pca.obj = pca.both.sp, 
                                            text = NULL,
                                            nb.pc = 3, 
                                            PCA.hull = T, tree = MCC,text_phylo = MCC$tip.label[match(sub(" ", "_", shape.data.both$infos.sp$Species),MCC$tip.label)],
                                            phylomorpho.hull = T,
                                            col.fill = factor(shape.data.both$infos.sp$Family), 
                                            name.fill = "Family",
                                            alpha = 0.3)

# Preliminary plot
PCA.plots.both.sp$PC1_VS_PC2 +theme(legend.position = "bottom")+scale_fill_manual(values=c("#681F79","#FECA04","#27BAE5","#71FFAA","#D5014B",rep("grey",5)))

### Size corrected shape ####

PCA.plots.both.sp.allo.rm <- morphospace.visualization(pca.obj = pca.both.sp.allo.rm, 
                                               text = NULL,
                                               nb.pc = 3, 
                                               PCA.hull = T, tree = MCC,text_phylo = MCC$tip.label[match(sub(" ", "_", shape.data.both$infos.sp$Species),MCC$tip.label)],
                                               phylomorpho.hull = T,
                                               col.fill = factor(shape.data.both$infos.sp$Family), 
                                               name.fill = "Family",
                                               alpha = 0.3)

# Preliminary plot
PCA.plots.both.sp.allo.rm$PC1_VS_PC2 +theme(legend.position = "bottom")+scale_shape_manual(values = c(21,7,8,9,22,23,10,11,24,25))+scale_fill_manual(values=c("#681F79","#FECA04","#27BAE5","#71FFAA","#D5014B",rep("grey",5)))

## Saving plots ####

ggsave(both_genus, filename = "phylomorph_size_nosize.pdf",path = "./output/PCA", height = 10, width = 14)

# Mantel tests ####

dist.br <- dist(pca.brain.sp$x)
dist.en <- dist(pca.endo.sp$x)

dist.br.allo.rm <- dist(pca.brain.sp.allo.rm$x)
dist.en.allo.rm <- dist(pca.endo.sp.allo.rm$x)

mantel(dist.br, dist.en, method = "pearson", permutations = 10000)
mantel(dist.br.allo.rm, dist.en.allo.rm, method = "pearson", permutations = 10000)
mantel(dist.br, dist.br.allo.rm, method = "pearson", permutations = 10000)
mantel(dist.en, dist.en.allo.rm, method = "pearson", permutations = 10000)


# Allometry

data.allom.br <- data.frame(PC1 = pca.brain.sp$x[,1], PC2 = pca.brain.sp$x[,2], CS = shape.data.both$sp.csize)

ggplot(data.allom.br, aes(x = PC1, y = PC2, fill = CS)) + geom_point(shape = 21, size = 3) + scale_fill_viridis_c()+ labs(x = "PC1: 34.6%", y = "PC2: 9.9%", fill = "Centroid size") + theme_bw() + theme(legend.position = "right")

data.allom.en <- data.frame(PC1 = pca.endo.sp$x[,1], PC2 = pca.endo.sp$x[,2], CS = shape.data.both$sp.csize)

ggplot(data.allom.en, aes(x = PC1, y = PC2, fill = CS)) + geom_point(shape = 21, size = 3) + scale_fill_viridis_c()+ labs(x = "PC1: 34.6%", y = "PC2: 9.9%", fill = "Centroid size") + theme_bw() + theme(legend.position = "right")

data.allom <- data.frame(PC1 = pca.both.sp$x[,1], PC2 = pca.both.sp$x[,2], CS = shape.data.both$sp.csize)

ggplot(data.allom, aes(x = PC1, y = PC2, fill = CS)) + geom_point(shape = 21, size = 3) + scale_fill_viridis_c()+ labs(x = "PC1: 34.6%", y = "PC2: 9.9%", fill = "Centroid size") + theme_bw() + theme(legend.position = "right")

# PACA ####

## Brain ####

### PCA objects ####

paca.brain.sp <- gm.prcomp(shape.data.both$symm.coords.sp[br.lands,,], align.to.phy = T, phy = MCC)

paca.brain.sp.allom.rm <- gm.prcomp(shape.data.both$shape.allom.rm.sp[br.lands,,], align.to.phy = T, phy = MCC)

### Symmetric shape ####

PACA.plots.brain.sp <- morphospace.visualization(pca.obj = paca.brain.sp, 
                                                       text = shape.data.both$infos.sp$acr,
                                                       tree = MCC,
                                                       text_phylo = sub(" ", "_", shape.data.both$infos.sp$Species),
                                                       nb.pc = 3, 
                                                       PCA.hull = T, 
                                                       phylomorpho.hull = T,
                                                       col.fill = shape.data.both$infos.sp$Genus, 
                                                       name.fill = "Genus")

# Preliminary plot
PACA.plots.brain.sp$PC1_VS_PC2 +theme(legend.position = "none")

# Data for plot

PACA.data.br <- data.frame(PC1 = PACA.plots.brain.sp$PC1_VS_PC2$plot_env$pca.data$PC1, PC2 = PACA.plots.brain.sp$PC1_VS_PC2$plot_env$pca.data$PC2, genus = PACA.plots.brain.sp$PC1_VS_PC2$plot_env$pca.data$FILL)

# Convex hull data

hull1 <- PACA.data.br %>%
  slice(chull(PC1, PC2))
hull1 <- hull1[NULL,]
for(j in 1:length(levels(PACA.plots.brain.sp$PC1_VS_PC2$plot_env$pca.data$FILL))){
  new <- filter(PACA.data.br, genus == levels(PACA.plots.brain.sp$PC1_VS_PC2$plot_env$pca.data$FILL)[j])
  
  if(length(new$PC1) >= 2){
    hull.new <- new %>%
      slice(chull(PC1, PC2))
    hull1 <- rbind(hull1,hull.new)
  }
  
}

# Text

text.genus <- c("Ambystoma",rep(NA,12),"Bolitoglossa",rep(NA,16),"Eurycea",rep(NA,5),"Hynobius",rep(NA,16),"Plethodon",rep(NA,5),"Pseudoeurycea",rep(NA,8),"Thorius",rep(NA,7))

# Phylogenetic tree

R.devices::suppressGraphics({
  lines.b <- phylomorphospace(MCC, matrix(c(PACA.data.br$PC1, PACA.data.br$PC2), ncol = 2, dimnames = list(MCC$tip.label[match(sub(" ", "_", shape.data.both$infos.sp$Species),MCC$tip.label)], c("PC1","PC2"))), nsteps = 200, ftype="off")
})

lines.PM <- data.frame(x = NULL,xend = NULL,y = NULL,yend = NULL)

for(h in 1:length(lines.b$edge[,1])){
  
  new <- data.frame(x = lines.b$xx[[lines.b$edge[h,1]]], xend = lines.b$xx[[lines.b$edge[h,2]]], y = lines.b$yy[[lines.b$edge[h,1]]], yend = lines.b$yy[[lines.b$edge[h,2]]])
  lines.PM <- rbind(lines.PM, new)
}

# Plot

neutral <- "grey"

PACA_brain_genus <- ggplot(PACA.data.br, aes(x = PC1, y = PC2, fill = genus)) + geom_polygon(data = hull1, aes(x = PC1, y = PC2, fill = genus), alpha = 0.5)+ geom_segment(data = lines.PM, aes(x=x,y=y,yend=yend,xend=xend),inherit.aes = F) + geom_point(size = 3, shape = 21) + scale_fill_manual(values = c("#444444",neutral,neutral,"#A95C2F",neutral,"#990099",neutral,"#109618",neutral,"#DC3912","#3366CC","#FF9900",rep(neutral,28)))+labs(title = "Brain", x = PACA.plots.brain.sp$PC1_VS_PC2$labels$x, y = PACA.plots.brain.sp$PC1_VS_PC2$labels$y) +geom_text_repel(aes(label = text.genus),inherit.aes = T) + theme_bw() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
PACA_brain_genus

### Size corrected shape ####

PACA.plots.brain.sp.allom.rm <- morphospace.visualization(pca.obj = paca.brain.sp.allom.rm, 
                                                text = shape.data.both$infos.sp$acr,
                                                tree = MCC,
                                                text_phylo = sub(" ", "_", shape.data.both$infos.sp$Species),
                                                nb.pc = 3, 
                                                PCA.hull = T, 
                                                phylomorpho.hull = T,
                                                col.fill = shape.data.both$infos.sp$Genus, 
                                                name.fill = "Genus")

# Preliminary plot
PACA.plots.brain.sp.allom.rm$PC1_VS_PC2 +theme(legend.position = "bottom")

# Data for plot

PACA.data.br.allom.rm <- data.frame(PC1 = PACA.plots.brain.sp.allom.rm$PC1_VS_PC2$plot_env$pca.data$PC1, PC2 = PACA.plots.brain.sp.allom.rm$PC1_VS_PC2$plot_env$pca.data$PC2, genus = PACA.plots.brain.sp.allom.rm$PC1_VS_PC2$plot_env$pca.data$FILL)

# Convex hull data

hull1 <- PACA.data.br.allom.rm %>%
  slice(chull(PC1, PC2))
hull1 <- hull1[NULL,]

for(j in 1:length(levels(PACA.plots.brain.sp.allom.rm$PC1_VS_PC2$plot_env$pca.data$FILL))){
  new <- filter(PACA.data.br.allom.rm, genus == levels(PACA.plots.brain.sp.allom.rm$PC1_VS_PC2$plot_env$pca.data$FILL)[j])
  
  if(length(new$PC1) >= 2){
    hull.new <- new %>%
      slice(chull(PC1, PC2))
    hull1 <- rbind(hull1,hull.new)
  }
  
}

# Text

text.genus <- c("Ambystoma",rep(NA,12),"Bolitoglossa",rep(NA,16),"Eurycea",rep(NA,5),"Hynobius",rep(NA,16),"Plethodon",rep(NA,5),"Pseudoeurycea",rep(NA,8),"Thorius",rep(NA,7))

# Phylogenetic tree

R.devices::suppressGraphics({
  lines.b <- phylomorphospace(MCC, matrix(c(PACA.data.br.allom.rm$PC1, PACA.data.br.allom.rm$PC2), ncol = 2, dimnames = list(MCC$tip.label[match(sub(" ", "_", shape.data.both$infos.sp$Species),MCC$tip.label)], c("PC1","PC2"))), nsteps = 200, ftype="off")
})

lines.PM <- data.frame(x = NULL,xend = NULL,y = NULL,yend = NULL)

for(h in 1:length(lines.b$edge[,1])){
  
  new <- data.frame(x = lines.b$xx[[lines.b$edge[h,1]]], xend = lines.b$xx[[lines.b$edge[h,2]]], y = lines.b$yy[[lines.b$edge[h,1]]], yend = lines.b$yy[[lines.b$edge[h,2]]])
  lines.PM <- rbind(lines.PM, new)
}

# Plot

neutral <- "grey"

PACA_brain_genus.allom.rm <- ggplot(PACA.data.br.allom.rm, aes(x = PC1, y = PC2, fill = genus)) + geom_polygon(data = hull1, aes(x = PC1, y = PC2, fill = genus), alpha = 0.5)+ geom_segment(data = lines.PM, aes(x=x,y=y,yend=yend,xend=xend),inherit.aes = F) + geom_point(size = 3, shape = 21) + scale_fill_manual(values = c("#444444",neutral,neutral,"#A95C2F",neutral,"#990099",neutral,"#109618",neutral,"#DC3912","#3366CC","#FF9900",rep(neutral,28)))+labs(title = "Brain (size correction)", x = PACA.plots.brain.sp.allom.rm$PC1_VS_PC2$labels$x, y = PACA.plots.brain.sp.allom.rm$PC1_VS_PC2$labels$y) +geom_text_repel(aes(label = text.genus),inherit.aes = T) + theme_bw() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
PACA_brain_genus.allom.rm


## Endocast ####

### PCA objects ####

paca.endo.sp <- gm.prcomp(shape.data.both$symm.coords.sp[en.lands,,], align.to.phy = T, phy = MCC)

paca.endo.sp.allom.rm <- gm.prcomp(shape.data.both$shape.allom.rm.sp[en.lands,,], align.to.phy = T, phy = MCC)

### Symmetric shape ####

PACA.plots.endo.sp <- morphospace.visualization(pca.obj = paca.endo.sp, 
                                                 text = shape.data.both$infos.sp$acr,
                                                 tree = MCC,
                                                 text_phylo = sub(" ", "_", shape.data.both$infos.sp$Species),
                                                 nb.pc = 3, 
                                                 PCA.hull = T, 
                                                 phylomorpho.hull = T,
                                                 col.fill = shape.data.both$infos.sp$Genus, 
                                                 name.fill = "Genus")

# Preliminary plot
PACA.plots.endo.sp$PC1_VS_PC2 +theme(legend.position = "none")

# Data for plot

PACA.data.en <- data.frame(PC1 = PACA.plots.endo.sp$PC1_VS_PC2$plot_env$pca.data$PC1, PC2 = PACA.plots.endo.sp$PC1_VS_PC2$plot_env$pca.data$PC2, genus = PACA.plots.endo.sp$PC1_VS_PC2$plot_env$pca.data$FILL)

# Convex hull data

hull1 <- PACA.data.en %>%
  slice(chull(PC1, PC2))
hull1 <- hull1[NULL,]
for(j in 1:length(levels(PACA.plots.endo.sp$PC1_VS_PC2$plot_env$pca.data$FILL))){
  new <- filter(PACA.data.en, genus == levels(PACA.plots.endo.sp$PC1_VS_PC2$plot_env$pca.data$FILL)[j])
  
  if(length(new$PC1) >= 2){
    hull.new <- new %>%
      slice(chull(PC1, PC2))
    hull1 <- rbind(hull1,hull.new)
  }
  
}

# Text

text.genus <- c("Ambystoma",rep(NA,12),"Bolitoglossa",rep(NA,16),"Eurycea",rep(NA,5),"Hynobius",rep(NA,16),"Plethodon",rep(NA,5),"Pseudoeurycea",rep(NA,8),"Thorius",rep(NA,7))

# Phylogenetic tree

R.devices::suppressGraphics({
  lines.b <- phylomorphospace(MCC, matrix(c(PACA.data.en$PC1, PACA.data.en$PC2), ncol = 2, dimnames = list(MCC$tip.label[match(sub(" ", "_", shape.data.both$infos.sp$Species),MCC$tip.label)], c("PC1","PC2"))), nsteps = 200, ftype="off")
})

lines.PM <- data.frame(x = NULL,xend = NULL,y = NULL,yend = NULL)

for(h in 1:length(lines.b$edge[,1])){
  
  new <- data.frame(x = lines.b$xx[[lines.b$edge[h,1]]], xend = lines.b$xx[[lines.b$edge[h,2]]], y = lines.b$yy[[lines.b$edge[h,1]]], yend = lines.b$yy[[lines.b$edge[h,2]]])
  lines.PM <- rbind(lines.PM, new)
}

# Plot

neutral <- "grey"

PACA_endo_genus <- ggplot(PACA.data.en, aes(x = PC1, y = PC2, fill = genus)) + geom_polygon(data = hull1, aes(x = PC1, y = PC2, fill = genus), alpha = 0.5)+ geom_segment(data = lines.PM, aes(x=x,y=y,yend=yend,xend=xend),inherit.aes = F) + geom_point(size = 3, shape = 21) + scale_fill_manual(values = c("#444444",neutral,neutral,"#A95C2F",neutral,"#990099",neutral,"#109618",neutral,"#DC3912","#3366CC","#FF9900",rep(neutral,28)))+labs(title = "Endocast", x = PACA.plots.endo.sp$PC1_VS_PC2$labels$x, y = PACA.plots.endo.sp$PC1_VS_PC2$labels$y) +geom_text_repel(aes(label = text.genus),inherit.aes = T) + theme_bw() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
PACA_endo_genus


all.genus <- ggarrange(brain_genus, endo_genus,PACA_brain_genus, PACA_endo_genus)
all.genus
ggsave(all.genus, filename = "PCA_PACA_br_en_genus.pdf",path = "./output/PCA", height = 10, width = 14)

### Size corrected shape ####

PACA.plots.endo.sp.allom.rm <- morphospace.visualization(pca.obj = paca.endo.sp.allom.rm, 
                                                          text = shape.data.both$infos.sp$acr,
                                                          tree = MCC,
                                                          text_phylo = sub(" ", "_", shape.data.both$infos.sp$Species),
                                                          nb.pc = 3, 
                                                          PCA.hull = T, 
                                                          phylomorpho.hull = T,
                                                          col.fill = shape.data.both$infos.sp$Genus, 
                                                          name.fill = "Genus")

# Preliminary plot
PACA.plots.endo.sp.allom.rm$PC1_VS_PC2 +theme(legend.position = "bottom")

# Data for plot

PACA.data.en.allom.rm <- data.frame(PC1 = PACA.plots.endo.sp.allom.rm$PC1_VS_PC2$plot_env$pca.data$PC1, PC2 = PACA.plots.endo.sp.allom.rm$PC1_VS_PC2$plot_env$pca.data$PC2, genus = PACA.plots.endo.sp.allom.rm$PC1_VS_PC2$plot_env$pca.data$FILL)

# Convex hull data

hull1 <- PACA.data.en.allom.rm %>%
  slice(chull(PC1, PC2))
hull1 <- hull1[NULL,]
for(j in 1:length(levels(PACA.plots.endo.sp.allom.rm$PC1_VS_PC2$plot_env$pca.data$FILL))){
  new <- filter(PACA.data.en.allom.rm, genus == levels(PACA.plots.endo.sp.allom.rm$PC1_VS_PC2$plot_env$pca.data$FILL)[j])
  
  if(length(new$PC1) >= 2){
    hull.new <- new %>%
      slice(chull(PC1, PC2))
    hull1 <- rbind(hull1,hull.new)
  }
  
}

# Text

text.genus <- c("Ambystoma",rep(NA,12),"Bolitoglossa",rep(NA,16),"Eurycea",rep(NA,5),"Hynobius",rep(NA,16),"Plethodon",rep(NA,5),"Pseudoeurycea",rep(NA,8),"Thorius",rep(NA,7))

# Phylogenetic tree

R.devices::suppressGraphics({
  lines.b <- phylomorphospace(MCC, matrix(c(PACA.data.en.allom.rm$PC1, PACA.data.en.allom.rm$PC2), ncol = 2, dimnames = list(MCC$tip.label[match(sub(" ", "_", shape.data.both$infos.sp$Species),MCC$tip.label)], c("PC1","PC2"))), nsteps = 200, ftype="off")
})

lines.PM <- data.frame(x = NULL,xend = NULL,y = NULL,yend = NULL)

for(h in 1:length(lines.b$edge[,1])){
  
  new <- data.frame(x = lines.b$xx[[lines.b$edge[h,1]]], xend = lines.b$xx[[lines.b$edge[h,2]]], y = lines.b$yy[[lines.b$edge[h,1]]], yend = lines.b$yy[[lines.b$edge[h,2]]])
  lines.PM <- rbind(lines.PM, new)
}

# Plot

neutral <- "grey"

PACA_endo_genus.allom.rm <- ggplot(PACA.data.en.allom.rm, aes(x = PC1, y = PC2, fill = genus)) + geom_polygon(data = hull1, aes(x = PC1, y = PC2, fill = genus), alpha = 0.5)+ geom_segment(data = lines.PM, aes(x=x,y=y,yend=yend,xend=xend),inherit.aes = F) + geom_point(size = 3, shape = 21) + scale_fill_manual(values = c("#444444",neutral,neutral,"#A95C2F",neutral,"#990099",neutral,"#109618",neutral,"#DC3912","#3366CC","#FF9900",rep(neutral,28)))+labs(title = "Endocast (size correction)", x = PACA.plots.endo.sp.allom.rm$PC1_VS_PC2$labels$x, y = PACA.plots.endo.sp.allom.rm$PC1_VS_PC2$labels$y) +geom_text_repel(aes(label = text.genus),inherit.aes = T) + theme_bw() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
PACA_endo_genus.allom.rm

PACA.genus <- ggarrange(PACA_brain_genus, PACA_endo_genus, PACA_brain_genus.allom.rm,PACA_endo_genus.allom.rm)
PACA.genus
ggsave(PACA.genus, filename = "PACA_size_nosize.pdf",path = "./output/PCA", height = 10, width = 14)


## Both ####

### PCA objects ####

paca.both.sp <- gm.prcomp(shape.data.both$symm.coords.sp, align.to.phy = T, phy = MCC)

paca.both.sp.allom.rm <- gm.prcomp(shape.data.both$shape.allom.rm.sp, align.to.phy = T, phy = MCC)

### Symmetric shape ####

PACA.plots.both.sp <- morphospace.visualization(pca.obj = paca.both.sp, 
                                                text = shape.data.both$infos.sp$acr,
                                                tree = MCC,
                                                text_phylo = sub(" ", "_", shape.data.both$infos.sp$Species),
                                                nb.pc = 3, 
                                                PCA.hull = T, 
                                                phylomorpho.hull = T,
                                                col.fill = shape.data.both$infos.sp$Genus, 
                                                name.fill = "Genus")

# Preliminary plot
PACA.plots.both.sp$PC1_VS_PC2 +theme(legend.position = "bottom")

### Size corrected shape ####

PACA.plots.both.sp.allom.rm <- morphospace.visualization(pca.obj = paca.both.sp.allom.rm, 
                                                         text = shape.data.both$infos.sp$acr,
                                                         tree = MCC,
                                                         text_phylo = sub(" ", "_", shape.data.both$infos.sp$Species),
                                                         nb.pc = 3, 
                                                         PCA.hull = T, 
                                                         phylomorpho.hull = T,
                                                         col.fill = shape.data.both$infos.sp$Genus, 
                                                         name.fill = "Genus")

# Preliminary plot
PACA.plots.both.sp.allom.rm$PC1_VS_PC2 +theme(legend.position = "bottom")

# 3D polygons ####

## New module definitions for visualization ####

# Brain

# Polygon colors

colors.region6 <- c("#115740","#0093DD","#AF8800","#E9052E","#D80BCF","#FD4513")
names(colors.region6) <- c("olfactory_bulb","telencephalon","thalamus","hypothalamus","optic_teg","rhombencephalon")

# New modules

modules_6_br <- land.d.brain$module_6

new.mod.br <- as.character(modules_6_br)
for(i in 1:length(modules_6_br)){
  
  if(modules_6_br[i] == "telencephalon"){
    new <- paste(modules_6_br[i], land.d.brain$Side[i], sep="_")
    new.mod.br[i] <- new
  } 
}
new.mod.br

# Endocast

# Polygon colors

colors.endo.bonesp2 <- c("#68BC28", "#4F2185","#FF9700","#A2A8AC","#00E7FF","#002143")
names(colors.endo.bonesp2) <- c("ant_parasphenoid", "frontal","occipito-otic","orbitoshpenoid","parietal","post_parasphenoid")

colors.region6 <- c("#115740","#0093DD","#AF8800","#E9052E","#D80BCF","#FD4513")
names(colors.region6) <- c("olfactory_bulb","telencephalon","thalamus","hypothalamus","optic_teg","rhombencephalon")

# New modules

modules_bones_p2 <- land.d.endo$mod_bones_para2
 
modules_mirror_6 <- land.d.endo$mirror_6

new.mod.en <- as.character(modules_mirror_6)
for(i in 1:length(modules_mirror_6)){
  
  if(modules_mirror_6[i] == "FF"){
    new <- paste(modules_mirror_6[i], land.d.endo$side[i], sep="_")
    new.mod.en[i] <- new
  } 
}
new.mod.en

## Brain ####

### Size-shape differences ####

# Preparing data

THPE <- as.data.frame(shape.data.both$symm.coords.sp[br.lands,,which(dimnames(shape.data.both$symm.coords.sp)[[3]] == "Thorius_pennatulus")])

CRAL <- as.data.frame(shape.data.both$symm.coords.sp[br.lands,,which(dimnames(shape.data.both$symm.coords.sp)[[3]] == "Cryptobranchus_alleganiensis")])
CRAL$V3 <- CRAL$V3+0.1

# Visualization

morphospace3d(data = THPE, 
              groups = new.mod.br, 
              show.axis = F, 
              alpha = 0.8, 
              names.group = F, 
              add.points = F, 
              color.pal = colors.region6[c(4,1,5,6,2,2,2,3)])

morphospace3d(data = CRAL, 
              groups = new.mod.br, 
              show.axis = F, 
              alpha = 0.8, 
              ini = F, 
              names.group = F, 
              color.pal = colors.region6[c(4,1,5,6,2,2,2,3)])

# Snapshot

snapshot3d(filename = "./output/Shape_deformation/3Dpol_brain_SS_front.png", height = 1000, width = 1200, webshot = F)


### Size corrected shape differences ####

#### PC1 ####

# Preparing data

MIN <- as.data.frame(pca.brain.sp.allo.rm$shapes$shapes.comp1$min)

MAX <- as.data.frame(pca.brain.sp.allo.rm$shapes$shapes.comp1$max)
MAX$`2` <- MAX$`2`+0.1

# Visualization

morphospace3d(data = MIN, 
              groups = new.mod.br, 
              show.axis = F, 
              alpha = 0.8, 
              names.group = F, 
              add.points = F, 
              color.pal = colors.region6[c(4,1,5,6,2,2,2,3)])

morphospace3d(data = MAX, 
              groups = new.mod.br, 
              show.axis = F, 
              alpha = 0.8, 
              ini = F, 
              names.group = F, 
              color.pal = colors.region6[c(4,1,5,6,2,2,2,3)])

# Snapshot

snapshot3d(filename = "./output/Shape_deformation/3Dpol_brain_allom_top.png", height = 1000, width = 1200, webshot = F)


#### PC2 ####

# Preparing data

MIN <- as.data.frame(pca.brain.sp.allo.rm$shapes$shapes.comp2$min)

MAX <- as.data.frame(pca.brain.sp.allo.rm$shapes$shapes.comp2$max)
MAX$`3` <- MAX$`3`+0.1

# Visualization

morphospace3d(data = MIN, 
              groups = new.mod.br, 
              show.axis = F, 
              alpha = 0.8, 
              names.group = F, 
              add.points = F, 
              color.pal = colors.region6[c(4,1,5,6,2,2,2,3)])

morphospace3d(data = MAX, 
              groups = new.mod.br, 
              show.axis = F, 
              alpha = 0.8, 
              ini = F, 
              names.group = F, 
              color.pal = colors.region6[c(4,1,5,6,2,2,2,3)])

# Snapshot

snapshot3d(filename = "./output/Shape_deformation/3Dpol_brain_allom2_front.png", height = 1000, width = 1200, webshot = F)

## Endocast ####

### Size-shape differences ####

# Preparing data

THPE <- as.data.frame(shape.data.both$symm.coords.sp[en.lands,,which(dimnames(shape.data.both$symm.coords.sp)[[3]] == "Thorius_pennatulus")])

CRAL <- as.data.frame(shape.data.both$symm.coords.sp[en.lands,,which(dimnames(shape.data.both$symm.coords.sp)[[3]] == "Cryptobranchus_alleganiensis")])
CRAL$V3 <- CRAL$V3+0.1

# Visualization

morphospace3d(data = THPE, 
              groups = new.mod.en, 
              show.axis = F, 
              alpha = 0.8, 
              names.group = F, 
              add.points = F, 
              color.pal = colors.region6[c(4,1,5,6,2,3)])

morphospace3d(data = CRAL, 
              groups = new.mod.en, 
              show.axis = F, 
              alpha = 0.8, 
              ini = F, 
              names.group = F, 
              color.pal = colors.region6[c(4,1,5,6,2,3)])

# Snapshot

snapshot3d(filename = "./output/Shape_deformation/3Dpol_endo_SS_front.png", height = 1000, width = 1200, webshot = F)


### Size corrected shape differences ####

#### PC1 ####

# Preparing data

MIN <- as.data.frame(pca.endo.sp.allo.rm$shapes$shapes.comp1$min)

MAX <- as.data.frame(pca.endo.sp.allo.rm$shapes$shapes.comp1$max)
MAX$`2` <- MAX$`2`+0.15

# Visualization

morphospace3d(data = MIN, 
              groups = new.mod.en, 
              show.axis = F, 
              alpha = 0.8, 
              names.group = F, 
              add.points = F, 
              color.pal = colors.region6[c(4,1,5,6,2,3)])

morphospace3d(data = MAX, 
              groups = new.mod.en, 
              show.axis = F, 
              alpha = 0.8, 
              ini = F, 
              names.group = F, 
              color.pal = colors.region6[c(4,1,5,6,2,3)])

# Snapshot

snapshot3d(filename = "./output/Shape_deformation/3Dpol_endo_allom_top.png", height = 1000, width = 1200, webshot = F)

#### PC2 ####

# Preparing data

MIN <- as.data.frame(pca.endo.sp.allo.rm$shapes$shapes.comp2$min)

MAX <- as.data.frame(pca.endo.sp.allo.rm$shapes$shapes.comp2$max)
MAX$`3` <- MAX$`3`+0.15

# Visualization

morphospace3d(data = MIN, 
              groups = new.mod.en, 
              show.axis = F, 
              alpha = 0.8, 
              names.group = F, 
              add.points = F, 
              color.pal = colors.region6[c(4,1,5,6,2,3)])

morphospace3d(data = MAX, 
              groups = new.mod.en, 
              show.axis = F, 
              alpha = 0.8, 
              ini = F, 
              names.group = F, 
              color.pal = colors.region6[c(4,1,5,6,2,3)])

# Snapshot

snapshot3d(filename = "./output/Shape_deformation/3Dpol_endo_allom2_front.png", height = 1000, width = 1200, webshot = F)


# Saving output ####

save(PHYPCA, tree.scaled, pca.brain.ind, pca.brain.sp, pca.brain.sp.allo.rm, pca.endo.ind, pca.endo.sp,pca.endo.sp.allo.rm, pca.both.ind, pca.both.sp, pca.both.sp.allo.rm, shape.data.both, MCC, nb.landmarks.both, nb.landmarks.brain, nb.landmarks.endo, file ="./output/RDA/PCA_shape_data.rda")

