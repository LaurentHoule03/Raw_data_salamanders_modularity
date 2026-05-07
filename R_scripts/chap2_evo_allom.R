
# Chapter 2

# Evolutionnary allometry - Step 4

source(paste0(getwd(),"/functions.R"))
load("./output/RDA/tree.rda")
load("./output/RDA/GPA_shape_data.rda")

shape.data.both$infos.sp$Species[which(shape.data.both$infos.sp$Species == "Isthmura bellii")] <- "Pseudoeurycea bellii"
shape.data.both$infos.sp$Species[which(shape.data.both$infos.sp$Species == "Aquiloeurycea cephalica")] <- "Pseudoeurycea cephalica"

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


# Creating a size corrected dataset ####

## Species ####

dimnames(shape.data.both$symm.coords.sp)[[3]] <- sub(" ", "_", shape.data.both$infos.sp$Species)

allometry.both.sp <- procD.lm(shape.data.both$symm.coords.sp ~ log(shape.data.both$sp.csize), iter = 9999)

summary(allometry.both.sp)


mean.s.sp <- mshape(shape.data.both$symm.coords.sp)

coords.res.sp <- arrayspecs(allometry.both.sp$residuals, p = nb.landmarks.both, k = 3)

allom.corr.coords.sp <- add.consensus.allom(resids = coords.res.sp, consensus = mean.s.sp, nb.land = nb.landmarks.both, nb.dims = 3)

shape.data.both$shape.allom.rm.sp <- allom.corr.coords.sp

allometry.phy.both.sp <- procD.pgls(shape.data.both$symm.coords.sp ~ log(shape.data.both$sp.csize), phy = tree.scaled)

summary(allometry.phy.both.sp)

plot.allom.sp <- plotAllometry(fit = allometry.phy.both.sp, size = shape.data.both$sp.csize, method = "RegScore")

plot.allom.sp <- plotAllometry(fit = allometry.phy.both.sp, size = shape.data.both$sp.csize, method = "PredLine")

plot.allom.sp <- plotAllometry(fit = allometry.phy.both.sp, size = shape.data.both$sp.csize, method = "CAC")

plot.allom.sp <- plotAllometry(fit = allometry.phy.both.sp, size = shape.data.both$sp.csize, method = "size.shape")


PCA.data.br <- data.frame(PC1 = plot.allom.sp$PC.points[,1], PC2 = plot.allom.sp$PC.points[,2], family = shape.data.both$infos.sp$Family)

hull1 <- PCA.data.br %>%
  slice(chull(PC1, PC2))
hull1 <- hull1[NULL,]
for(j in 1:length(levels(factor(PCA.data.br$family)))){
  new <- dplyr::filter(PCA.data.br, family == levels(factor(PCA.data.br$family))[j])
  
  if(length(new$PC1) >= 2){
    hull.new <- new %>%
      slice(chull(PC1, PC2))
    hull1 <- rbind(hull1,hull.new)
  }
  
}

R.devices::suppressGraphics({
  lines.b <- phylomorphospace(MCC, matrix(c(PCA.data.br$PC1, PCA.data.br$PC2), ncol = 2, dimnames = list(MCC$tip.label[match(sub(" ", "_", shape.data.both$infos.sp$Species),MCC$tip.label)], c("PC1","PC2"))), nsteps = 200, ftype="off")
})

lines.PM <- data.frame(x = NULL,xend = NULL,y = NULL,yend = NULL)

for(h in 1:length(lines.b$edge[,1])){
  
  new <- data.frame(x = lines.b$xx[[lines.b$edge[h,1]]], xend = lines.b$xx[[lines.b$edge[h,2]]], y = lines.b$yy[[lines.b$edge[h,1]]], yend = lines.b$yy[[lines.b$edge[h,2]]])
  lines.PM <- rbind(lines.PM, new)
}

neutral <- "grey"

both_fam <- ggplot(PCA.data.br, aes(x = PC1, y = PC2, fill = family)) + geom_polygon(data = hull1, aes(x = PC1, y = PC2, fill = family), alpha = 0.5)+ geom_segment(data = lines.PM, aes(x=x,y=y,yend=yend,xend=xend),inherit.aes = F) + geom_point(size = 3, shape = 21) +labs(title = "Brain", x = plot.allom.sp$plot.args$xlab, y = plot.allom.sp$plot.args$ylab) + theme_bw() + theme(legend.position = "bottom")
both_fam

plot.allom.sp <- plotAllometry(fit = allometry.phy.both.sp, size = shape.data.both$sp.csize, method = "RegScore")

data.plot.allom.sp <- data.frame(x = plot.allom.sp$plot_args$x, y = plot.allom.sp$plot_args$y, Family = shape.data.both$infos.sp$Family)

hull.tab.sp <- convex.hulls(data = data.plot.allom.sp)

Allometry.plot.sp <- ggplot(data.plot.allom.sp, aes(x = x, y = y, fill = Family, color = Family)) + geom_polygon(data = hull.tab.sp, alpha = 0.2, show.legend = FALSE)  + geom_point(size = 3, shape = 21, color = "black") + labs(x = "log(Centroid size)", y = "Regression score", fill = "Family") + theme_bw() + theme(legend.position = "bottom")
Allometry.plot.sp

plot.allom.sp <- plotAllometry(fit = allometry.phy.both.sp, size = shape.data.both$sp.csize, method = "RegScore")

data.plot.allom.sp <- data.frame(x = plot.allom.sp$plot_args$x, y = plot.allom.sp$plot_args$y, Family = shape.data.both$infos.sp$Family)

hull.tab.sp <- convex.hulls(data = data.plot.allom.sp)

Allometry.plot.sp.both <- ggplot(data.plot.allom.sp, aes(x = x, y = y, fill = Family, shape = Family)) + geom_polygon(data = hull.tab.sp, alpha = 0.2)  + scale_fill_manual(values = c(viridis(100)[1], "grey","grey","grey",viridis(100)[25], viridis(100)[50], "grey","grey",viridis(100)[75], viridis(100)[100])) + geom_point(size = 3, color = "black")+ scale_shape_manual(values = c(21,10,11,9,22,23,8,7,24,25)) + labs(x = "log(Centroid size)", y = "Regression score", fill = "Family", shape = "Family", title = "Brain and endocast") + theme_bw() + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
Allometry.plot.sp.both

# Brain 

allometry.phy.br.sp <- procD.pgls(shape.data.both$symm.coords.sp[br.lands,,] ~ log(shape.data.both$sp.csize), phy = tree.scaled.br)

summary(allometry.phy.br.sp)

plot.allom.sp <- plotAllometry(fit = allometry.phy.br.sp, size = shape.data.both$sp.csize, method = "RegScore")

plot.allom.sp <- plotAllometry(fit = allometry.phy.br.sp, size = shape.data.both$sp.csize, method = "PredLine")

plot.allom.sp <- plotAllometry(fit = allometry.phy.br.sp, size = shape.data.both$sp.csize, method = "CAC")

plot.allom.sp <- plotAllometry(fit = allometry.phy.br.sp, size = shape.data.both$sp.csize, method = "size.shape")


PCA.data.br <- data.frame(PC1 = plot.allom.sp$PC.points[,1], PC2 = plot.allom.sp$PC.points[,2], family = shape.data.both$infos.sp$Family)

hull1 <- PCA.data.br %>%
  slice(chull(PC1, PC2))
hull1 <- hull1[NULL,]
for(j in 1:length(levels(factor(PCA.data.br$family)))){
  new <- dplyr::filter(PCA.data.br, family == levels(factor(PCA.data.br$family))[j])
  
  if(length(new$PC1) >= 2){
    hull.new <- new %>%
      slice(chull(PC1, PC2))
    hull1 <- rbind(hull1,hull.new)
  }
  
}

R.devices::suppressGraphics({
  lines.b <- phylomorphospace(MCC, matrix(c(PCA.data.br$PC1, PCA.data.br$PC2), ncol = 2, dimnames = list(MCC$tip.label[match(sub(" ", "_", shape.data.both$infos.sp$Species),MCC$tip.label)], c("PC1","PC2"))), nsteps = 200, ftype="off")
})

lines.PM <- data.frame(x = NULL,xend = NULL,y = NULL,yend = NULL)

for(h in 1:length(lines.b$edge[,1])){
  
  new <- data.frame(x = lines.b$xx[[lines.b$edge[h,1]]], xend = lines.b$xx[[lines.b$edge[h,2]]], y = lines.b$yy[[lines.b$edge[h,1]]], yend = lines.b$yy[[lines.b$edge[h,2]]])
  lines.PM <- rbind(lines.PM, new)
}

neutral <- "grey"

both_fam <- ggplot(PCA.data.br, aes(x = PC1, y = PC2, fill = family)) + geom_polygon(data = hull1, aes(x = PC1, y = PC2, fill = family), alpha = 0.5)+ geom_segment(data = lines.PM, aes(x=x,y=y,yend=yend,xend=xend),inherit.aes = F) + geom_point(size = 3, shape = 21) +labs(title = "Brain", x = plot.allom.sp$plot.args$xlab, y = plot.allom.sp$plot.args$ylab) + theme_bw() + theme(legend.position = "bottom")
both_fam

plot.allom.sp <- plotAllometry(fit = allometry.phy.br.sp, size = shape.data.both$sp.csize, method = "RegScore")

data.plot.allom.sp <- data.frame(x = plot.allom.sp$plot_args$x, y = plot.allom.sp$plot_args$y, Family = shape.data.both$infos.sp$Family)

hull.tab.sp <- convex.hulls(data = data.plot.allom.sp)

Allometry.plot.sp <- ggplot(data.plot.allom.sp, aes(x = x, y = y, fill = Family, color = Family)) + geom_polygon(data = hull.tab.sp, alpha = 0.2, show.legend = FALSE)  + geom_point(size = 3, shape = 21, color = "black") + labs(x = "log(Centroid size)", y = "Regression score", fill = "Family") + theme_bw() + theme(legend.position = "bottom")
Allometry.plot.sp

plot.allom.sp <- plotAllometry(fit = allometry.phy.br.sp, size = shape.data.both$sp.csize, method = "RegScore")

data.plot.allom.sp <- data.frame(x = plot.allom.sp$plot_args$x, y = plot.allom.sp$plot_args$y, Family = shape.data.both$infos.sp$Family)

hull.tab.sp <- convex.hulls(data = data.plot.allom.sp)

Allometry.plot.sp.br <- ggplot(data.plot.allom.sp, aes(x = x, y = y, fill = Family, shape = Family)) + geom_polygon(data = hull.tab.sp, alpha = 0.2)  + scale_fill_manual(values = c(viridis(100)[1], "grey","grey","grey",viridis(100)[25], viridis(100)[50], "grey","grey",viridis(100)[75], viridis(100)[100])) + geom_point(size = 3, color = "black")+ scale_shape_manual(values = c(21,10,11,9,22,23,8,7,24,25)) + labs(x = "log(Centroid size)", y = "Regression score", fill = "Family", shape = "Family", title = "Brain") + theme_bw() + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
Allometry.plot.sp.br

# Endocast

allometry.phy.en.sp <- procD.pgls(shape.data.both$symm.coords.sp[en.lands,,] ~ log(shape.data.both$sp.csize), phy = tree.scaled.en)

summary(allometry.phy.en.sp)

plot.allom.sp <- plotAllometry(fit = allometry.phy.en.sp, size = shape.data.both$sp.csize, method = "RegScore")

plot.allom.sp <- plotAllometry(fit = allometry.phy.en.sp, size = shape.data.both$sp.csize, method = "PredLine")

plot.allom.sp <- plotAllometry(fit = allometry.phy.en.sp, size = shape.data.both$sp.csize, method = "CAC")

plot.allom.sp <- plotAllometry(fit = allometry.phy.en.sp, size = shape.data.both$sp.csize, method = "size.shape")


PCA.data.br <- data.frame(PC1 = plot.allom.sp$PC.points[,1], PC2 = plot.allom.sp$PC.points[,2], family = shape.data.both$infos.sp$Family)

hull1 <- PCA.data.br %>%
  slice(chull(PC1, PC2))
hull1 <- hull1[NULL,]
for(j in 1:length(levels(factor(PCA.data.br$family)))){
  new <- dplyr::filter(PCA.data.br, family == levels(factor(PCA.data.br$family))[j])
  
  if(length(new$PC1) >= 2){
    hull.new <- new %>%
      slice(chull(PC1, PC2))
    hull1 <- rbind(hull1,hull.new)
  }
  
}

R.devices::suppressGraphics({
  lines.b <- phylomorphospace(MCC, matrix(c(PCA.data.br$PC1, PCA.data.br$PC2), ncol = 2, dimnames = list(MCC$tip.label[match(sub(" ", "_", shape.data.both$infos.sp$Species),MCC$tip.label)], c("PC1","PC2"))), nsteps = 200, ftype="off")
})

lines.PM <- data.frame(x = NULL,xend = NULL,y = NULL,yend = NULL)

for(h in 1:length(lines.b$edge[,1])){
  
  new <- data.frame(x = lines.b$xx[[lines.b$edge[h,1]]], xend = lines.b$xx[[lines.b$edge[h,2]]], y = lines.b$yy[[lines.b$edge[h,1]]], yend = lines.b$yy[[lines.b$edge[h,2]]])
  lines.PM <- rbind(lines.PM, new)
}

neutral <- "grey"

both_fam <- ggplot(PCA.data.br, aes(x = PC1, y = PC2, fill = family)) + geom_polygon(data = hull1, aes(x = PC1, y = PC2, fill = family), alpha = 0.5)+ geom_segment(data = lines.PM, aes(x=x,y=y,yend=yend,xend=xend),inherit.aes = F) + geom_point(size = 3, shape = 21) +labs(title = "Brain", x = plot.allom.sp$plot.args$xlab, y = plot.allom.sp$plot.args$ylab) + geom_text(aes(label = shape.data.both$infos.sp$acr)) + theme_bw() + theme(legend.position = "bottom")
both_fam

plot.allom.sp <- plotAllometry(fit = allometry.phy.en.sp, size = shape.data.both$sp.csize, method = "RegScore")

data.plot.allom.sp <- data.frame(x = plot.allom.sp$plot_args$x, y = plot.allom.sp$plot_args$y, Family = shape.data.both$infos.sp$Family)

hull.tab.sp <- convex.hulls(data = data.plot.allom.sp)

Allometry.plot.sp <- ggplot(data.plot.allom.sp, aes(x = x, y = y, fill = Family, color = Family)) + geom_polygon(data = hull.tab.sp, alpha = 0.2, show.legend = FALSE)  + geom_point(size = 3, shape = 21, color = "black") + labs(x = "log(Centroid size)", y = "Regression score", fill = "Family") + theme_bw() + theme(legend.position = "bottom")
Allometry.plot.sp

plot.allom.sp.en <- plotAllometry(fit = allometry.phy.en.sp, size = shape.data.both$sp.csize, method = "RegScore")
plot.allom.sp$plot.args$y
data.plot.allom.sp.en <- data.frame(x = plot.allom.sp.en$plot_args$x, y = plot.allom.sp.en$plot_args$y, Family = shape.data.both$infos.sp$Family)

hull.tab.sp.en <- convex.hulls(data = data.plot.allom.sp.en)

Allometry.plot.sp.en <- ggplot(data.plot.allom.sp.en, aes(x = x, y = y, shape = Family, fill = Family)) + geom_polygon(data = hull.tab.sp.en, alpha = 0.2)  + geom_point(size = 3) + scale_shape_manual(values = c(21,10,11,9,22,23,8,7,24,25)) + labs(x = "log(Centroid size)", y = "Regression score", fill = "Family", shape = "Family", title = "Endocast") + scale_fill_manual(values = c(viridis(100)[1], "grey","grey","grey",viridis(100)[25], viridis(100)[50], "grey","grey",viridis(100)[75], viridis(100)[100])) + theme_bw() + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
Allometry.plot.sp.en

both.CAC <- ggarrange(Allometry.plot.sp.br, Allometry.plot.sp.en, Allometry.plot.sp.both, common.legend = T, legend = "bottom",ncol= 1, nrow = 3)
both.CAC
ggsave(both.CAC, filename = "RS_br_en.pdf", path = "./output/Allometry", height = 15, width = 7)

# Saving size corrected data ####

save(shape.data.both, lm.pairs.both, lm.pairs.brain, lm.pairs.endo, nb.landmarks.both, nb.landmarks.brain, nb.landmarks.endo, br.lands, en.lands, file ="./output/RDA/allom_rm_data.rda")


