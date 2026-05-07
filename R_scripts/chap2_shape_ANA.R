
# Chapter 2

# Phylogenetic signal - Step 7

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
library(writexl)
library(remotes)

# Importation ####
load("F:/MicroCT_Scan/xx_landmark_brain_salamanders_xx/Chapiters/doc-chap2/output/RDA/tree.rda")
load("F:/MicroCT_Scan/xx_landmark_brain_salamanders_xx/Chapiters/doc-chap2/output/RDA/modularity.rda")
load("F:/MicroCT_Scan/xx_landmark_brain_salamanders_xx/Chapiters/doc-chap2/output/RDA/PCA_shape_data.rda")

# Phylogenetic signal ####

## Shape ####

dimnames(shape.data.both$symm.coords.sp)[[3]] <- sub(" ", "_", shape.data.both$infos.sp$Species)

phy.shape.brain <- physignal.z(A = shape.data.both$symm.coords.sp[br.lands,,], phy = tree.scaled.br, lambda = "burn", PAC.no = 5)
plot(phy.shape.brain$PACA$sdev)
phy.shape.brain$PACA
summary(phy.shape.brain)

phy.shape.brain.allo.rm <- physignal.z(A = shape.data.both$shape.allom.rm.sp[br.lands,,], phy = tree.scaled.br, lambda = "burn", PAC.no = 7)
phy.shape.brain.allo.rm$PACA
plot(phy.shape.brain.allo.rm$PACA$sdev)
summary(phy.shape.brain.allo.rm)

phy.shape.endo <- physignal.z(A = shape.data.both$symm.coords.sp[en.lands,,], phy = tree.scaled.en, lambda = "burn", PAC.no = 5)
phy.shape.endo$PACA
plot(phy.shape.endo$PACA$sdev)
summary(phy.shape.endo)

phy.shape.endo.allo.rm <- physignal.z(A = shape.data.both$shape.allom.rm.sp[en.lands,,], phy = tree.scaled.en, lambda = "burn", PAC.no = 8)
phy.shape.endo.allo.rm$PACA
plot(phy.shape.endo.allo.rm$PACA$sdev)
summary(phy.shape.endo.allo.rm)

phy.shape.both <- physignal.z(A = shape.data.both$symm.coords.sp, phy = tree.scaled, PAC.no = 6)
phy.shape.both$PACA
plot(phy.shape.both$PACA$sdev)
phy.shape.both


phy.shape.both.allo.rm <- physignal.z(A = shape.data.both$shape.allom.rm.sp, phy = tree.scaled, lambda = "burn", PAC.no = 8)
phy.shape.both.allo.rm$PACA
plot(phy.shape.both.allo.rm$PACA$sdev)
summary(phy.shape.both.allo.rm)


cmp.list <- list(brain = phy.shape.brain, brain_allom = phy.shape.brain.allo.rm, endo = phy.shape.endo, endo_allom = phy.shape.endo.allo.rm)

comp.j <- compare.physignal.z(cmp.list)
summary(comp.j)
comp.j$pairwise.z
comp.j$pairwise.P

cmp.list <- list(symm = phy.shape.both, allom = phy.shape.both.allo.rm)

comp.j <- compare.physignal.z(cmp.list)
summary(comp.j)
comp.j$pairwise.z
comp.j$pairwise.P

## Centroid size ####

names(shape.data.both$sp.csize) <- sub(" ", "_", shape.data.both$infos.sp$Species)

phy.size.both <- physignal.z(A = shape.data.both$sp.csize, phy = tree.scaled) # geomorph 4.0.8

summary(phy.size.both)

## Phylogenetic signal per module ####
# HYP 20

PZ <- physignal.mult(modules = mod.hyp.both$`hyp_6_br mirror_6`, coords = shape.data.both$shape.allom.rm.sp, tree = tree.scaled, PAC.no = "find", threshold = 0.99)

PZ$Results.data.frame
PZ$Comparison

Ps <- as.data.frame(PZ$Comparison$pairwise.P)
Ps$modules <- rownames(Ps)
Zs <- as.data.frame(PZ$Comparison$pairwise.z)
Zs$modules <- rownames(Zs)

#write_xlsx(Ps, "./output/Phylogenetic_signal/phy_pairw_p.xlsx")
#write_xlsx(Zs, "./output/Phylogenetic_signal/phy_pairw_Z.xlsx")
#write_xlsx(PZ$Results.data.frame, "./output/Phylogenetic_signal/phy_res_mod.xlsx")
