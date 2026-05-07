
# Chapter 2

# Phylogenetic analyses - Step 3

# Laurent Houle

# Importing functions ####

source(paste0(getwd(),"/functions.R"))
source(paste0(getwd(),"/Phylogeny/btrtools.R"))

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
library(pbapply)
library(coda)
library(ggpubr)
library(ggridges)
library(scales)
library(writexl)

# Importation ####

load("./output/RDA/GPA_shape_data.rda")

##filelist <- list.files(path = "F:/MicroCT_Scan/xx_landmark_brain_salamanders_xx/Chapiters/doc-chap2/Phylogeny", pattern = ".trees", recursive = F, full.names = T)

##Ptrees <- read_all_trees(filenames = filelist, keep = F)

# Write nexus tree for TreeAnnotator ####

##write.nexus(Ptrees[[1]], file = "Post_JP.trees")

# Post_JP.trees must be imported in TreeAnnotator.exe from the BEAST series
# MCC tree must be produced using Maximum clade credibility tree and common ancestor heights

# Import MCC tree from TreeAnnotator ####

MCC_all <- read.nexus(file = "./MCC.nexus")

# Pruning the tree ####

find.species.trees(tree = MCC_all, find = "cephalica")
find.species.trees(tree = MCC_all, find = "bellii")

shape.data.both$infos.sp$Species[which(shape.data.both$infos.sp$Species == "Isthmura bellii")] <- "Pseudoeurycea bellii"
shape.data.both$infos.sp$Species[which(shape.data.both$infos.sp$Species == "Aquiloeurycea cephalica")] <- "Pseudoeurycea cephalica"

tips.both <- check.tips(tree = MCC_all, sub(" ", "_", shape.data.both$infos.sp$Species))

MCC <- keep.tip(MCC_all, tip = tips.both)

#write.nexus(MCC, file = "MCC_pruned.nexus")

# Phylogenetic PCA ####

dimnames(shape.data.both$symm.coords.sp)[[3]] <- sub(" ", "_", shape.data.both$infos.sp$Species)

PHYPCA <- gm.prcomp(shape.data.both$symm.coords.sp, align.to.phy = F, phy = MCC, GLS = T)
PHYPCA

PHYPCA.plot <- morphospace.visualization(pca.obj = PHYPCA, 
                                                 text = shape.data.both$infos.sp$acr,
                                                 nb.pc = 3, 
                                                 PCA.hull = T, 
                                                 phylomorpho.hull = F,
                                                 col.fill = shape.data.both$infos.sp$Life_cycle, 
                                                 name.fill = "Species")
PHYPCA.plot$PC1_VS_PC2 +theme(legend.position = "none")

PHYPCA.br <- gm.prcomp(shape.data.both$symm.coords.sp[1:181,,], align.to.phy = F, phy = MCC, GLS = T)
PHYPCA.br

PHYPCA.plot.br <- morphospace.visualization(pca.obj = PHYPCA.br, 
                                         text = shape.data.both$infos.sp$acr,
                                         nb.pc = 3, 
                                         PCA.hull = T, 
                                         phylomorpho.hull = F,
                                         col.fill = shape.data.both$infos.sp$Life_cycle, 
                                         name.fill = "Species")
PHYPCA.plot.br$PC1_VS_PC2 +theme(legend.position = "none")

PHYPCA.en <- gm.prcomp(shape.data.both$symm.coords.sp[182:258,,], align.to.phy = F, phy = MCC, GLS = T)
PHYPCA.en

PHYPCA.plot.en <- morphospace.visualization(pca.obj = PHYPCA.en, 
                                            text = shape.data.both$infos.sp$acr,
                                            nb.pc = 3, 
                                            PCA.hull = T, 
                                            phylomorpho.hull = F,
                                            col.fill = shape.data.both$infos.sp$Life_cycle, 
                                            name.fill = "Species")
PHYPCA.plot.en$PC1_VS_PC2 +theme(legend.position = "none")

# Creating the traits data file for BayesTraits ####

traits <- PHYPCA$x[,1:29]
fc <- sub(" ", "_", shape.data.both$infos.sp$Species)
traits <- cbind(fc, traits)
#write.table(traits, file = "PCS.txt", append = FALSE, sep = " ", dec = ".",row.names = F, col.names = F, quote = F)

traits.br <- PHYPCA.br$x[,1:27]
fc.br <- sub(" ", "_", shape.data.both$infos.sp$Species)
traits.br <- cbind(fc.br, traits.br)
#write.table(traits.br, file = "PCS_br.txt", append = FALSE, sep = " ", dec = ".",row.names = F, col.names = F, quote = F)

traits.en <- PHYPCA.en$x[,1:19]
fc.en <- sub(" ", "_", shape.data.both$infos.sp$Species)
traits.en <- cbind(fc.en, traits.en)
#write.table(traits.en, file = "PCS_en.txt", append = FALSE, sep = " ", dec = ".",row.names = F, col.names = F, quote = F)

# Four rj MCMC chains must be run from the BayesTraitsV4 software using PCS.txt, MCC_pruned.nexus, and Script_BM_rj_MCMC_bayestrait.txt files

# Calculating branch specific rates ###

## Both ####

results1 <- rjpp(rjlog = "./Phylogeny/CHAIN1_BT/PCS.txt.VarRates.txt", rjtrees = "./Phylogeny/CHAIN1_BT/PCS.txt.Output.trees", tree = "./Phylogeny/CHAIN1_BT/MCC_pruned.nexus")
results2 <- rjpp(rjlog = "./Phylogeny/CHAIN2_BT/PCS.txt.VarRates.txt", rjtrees = "./Phylogeny/CHAIN2_BT/PCS.txt.Output.trees", tree = "./Phylogeny/CHAIN2_BT/MCC_pruned.nexus")
results3 <- rjpp(rjlog = "./Phylogeny/CHAIN3_BT/PCS.txt.VarRates.txt", rjtrees = "./Phylogeny/CHAIN3_BT/PCS.txt.Output.trees", tree = "./Phylogeny/CHAIN3_BT/MCC_pruned.nexus")
results4 <- rjpp(rjlog = "./Phylogeny/CHAIN4_BT/PCS.txt.VarRates.txt", rjtrees = "./Phylogeny/CHAIN4_BT/PCS.txt.Output.trees", tree = "./Phylogeny/CHAIN4_BT/MCC_pruned.nexus")

EB1 <- results1$meantree$edge.length
EB2 <- results2$meantree$edge.length
EB3 <- results3$meantree$edge.length
EB4 <- results4$meantree$edge.length
EBS <- data.frame(V1 = EB1,V2 = EB2,V3 = EB3,V4 = EB4)
EB_m <- apply(EBS, 1, mean)
tree.scaled <- results1$meantree
tree.scaled$edge.length <- EB_m
plot(tree.scaled)

## Brain ####

results1.br <- rjpp(rjlog = "./Phylogeny/Brain/CHAIN1_BT/PCS_br.txt.VarRates.txt", rjtrees = "./Phylogeny/Brain/CHAIN1_BT/PCS_br.txt.Output.trees", tree = "./Phylogeny/Brain/CHAIN1_BT/MCC_pruned.nexus")
results2.br <- rjpp(rjlog = "./Phylogeny/CHAIN2_BT/PCS.txt.VarRates.txt", rjtrees = "./Phylogeny/CHAIN2_BT/PCS.txt.Output.trees", tree = "./Phylogeny/CHAIN2_BT/MCC_pruned.nexus")
results3.br <- rjpp(rjlog = "./Phylogeny/CHAIN3_BT/PCS.txt.VarRates.txt", rjtrees = "./Phylogeny/CHAIN3_BT/PCS.txt.Output.trees", tree = "./Phylogeny/CHAIN3_BT/MCC_pruned.nexus")
results4.br <- rjpp(rjlog = "./Phylogeny/CHAIN4_BT/PCS.txt.VarRates.txt", rjtrees = "./Phylogeny/CHAIN4_BT/PCS.txt.Output.trees", tree = "./Phylogeny/CHAIN4_BT/MCC_pruned.nexus")

EB1.br <- results1.br$meantree$edge.length
EB2.br <- results2.br$meantree$edge.length
EB3.br <- results3.br$meantree$edge.length
EB4.br <- results4.br$meantree$edge.length
EBS.br <- data.frame(V1 = EB1.br,V2 = EB2.br,V3 = EB3.br,V4 = EB4.br)
EB_m.br <- apply(EBS.br, 1, mean)
tree.scaled.br <- results1$meantree
tree.scaled.br$edge.length <- EB_m.br
plot(tree.scaled.br)

## Endocast ####

results1.en <- rjpp(rjlog = "./Phylogeny/Endocast/CHAIN1_BT/PCS_en.txt.VarRates.txt", rjtrees = "./Phylogeny/Endocast/CHAIN1_BT/PCS_en.txt.Output.trees", tree = "./Phylogeny/Endocast/CHAIN1_BT/MCC_pruned.nexus")
results2.en <- rjpp(rjlog = "./Phylogeny/Endocast/CHAIN2_BT/PCS_en.txt.VarRates.txt", rjtrees = "./Phylogeny/Endocast/CHAIN2_BT/PCS_en.txt.Output.trees", tree = "./Phylogeny/Endocast/CHAIN2_BT/MCC_pruned.nexus")
results3.en <- rjpp(rjlog = "./Phylogeny/Endocast/CHAIN3_BT/PCS_en.txt.VarRates.txt", rjtrees = "./Phylogeny/Endocast/CHAIN3_BT/PCS_en.txt.Output.trees", tree = "./Phylogeny/Endocast/CHAIN3_BT/MCC_pruned.nexus")
results4.en <- rjpp(rjlog = "./Phylogeny/Endocast/CHAIN4_BT/PCS_en.txt.VarRates.txt", rjtrees = "./Phylogeny/Endocast/CHAIN4_BT/PCS_en.txt.Output.trees", tree = "./Phylogeny/Endocast/CHAIN4_BT/MCC_pruned.nexus")

EB1.en <- results1.en$meantree$edge.length
EB2.en <- results2.en$meantree$edge.length
EB3.en <- results3.en$meantree$edge.length
EB4.en <- results4.en$meantree$edge.length
EBS.en <- data.frame(V1 = EB1.en,V2 = EB2.en,V3 = EB3.en,V4 = EB4.en)
EB_m.en <- apply(EBS.en, 1, mean)
tree.scaled.en <- results1.en$meantree
tree.scaled.en$edge.length <- EB_m.en
plot(tree.scaled.en)


# MCMC chains validation ####

## Creating mcmc lists ####

### Both ####

mcmc1 <- btmcmc("./Phylogeny/CHAIN1_BT/PCS.txt.Log.txt")[,-c(3)]
mcmc2 <- btmcmc("./Phylogeny/CHAIN2_BT/PCS.txt.Log.txt")[,-c(3)]
mcmc3 <- btmcmc("./Phylogeny/CHAIN3_BT/PCS.txt.Log.txt")[,-c(3)]
mcmc4 <- btmcmc("./Phylogeny/CHAIN4_BT/PCS.txt.Log.txt")[,-c(3)]

mcmc1 <- as.mcmc(mcmc1, start = 55000000, thin = 50000)
mcmc2 <- as.mcmc(mcmc2, start = 55000000, thin = 50000)
mcmc3 <- as.mcmc(mcmc3, start = 55000000, thin = 50000)
mcmc4 <- as.mcmc(mcmc4, start = 55000000, thin = 50000)

mcmcs.list <- as.mcmc.list(list(mcmc1, mcmc2, mcmc3, mcmc4))

### Brain ####

mcmc1.br <- btmcmc("./Phylogeny/Brain/CHAIN1_BT/PCS_br.txt.Log.txt")[,-c(3)]
mcmc2.br <- btmcmc("./Phylogeny/Brain/CHAIN2_BT/PCS_br.txt.Log.txt")[,-c(3)]
mcmc3.br <- btmcmc("./Phylogeny/Brain/CHAIN3_BT/PCS_br.txt.Log.txt")[,-c(3)]
mcmc4.br <- btmcmc("./Phylogeny/Brain/CHAIN4_BT/PCS_br.txt.Log.txt")[,-c(3)]

mcmc1.br <- as.mcmc(mcmc1.br, start = 55000000, thin = 50000)
mcmc2.br <- as.mcmc(mcmc2.br, start = 55000000, thin = 50000)
mcmc3.br <- as.mcmc(mcmc3.br, start = 55000000, thin = 50000)
mcmc4.br <- as.mcmc(mcmc4.br, start = 55000000, thin = 50000)

mcmcs.list.br <- as.mcmc.list(list(mcmc1.br, mcmc2.br, mcmc3.br, mcmc4.br))

### Endocast ####

mcmc1.en <- btmcmc("./Phylogeny/Endocast/CHAIN1_BT/PCS_en.txt.Log.txt")[,-c(3)]
mcmc2.en <- btmcmc("./Phylogeny/Endocast/CHAIN2_BT/PCS_en.txt.Log.txt")[,-c(3)]
mcmc3.en <- btmcmc("./Phylogeny/Endocast/CHAIN3_BT/PCS_en.txt.Log.txt")[,-c(3)]
mcmc4.en <- btmcmc("./Phylogeny/Endocast/CHAIN4_BT/PCS_en.txt.Log.txt")[,-c(3)]

mcmc1.en <- as.mcmc(mcmc1.en, start = 55000000, thin = 50000)
mcmc2.en <- as.mcmc(mcmc2.en, start = 55000000, thin = 50000)
mcmc3.en <- as.mcmc(mcmc3.en, start = 55000000, thin = 50000)
mcmc4.en <- as.mcmc(mcmc4.en, start = 55000000, thin = 50000)

mcmcs.list.en <- as.mcmc.list(list(mcmc1.en, mcmc2.en, mcmc3.en, mcmc4.en))

## Effective sample size ####

### Both ####

(ES1 <- effectiveSize(mcmc1))
(ES2 <- effectiveSize(mcmc2))
(ES3 <- effectiveSize(mcmc3))
(ES4 <- effectiveSize(mcmc4))

ES_table <- data.frame(Chain1 = ES1, Chain2 = ES2, Chain3 = ES3, Chain4 = ES4)[-1,]
ES_table$parameter <- rownames(ES_table)
write_xlsx(ES_table, "./output/PHYLO_effective_size.xlsx")

### Brain ####

(ES1.br <- effectiveSize(mcmc1.br))
(ES2.br <- effectiveSize(mcmc2.br))
(ES3.br <- effectiveSize(mcmc3.br))
(ES4.br <- effectiveSize(mcmc4.br))

ES_table.br <- data.frame(Chain1 = ES1.br, Chain2 = ES2.br, Chain3 = ES3.br, Chain4 = ES4.br)[-1,]
ES_table.br$parameter <- rownames(ES_table.br)
#write_xlsx(ES_table.br, "./output/PHYLO_effective_size_br.xlsx")

### Endocast ####

(ES1.en <- effectiveSize(mcmc1.en))
(ES2.en <- effectiveSize(mcmc2.en))
(ES3.en <- effectiveSize(mcmc3.en))
(ES4.en <- effectiveSize(mcmc4.en))

ES_table.en <- data.frame(Chain1 = ES1.en, Chain2 = ES2.en, Chain3 = ES3.en, Chain4 = ES4.en)[-1,]
ES_table.en$parameter <- rownames(ES_table.en)
#write_xlsx(ES_table.en, "./output/PHYLO_effective_size_en.xlsx")

## Gelman diagnostic ####

### Both ####

GD <- gelman.diag(mcmcs.list)
GD
GD_table <- data.frame(point_est = GD$psrf[,1], upper_CI = GD$psrf[,2])[-1,]
GD_table$parameter <- rownames(GD_table)
write_xlsx(GD_table, "./output/PHYLO_GelmanDiag.xlsx")

### Brain ####

GD.br <- gelman.diag(mcmcs.list.br)
GD.br
GD_table.br <- data.frame(point_est = GD.br$psrf[,1], upper_CI = GD.br$psrf[,2])[-1,]
GD_table.br$parameter <- rownames(GD_table.br)
#write_xlsx(GD_table.br, "./output/PHYLO_GelmanDiag_br.xlsx")

### Endocast ####

GD.en <- gelman.diag(mcmcs.list.en)
GD.en
GD_table.en <- data.frame(point_est = GD.en$psrf[,1], upper_CI = GD.en$psrf[,2])[-1,]
GD_table.en$parameter <- rownames(GD_table.en)
#write_xlsx(GD_table.en, "./output/PHYLO_GelmanDiag_en.xlsx")

## Trace plots ####

### Both ####

mcmc_df_tr <- do.call(rbind, lapply(seq_along(mcmcs.list), function(chain_id) {
  chain_df <- as.data.frame(as.matrix(mcmcs.list[[chain_id]]))
  chain_df$Iteration <- 1:nrow(chain_df)
  chain_df$Chain <- paste0("Chain", chain_id)
  chain_df
}))

mcmc_long_tr <- pivot_longer(mcmc_df_tr, 
                             cols = -c(Iteration, Chain), 
                             names_to = "Parameter", 
                             values_to = "Value")


mcmc_long_tr$Parameter <- factor(mcmc_long_tr$Parameter, levels = unique(mcmc_long_tr$Parameter))

trace.plots <- ggplot(mcmc_long_tr, aes(x = Iteration, y = Value, color = Chain)) +
  geom_line(alpha = 0.7) +
  scale_color_viridis_d()+
  facet_wrap(~ Parameter, scales = "free_y") +
  theme_minimal() +
  labs(x = "Iteration",
       y = "Parameter Value",
       color = "Chain") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 45))
trace.plots

ggsave(trace.plots, filename = "./output/PHYLO_trace_plots.bmp", height = 10, width = 12, dpi = 500)

### Brain ####

mcmc_df_tr_br <- do.call(rbind, lapply(seq_along(mcmcs.list.br), function(chain_id) {
  chain_df <- as.data.frame(as.matrix(mcmcs.list.br[[chain_id]]))
  chain_df$Iteration <- 1:nrow(chain_df)
  chain_df$Chain <- paste0("Chain", chain_id)
  chain_df
}))

mcmc_long_tr_br <- pivot_longer(mcmc_df_tr_br, 
                             cols = -c(Iteration, Chain), 
                             names_to = "Parameter", 
                             values_to = "Value")


mcmc_long_tr_br$Parameter <- factor(mcmc_long_tr_br$Parameter, levels = mcmc_long_tr_br$Parameter[1:57])

ggplot(mcmc_long_tr_br, aes(x = Iteration, y = Value, color = Chain)) +
  geom_line(alpha = 0.7) +
  scale_color_viridis_d()+
  facet_wrap(~ Parameter, scales = "free_y") +
  theme_minimal() +
  labs(x = "Iteration",
       y = "Parameter Value",
       color = "Chain") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 45))
ggsave(filename = "./output/PHYLO_trace_plots_br.bmp", height = 10, width = 12, dpi = 500)

### Endocast ####

mcmc_df_tr_en <- do.call(rbind, lapply(seq_along(mcmcs.list.en), function(chain_id) {
  chain_df <- as.data.frame(as.matrix(mcmcs.list.en[[chain_id]]))
  chain_df$Iteration <- 1:nrow(chain_df)
  chain_df$Chain <- paste0("Chain", chain_id)
  chain_df
}))

mcmc_long_tr_en <- pivot_longer(mcmc_df_tr_en, 
                                cols = -c(Iteration, Chain), 
                                names_to = "Parameter", 
                                values_to = "Value")


mcmc_long_tr_en$Parameter <- factor(mcmc_long_tr_en$Parameter, levels = mcmc_long_tr_en$Parameter[1:41])

ggplot(mcmc_long_tr_en, aes(x = Iteration, y = Value, color = Chain)) +
  geom_line(alpha = 0.7) +
  scale_color_viridis_d()+
  facet_wrap(~ Parameter, scales = "free_y") +
  theme_minimal() +
  labs(x = "Iteration",
       y = "Parameter Value",
       color = "Chain") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 45))
#ggsave(filename = "./output/PHYLO_trace_plots_en.bmp", height = 10, width = 12, dpi = 500)


## Density of posterior values plots ####

### Both ####

mcmc_df_ds <- do.call(rbind, lapply(seq_along(mcmcs.list), function(chain_id) {
  chain_df <- as.data.frame(as.matrix(mcmcs.list[[chain_id]]))
  chain_df$Chain <- paste0("Chain", chain_id)
  chain_df
}))
mcmc_df_ds <- mcmc_df_ds[,-1]
mcmc_long_ds <- pivot_longer(mcmc_df_ds, 
                             cols = -Chain, 
                             names_to = "Parameter", 
                             values_to = "Value")

mcmc_long_ds$Parameter <- factor(mcmc_long_ds$Parameter, levels = unique(mcmc_long_ds$Parameter))

density.plots <- ggplot(mcmc_long_ds, aes(x = Value, fill = Chain, color = Chain)) +
  geom_density(alpha = 0.4) +
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  facet_wrap(~ Parameter, scales = "free") +
  theme_minimal() +
  labs(x = "Parameter Value",
       y = "Density",
       fill = "Chain",
       color = "Chain") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 45, vjust = 0.8))

ggsave(density.plots, filename = "./output/PHYLO_density_plots.bmp", height = 10, width = 12, dpi = 500)

### Brain ####

mcmc_df_ds_br <- do.call(rbind, lapply(seq_along(mcmcs.list.br), function(chain_id) {
  chain_df <- as.data.frame(as.matrix(mcmcs.list.br[[chain_id]]))
  chain_df$Chain <- paste0("Chain", chain_id)
  chain_df
}))
mcmc_df_ds_br <- mcmc_df_ds_br[,-1]
mcmc_long_ds_br <- pivot_longer(mcmc_df_ds_br, 
                             cols = -Chain, 
                             names_to = "Parameter", 
                             values_to = "Value")

mcmc_long_ds_br$Parameter <- factor(mcmc_long_ds_br$Parameter, levels = mcmc_long_ds_br$Parameter[1:57])

ggplot(mcmc_long_ds_br, aes(x = Value, fill = Chain, color = Chain)) +
  geom_density(alpha = 0.4) +
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  facet_wrap(~ Parameter, scales = "free") +
  theme_minimal() +
  labs(x = "Parameter Value",
       y = "Density",
       fill = "Chain",
       color = "Chain") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 45, vjust = 0.8))

#ggsave(filename = "./output/PHYLO_density_plots_br.bmp", height = 10, width = 12, dpi = 500)

### Endocast ####

mcmc_df_ds_en <- do.call(rbind, lapply(seq_along(mcmcs.list.en), function(chain_id) {
  chain_df <- as.data.frame(as.matrix(mcmcs.list.en[[chain_id]]))
  chain_df$Chain <- paste0("Chain", chain_id)
  chain_df
}))
mcmc_df_ds_en <- mcmc_df_ds_en[,-1]
mcmc_long_ds_en <- pivot_longer(mcmc_df_ds_en, 
                                cols = -Chain, 
                                names_to = "Parameter", 
                                values_to = "Value")

mcmc_long_ds_en$Parameter <- factor(mcmc_long_ds_en$Parameter, levels = mcmc_long_ds_en$Parameter[1:41])

ggplot(mcmc_long_ds_en, aes(x = Value, fill = Chain, color = Chain)) +
  geom_density(alpha = 0.4) +
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  facet_wrap(~ Parameter, scales = "free") +
  theme_minimal() +
  labs(x = "Parameter Value",
       y = "Density",
       fill = "Chain",
       color = "Chain") +
  theme(legend.position = "right", axis.text.x = element_text(angle = 45, vjust = 0.8))
#ggsave(filename = "./output/PHYLO_density_plots_en.bmp", height = 10, width = 12, dpi = 500)

# Shifts plot ####

shifts.plot <- plotShifts(results1, scalar = "node", excludeones = T, tips = T, nodecex = 2, shp = 24,typ = "fan")

shifts.plot.br <- plotShifts(results1.br, scalar = "node", excludeones = T, tips = T, nodecex = 2, shp = 24,typ = "fan")

shifts.plot.en <- plotShifts(results1.en, scalar = "node", excludeones = T, tips = T, nodecex = 2, shp = 24, typ = "fan")


# Discrete character evolution (life cycle and microhabitat) ####

LC <- shape.data.both$infos.sp$Life_cycle
names(LC) <- sub(" ", "_", shape.data.both$infos.sp$Species)

LC[which(names(LC) == "Lyciasalamandra_luschani" | names(LC) == "Lissotriton_montandoni" | names(LC) ==  "Salamandra_salamandra")] <- "Complex"

MH <- shape.data.both$infos.sp$Habitat
names(MH) <- sub(" ", "_", shape.data.both$infos.sp$Species)


ER.model <- fitMk(tree.scaled.br, LC, model = "ER")
AIC(ER.model)
SYM.model <- fitMk(tree.scaled.br, LC, model = "SYM")
AIC(SYM.model)
ARD.model <- fitMk(tree.scaled.br, LC, model = "ARD")
AIC(ARD.model)

ER.model <- fitMk(tree.scaled.en, LC, model = "ER")
AIC(ER.model)
SYM.model <- fitMk(tree.scaled.en, LC, model = "SYM")
AIC(SYM.model)
ARD.model <- fitMk(tree.scaled.en, LC, model = "ARD")
AIC(ARD.model)

ER.model <- fitMk(tree.scaled.br, MH, model = "ER")
AIC(ER.model)
SYM.model <- fitMk(tree.scaled.br, MH, model = "SYM")
AIC(SYM.model)
ARD.model <- fitMk(tree.scaled.br, MH, model = "ARD")
AIC(ARD.model)

ER.model <- fitMk(tree.scaled.en, MH, model = "ER")
AIC(ER.model)
SYM.model <- fitMk(tree.scaled.en, MH, model = "SYM")
AIC(SYM.model)
ARD.model <- fitMk(tree.scaled.en, MH, model = "ARD")
AIC(ARD.model)

mappi<-make.simmap(MCC.lad,MH, model = "ER")
mm <- simmap(mappi)
plot(summary(mm))
onj <- describe.simmap(mm)
plot(onj,type = "fan",fsize = 0.5)
plotSimmap(mappi,type = "fan")

gg <- ancr(fitMk(MCC.lad, LC, model = "ER"))
plot(gg,type = "fan", args.nodelabels = list(piecol = hue_pal()(5)))

gg <- ancr(fitMk(MCC.lad, MH, model = "ER"))
plot(gg,type = "fan", args.nodelabels = list(piecol = c(viridis(100)[30],viridis(100)[85],viridis(100)[1], "#56976F","#D2B48C")))


# Evolutionnary rates and rate shifts plot ####

## Both ####

tree.scaled.lad <- ladderize(tree.scaled)
MCC.lad <- ladderize(MCC)
new.cex <- shifts.plot$nodecex/2
new.cex[which(new.cex < 0.1)] <- 0
names(new.cex) <- as.character(shifts.plot$nodes)
cex.df <- data.frame(cex = new.cex, node = shifts.plot$nodes)
cex.df <- cex.df[order(cex.df$node),]

EVRnRS.plot(color_palette = "viridis", tree = MCC.lad, node_cex = c(0,(new.cex*2)), edge_variable = log(tree.scaled.lad$edge.length), export = T, export.params = list(file = "./output/Phylogeny_EVO.pdf", width = 10, height = 10))

hist.rates <- density.rate.plot(values = log(tree.scaled.lad$edge.length), ylim1 = c(0,120), export = T, export.params = list(path = "./output",filename = "Phylogeny_rate_PP.pdf", width = 20, height = 10))

PP.node.plot(tree = MCC.lad, PP.node = cex.df$cex, export = T, export.params = list(file = "./output/Phylogeny_PP_node.pdf", width = 10, height = 10))

## Brain ####

tree.scaled.lad.br <- ladderize(tree.scaled.br)
MCC.lad.br <- ladderize(MCC)
new.cex.br <- shifts.plot.br$nodecex/2
new.cex.br[which(new.cex.br < 0.1)] <- 0
cex.df <- data.frame(cex = new.cex.br, node = shifts.plot.br$nodes)
cex.df <- cex.df[order(cex.df$node),]

EVRnRS.plot(color_palette = "viridis", tree = MCC.lad.br, node_cex = c(0,cex.df$cex*2), edge_variable = log(tree.scaled.lad.br$edge.length), export = T, export.params = list(file = "./output/Phylogeny_EVO_br.pdf", width = 10, height = 10))

hist.rates.br <- density.rate.plot(values = log(tree.scaled.lad.br$edge.length), ylim1 = c(0,120), export = T, export.params = list(path = "./output",filename = "Phylogeny_rate_PP_br.pdf", width = 20, height = 10))

PP.node.plot(tree = MCC.lad.br, PP.node = cex.df$cex, export = T, export.params = list(file = "./output/Phylogeny_PP_node_br.pdf", width = 10, height = 10))

## Endocast ####

tree.scaled.lad.en <- ladderize(tree.scaled.en)
MCC.lad.en <- ladderize(MCC)
new.cex.en <- shifts.plot.en$nodecex/2
new.cex.en[which(new.cex.en < 0.1)] <- 0
names(new.cex.en) <- as.character(shifts.plot.en$nodes)
cex.df <- data.frame(cex = new.cex.en, node = shifts.plot.en$nodes)
cex.df <- cex.df[order(cex.df$node),]

EVRnRS.plot(color_palette = "viridis", tree = MCC.lad.en, node_cex = c(0,cex.df$cex*2), edge_variable = log(tree.scaled.lad.en$edge.length), export = T, export.params = list(file = "./output/Phylogeny_EVO_en.pdf", width = 10, height = 10))

hist.rates.en <- density.rate.plot(values = log(tree.scaled.lad.en$edge.length), ylim1 = c(0,120), export = T, export.params = list(path = "./output",filename = "Phylogeny_rate_PP_en.pdf", width = 20, height = 10))

PP.node.plot(tree = MCC.lad.en, PP.node = cex.df$cex, export = T, export.params = list(file = "./output/Phylogeny_PP_node_en.pdf", width = 10, height = 10))

# Saving tree ####

save(tree.scaled, tree.scaled.br, tree.scaled.en, MCC, PHYPCA, shape.data.both, lm.pairs.both, lm.pairs.brain, lm.pairs.endo, nb.landmarks.both, nb.landmarks.brain, nb.landmarks.endo, br.lands, en.lands, file ="./output/RDA/tree.rda")

