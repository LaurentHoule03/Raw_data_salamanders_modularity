
# Chapter 2

# Modularity and morphological integration - Step 6

# Laurent Houle

# Importing functions ####

source(paste0("F:/MicroCT_Scan/xx_landmark_brain_salamanders_xx/Chapiters/doc-chap2/functions.R"))

# Librairies ####

library(EMMLi)
library(stringr)
library(readxl)
library(geomorph)
library(Morpho)
library(tidyverse)
library(ape)
library(phytools)
library(viridis)
library(ggplot2)
library(qgraph)
library(writexl)
library(evolqg)

# Importation ####

load("./output/RDA/tree.rda")
load("./output/RDA/PCA_shape_data.rda")


lm.d.brain <- read_excel(path = "F:/MicroCT_Scan/xx_landmark_brain_salamanders_xx/Chapiters/doc-chap2/Patches/description_land_brain.xlsx",1)
lm.d.endo <- read_excel("F:/MicroCT_Scan/xx_landmark_brain_salamanders_xx/Chapiters/doc-chap2/Patches/description_land_endo.xlsx",1)

# Module definition

## Brain ####

modules_3_br <- factor(lm.d.brain$module_3)

modules_extremities <- lm.d.brain$module_3
modules_extremities[which(modules_extremities == "forebrain")] <- "forebrain_and_hindbrain"
modules_extremities[which(modules_extremities == "hindbrain")] <- "forebrain_and_hindbrain"
modules_extremities <- factor(modules_extremities)

modules_4_br <- factor(lm.d.brain$module_4)

modules_DV_br <- factor(lm.d.brain$dorso_ventral)

modules_5_br <- factor(lm.d.brain$module_5)

modules_5_tel_br <- factor(lm.d.brain$module_5_tel)

modules_6_br <- factor(lm.d.brain$module_6)

## Endocast ####

modules_ext <- lm.d.endo$module_3
modules_ext[which(modules_ext == "forebrain")] <- "forebrain_and_hindbrain"
modules_ext[which(modules_ext == "hindbrain")] <- "forebrain_and_hindbrain"
modules_ext <- factor(modules_ext)

modules_CNC <- factor(lm.d.endo$module_CNC)

modules_UD <- factor(lm.d.endo$module_UPDOWN)

modules_otic <- factor(lm.d.endo$module_otic)

modules_bones <- factor(lm.d.endo$module_bones)

modules_bones_para2 <- factor(lm.d.endo$mod_bones_para2)

#mirror

modules_6_en <- factor(lm.d.endo$mirror_6)

modules_5_tel_en <- factor(lm.d.endo$mirror_5_tel)

modules_5_dien_en <- factor(lm.d.endo$mirror_5_dien)

modules_4_en <- factor(lm.d.endo$mirror_4)

modules_3_en <- factor(lm.d.endo$mirror_3)

# Data cleaning ####

shape.data.sp <- list(shape = shape.data.both$symm.coords.sp, shape.allom.rm = shape.data.both$shape.allom.rm.sp,csize = shape.data.both$sp.csize, infos = shape.data.both$infos.sp, module = list(brain = list(hyp_m3_br = modules_3_br, 
                                                                                                                                                                                                                hyp_5_br = modules_5_br,
                                                                                                                                                                                                                hyp_5_tel_br = modules_5_tel_br,
                                                                                                                                                                                                                hyp_6_br = modules_6_br,
                                                                                                                                                                                                                hyp_ext_br = modules_extremities, hyp_4_br = modules_4_br, hyp_DV_br = modules_DV_br), 
                                                                                                                                                                                                   endocast = list(hyp_UD = modules_UD,
                                                                                                                                                                                                                   hyp_otic = modules_otic,
                                                                                                                                                                                                                   bones = modules_bones,
                                                                                                                                                                                                                   hyp_ext_en = modules_ext, hyp_CNC = modules_CNC, hyp_bones_p2 = modules_bones_para2,
                                                                                                                                                                                                                   mirror_6 = modules_6_en,
                                                                                                                                                                                                                   mirror_5_tel = modules_5_tel_en,
                                                                                                                                                                                                                   mirror_5_dien = modules_5_dien_en, 
                                                                                                                                                                                                                   mirror_4 = modules_4_en,
                                                                                                                                                                                                                   mirror_3 = modules_3_en)))


dimnames(shape.data.sp$shape)[[3]] <- sub(" ", "_", shape.data.sp$infos$Species)
dimnames(shape.data.sp$shape.allom.rm)[[3]] <- sub(" ", "_", shape.data.sp$infos$Species)

# Modularity and integration tests ####

## Brain ####

### A priori hypotheses ####

#### CR method ####

#results.brain.or <- pairwise.mod.analysesV3(shape = shape.data.sp$shape.allom.rm[br.lands,,], phy = tree.scaled.br, mod.hyp = shape.data.sp$module$brain, CI = F, comp = T)

#save(results.brain.or, file = "./output/RDA/results_mod_brain_or.rda")
load("./output/RDA/results_mod_brain_or.rda")

tab.br <- data.frame(CI = results.brain.or[["Modularity"]][["comp.CR"]][["sample.se"]]*1.96, nam = rownames(results.brain.or[["Modularity"]][["comp.CR"]][["pairwise.pooled.se"]]), Z = results.brain.or[["Modularity"]][["comp.CR"]]$sample.z)

tab.br <- tab.br[order(tab.br$Z),]
tab.br$nam <- factor(tab.br$nam, levels = tab.br$nam)
ggplot(tab.br, aes(y = nam, x = Z)) + geom_point(size = 3) +geom_errorbarh(aes(xmin = Z - CI, xmax = Z + CI), height = 0.1) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.8, hjust = 0.8))

#### EMMLi method ####

EMMLi.results.br.or <- EMMLi.mod(shape = shape.data.sp$shape.allom.rm[br.lands,,], 
                              landmark.names = lm.d.brain$Name, 
                              N_sample = length(shape.data.sp$infos[[1]]),
                              hypotheses = shape.data.sp$module$brain,
                              show.comparison = T, EMMLiv2 = T, phylo = tree.scaled.br)


vec.names <- names(shape.data.sp$module$brain)
vec.num <- 1:length(names(shape.data.sp$module$brain))
df.names <- data.frame(num = vec.num, name = vec.names)

EMMLi.results.br.or$model.comparison$Model.name <- rownames(EMMLi.results.br.or$model.comparison)
n.step1 <- str_split_i(EMMLi.results.br.or$model.comparison$Model.name, pattern = "mod", i = 2)
n.step2 <- as.numeric(str_split_i(n.step1, pattern = ".s", i = 1))
vec.names2 <- df.names$name[n.step2]

EMMLi.results.br.or$model.comparison$Model.name2 <- vec.names2

# Hypothesis 7 and 6 retained using CR method 
# Hypothesis 7 retained using EMMLi

### A posteriori ####

mod.hyp.br <- shape.data.sp$module$brain

mod.hyp.br$hyp_8_br <- as.character(shape.data.sp$module$brain$hyp_6_br)
mod.hyp.br$hyp_8_br[which(mod.hyp.br$hyp_8_br == "thalamus")] <- "thalamus_optic"
mod.hyp.br$hyp_8_br[which(mod.hyp.br$hyp_8_br == "optic_teg")] <- "thalamus_optic"
mod.hyp.br$hyp_8_br <- factor(mod.hyp.br$hyp_8_br)

mod.hyp.br$hyp_9_br <- as.character(shape.data.sp$module$brain$hyp_5_tel_br)
mod.hyp.br$hyp_9_br[which(mod.hyp.br$hyp_9_br == "thalamus")] <- "thalamus_optic"
mod.hyp.br$hyp_9_br[which(mod.hyp.br$hyp_9_br == "optic_teg")] <- "thalamus_optic"
mod.hyp.br$hyp_9_br <- factor(mod.hyp.br$hyp_9_br)

#results.brain <- pairwise.mod.analysesV3(shape = shape.data.sp$shape.allom.rm[br.lands,,], phy = tree.scaled.br, mod.hyp = mod.hyp.br, CI = F, comp = T)

#save(results.brain, file = "./output/RDA/results_mod_brain.rda")
load("./output/RDA/results_mod_brain.rda")

tab.br <- data.frame(CI = results.brain[["Modularity"]][["comp.CR"]][["sample.se"]]*1.96, nam = rownames(results.brain[["Modularity"]][["comp.CR"]][["pairwise.pooled.se"]]), Z = results.brain[["Modularity"]][["comp.CR"]]$sample.z)

tab.br <- tab.br[order(tab.br$Z),]
tab.br$nam <- factor(tab.br$nam, levels = tab.br$nam)
ggplot(tab.br, aes(y = nam, x = Z)) + geom_point(size = 3) +geom_errorbarh(aes(xmin = Z - CI, xmax = Z + CI), height = 0.1) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.8, hjust = 0.8))

EMMLi.results.br <- EMMLi.mod(shape = shape.data.sp$shape.allom.rm[br.lands,,], 
                                 landmark.names = lm.d.brain$Name, 
                                 N_sample = length(shape.data.sp$infos[[1]]),
                                 hypotheses = mod.hyp.br,
                                 show.comparison = T, EMMLiv2 = T, phylo = tree.scaled.br)


vec.names <- names(mod.hyp.br)
vec.num <- 1:length(names(mod.hyp.br))
df.names <- data.frame(num = vec.num, name = vec.names)

EMMLi.results.br$model.comparison$Model.name <- rownames(EMMLi.results.br$model.comparison)
n.step1 <- str_split_i(EMMLi.results.br$model.comparison$Model.name, pattern = "mod", i = 2)
n.step2 <- as.numeric(str_split_i(n.step1, pattern = ".s", i = 1))
vec.names2 <- df.names$name[n.step2]

EMMLi.results.br$model.comparison$Model.name2 <- vec.names2

# Hypothesis 5, 6, 7 and 8 are selected using CR
# Hypothesis 7 selected using EMMLi

## Endocast ####

### A priori hypotheses ####

#### CR method ####

results.endocast.or <- pairwise.mod.analysesV3(shape = shape.data.sp$shape.allom.rm[en.lands,,], phy = tree.scaled.en, mod.hyp = shape.data.sp$module$endocast[1:4], CI = F, comp = T)

save(results.endocast.or, file = "./output/RDA/results_mod_endocast_or.rda")
load("./output/RDA/results_mod_endocast_or.rda")

tab.en <- data.frame(CI = results.endocast.or[["Modularity"]][["comp.CR"]][["sample.se"]]*1.96, nam = rownames(results.endocast.or[["Modularity"]][["comp.CR"]][["pairwise.pooled.se"]]), Z = results.endocast.or[["Modularity"]][["comp.CR"]]$sample.z)

tab.en <- tab.en[order(tab.en$Z),]
tab.en$nam <- factor(tab.en$nam, levels = tab.en$nam)
ggplot(tab.en, aes(x = nam, y = Z)) + geom_point(size = 3) +geom_errorbar(aes(ymin = Z - CI, ymax = Z + CI), width = 0.2) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.8, hjust = 0.8))

#### EMMLi method ####

EMMLi.results.en.or <- EMMLi.mod(shape = shape.data.sp$shape.allom.rm[en.lands,,], 
                              landmark.names = lm.d.endo$desc, 
                              N_sample = length(shape.data.sp$infos[[1]]),
                              hypotheses = shape.data.sp$module$endocast,
                              show.comparison = T, EMMLiv2 = T, phylo = tree.scaled.en)



vec.names <- names(shape.data.sp$module$endocast)
vec.num <- 1:length(names(shape.data.sp$module$endocast))
df.names <- data.frame(num = vec.num, name = vec.names)

EMMLi.results.en.or$model.comparison$Model.name <- rownames(EMMLi.results.en.or$model.comparison)
n.step1 <- str_split_i(EMMLi.results.en.or$model.comparison$Model.name, pattern = "mod", i = 2)
n.step2 <- as.numeric(str_split_i(n.step1, pattern = ".s", i = 1))
vec.names2 <- df.names$name[n.step2]

EMMLi.results.en.or$model.comparison$Model.name2 <- vec.names2

# Hypothesis bones_p2 and bones retained from CR method
# Hypothesis bones_p2 retained from EMMLi

### A posteriori hypotheses ####

#### CR method ####

mod.hyp.endocast <- generate.comb.mod.hyp(mod.hyp = shape.data.sp$module$endocast)

results.endocast <- pairwise.mod.analysesV3(shape = shape.data.sp$shape.allom.rm[en.lands,,], phy = tree.scaled.en, mod.hyp = mod.hyp.endocast[c(names(shape.data.sp$module$endocast), "hyp_bones_p2_with_occipito-otic_and_parietal", "hyp_bones_p2_with_occipito-otic_and_parietal_and_post_parasphenoid", "hyp_bones_p2_with_frontal_and_parietal")], CI = F, comp = T)

save(results.endocast, file = "./output/RDA/results_mod_endocast.rda")
load("./output/RDA/results_mod_endocast.rda")

tab.en <- data.frame(CI = results.endocast[["Modularity"]][["comp.CR"]][["sample.se"]]*1.96, nam = rownames(results.endocast[["Modularity"]][["comp.CR"]][["pairwise.pooled.se"]]), Z = results.endocast[["Modularity"]][["comp.CR"]]$sample.z)

tab.en <- tab.en[order(tab.en$Z),]
tab.en$nam <- factor(tab.en$nam, levels = tab.en$nam)
ggplot(tab.en, aes(x = nam, y = Z)) + geom_point(size = 3) +geom_errorbar(aes(ymin = Z - CI, ymax = Z + CI), width = 0.2) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.8, hjust = 0.8))

#### EMMLi method ####

EMMLi.results.en <- EMMLi.mod(shape = shape.data.sp$shape.allom.rm[en.lands,,], 
                              landmark.names = lm.d.endo$desc, 
                              N_sample = length(shape.data.sp$infos[[1]]),
                              hypotheses = mod.hyp.endocast[c(names(shape.data.sp$module$endocast), "hyp_bones_p2_with_occipito-otic_and_parietal", "hyp_bones_p2_with_occipito-otic_and_parietal_and_post_parasphenoid", "hyp_bones_p2_with_frontal_and_parietal")],
                              show.comparison = T,
                              EMMLiv2 = T,
                              phylo = tree.scaled.en)

vec.names <- c(names(shape.data.sp$module$endocast), "hyp_bones_p2_with_occipito-otic_and_parietal", "hyp_bones_p2_with_occipito-otic_and_parietal_and_post_parasphenoid","hyp_bones_p2_with_frontal_and_parietal")
vec.num <- 1:length(vec.names)
df.names <- data.frame(num = vec.num, name = vec.names)

EMMLi.results.en$model.comparison$Model.name <- rownames(EMMLi.results.en$model.comparison)
n.step1 <- str_split_i(EMMLi.results.en$model.comparison$Model.name, pattern = "mod", i = 2)
n.step2 <- as.numeric(str_split_i(n.step1, pattern = ".s", i = 1))
vec.names2 <- df.names$name[n.step2]

EMMLi.results.en$model.comparison$Model.name2 <- vec.names2
EMMLi.results.en$model.comparison$Model.name2[which(is.na(EMMLi.results.en$model.comparison$Model.name2) == T)] <- " "

# Hypothesis bones_p2 and bones_p2_OP retained from CR method
# Hypothesis bones_p2 retained from CR and EMMLi method



## Endocast mirror ####

### A priori hypotheses ####

#### CR method ####

results.endocast.or.mirror <- pairwise.mod.analysesV3(shape = shape.data.sp$shape.allom.rm[en.lands,,], phy = tree.scaled.en, mod.hyp = shape.data.sp$module$endocast, CI = F, comp = T)

save(results.endocast.or.mirror, file = "./output/RDA/results_mod_endocast_or_mirror.rda")
load("./output/RDA/results_mod_endocast_or_mirror.rda")

tab.en <- data.frame(CI = results.endocast.or.mirror[["Modularity"]][["comp.CR"]][["sample.se"]]*1.96, nam = rownames(results.endocast.or.mirror[["Modularity"]][["comp.CR"]][["pairwise.pooled.se"]]), Z = results.endocast.or.mirror[["Modularity"]][["comp.CR"]]$sample.z)

tab.en <- tab.en[order(tab.en$Z),]
tab.en$nam <- factor(tab.en$nam, levels = tab.en$nam)
ggplot(tab.en, aes(x = nam, y = Z)) + geom_point(size = 3) +geom_errorbar(aes(ymin = Z - CI, ymax = Z + CI), width = 0.2) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.8, hjust = 0.8))

#### EMMLi method ####

EMMLi.results.en.or.mirror <- EMMLi.mod(shape = shape.data.sp$shape.allom.rm[en.lands,,], 
                                 landmark.names = lm.d.endo$desc, 
                                 N_sample = length(shape.data.sp$infos[[1]]),
                                 hypotheses = shape.data.sp$module$endocast,
                                 show.comparison = T, EMMLiv2 = T, phylo = tree.scaled.en)



vec.names <- names(shape.data.sp$module$endocast)
vec.num <- 1:length(names(shape.data.sp$module$endocast))
df.names <- data.frame(num = vec.num, name = vec.names)

EMMLi.results.en.or.mirror$model.comparison$Model.name <- rownames(EMMLi.results.en.or.mirror$model.comparison)
n.step1 <- str_split_i(EMMLi.results.en.or.mirror$model.comparison$Model.name, pattern = "mod", i = 2)
n.step2 <- as.numeric(str_split_i(n.step1, pattern = ".s", i = 1))
vec.names2 <- df.names$name[n.step2]

EMMLi.results.en.or.mirror$model.comparison$Model.name2 <- vec.names2


## Both ####

### Best brain and endocast hypotheses ####

mirror_6 <- paste0(shape.data.sp$module$endocast$mirror_6,"_en")
mirror_5_tel <- paste0(shape.data.sp$module$endocast$mirror_5_tel,"_en")

shape.data.sp$module[["both"]] <- list(brain = list(hyp_6_br = modules_6_br), endocast = list(hyp_bones_p2 = mod.hyp.endocast$hyp_bones_p2, mirror_6 = mirror_6, mirror_5_tel = mirror_5_tel))

#### CR method ####

mod.hyp.both <- generate.mod.hyp.comp(mod.hyp.1 = shape.data.sp$module$both$brain, mod.hyp.2 = shape.data.sp$module$both$endocast)

results.both <- pairwise.mod.analysesV3(shape = shape.data.sp$shape.allom.rm, phy = tree.scaled, mod.hyp = mod.hyp.both[c("hyp_6_br mirror_6", "hyp_6_br mirror_5_tel", "hyp_6_br hyp_bones_p2")], CI = F, comp = T)

save(results.both, file = "./output/RDA/results_mod_both.rda")
load("./output/RDA/results_mod_both.rda")

tab.br <- data.frame(CI = results.both[["Modularity"]][["comp.CR"]][["sample.se"]]*1.96, nam = rownames(results.both[["Modularity"]][["comp.CR"]][["pairwise.pooled.se"]]), Z = results.both[["Modularity"]][["comp.CR"]]$sample.z)

tab.br <- tab.br[order(tab.br$Z),]
tab.br$nam <- factor(tab.br$nam, levels = tab.br$nam)
ggplot(tab.br, aes(y = nam, x = Z)) + geom_point(size = 3) +geom_errorbarh(aes(xmin = Z - CI, xmax = Z + CI), height = 0.1) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.8, hjust = 0.8))

#### EMMLi method ####

EMMLi.results.both <- EMMLi.mod(shape = shape.data.sp$shape.allom.rm, 
                              landmark.names = c(lm.d.brain$Name, lm.d.endo$desc), 
                              N_sample = length(shape.data.sp$infos[[1]]),
                              hypotheses = mod.hyp.both[c("hyp_6_br mirror_6", "hyp_6_br mirror_5_tel", "hyp_6_br hyp_bones_p2")],
                              show.comparison = T,
                              EMMLiv2 = T,
                              phylo = MCC)

vec.names <- c("hyp_6_br mirror_6", "hyp_6_br mirror_5_tel", "hyp_6_br hyp_bones_p2")
vec.num <- 1:length(vec.names)
df.names <- data.frame(num = vec.num, name = vec.names)

EMMLi.results.both$model.comparison$Model.name <- rownames(EMMLi.results.both$model.comparison)
n.step1 <- str_split_i(EMMLi.results.both$model.comparison$Model.name, pattern = "mod", i = 2)
n.step2 <- as.numeric(str_split_i(n.step1, pattern = ".s", i = 1))
vec.names2 <- df.names$name[n.step2]

EMMLi.results.both$model.comparison$Model.name2 <- vec.names2
EMMLi.results.both$model.comparison$Model.name2[which(is.na(EMMLi.results.both$model.comparison$Model.name2) == T)] <- " "

# Hypothesis 7 + bones_p2 retained from CR method
# Hypothesis 7 + bones_p2 retained from EMMLi

# Export results ####

## Merging the CR and EMMLi results ####

tab.res.merged.brain.or <- merge.results.CR.EMMLi(CR = results.brain.or, EMMLi = EMMLi.results.br.or, names.hyp = shape.data.sp$module$brain)

tab.res.merged.brain <- merge.results.CR.EMMLi(CR = results.brain, EMMLi = EMMLi.results.br, names.hyp = mod.hyp.br)

tab.res.merged.endocast.or.mirror <- merge.results.CR.EMMLi(CR = results.endocast.or.mirror, EMMLi = EMMLi.results.en.or.mirror, names.hyp = shape.data.sp$module$endocast)

tab.res.merged.endocast <- merge.results.CR.EMMLi(CR = results.endocast, EMMLi = EMMLi.results.en, names.hyp = mod.hyp.endocast[c(names(shape.data.sp$module$endocast), "hyp_bones_p2_with_occipito-otic_and_parietal", "hyp_bones_p2_with_occipito-otic_and_parietal_and_post_parasphenoid", "hyp_bones_p2_with_frontal_and_parietal", "hyp_6_br mirror_6")])

tab.res.merged.both <- merge.results.CR.EMMLi(CR = results.both, EMMLi = EMMLi.results.both, names.hyp = mod.hyp.both[c("hyp_6_br mirror_6", "hyp_6_br mirror_5_tel")])

## Exporting ####

write_xlsx(tab.res.merged.brain.or, "./output/Modularity/Modularity_brain_apriori_results.xlsx")
write_xlsx(tab.res.merged.brain, "./output/Modularity/Modularity_brain_apost_results.xlsx")
write_xlsx(tab.res.merged.endocast.or.mirror, "./output/Modularity/Modularity_endo_apriori_mirror_results.xlsx")
write_xlsx(tab.res.merged.endocast, "./output/Modularity/Modularity_endo_apost_results.xlsx")

write_xlsx(tab.res.merged.both, "./output/Modularity/Modularity_both_apost_results.xlsx")

# CR pairwise comparison 

list.p.pairwise <- list(brain.apriori = as.data.frame(results.brain.or$Modularity$comp.CR$pairwise.P),
                        brain.apost = as.data.frame(results.brain$Modularity$comp.CR$pairwise.P),
                        endocast.apriori = as.data.frame(results.endocast.or.mirror$Modularity$comp.CR$pairwise.P),
                        endocast.apost = as.data.frame(results.endocast$Modularity$comp.CR$pairwise.P),
                        both = as.data.frame(results.both$Modularity$comp.CR$pairwise.P))
write_xlsx(list.p.pairwise, "./output/Modularity/Modularity_p_pairwise_hyps.xlsx")


list.Z.pairwise <- list(brain.apriori = as.data.frame(results.brain.or$Modularity$comp.CR$pairwise.z),
                        brain.apost = as.data.frame(results.brain$Modularity$comp.CR$pairwise.z),
                        endocast.apriori = as.data.frame(results.endocast.or.mirror$Modularity$comp.CR$pairwise.z),
                        endocast.apost = as.data.frame(results.endocast$Modularity$comp.CR$pairwise.z),
                        both = as.data.frame(results.both$Modularity$comp.CR$pairwise.z))
write_xlsx(list.Z.pairwise, "./output/Modularity/Modularity_Z_pairwise_hyps.xlsx")

# PLS pairwise comp (integration)

list.p.pairwise.int <- list(brain.apriori = as.data.frame(results.brain.or$Integration$comp.pls$pairwise.P),
                        brain.apost = as.data.frame(results.brain$Integration$comp.pls$pairwise.P),
                        endocast.apriori = as.data.frame(results.endocast.or.mirror$Integration$comp.pls$pairwise.P),
                        endocast.apost = as.data.frame(results.endocast$Integration$comp.pls$pairwise.P),
                        both = as.data.frame(results.both$Integration$comp.pls$pairwise.P))
write_xlsx(list.p.pairwise.int, "./output/Modularity/Integration_p_pairwise_hyps.xlsx")


list.Z.pairwise.int <- list(brain.apriori = as.data.frame(results.brain.or$Integration$comp.pls$pairwise.z),
                        brain.apost = as.data.frame(results.brain$Integration$comp.pls$pairwise.z),
                        endocast.apriori = as.data.frame(results.endocast.or.mirror$Integration$comp.pls$pairwise.z),
                        endocast.apost = as.data.frame(results.endocast$Integration$comp.pls$pairwise.z),
                        both = as.data.frame(results.both$Integration$comp.pls$pairwise.z))
write_xlsx(list.Z.pairwise.int, "./output/Modularity/Integration_Z_pairwise_hyps.xlsx")

# Sample Z

list.sample.Z.int <- list(brain.apriori = cbind(as.data.frame(results.brain.or$Integration$comp.pls$sample.z), names(results.brain.or$Integration$comp.pls$sample.z)),
                            brain.apost = cbind(as.data.frame(results.brain$Integration$comp.pls$sample.z), names(results.brain$Integration$comp.pls$sample.z)),
                          endocast.apriori = cbind(as.data.frame(results.endocast.or$Integration$comp.pls$sample.z), names(results.endocast.or$Integration$comp.pls$sample.z)),
                          endocast.apost = cbind(as.data.frame(results.endocast$Integration$comp.pls$sample.z), names(results.endocast$Integration$comp.pls$sample.z)),
                          both = cbind(as.data.frame(results.both$Integration$comp.pls$sample.z), names(results.both$Integration$comp.pls$sample.z)))

write_xlsx(list.sample.Z.int, "./output/Modularity/Integration_sampleZ_hyps.xlsx")

# Global integration

brain.GB <- globalIntegration(A = shape.data.both$shape.allom.rm.sp[br.lands,,])
endo.GB <- globalIntegration(A = shape.data.both$shape.allom.rm.sp[en.lands,,])
both.GB <- globalIntegration(A = shape.data.both$shape.allom.rm.sp)

# Qgraphs ####

colors.region7 <- c("#115740", "#0093DD","#AF8800","#E9052E","#D80BCF","#FD4513")
names(colors.region7) <- c("olfactory_bulb", "telencephalon","thalamus","hypothalamus","optic_teg","rhombencephalon")

colors.region7_en <- c("#115740", "#0093DD","#AF8800","#E9052E","#D80BCF","#FD4513")
names(colors.region7) <- c("olfactory_bulb_en", "telencephalon_en","thalamus_en","hypothalamus_en","optic_teg_en","rhombencephalon")

colors.region5_tel_en <- c("#00DB82","#AF8800","#E9052E","#D80BCF","#FD4513")
names(colors.region7) <- c("telencephalon_en","thalamus_en","hypothalamus_en","optic_teg_en","rhombencephalon_en")

col.list.br <- list(hyp_6_br = colors.region7)

colors.endo.bonesp2 <- c("#68BC28", "#4F2185","#FF9700","#A2A8AC","#00E7FF","#002143")
names(colors.endo.bonesp2) <- c("ant_parasphenoid", "frontal","occipito-otic","orbitoshpenoid","parietal","post_parasphenoid")

colors.endo.bonesp2_OP <- c("#68BC28", "#4F2185","#FFC704","#A2A8AC","#002143")
names(colors.endo.bonesp2_OP) <- c("ant_parasphenoid", "frontal","occipito-otic_and_parietal","orbitoshpenoid","post_parasphenoid")

colors.endo.bonesp2_FP <- c("#68BC28", "#737F49","#FF9700","#A2A8AC","#002143")
names(colors.endo.bonesp2_FP) <- c("ant_parasphenoid", "frontal_and_parietal","occipito-otic","orbitoshpenoid","post_parasphenoid")

col.list.en <- list(hyp_bones_p2 = colors.endo.bonesp2,
                    "hyp_bones_p2_with_occipito-otic_and_parietal" = colors.endo.bonesp2_OP,
                    hyp_bones_p2_with_frontal_and_parietal = colors.endo.bonesp2_FP,
                    mirror_6 = colors.region7_en,
                    mirror_5_tel = colors.region5_tel_en)

col.list.both <- list("hyp_6_br hyp_bones_p2_with_occipito_otic_and_parietal" = c(colors.region7, colors.endo.bonesp2_OP),
                      "hyp_6_br hyp_bones_p2" = c(colors.region7, colors.endo.bonesp2),
                      "hyp_6_br mirror_6" = c(colors.region7, colors.region7_en),
                      "hyp_6_br mirror_5_tel" = c(colors.region7, colors.region5_tel_en))

## Brain ####

### Modularity ####

#### CR method ####

# A priori
qgraph.mod.mult(results.brain.or$Modularity, shape.data.sp$module$brain[c("hyp_6_br")], hyp.selected = c("hyp_6_br"), prefix = "./output/Modularity/qgraphs/Modularity/CR/MOD_apriori_", treshold = 0.25, color = col.list.br)

#### EMMLi ####

# A priori
EMMLi.results.br.or.mult <- EMMLi.mod(shape = shape.data.sp$shape.allom.rm[br.lands,,], 
                              landmark.names = lm.d.brain$Name, 
                              N_sample = length(shape.data.sp$infos[[1]]),
                              hypotheses = shape.data.sp$module$brain[c("hyp_6_br")], mult = T,EMMLiv2 = T,phylo = tree.scaled.br)

qgraph.mod.mult(res.list = EMMLi.results.br.or.mult, all.levels = shape.data.sp$module$brain[c("hyp_6_br")], hyp.selected = c("hyp_6_br"), EMMLi = T, prefix = "./output/Modularity/qgraphs/Modularity/EMMLi/MOD_apriori_EM_", vsize = EMMLi.results.br.or.mult[c("hyp_6_br")], treshold= 0.25, color = col.list.br[c("hyp_6_br")])


#### Integration ####

# A priori
qgraph.int.mult(results.brain.or$Integration, shape.data.sp$module$brain, hyp.selected = c("hyp_6_br"), prefix = "./output/Modularity/qgraphs/Integration/INT_apriori_", treshold = 0.25, color = col.list.br[c("hyp_6_br")])

## Endocast ####

### Modularity ####

#### CR method ####

# A priori
qgraph.mod.mult(results.endocast.or.mirror$Modularity, shape.data.sp$module$endocast[c("hyp_bones_p2","mirror_6","mirror_5_tel")], hyp.selected = c("hyp_bones_p2" ,"mirror_6","mirror_5_tel"), prefix = "./output/Modularity/qgraphs/Modularity/CR/MOD_apriori_", treshold = 0.25, color = col.list.en[c("hyp_bones_p2","mirror_6","mirror_5_tel")])

#### EMMLi ####

# A priori
EMMLi.results.en.or.mult <- EMMLi.mod(shape = shape.data.sp$shape.allom.rm[en.lands,,], 
                                      landmark.names = lm.d.endo$desc, 
                                      N_sample = length(shape.data.sp$infos[[1]]),
                                      hypotheses = shape.data.sp$module$endocast[c("hyp_bones_p2","mirror_6","mirror_5_tel")], mult = T,EMMLiv2 = T,phylo = tree.scaled.en)

qgraph.mod.mult(res.list = EMMLi.results.en.or.mult, all.levels = shape.data.sp$module$endocast, hyp.selected = names(EMMLi.results.en.or.mult), EMMLi = T, prefix = "./output/Modularity/qgraphs/Modularity/EMMLi/MOD_apriori_EM_",vsize = EMMLi.results.en.or.mult, treshold = 0.25, color = col.list.en[c("hyp_bones_p2","mirror_6","mirror_5_tel")])

#### Integration ####

# A priori
qgraph.int.mult(results.endocast.or.mirror$Integration, shape.data.sp$module$endocast[c("mirror_6","mirror_5_tel")], hyp.selected = c("mirror_6","mirror_5_tel"), prefix = "./output/Modularity/qgraphs/Integration/INT_apriori_", treshold = 0.25, color = col.list.en[c("mirror_6","mirror_5_tel")])


## Both ####

### Modularity ####

#### CR method ####

qgraph.mod.mult(results.both$Modularity, mod.hyp.both[c("hyp_6_br mirror_6", "hyp_6_br mirror_5_tel")], hyp.selected = c("hyp_6_br mirror_6", "hyp_6_br mirror_5_tel"), prefix = "./output/Modularity/qgraphs/Modularity/CR/MOD_both_", treshold = 0.5, color = col.list.both)


#### EMMLi ####

EMMLi.results.both.mult <- EMMLi.mod(shape = shape.data.sp$shape.allom.rm, 
                                      landmark.names = c(lm.d.brain$Name,lm.d.endo$desc), 
                                      N_sample = length(shape.data.sp$infos[[1]]),
                                      hypotheses = mod.hyp.both[c("hyp_6_br mirror_6", "hyp_6_br mirror_5_tel")], mult = T,EMMLiv2 = T,phylo = MCC)

qgraph.mod.mult(res.list = EMMLi.results.both.mult[c("hyp_6_br mirror_6", "hyp_6_br mirror_5_tel")], all.levels = mod.hyp.both[c("hyp_6_br mirror_6", "hyp_6_br mirror_5_tel")], hyp.selected = c("hyp_6_br mirror_6", "hyp_6_br mirror_5_tel"), EMMLi = T, prefix = "./output/Modularity/qgraphs/Modularity/EMMLi/MOD_both_EM_", vsize = EMMLi.results.both.mult[c("hyp_6_br mirror_6", "hyp_6_br mirror_5_tel")], treshold = 0.5, color = col.list.both[c("hyp_6_br mirror_6", "hyp_6_br mirror_5_tel")])

#### Integration ####

qgraph.int.mult(results.both$Integration, mod.hyp.both[c("hyp_6_br hyp_bones_p2", "hyp_6_br mirror_6", "hyp_6_br mirror_5_tel")], hyp.selected = c("hyp_6_br hyp_bones_p2", "hyp_6_br mirror_6", "hyp_6_br mirror_5_tel"), prefix = "./output/Modularity/qgraphs/Integration/INT_both_", treshold = 0.5, color = col.list.both[c("hyp_6_br hyp_bones_p2", "hyp_6_br mirror_6", "hyp_6_br mirror_5_tel")])

# Saving output ####

save(shape.data.sp, results.brain, results.brain.or, results.endocast, results.endocast.or, results.both, tree.scaled, MCC, mod.hyp.both, mod.hyp.endocast, file ="./output/RDA/modularity.rda")


lm.d.brain$module_6
mod.hyp.both$`hyp_6_br mirror_6`[1:181]
