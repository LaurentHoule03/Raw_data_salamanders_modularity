# Chapter 2

# Import shape - Step 1

source(paste0(getwd(),"/functions.R"))

# Librairies ####

library(stringr)
library(readxl)
library(geomorph)
library(Morpho)
library(tidyverse)


# Importation ####

wd_base <- getwd()

sample.data <- read_excel(paste(wd_base, "/DATA/Sample_salamanders.xlsx", sep = ""),1)

filelist <- list.files(path = paste(wd_base, "/landmarks_coords", sep = ""), pattern = ".nts", recursive = T, full.names = T)
filenames <- list.files(path = paste(wd_base, "/landmarks_coords", sep = ""), pattern = ".nts", recursive = T, full.names = F)

Data.temp <- Import.nts.spec(filelist = filelist, filenames = filenames, sep = "_", words = c("Brain", "Endocast"), remove = ".nts")

Data.2d.brain <- two.d.array(Data.temp[[1]])
Data.2d.endo <- two.d.array(Data.temp[[2]])
raw.land.brain <- arrayspecs(Data.2d.brain, dim(Data.temp[[1]])[1], dim(Data.temp[[1]])[2])
raw.land.endo <- arrayspecs(Data.2d.endo, dim(Data.temp[[2]])[1], dim(Data.temp[[2]])[2])
rm(Data.temp, Data.2d.brain, Data.2d.endo)

# Landmark types ####

## Brain ####

type.land.brain <- read_excel(path = paste(wd_base, "/patches/type_raw_land_patch_brain.xlsx", sep = ""),1)
delete.land.brain <- filter(type.land.brain, type == "DELETE")
surf.land.brain <- filter(type.land.brain, type == "surfaceSM")
fixed.land.brain <- filter(type.land.brain, type == "fixed" | type == "curveSM")
keep.fixed.brain <- filter(type.land.brain, type == "fixed")

## Endocast ####

type.land.endo <- read_excel("./patches/type_raw_land_patch_endocast.xlsx",1)
delete.land.endo <- filter(type.land.endo, type == "DELETE")
surf.land.endo <- filter(type.land.endo, type == "surfaceSM")
fixed.land.endo <- filter(type.land.endo, type == "fixed")

# Adding automatic surface semilandmarks ####

## Brain ####

filelist <- list.files(path = paste(wd_base, "/Patches", sep = ""), pattern = ".nts", full.names = T)
Data.temp <- readmulti.nts(filelist[2])
Data.2d <- two.d.array(Data.temp)
raw.land.atlas <- arrayspecs(Data.2d, dim(Data.temp)[1], dim(Data.temp)[2])

rm(Data.temp, Data.2d)

rwl <- raw.land.atlas[-delete.land.brain$no,,]


lgt.land.wdel <- length(type.land.brain[[1]]) - length(delete.land.brain[[1]])

rwl <- array(rwl, dim = c(lgt.land.wdel, 3, 1), dimnames = list(NULL, NULL, "Brain_HYGE_rep1.ply"))

# Creating the atlas

Mean_shape_atlas <-ply2mesh(filename = paste(wd_base, "/PLYS/Brain_HYGE_rep1.ply", sep = ""), silent  = T)

Atlas <- createAtlas(Mean_shape_atlas, landmarks =  raw.land.atlas[fixed.land.brain$no,,]/1000, patch = raw.land.atlas[surf.land.brain$no,,]/1000, corrCurves = list(c(14:15), c(17:18), c(21:22), c(24:25), c(28:29), c(32:33), c(36:40), c(42:46), c(49:50), c(53:52)), keep.fix = as.integer(keep.fixed.brain$no))

plotAtlas(Atlas, meshcol = "white", alpha = 1, legend = T, cols = c("red", "blue", "cyan", "white"))

#snapshot3d(filename = "./output/atlas.png", width = 1000, height = 1000)
#snapshot3d(filename = "./output/atlas_side.png", width = 1000, height = 1000)
#snapshot3d(filename = "./output/atlas_ventral.png", width = 1000, height = 1000)


# Changing the name of each specimen with the filename of the mesh

filelistply <- list.files(path = paste(wd_base, "/PLYS", sep = ""), pattern = ".ply")

raw.land.brain <- change_coord_names_w_file_names(file.names = filelistply, coords = raw.land.brain, exclude.extension = T, extension = ".ply")

# Placing surface semilandmarks

raw.land.brain <- placePatch(Atlas, raw.land.brain/1000, path = paste(wd_base, "/PLYS", sep = ""), mc.cores = 4)

raw.data.brain <- geomorph.data.frame(raw = raw.land.brain,
                                     species = str_split_i(dimnames(raw.land.brain)[[3]], i = 2, pattern = "_"), 
                                     reps = str_split_i(dimnames(raw.land.brain)[[3]], i = 3, pattern = "_"), 
                                     ID = str_split_i(dimnames(raw.land.brain)[[3]], i = 4, pattern = "_"))

rm(raw.land.brain, rwl,raw.land.atlas, lgt.land.wdel)


## Endocast ####

# Importing atlas landmarks

filelist <- list.files(path = paste(wd_base, "/Patches", sep = ""), pattern = ".nts", full.names = T)
Data.temp <- readmulti.nts(filelist[4])

Data.2d <- two.d.array(Data.temp)
raw.land.atlas <- arrayspecs(Data.2d, dim(Data.temp)[1], dim(Data.temp)[2])

rm(Data.temp, Data.2d)

rwl <- raw.land.atlas[-delete.land.endo$no,,]

lgt.land.wdel <- length(type.land.endo[[1]]) - length(delete.land.endo[[1]])

rwl <- array(rwl, dim = c(lgt.land.wdel, 3, 1), dimnames = list(NULL, NULL, "Endocast_PSLE_rep1.ply"))

# Creating the atlas

Mean_shape_atlas <- ply2mesh(filename = paste(wd_base, "/PLYS/Endocast_PSLE_rep1.ply", sep = ""), silent  = T)

Atlas <- createAtlas(Mean_shape_atlas, landmarks =  raw.land.atlas[fixed.land.endo$no,,]/1000, patch = raw.land.atlas[surf.land.endo$no,,]/1000)

plotAtlas(Atlas, meshcol = "white", alpha = 1, legend = T, cols = c("red", "blue", "cyan", "white"))

#snapshot3d(filename = "./output/atlas_endo.png", width = 1000, height = 1000)
#snapshot3d(filename = "./output/atlas_side.png", width = 1000, height = 1000)
#snapshot3d(filename = "./output/atlas_ventral.png", width = 1000, height = 1000)

# Changing the name of each specimen with the filename of the mesh

filelistply <- list.files(path = paste(wd_base, "/PLYS", sep = ""), pattern = ".ply")

raw.land.endo <- change_coord_names_w_file_names(file.names = filelistply, coords = raw.land.endo, exclude.extension = T, extension = ".ply")

# Placing surface semilandmarks

raw.land.endo <- placePatch(Atlas, raw.land.endo/1000, path = paste(wd_base, "/PLYS", sep = ""), mc.cores = 4)

# Deleting landmark 4 and 5 ####

raw.land.endo <- raw.land.endo[-c(4,5),,]

raw.data.endo <- geomorph.data.frame(raw = raw.land.endo,
                                     species = str_split_i(dimnames(raw.land.endo)[[3]], i = 2, pattern = "_"), 
                                     reps = str_split_i(dimnames(raw.land.endo)[[3]], i = 3, pattern = "_"), 
                                     ID = str_split_i(dimnames(raw.land.endo)[[3]], i = 4, pattern = "_"))

rm(raw.land.endo, rwl,raw.land.atlas, lgt.land.wdel)




# Combining brain and endocast landmarks ####

both.raw <- combine.lands(coords.1 = raw.data.brain$raw, coords.2 = raw.data.endo$raw, nb.lands.total = 374)

raw.data.both <- geomorph.data.frame(raw = both.raw,
                                     species = str_split_i(dimnames(both.raw)[[3]], i = 2, pattern = "_"), 
                                     reps = str_split_i(dimnames(both.raw)[[3]], i = 3, pattern = "_"), 
                                     ID = str_split_i(dimnames(both.raw)[[3]], i = 4, pattern = "_"))

rm(both.raw)


# Paired landmarks ####

nb.landmarks.both <- 374
nb.landmarks.brain <- 181
nb.landmarks.endo <- 193


lm.pairs.brain <- matrix(nrow = 80, ncol = 2)
lm.pairs.brain[,1] <- c(5,8,11,13:19,27:30,35:40,48:50,55,56,65:69,70,71,72,73,74:77,84,85,91:95,96:99,104:108,114,115,118,119,122,123,126,128:131,136:140,155:159,167,168,171:173) # Left
lm.pairs.brain[,2] <- c(6,9,12,20:26,31:34,47:42,54:52,59,58,60:64,151,152,153,154,78:81,82,83,86:90,100:103,109:113,116,117,120,121,124,125,127,132:135,141:145,162:166,169,170,177:175) # Right

lm.pairs.endo <- matrix(nrow = 81, ncol = 2)
lm.pairs.endo[,1] <- c(2,4,7,9,12,16,17,22:30,59:50,60:69,80,82,88:91,92,93,98:104,112:118,126,127,136:140,157:161,167:173,181,182,185,186) # Left
lm.pairs.endo[,2] <- c(3,5,8,10,13,18,19,39:31,40:49,79:70,81,83,84:87,95,94,105:111,119:125,130,129,141:145,162:166,174:180,183,184,187,188) # Right

lm.pairs.both <- matrix(nrow = length(lm.pairs.brain[,1]) + length(lm.pairs.endo[,1]), ncol = 2)
lm.pairs.both[,1] <- c(lm.pairs.brain[,1], lm.pairs.endo[,1] + nb.landmarks.brain)
lm.pairs.both[,2] <- c(lm.pairs.brain[,2], lm.pairs.endo[,2] + nb.landmarks.brain)




# Other useful informartions ####

nb.ind.both <- length(raw.data.both$raw[1,1,])
nb.ind.brain <- length(raw.data.brain$raw[1,1,])
nb.ind.endo <- length(raw.data.endo$raw[1,1,])

nb.dims <- 3

curves.s.brain <- read_excel("./DATA/curves_semi_brain.xlsx",1)

raw.data.brain$infos <- find.infos.acr(sample.data[,c(1:6,8,10,13,14,15,16,18:20,22,23,56)], dimnames(raw.data.brain$raw)[[3]])

raw.data.endo$infos <- find.infos.acr(sample.data[,c(1:6,8,10,13,14,15,16,18:20,22,23,56)], dimnames(raw.data.endo$raw)[[3]])

raw.data.both$infos <- find.infos.acr(sample.data[,c(1:6,8,10,13,14,15,16,18:20,22,23,56)], dimnames(raw.data.both$raw)[[3]])

# Landmark validation ####

## Brain ####

names.split.brain <- str_split(dimnames(raw.data.brain$raw)[[3]], pattern = "_")
reps.order.brain <- find_reps_num(names.split.brain, string.rep1 = "rep1", string.rep2 = "rep2", position.reps = 3, position.name = 2)
rep1.brain <- raw.data.brain$raw[,,reps.order.brain$rep1]
rep2.brain <- raw.data.brain$raw[,,reps.order.brain$rep2]

super.rep1.brain <- gpagen(rep1.brain,verbose = T,print.progress = F)
super.rep2.brain <- gpagen(rep2.brain,verbose = T,print.progress = F)

plot.land.rep1.brain <- coords.hull.plot(coords = super.rep1.brain$coords, land.id.in.plot = F)
plot.land.rep1.brain

plot.land.rep2.brain <- coords.hull.plot(coords = super.rep2.brain$coords, land.id.in.plot = F)
plot.land.rep2.brain

## Endocast ####

names.split.endo <- str_split(dimnames(raw.data.endo$raw)[[3]], pattern = "_")
reps.order.endo <- find_reps_num(names.split.endo, string.rep1 = "rep1", string.rep2 = "rep2", position.reps = 3, position.name = 2)
rep1.endo <- raw.data.endo$raw[,,reps.order.endo$rep1]
rep2.endo <- raw.data.endo$raw[,,reps.order.endo$rep2]

super.rep1.endo <- gpagen(rep1.endo,verbose = T,print.progress = F)
super.rep2.endo <- gpagen(rep2.endo,verbose = T,print.progress = F)

plot.land.rep1.endo <- coords.hull.plot(coords = super.rep1.endo$coords, land.id.in.plot = F)
plot.land.rep1.endo

plot.land.rep2.endo <- coords.hull.plot(coords = super.rep2.endo$coords, land.id.in.plot = F)
plot.land.rep2.endo


# Saving output ####

save(raw.data.both, raw.data.brain, raw.data.endo, type.land.brain, type.land.endo, curves.s.brain, nb.landmarks.both, nb.landmarks.brain, nb.landmarks.endo, nb.dims, nb.ind.both, nb.ind.brain, nb.ind.endo, lm.pairs.both, lm.pairs.brain, lm.pairs.endo, file ="./output/RDA/raw_shape_data.rda")


land.id.brain <- read_excel("./patches/description_land_brain.xlsx")
land.id.endo <- read_excel("./patches/description_land_endo.xlsx")


polygone.br <- ply2mesh(filename = "./PLYS/Brain_AMME_rep1.ply", 
                          silent  = T)

shade3d(polygone.br, col = "white")
spheres3d(raw.data.both$raw[1:nb.landmarks.brain,,"Brain_AMME_rep1"], 
          radius = 0.00007, color = as.numeric(factor(land.id.brain$type, levels = c("fixed","surfaceSM","curveSM"))))

#snapshot3d(filename = "./output/landmarks_brain_side.png", height = 1000, width = 1200, webshot = F)

polygone.en <- ply2mesh(filename = "./PLYS/Endocast_PSLE_rep1.ply", 
                     silent  = T)

shade3d(polygone.en, col = "white")
spheres3d(raw.data.endo$raw[1:nb.landmarks.endo,,"Endocast_PSLE_rep1"], 
          radius = 0.00004, color = as.numeric(factor(land.id.endo$type)))

#snapshot3d(filename = "./output/landmarks_endo_side.png", height = 1000, width = 1200, webshot = F)

