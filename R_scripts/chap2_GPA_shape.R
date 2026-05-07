# Chapter 2

# GPA shape - Step 2

# Laurent Houle

# Importing functions ####

source(paste0(getwd(),"/functions.R"))

# Librairies ####

library(stringr)
library(readxl)
library(geomorph)
library(Morpho)
library(tidyverse)

# Importation ####

load("./output/RDA/raw_shape_data.rda")

# Get landmark no of surface semilandmarks ####

## Brain ####

type.land.ND.brain <- type.land.brain[-which(type.land.brain$type == "DELETE"),]
type.land.ND.brain$no2 <- 1:length(type.land.ND.brain$type)
surface.land.brain <- filter(type.land.ND.brain, type == "surfaceSM")
raw.data.brain$surface <- surface.land.brain$no2

## Endocast ####

type.land.ND.endo <- type.land.endo[-which(type.land.endo$type == "DELETE"),]
type.land.ND.endo <- type.land.ND.endo[-c(4,5),]
type.land.ND.endo$no2 <- 1:length(type.land.ND.endo$type)
surface.land.endo <- filter(type.land.ND.endo, type == "surfaceSM")
raw.data.endo$surface <- surface.land.endo$no2

## Both ####

raw.data.both$surface <- c(surface.land.brain$no2, surface.land.endo$no2 + nb.landmarks.brain)

rm(surface.land.brain, surface.land.endo, type.land.ND.brain, type.land.ND.endo)

# GPA ####

super.both <- gpagen(raw.data.both$raw, verbose = T, curves = curves.s.brain, surfaces = raw.data.both$surface, print.progress = F, ProcD = TRUE)

dimnames(super.both$coords)[[3]] <- dimnames(raw.data.both$raw)[[3]]

# Structure shape data ####

individuals.both <- find.ind(super.both$coords)

shape.data.both <- list(shape = super.both$coords, size = super.both$Csize, consensus = super.both$consensus, species.acr = raw.data.both$species, reps = raw.data.both$reps, ID = raw.data.both$ID, surface = raw.data.both$surface, lm.pairs = lm.pairs.both, ind = individuals.both$ind, infos = raw.data.both$infos)

## Averaging Csize by individual ####

shape.data.both$ind.csize <- mean_centroid_specimens(super = super.both, ind = shape.data.both$ind, reps = shape.data.both$reps)

# Removing asymmetry ####

new.coords.symm.both <- bilat.symmetry(A = shape.data.both$shape, ind = shape.data.both$ind, land.pairs = shape.data.both$lm.pairs, replicate = shape.data.both$reps, object.sym = TRUE, print.progress = FALSE)

summary(new.coords.symm.both)
shape.data.both$symm.coords.ind <- new.coords.symm.both$symm.shape

dimnames(shape.data.both$symm.coords.ind)[[3]] <- sub(" ", "_", shape.data.both$infos$Species)
rm(new.coords.symm.both)

# Averaging shape by species ####

shape.data.both$symm.coords.sp <- (aggregate(two.d.array(shape.data.both$symm.coords.ind) ~ shape.data.both$infos$acr, FUN=mean))[,-1]
rownames(shape.data.both$symm.coords.sp) <- unique(shape.data.both$infos$acr)
shape.data.both$symm.coords.sp <- arrayspecs(shape.data.both$symm.coords.sp, p = dim(shape.data.both$symm.coords.ind)[1], k = 3)

# Averaging centroid size by species ####

shape.data.both$sp.csize <- as.vector(aggregate(as.matrix(shape.data.both$ind.csize) ~ shape.data.both$infos$acr, FUN = mean)[,-1])

# Averaging infos by species ####

shape.data.both$infos.sp <- averaging.infos(infos = shape.data.both$infos[,-c(9:12,16,17)], avg.column = "Species")


# Arranging the order of the species names ####

shape.data.both$infos.sp <- shape.data.both$infos.sp[order(shape.data.both$infos.sp$acr),]

# Landmarks identification in entire shape dataset ####

br.lands <- 1:nb.landmarks.brain
en.lands <- (nb.landmarks.brain + 1):nb.landmarks.both

# Saving output ####

save(shape.data.both, nb.landmarks.both, nb.landmarks.brain, nb.landmarks.endo, lm.pairs.both, lm.pairs.brain, lm.pairs.endo, br.lands, en.lands, file ="./output/RDA/GPA_shape_data.rda")

