
# Chapter 2

# Fonctions

# Laurent Houle

source(paste0("F:/MicroCT_Scan/xx_landmark_brain_salamanders_xx/Chapiters/doc-chap2/EMMLiv2.R"))

# Import.nts.spec() ####p

Import.nts.spec <- function(filelist, filenames, words, sep, remove = NULL){
  
  require(geomorph)
  require(stringr)
  
  sep.names <- str_split(filenames, pattern = sep)
  
  list.multnts <- list()
  
  for(k in 1:length(words)){
    vec <- c()
    
    wordk <- words[k]
    
    for(i in 1:length(sep.names)){
    
      for(j in 1:length(sep.names[[i]])){
      
        if(sep.names[[i]][j] == wordk) vec <- c(vec,i)
      
      }
    
    }
    
    new <- readmulti.nts(filelist = filelist[vec])
    
    if(is.null(remove) == FALSE)filenames <- str_split_i(filenames, i = 1, pattern = remove)
    
    dimnames(new)[[3]] <- filenames[vec]
    
    list.multnts[[k]] <- new
    
  }
  
  names(list.multnts) <- words
  
  return(list.multnts)
}

# combine.lands ####

combine.lands <- function(coords.1, coords.2, using = "dimnames", name.position = 2, rep.position = 3, ID.position = 4, pattern = "_", nb.lands.total = 75){
  
  require(stringr)
  
  if(using == "dimnames"){
    
    split.c1 <- str_split(dimnames(coords.1)[[3]], pattern = pattern)
    split.c2 <- str_split(dimnames(coords.2)[[3]], pattern = pattern)
    
    name.c1 <- str_split_i(dimnames(coords.1)[[3]], pattern = pattern, i = name.position)
    name.c2 <- str_split_i(dimnames(coords.2)[[3]], pattern = pattern, i = name.position)
    
    rep.c1 <- str_split_i(dimnames(coords.1)[[3]], pattern = pattern, i = rep.position)
    rep.c2 <- str_split_i(dimnames(coords.2)[[3]], pattern = pattern, i = rep.position)
    
    ID.c1 <- str_split_i(dimnames(coords.1)[[3]], pattern = pattern, i = ID.position)
    ID.c2 <- str_split_i(dimnames(coords.2)[[3]], pattern = pattern, i = ID.position)
  }
  
  mats <- as.data.frame(coords.1[,,1])[NULL,]

  ind.names <- c()
  
  for(i in 1:length(coords.1[1,1,])){
    
    n.i <- dimnames(coords.1)[[3]][i]

    for(j in 1:length(coords.2[1,1,])){
      
      if(is.na(ID.c1[i]) == T){
        
        if(rep.c1[i] == rep.c2[j]){
          
          if(name.c1[i] == name.c2[j]){
            # ID NA, rep + name equal
            
            n1 <- as.data.frame(coords.1[,,i])
            n2 <- as.data.frame(coords.2[,,j])
            
            new <- rbind(n1,n2)
            mats <- rbind(mats, new)
            ind.names <- c(ind.names, n.i)
            
          }
          
        }
        
      }else{
        if(is.na(ID.c2[j]) == F){
          
          if(ID.c1[i] == ID.c2[j]){
            
            if(rep.c1[i] == rep.c2[j]){
              
              if(name.c1[i] == name.c2[j]){
                #ID + rep + name are equal
                
                n1 <- as.data.frame(coords.1[,,i])
                n2 <- as.data.frame(coords.2[,,j])
                
                new <- rbind(n1,n2)
                mats <- rbind(mats, new)
                ind.names <- c(ind.names, n.i)
                
              }
              
            }
            
          }
          
        }
        
        
        
      }
      
    }
  }
  nmats <- arrayspecs(as.matrix(mats), p = nb.lands.total, k = length(coords.1[1,,1]))
  dimnames(nmats)[[3]] <- ind.names
  
  return(nmats)
}

# change_coord_names_w_file_names() ####


change_coord_names_w_file_names <- function(file.names, coords, pattern.coords = ".nts", position.overlap.coords = 1, exclude.extension = F, extension = NULL){
  
  require(stringr)
  filename.or <- file.names
  if(exclude.extension == T) file.names <- str_split_i(file.names, pattern = extension, i = 1)
  
  file.ov <- file.names
  coords.ov <- str_split_i(dimnames(coords)[[3]], pattern = pattern.coords, i = position.overlap.coords)

  new.coords.names <- c()
  
  for(i in 1:length(file.ov)){
    str.ckpt <- str_split_i(filename.or[i], i = 2, pattern = ".ply")
    for(j in 1:length(coords.ov)){
      
      
      if(file.ov[i] == coords.ov[j] & str.ckpt[[1]] == ""){
        
        new.coords.names <- c(new.coords.names, (file.names[i]))
        
        break
      }
      
    }
    
  }
  
  if(length(dimnames(coords)[[3]]) == length(new.coords.names)){
    
    dimnames(coords)[[3]] <- new.coords.names
  return(coords)
  }else{
    cat("The following filenames are incorrect:", "\n")
    print(setdiff(dimnames(coords)[[3]], new.coords.names))
    
    
  }
  
  
}

# find.ind() ####

find.ind <- function(coords, pattern = "_"){
  
  require(stringr)
  
  struct <- str_split_i(dimnames(coords)[[3]], pattern = pattern, i = 1)
  acr <- str_split_i(dimnames(coords)[[3]], pattern = pattern, i = 2)
  reps <- str_split_i(dimnames(coords)[[3]], pattern = pattern, i = 3)
  ID <- str_split_i(dimnames(coords)[[3]], pattern = pattern, i = 4)
  
  ind <- c()
  
  for(i in 1:length(dimnames(coords)[[3]])){
    
    new <- paste0(struct[[i]], acr[[i]], ID[[i]])
    
    ind <- c(ind, new)
    
  }
  ind <- as.numeric(factor(ind))
  return(data.frame(ind = ind, ID = dimnames(coords)[[3]]))
}

# morphospace.visualization() ####

morphospace.visualization <- function(pca.obj, tree = NULL, PCA = T, PCA.hull = F, phylomorpho = F, phylomorpho.hull = F, text = NULL, text_phylo = NULL,  nb.pc = 3, col.fill = NULL, name.fill = NULL,  theme_custom = theme_bw(), change.levels = NULL, alpha = 0.5,point.size = 3){
  
  require(combinat)
  require(ggplot2)
  require(ggrepel)
  require(R.devices)
  
  if(class(pca.obj)[1] == "princomp"){
    
    npcs <- length(pca.obj$sdev)
    tot <- sum(pca.obj$sdev)
    pour <- pca.obj$sdev/tot*100
    pour <- as.numeric(round(pour, 1))
    
  }else if(class(pca.obj)[1] == "gm.prcomp"){
    
    npcs <- length(pca.obj$d)
    pour <- pca.obj$d/sum(pca.obj$d)*100
    pour <- as.numeric(round(pour, 1))
  }
  
  if(npcs < nb.pc)stop("nb.pc -> number of pcs to retain higher than actual number of pcs in pca.obj")
  
  if(class(pca.obj)[1] == "princomp"){
    
    
    pca.data <- as.data.frame(pca.obj$scores[,1:nb.pc])
    
  }else if(class(pca.obj)[1] == "gm.prcomp"){
    pca.data <- as.data.frame(pca.obj$x[,1:nb.pc])
    
  }
  
  pcs.nam <- c()
  pca.labs <- c()
  
  for(i in 1:nb.pc){
    
    pcs.nam <- c(pcs.nam, paste("PC", i, sep = ""))
    
    pca.labs <- c(pca.labs, paste(paste("PC", i, sep = ""),": ", pour[i],"%", sep = ""))
    
  } 
  
  colnames(pca.data) <- pcs.nam
  
  comb <- as.data.frame(combn(1:nb.pc, m = 2))
  
  PCA.list <- list()
  
  for(i in 1:length(comb)){
    
    x.val <- comb[1,i]
    y.val <- comb[2,i]
    pca.data$X <- pca.data[[comb[1,i]]]
    pca.data$Y <- pca.data[[comb[2,i]]]
    
    PCA.name.i <- paste("PC", x.val, "_VS_PC", y.val, sep ="")
    
    if(phylomorpho == T || phylomorpho.hull == T){
      
      if(is.null(tree) == T) stop("tree -> No phylo object provided")
      if(is.null(text_phylo) == T) stop("text -> No names provided with phylogenetic tree")
      
      R.devices::suppressGraphics({
        lines.b <- phylomorphospace(tree, matrix(c(pca.data$X, pca.data$Y), ncol = 2, dimnames = list(text_phylo, c("X","Y"))), nsteps = 200, ftype="off")
      })
      
      lines.PM <- data.frame(x = NULL,xend = NULL,y = NULL,yend = NULL)
      
      for(h in 1:length(lines.b$edge[,1])){
        
        new <- data.frame(x = lines.b$xx[[lines.b$edge[h,1]]], xend = lines.b$xx[[lines.b$edge[h,2]]], y = lines.b$yy[[lines.b$edge[h,1]]], yend = lines.b$yy[[lines.b$edge[h,2]]])
        lines.PM <- rbind(lines.PM, new)
      }
      
      
    }
    
    if(is.null(col.fill) == T){
      
      PCA <- ggplot(data = pca.data, aes(x = X, y = Y)) + geom_point(size = 2) + labs(x = pca.labs[x.val], y = pca.labs[y.val]) + theme_custom
      
      
      if(is.null(text) == F) PCA <- PCA + geom_text_repel(label = text, fontface = "italic")
      if(phylomorpho == T) PCA <- PCA + geom_segment(data = lines.PM,aes(x = x, y = y, xend = xend, yend = yend), color = "black")
      
      if(PCA.hull == T || phylomorpho.hull == T) stop("col.fill is NULL -> must contain a character vector or factor to produce hulls")
      
    }else{

      pca.data$FILL <- factor(col.fill)
      
      if(is.null(change.levels) == F)pca.data$FILL <- factor(pca.data$FILL, levels = change.levels)

      PCA <- ggplot(data = pca.data, aes(x = X, y = Y, fill = FILL)) + geom_point(size = point.size, shape = 21) + labs(x = pca.labs[x.val], y = pca.labs[y.val], fill = name.fill) + theme_custom
      
      if(is.null(text) == F) PCA <- PCA + geom_text_repel(label = text, fontface = "italic")
      if(phylomorpho == T) PCA <- PCA + geom_segment(data = lines.PM,aes(x = x, y = y, xend = xend, yend = yend), color = "black")
      
      if(PCA.hull == T){
        
        hull1 <- pca.data %>%
          slice(chull(X, Y))
        hull1 <- hull1[NULL,]
        for(j in 1:length(levels(pca.data$FILL))){
          new <- filter(pca.data, FILL == levels(pca.data$FILL)[j])
          
          if(length(new$X) >= 2){
            hull.new <- new %>%
              slice(chull(X, Y))
            hull1 <- rbind(hull1,hull.new)
          }
          
        }
        
        PCA <- ggplot(data = pca.data, aes(x = X, y = Y, fill = FILL))  + labs(x = pca.labs[x.val], y = pca.labs[y.val], fill = name.fill) + geom_polygon(data = hull1, alpha = alpha, color = "black")+ geom_point(size = point.size, shape = 21, aes(x = X, y = Y, fill = FILL))+ theme_custom
        
        if(is.null(text) == F) PCA <- PCA + geom_text_repel(label = text, fontface = "italic")
        
        if(phylomorpho.hull == T){
          
          PCA <- ggplot()+ geom_segment(data = lines.PM,aes(x = x, y = y, xend = xend, yend = yend))+ geom_polygon(data = hull1, aes(x = X, y = Y, fill = FILL), alpha = alpha)+geom_point(data = pca.data, aes(x = X, y = Y, fill = FILL), size = point.size, inherit.aes = F, shape = 21) + labs(x = pca.labs[x.val], y = pca.labs[y.val], fill = name.fill,shape = name.fill) +  theme_custom 
          
          if(is.null(text) == F) PCA <- PCA + geom_text_repel(data = pca.data, aes(x = X, y = Y), label = text, fontface = "italic")
        }
        
      }
      
      
    }
    
    PCA.list[[i]] <- PCA
    names(PCA.list)[i] <- PCA.name.i
  }
  return(PCA.list)
}

convex.hulls <- function(data, name.x = colnames(data)[1], name.y = colnames(data)[2], name.fill = colnames(data)[3]){
  
  require(tidyverse)
  require(grDevices)
  
  data$X <- data[[name.x]]
  data$Y <- data[[name.y]]
  data$FILL <- factor(data[[name.fill]])
  
  hull1 <- data %>%
    slice(chull(X, Y))
  hull1 <- hull1[NULL,]
  for(j in 1:length(levels(data$FILL))){
    new <- filter(data, FILL == levels(data$FILL)[j])
    
    if(length(new$X) > 2){
      hull.new <- new %>%
        slice(chull(X, Y))
      hull1 <- rbind(hull1,hull.new)
    }
    
  }
  return(hull1)
}

# find.infos.acr() ####

find.infos.acr <- function(sample.data, filenames, pattern = "_",position.acr.fn = 2, acr.SD = "acr", id.SD = "ID2", position.acr.id = 4, position.acr.reps = 3){
  
  require(stringr)
  require(tidyverse)
  
  fn <- str_split_i(filenames, pattern = pattern, i = position.acr.fn)
  id <- str_split_i(filenames, pattern = pattern, i = position.acr.id)
  reps <- str_split_i(filenames, pattern = pattern, i = position.acr.reps)
  
  out <- sample.data[NULL,]
  
  for(i in 1:length(fn)){
    
    if(reps[i] == "rep2")next
    
    sd <- sample.data[sample.data[[acr.SD]] %in% fn[i],]
    
    if(length(sd[[1]]) > 1){
      sd <- sd[sd[[id.SD]] %in% id[i],]
      
    }
    
    out <- rbind(out, sd)
  }
  
  return(out)
  
  
}

mean_centroid_specimens <- function(super, ind, reps){
  
  m.csize <- c()
  
  N <- length(unique(ind))
  out <- c()
  
  for(i in 1:N){
    
    indi <- as.numeric(levels(factor(ind)))[i]
    rep1.cs <- c()
    rep2.cs <- c()
    for(j in 1:length(ind)){

      if(indi == ind[j] & reps[j] == "rep1"){

        
        rep1.cs <- super$Csize[j]
        
      }else if(indi == ind[j] & reps[j] == "rep2"){

        rep2.cs <- super$Csize[j]
        
      }
      
      
    }
    mcs <- mean(c(as.numeric(rep1.cs), as.numeric(rep2.cs)))

    out <- c(out, mcs)
    
  }
  
  return(out)
}

# averaging.infos() ####

averaging.infos <- function(infos, avg.column){
  
  require(tidyverse)
  
  out <- infos[NULL,]
  
  infos$fact <- factor(infos[[avg.column]])
  sls <- levels(infos$fact)
  
  for(i in 1:length(sls)){

    new <- filter(infos, fact == sls[i])
    
    if(length(new[[1]]) > 1){
      
      linei <- c()
      
      for(j in 1:length(new)){
        
        cl <- class(new[[j]])
        
        if(cl == "numeric"){
          
          linei <- c(linei, mean(new[[j]], na.rm = T))
          
        }else{
          
          linei <- c(linei, new[1,j])
          
        }
        
      }
      temp <- new[NULL,]
      
      new <- rbind(temp,linei)
      colnames(new) <- colnames(temp)
    }
    
    out <- rbind(out,new)
  }
  return(out)
}

# check.tips()####

check.tips <- function(tree, sample.names){
  
  all.names <- tree$tip.label
  
  tips <- c()
  for(i in 1:length(sample.names)){
    
    ind = "OUI"
    
    for(j in 1:length(all.names)){
      
      if(sample.names[i] == all.names[j]){
        ind = "NON"
        tips <- c(tips, sample.names[i])
      }
      
    }
    
    if(ind == "OUI"){
      
      cat("NOT in tree:", sample.names[i],"\n")
      
    }
    
  }
  
  return(tips)
}

# find.species.trees() ####

find.species.trees <- function(tree, pattern = "_", find, position = 2){
  
  require(stringr)
  
  nn <- str_split(tree$tip.label, pattern = pattern)
  
  for(i in 1:length(tree$tip.label)){
    
    if(nn[[i]][position] == find){
      
      print(tree$tip.label[i])
      
    }
    
  }
  
}

# filter.species() ####

filter.species <- function(data, suppress = T, sup.names = NULL){
  
  if(suppress == T){
    
    for(i in 1:length(data)){
      
      if(class(data[[i]]) == "array"){
        
        num.sup.n <- c()
        
        cat("Removed ")
        
        for(j in 1:length(sup.names)){
          num.sup.j <- c(1:length(dimnames(data[[i]])[[3]]))[which(dimnames(data[[i]])[[3]] == sup.names[j])]
        
          num.sup.n <- c(num.sup.n, num.sup.j)
          cat(paste0(sup.names[j], " "))
        }
        
        data[[i]] <- data[[i]][,,-num.sup.n]
        cat("from array.","\n")
        if(length(num.sup.n) != length(sup.names)) stop("ERROR with sup.names")
      }
      
      if(class(data[[i]]) == "character" | class(data[[i]]) == "factor" | class(data[[i]]) == "numeric"){
        l.b <- length(data[[i]])
        data[[i]] <- data[[i]][-num.sup.n]
        l.a <- length(data[[i]])
        
        dif <- l.b - l.a
        
        if(dif != length(sup.names))stop("ERROR with sup.names")
        
        cat("Removed", sup.names, "from vector.","\n")
        
      }
      
      if(class(data[[i]]) == "phylo"){
        require(ape)
        
        tips.ini <- data[[i]]$tip.label
        
        l.ti <- length(tips.ini)
        
        new.tips <- setdiff(tips.ini, sup.names)
        l.nt <- length(new.tips)
        
        dif <- l.ti - l.nt
        
        if(dif != length(sup.names))stop("ERROR with sup.names")
        
        ntree <- keep.tip(data[[i]], new.tips)
        
        data[[i]] <- ntree
        cat("Removed", sup.names, "from phylo.","\n")
      }
    }
    
  }
  
  return(data)
}

# filter.ind() ####

filter.ind <- function(data, suppress = T, sup.names = NULL){
  
  if(suppress == T){
    
    for(i in 1:length(data)){
      
      if(class(data[[i]]) == "array"){
        
        num.sup.n <- c()
        
        cat("Removed ")
        
        for(j in 1:length(sup.names)){
          num.sup.j <- c(1:length(dimnames(data[[i]])[[3]]))[which(dimnames(data[[i]])[[3]] == sup.names[j])]

          num.sup.n <- c(num.sup.n, num.sup.j)
          cat(paste0(sup.names[j], " "))
        }
        
        data[[i]] <- data[[i]][,,-num.sup.n]
        cat("from array.","\n")

        if(length(num.sup.n) < length(sup.names)) stop("ERROR with sup.names")
      }
      
      if(class(data[[i]]) == "character" | class(data[[i]]) == "factor" | class(data[[i]]) == "numeric"){
        l.b <- length(data[[i]])
        data[[i]] <- data[[i]][-num.sup.n]
        l.a <- length(data[[i]])
        
        dif <- l.b - l.a
        
        if(dif < length(sup.names))stop("ERROR with sup.names")
        
        cat("Removed", sup.names, "from vector.","\n")
        
      }
      
      if(class(data[[i]]) == "phylo"){
        require(ape)
        
        tips.ini <- data[[i]]$tip.label
        
        l.ti <- length(tips.ini)
        
        new.tips <- setdiff(tips.ini, sup.names)
        l.nt <- length(new.tips)
        
        dif <- l.ti - l.nt
        
        if(dif != length(sup.names))stop("ERROR with sup.names")
        
        ntree <- keep.tip(data[[i]], new.tips)
        
        data[[i]] <- ntree
        cat("Removed", sup.names, "from phylo.","\n")
      }
    }
    
  }
  
  return(data)
}

# coords.hull.plot() ####


coords.hull.plot <- function(coords, hull = T, fact = 0.01, points.in.hull = F, land.id.in.plot = F){
  
  require(grDevices)
  require(Irescale)
  require(R.devices)
  
  land.no <- as.character(1:length(coords[,1,1]))
  dat.temp <- data.frame(x= NULL,y= NULL,z= NULL, land = NULL)
  
  max.x <- max(coords[,1,])
  max.y <- max(coords[,2,])
  min.x <- min(coords[,1,])
  min.y <- min(coords[,2,])
  
  if(land.id.in.plot == T) centroids <- data.frame(x = NULL, y = NULL, land.no = NULL)
  
  for(j in 1:length(coords[1,1,])){
    
    new <- data.frame(x = as.numeric(coords[,1,j]),y= as.numeric(coords[,2,j]),z= as.numeric(coords[,3,j]), land = land.no)
    
    dat.temp <- rbind(dat.temp, new)
    
  }
  
  dat.temp$land <- factor(dat.temp$land, levels = land.no)
  
  rownames(dat.temp) <- NULL
  
  gr <- ggplot(data = dat.temp, aes(x=x,y=y, color = land))+
    ylim(min.y-fact, max.y+fact)+
    xlim(min.x-fact, max.x+fact)+
    geom_point()+
    theme_bw()+
    theme(legend.position = "none")
  
  
  if(hull == T){
    
    hull1 <- dat.temp %>%
      slice(chull(x, y))
    hull1 <- hull1[NULL,]
    for(i in 1:length(levels(dat.temp$land))){
      new <- filter(dat.temp, land == levels(dat.temp$land)[i])
      
      if(length(new$x) > 2){
        hull.new <- new %>%
          slice(chull(x, y))
        
        if(land.id.in.plot == T){
          
          suppressWarnings(suppressGraphics(area_centroid<-convexHull(hull.new[,1:2], hull.new$z)))
          
          new.centroid <- data.frame(x = area_centroid$centroid[1], y = area_centroid$centroid[2], land.no = as.character(i))
          
          centroids <- rbind(centroids, new.centroid)
        }
        
        
        hull1 <- rbind(hull1,hull.new)
      }
      
    }
    
    if(points.in.hull == T){
      
      gr <- ggplot(data = dat.temp, aes(x = x, y = y))+
        geom_point(size = 2, shape = 21, color = "black",aes(fill = dat.temp$land))+
        theme_bw()+
        geom_polygon(data = hull1, alpha = 0.2,color = "black",aes(fill = land))+
        labs(fill = "", x = "X", y = "Y")+
        theme(legend.position = "none")
      
      if(land.id.in.plot == T)gr <- gr + geom_text(data = centroids, aes(x = x, y = y), label = land.no, inherit.aes = F)
      
      
    }else{
      
      gr <- ggplot(data = dat.temp, aes(x = x, y = y))+
        theme_bw()+
        geom_polygon(data = hull1, alpha = 0.2,color = "black",aes(fill = land))+
        labs(fill = "", x = "X", y = "Y")+
        theme(legend.position = "none")
      
      if(land.id.in.plot == T)gr <- gr + geom_text(data = centroids, aes(x = x, y = y), label = land.no, inherit.aes = F)
      
    }
    
    
  }
  return(gr)
  
}

# find_reps_num() ####

find_reps_num <- function(names.split, string.rep1, string.rep2, position.reps, position.name){
  
  type.spec <- list()
  rep1.num <- c()
  rep2.num <- c()
  for(i in 1:length(names.split)){
    
    type.spec <- c(type.spec, str_split(names.split[[i]][position.reps], pattern = ".nts")[[1]][1])
    names(type.spec)[[i]] <- names.split[[i]][position.name]
    
    if(type.spec[[i]][1] == string.rep1){
      
      rep1.num <- c(rep1.num, i)
      
    }else if(type.spec[[i]][1] == string.rep2){
      
      rep2.num <- c(rep2.num, i)
      
    }else{
      
      type.spec[[i]] <- str_split(names.split[[i]][position.reps+1], pattern = ".nts")[[1]][1]
      names(type.spec)[[i]] <- names.split[[i]][position.name]
      
      if(type.spec[[i]][1] == string.rep1){
        
        rep1.num <- c(rep1.num, i)
      }else if(type.spec[[i]][1] == string.rep2){
        rep2.num <- c(rep2.num, i)
      }else{
        stop("position.rep does not match with string.rep1 and 2")
      }
      
    }
  }
  names.species <- names(type.spec)
  
  out <- data.frame(rep1 = rep1.num, rep2 = rep2.num, names = names.species)
  
  return(out)
}

# morphospace3d() ####

morphospace3d <- function(data, groups, color.pal = viridis(10), lims = list(c(-0.2,0.2), c(-0.2,0.2), c(-0.2,0.2)), labs = list("", "", ""), add.points = F, names = NULL, show.axis = T, names.group = T, alpha = 0.5, shape = T, color.shape = "grey", ini = T, alpha.shape = 0.1, background = "white"){
  
  require(viridis)
  require(tidyverse)
  require(gMOIP)
  
  if(length(data) < 3)stop("Number of PCs inferior to 3")
  
  # loop to know how much convex hull 3d to produce
  
  data$group <- factor(groups)
  data$colors <- NA
  k <- c()
  new.data <- data[NULL,]
  centroids <- data.frame(x = NULL, y = NULL, z = NULL, group.text = NULL)
  
  for(i in 1:length(levels(data$group))){
    new <- filter(data, group == levels(data$group)[i])
    mx <- mean(new[,1], na.rm = T)
    my <- mean(new[,2], na.rm = T)
    mz <- mean(new[,3], na.rm = T)
    dat <- data.frame(x = mx, y = my, z = mz, group.text = levels(data$group)[i])
    if(add.points == T){
      if(length(new[[1]]) > 1)centroids <- rbind(centroids, dat)
    }else{
      
      centroids <- rbind(centroids, dat)
      
    }
    
    
    
    
    new$colors <- color.pal[i]
    new.data <- rbind(new.data, new)
    if(length(new[[1]]) > 2){
      new.n <- levels(data$group)[i]
      k <- c(k,new.n)
    }
  }
  data <- new.data
  # the plot
  
  print(centroids)
  
  if(ini == T){
    
    ini3D(argsPlot3d = list(xlim = lims[[1]],
                            ylim = lims[[2]],
                            zlim = lims[[3]],
                            xlab = labs[[1]],
                            ylab = labs[[2]],
                            zlab = labs[[3]],
                            bg = background))
    bg3d(background)
    
  }
  
  
  
  if(add.points == T){
    
    rownames(data) <- names
    
    plotPoints3D(data[,1:3],argsPlot3d = list(col = data$colors,type ="s",size = 0.5), addText = names,argsText3d = list(adj = 1))
    
    if(names.group == T)text3d(centroids[,1:3], texts = centroids$group.text)
    
  }else{
    
    if(names.group == T)text3d(centroids[,1:3], texts = centroids$group.text)
  }
  
  
  for(i in 1:length(k)){
    
    new <- filter(data,group == k[i])
    colori <- new$colors[1]
    
    lst <- plotHull3D(as.matrix(new[,1:3]), argsPolygon3d = list(alpha=alpha, color = colori), argsSegments3d = list(color=colori))
  }
  
  if(shape == T)plotHull3D(as.matrix(data[,1:3]), argsPolygon3d = list(alpha=alpha.shape, color = color.shape), argsSegments3d = list(color=color.shape))
  
  if(show.axis == T)finalize3D(argsTitle3d = list(xlab = "",ylab = "",zlab = ""))
  
  
}

# ProcD.comp() ####

procD.comp <- function(spec.shape, groups, plot = T, category.names = c("Overall", "Within groups", "Among groups"), title = NULL){
  
  D.overall <- dist(two.d.array(spec.shape))
  
  D.W.ind <- NULL
  for(i in unique(groups)){
    tmp <- which(groups == i)
    if(length(tmp) > 1){
      
      D.W.ind <- c(D.W.ind, dist(two.d.array(spec.shape[,,tmp]))) }
  }
  
  group.shape <- (aggregate(two.d.array(spec.shape) ~ groups, FUN=mean))[,-1]
  rownames(group.shape) <- unique(groups)
  group.shape <- arrayspecs(group.shape,p=dim(spec.shape)[1],k=dim(spec.shape)[2])
  
  
  D.among.groups <- dist(two.d.array(group.shape))
  
  out <- list(Overall.distance = D.overall, 
              Within.spec.distance = D.W.ind,
              Among.groups.distance = D.among.groups)
  
  if(plot == T){
    
    dat <- data.frame(x = c(D.overall, D.W.ind, D.among.groups), FILL = c(rep(category.names[1],length(D.overall)), rep(category.names[2], length(D.W.ind)), rep(category.names[3], length(D.among.groups))))
    
    gr <- ggplot(dat, aes(x=x, fill = FILL))+geom_density(alpha = 0.5, outline.type = "full")+labs(x = "Procrustes distance", y = "Density", fill= "", title = title)+theme_bw()
    
    out[[4]] <- gr
    out[[5]] <- dat
    names(out)[4] <- "Plot"
    names(out)[5] <- "Density.data"
  }
  
  
  return(out)
  
}

# pairwise.mod.analyses() ####

pairwise.mod.analyses <- function(shape, phy, mod.hyp, CI = T, comp = T, integration = T, extract = T){
  library(knitr)
  library(geomorph)
  library(progress)
  
  new.mod <- mod.hyp
  nb.hyps <- length(new.mod)
  
  for(i in 1:length(mod.hyp)){
    
    hypi <- mod.hyp[[i]]
    n.hypi <- names(mod.hyp)[i]
    lgt.levs <- length(levels(factor(hypi)))
    
    for(j in 1:lgt.levs){
      
      levj <- levels(factor(hypi))[j]
      if(j < lgt.levs){
        
        for(k in j+1:lgt.levs){
          levk <- levels(factor(hypi))[k]
          levkp1 <- levels(factor(hypi))[k+1]
          levkp2 <- levels(factor(hypi))[k+2]
          levkp3 <- levels(factor(hypi))[k+3]
          levkp4 <- levels(factor(hypi))[k+4]
          levkp5 <- levels(factor(hypi))[k+5]
          levkp6 <- levels(factor(hypi))[k+6]
          levkp7 <- levels(factor(hypi))[k+7]
          levkp8 <- levels(factor(hypi))[k+8]
          levkp9 <- levels(factor(hypi))[k+9]
          levkp10 <- levels(factor(hypi))[k+10]
          levkp11 <- levels(factor(hypi))[k+11]
          
          if(is.na(levk) == F){
            nb.hyps <- nb.hyps + 1
            new <- paste(levj,levk, sep = " + ")
            hyp.new <- as.character(hypi)
            hyp.new[which(hyp.new == levj)] <- new
            hyp.new[which(hyp.new == levk)] <- new
            new.mod[[nb.hyps]] <- factor(hyp.new)
            names(new.mod)[nb.hyps] <- paste(n.hypi, new, sep = "_with_")
          }
          
          if(is.na(levkp1) == F){
            nb.hyps <- nb.hyps + 1
            new <- paste(levj,levk,levkp1, sep = " + ")
            hyp.new <- as.character(hypi)
            hyp.new[which(hyp.new == levj)] <- new
            hyp.new[which(hyp.new == levk)] <- new
            hyp.new[which(hyp.new == levkp1)] <- new
            new.mod[[nb.hyps]] <- factor(hyp.new)
            names(new.mod)[nb.hyps] <- paste(n.hypi, new, sep = "_with_")
          }
          
          if(is.na(levkp2) == F){
            nb.hyps <- nb.hyps + 1
            new <- paste(levj,levk,levkp1, levkp2, sep = " + ")
            hyp.new <- as.character(hypi)
            hyp.new[which(hyp.new == levj)] <- new
            hyp.new[which(hyp.new == levk)] <- new
            hyp.new[which(hyp.new == levkp1)] <- new
            hyp.new[which(hyp.new == levkp2)] <- new
            new.mod[[nb.hyps]] <- factor(hyp.new)
            names(new.mod)[nb.hyps] <- paste(n.hypi, new, sep = "_with_")
          }
          
          if(is.na(levkp3) == F){
            nb.hyps <- nb.hyps + 1
            new <- paste(levj,levk,levkp1, levkp2, levkp3, sep = " + ")
            hyp.new <- as.character(hypi)
            hyp.new[which(hyp.new == levj)] <- new
            hyp.new[which(hyp.new == levk)] <- new
            hyp.new[which(hyp.new == levkp1)] <- new
            hyp.new[which(hyp.new == levkp2)] <- new
            hyp.new[which(hyp.new == levkp3)] <- new
            new.mod[[nb.hyps]] <- factor(hyp.new)
            names(new.mod)[nb.hyps] <- paste(n.hypi, new, sep = "_with_")
          }
          
          if(is.na(levkp4) == F){
            nb.hyps <- nb.hyps + 1
            new <- paste(levj,levk,levkp1, levkp2, levkp3, levkp4, sep = " + ")
            hyp.new <- as.character(hypi)
            hyp.new[which(hyp.new == levj)] <- new
            hyp.new[which(hyp.new == levk)] <- new
            hyp.new[which(hyp.new == levkp1)] <- new
            hyp.new[which(hyp.new == levkp2)] <- new
            hyp.new[which(hyp.new == levkp3)] <- new
            hyp.new[which(hyp.new == levkp4)] <- new
            new.mod[[nb.hyps]] <- factor(hyp.new)
            names(new.mod)[nb.hyps] <- paste(n.hypi, new, sep = "_with_")
          }
          
          if(is.na(levkp5) == F){
            nb.hyps <- nb.hyps + 1
            new <- paste(levj,levk,levkp1, levkp2, levkp3, levkp4, levkp5, sep = " + ")
            hyp.new <- as.character(hypi)
            hyp.new[which(hyp.new == levj)] <- new
            hyp.new[which(hyp.new == levk)] <- new
            hyp.new[which(hyp.new == levkp1)] <- new
            hyp.new[which(hyp.new == levkp2)] <- new
            hyp.new[which(hyp.new == levkp3)] <- new
            hyp.new[which(hyp.new == levkp4)] <- new
            hyp.new[which(hyp.new == levkp5)] <- new
            new.mod[[nb.hyps]] <- factor(hyp.new)
            names(new.mod)[nb.hyps] <- paste(n.hypi, new, sep = "_with_")
          }
          
          if(is.na(levkp6) == F){
            nb.hyps <- nb.hyps + 1
            new <- paste(levj,levk,levkp1, levkp2, levkp3, levkp4, levkp5, levkp6, sep = " + ")
            hyp.new <- as.character(hypi)
            hyp.new[which(hyp.new == levj)] <- new
            hyp.new[which(hyp.new == levk)] <- new
            hyp.new[which(hyp.new == levkp1)] <- new
            hyp.new[which(hyp.new == levkp2)] <- new
            hyp.new[which(hyp.new == levkp3)] <- new
            hyp.new[which(hyp.new == levkp4)] <- new
            hyp.new[which(hyp.new == levkp5)] <- new
            hyp.new[which(hyp.new == levkp6)] <- new
            new.mod[[nb.hyps]] <- factor(hyp.new)
            names(new.mod)[nb.hyps] <- paste(n.hypi, new, sep = "_with_")
          }
          
          if(is.na(levkp7) == F){
            nb.hyps <- nb.hyps + 1
            new <- paste(levj,levk,levkp1, levkp2, levkp3, levkp4, levkp5, levkp6, levkp7, sep = " + ")
            hyp.new <- as.character(hypi)
            hyp.new[which(hyp.new == levj)] <- new
            hyp.new[which(hyp.new == levk)] <- new
            hyp.new[which(hyp.new == levkp1)] <- new
            hyp.new[which(hyp.new == levkp2)] <- new
            hyp.new[which(hyp.new == levkp3)] <- new
            hyp.new[which(hyp.new == levkp4)] <- new
            hyp.new[which(hyp.new == levkp5)] <- new
            hyp.new[which(hyp.new == levkp6)] <- new
            hyp.new[which(hyp.new == levkp7)] <- new
            new.mod[[nb.hyps]] <- factor(hyp.new)
            names(new.mod)[nb.hyps] <- paste(n.hypi, new, sep = "_with_")
          }
          
          if(is.na(levkp8) == F){
            nb.hyps <- nb.hyps + 1
            new <- paste(levj,levk,levkp1, levkp2, levkp3, levkp4, levkp5, levkp6, levkp7, levkp8, sep = " + ")
            hyp.new <- as.character(hypi)
            hyp.new[which(hyp.new == levj)] <- new
            hyp.new[which(hyp.new == levk)] <- new
            hyp.new[which(hyp.new == levkp1)] <- new
            hyp.new[which(hyp.new == levkp2)] <- new
            hyp.new[which(hyp.new == levkp3)] <- new
            hyp.new[which(hyp.new == levkp4)] <- new
            hyp.new[which(hyp.new == levkp5)] <- new
            hyp.new[which(hyp.new == levkp6)] <- new
            hyp.new[which(hyp.new == levkp7)] <- new
            hyp.new[which(hyp.new == levkp8)] <- new
            new.mod[[nb.hyps]] <- factor(hyp.new)
            names(new.mod)[nb.hyps] <- paste(n.hypi, new, sep = "_with_")
          }
          
          if(is.na(levkp9) == F){
            nb.hyps <- nb.hyps + 1
            new <- paste(levj,levk,levkp1, levkp2, levkp3, levkp4, levkp5, levkp6, levkp7, levkp8, levkp9, sep = " + ")
            hyp.new <- as.character(hypi)
            hyp.new[which(hyp.new == levj)] <- new
            hyp.new[which(hyp.new == levk)] <- new
            hyp.new[which(hyp.new == levkp1)] <- new
            hyp.new[which(hyp.new == levkp2)] <- new
            hyp.new[which(hyp.new == levkp3)] <- new
            hyp.new[which(hyp.new == levkp4)] <- new
            hyp.new[which(hyp.new == levkp5)] <- new
            hyp.new[which(hyp.new == levkp6)] <- new
            hyp.new[which(hyp.new == levkp7)] <- new
            hyp.new[which(hyp.new == levkp8)] <- new
            hyp.new[which(hyp.new == levkp9)] <- new
            new.mod[[nb.hyps]] <- factor(hyp.new)
            names(new.mod)[nb.hyps] <- paste(n.hypi, new, sep = "_with_")
          }
          
          if(is.na(levkp10) == F){
            nb.hyps <- nb.hyps + 1
            new <- paste(levj,levk,levkp1, levkp2, levkp3, levkp4, levkp5, levkp6, levkp7, levkp8, levkp9, levkp10, sep = " + ")
            hyp.new <- as.character(hypi)
            hyp.new[which(hyp.new == levj)] <- new
            hyp.new[which(hyp.new == levk)] <- new
            hyp.new[which(hyp.new == levkp1)] <- new
            hyp.new[which(hyp.new == levkp2)] <- new
            hyp.new[which(hyp.new == levkp3)] <- new
            hyp.new[which(hyp.new == levkp4)] <- new
            hyp.new[which(hyp.new == levkp5)] <- new
            hyp.new[which(hyp.new == levkp6)] <- new
            hyp.new[which(hyp.new == levkp7)] <- new
            hyp.new[which(hyp.new == levkp8)] <- new
            hyp.new[which(hyp.new == levkp9)] <- new
            hyp.new[which(hyp.new == levkp10)] <- new
            new.mod[[nb.hyps]] <- factor(hyp.new)
            names(new.mod)[nb.hyps] <- paste(n.hypi, new, sep = "_with_")
          }
          
          if(is.na(levkp11) == F){
            nb.hyps <- nb.hyps + 1
            new <- paste(levj,levk,levkp1, levkp2, levkp3, levkp4, levkp5, levkp6, levkp7, levkp8, levkp9, levkp10, levkp11, sep = " + ")
            hyp.new <- as.character(hypi)
            hyp.new[which(hyp.new == levj)] <- new
            hyp.new[which(hyp.new == levk)] <- new
            hyp.new[which(hyp.new == levkp1)] <- new
            hyp.new[which(hyp.new == levkp2)] <- new
            hyp.new[which(hyp.new == levkp3)] <- new
            hyp.new[which(hyp.new == levkp4)] <- new
            hyp.new[which(hyp.new == levkp5)] <- new
            hyp.new[which(hyp.new == levkp6)] <- new
            hyp.new[which(hyp.new == levkp7)] <- new
            hyp.new[which(hyp.new == levkp8)] <- new
            hyp.new[which(hyp.new == levkp9)] <- new
            hyp.new[which(hyp.new == levkp10)] <- new
            hyp.new[which(hyp.new == levkp11)] <- new
            new.mod[[nb.hyps]] <- factor(hyp.new)
            names(new.mod)[nb.hyps] <- paste(n.hypi, new, sep = "_with_")
          }
          
          }
          
        
      }
      
      
    }
    
  }
  mod.hyp <- new.mod
  print(names(new.mod))
  
  if(extract == T){
    
    out.gen <- data.frame(HYP = NULL, CR = NULL, Z = NULL, P.Z = NULL, pls = NULL, P.pls = NULL, ES = NULL)
    
    
  }
  
  CR.list<- list()
  
  I.list <- list()
    
    pb = txtProgressBar(min = 0, max = length(mod.hyp), initial = 0)
    
    # Brain
    ind = 0
    nb.null.int <- c()
    for(a in 1:length(mod.hyp)){
      if(length(levels(mod.hyp[[a]])) < 2){
        
        setTxtProgressBar(pb,a)
        next
      }
      
      ind <- ind + 1
      new <- phylo.modularity(A = shape, partition.gp = mod.hyp[[a]],CI = CI, phy = phy, print.progress = F)
      
      if(integration == T){
        
        indi <- "good"
        
        trying <- try(new.i <- phylo.integration(A = shape, partition.gp = mod.hyp[[a]], phy = phy, print.progress = F), silent = T)
        
        if(class(trying) == "try-error")indi <- "bad" else new.i <- phylo.integration(A = shape, partition.gp = mod.hyp[[a]], phy = phy, print.progress = F)
        
      }
      
      setTxtProgressBar(pb,a)
      
      CR.list[[ind]] <- new
      names(CR.list)[ind] <- names(mod.hyp)[a]
      
      if(integration == T & extract == T & indi == "good"){
        
        n.out <- data.frame(HYP = names(mod.hyp)[a], CR = new$CR, Z = new$Z, P.Z = new$P.value[1], pls = new.i$r.pls, P.pls = new.i$P.value[1], ES = new.i$Z)
        
      }else if(integration == T & extract == T & indi == "bad"){
        n.out <- data.frame(HYP = names(mod.hyp)[a], CR = new$CR, Z = new$Z, P.Z = new$P.value[1], pls = NA, P.pls = NA, ES = NA)
        
      }else if(integration == F & extract == T){
        
        n.out <- data.frame(HYP = names(mod.hyp)[a], CR = new$CR, Z = new$Z, P.Z = new$P.value[1], pls = NA, P.pls = NA, ES = NA)
      }
      
      if(extract == T)out.gen <- rbind(out.gen, n.out)
      
      if(integration == T & indi == "good")I.list[[ind]] <- new.i
      if(integration == T & indi == "bad"){
        I.list[[ind]] <- NA
        nb.null.int <- c(nb.null.int, ind)
        }
      if(integration == T)names(I.list)[ind] <- names(mod.hyp)[a]
    }
    
    if(comp == T){
      
      if(length(CR.list) > 1){
        
        comp.f <- compare.CR(CR.list, two.tailed = F)
        if(integration == T){
          
          if(length(nb.null.int) > 0){
            
            I.list <- I.list[-nb.null.int]
            
          }
          
          comp.i <- compare.pls(I.list, two.tailed = F)
          }
        
        CR.list$comp.CR <- comp.f
        if(integration == T)I.list$comp.pls <- comp.i
          
        
      }
      
    }
    close(pb)
    
    
    
  
  if(integration == T){
    
    out.list <- list(Modularity = CR.list, Integration = I.list)
  }else{
    
    out.list <- list(Modularity = CR.list)
  }
  
  if(extract == T)out.list$df <- out.gen
  out.list$new.hyps <- new.mod
  return(out.list)
}

# pairwise.mod.analysesV2() ####

pairwise.mod.analysesV2 <- function(shape, phy, mod.hyp, CI = T, comp = T, integration = T, extract = T){
  
  library(knitr)
  library(geomorph)
  library(progress)
  library(stringr)
  
  nb.hyps <- length(mod.hyp)
  
  all_combinations <- list()
  names_comb <- list()
  new.hyps <- list()
  num <- 0
  
  cat("","\n")
  cat("Number of hypotheses:",nb.hyps,"\n")
    
    for(i in 1:nb.hyps){
    
    cat("","\n")
    cat("Hypothesis:",names(mod.hyp)[i],"\n")
    cat("","\n")
    cat("Number of modules:",length(levels(mod.hyp[[i]])),"\n")
    if(length(levels(mod.hyp[[i]])) == 1)next
    nb.comb <- 0
    
    pb = txtProgressBar(min = 0, max = length(levels(mod.hyp[[i]])), initial = 0)
      
      for(k in 1:length(levels(mod.hyp[[i]]))){
      
      if(k == 1 | k == length(levels(mod.hyp[[i]])))next

      combk <- combn(as.character(levels(mod.hyp[[i]])), k, simplify = FALSE)
      nb.comb <- nb.comb + length(combk)

      for(j in 1:length(combk)){
        num <- num + 1
        all_combinations[[num]] <- combk[[j]]
        names_comb[[num]] <- paste0(combk[[j]], collapse  = "_and_")
        new <- as.character(mod.hyp[[i]])

        for(a in 1:length(combk[[j]])){
          
          new[which(new == combk[[j]][a])] <- names_comb[[num]]
          
        }
        
        new.hyps[[num]] <- factor(new)
        names(new.hyps)[num] <- paste(names(mod.hyp)[i],names_comb[[num]],sep = "_with_")

        if(length(levels(new.hyps[[num]])) > 2){
          
          nn <- as.character(levels(new.hyps[[num]]))
          del <- c()

          for(h in 1:length(nn)){

            str_h <- str_split_1(nn[h], pattern = "_and_")

            if(length(str_h) > 1){
              
              del <- c(del,h)
              
            }
            
          }
          
          if(length(del) > 0) nn <- nn[-del]
          if(length(nn) > 2){
            
            HYPO1 <- factor(nn)
            name.prior <- names(new.hyps)[num]

            #print("lev1")
          for(r in 1:length(HYPO1)){
            if(r == 1)next

            combk1 <- combn(as.character(HYPO1), r, simplify = FALSE)
            combk1 <- combk1

            nb.comb <- nb.comb + length(combk1)
            
            for(q in 1:length(combk1)){
              num <- num + 1
              all_combinations[[num]] <- combk1[[q]]
              names_comb[[num]] <- paste0(combk1[[q]], collapse  = "_and_")
              new1 <- as.character(new.hyps[[num - 1]])
              
              for(b in 1:length(combk1[[q]])){
                
                new1[which(new1 == combk1[[q]][b])] <- names_comb[[num]]
                
              }
              
              new.hyps[[num]] <- factor(new1)
              names(new.hyps)[num] <- paste(name.prior, names_comb[[num]], sep = "_with_")

              if(length(levels(new.hyps[[num]])) > 2){
                
                nn <- as.character(levels(new.hyps[[num]]))
                del <- c()
                
                for(h in 1:length(nn)){
                  
                  str_h <- str_split_1(nn[h], pattern = "_and_")
                  
                  if(length(str_h) > 1){
                    
                    del <- c(del,h)
                    
                  }
                  
                }
                
                if(length(del) > 0) nn <- nn[-del]
                if(length(nn) > 2){
                  
                  HYPO2 <- factor(nn)
                  name.prior <- names(new.hyps)[num]
                  #print("lev2")
                  for(f in 1:length(HYPO2)){
                    if(f == 1)next
                    
                    combk2 <- combn(as.character(HYPO2), f, simplify = FALSE)
                    nb.comb <- nb.comb + length(combk2)
                    
                    for(t in 1:length(combk2)){
                      num <- num + 1
                      all_combinations[[num]] <- combk2[[t]]
                      names_comb[[num]] <- paste0(combk2[[t]], collapse  = "_and_")
                      new2 <- as.character(new.hyps[[num - 1]])
                      
                      for(c in 1:length(combk2[[t]])){
                        
                        new2[which(new2 == combk2[[t]][c])] <- names_comb[[num]]
                        
                      }
                      
                      new.hyps[[num]] <- factor(new2)
                      names(new.hyps)[num] <- paste(name.prior, names_comb[[num]], sep = "_with_")
                      
                      if(length(levels(new.hyps[[num]])) > 2){
                        
                        nn <- as.character(levels(new.hyps[[num]]))
                        del <- c()
                        
                        for(h in 1:length(nn)){
                          
                          str_h <- str_split_1(nn[h], pattern = "_and_")
                          
                          if(length(str_h) > 1){
                            
                            del <- c(del,h)
                            
                          }
                          
                        }
                        
                        if(length(del) > 0) nn <- nn[-del]
                        if(length(nn) > 2){
                          
                          HYPO3 <- factor(nn)
                          name.prior <- names(new.hyps)[num]
                          #print("lev3")
                          for(y in 1:length(HYPO3)){
                            if(y == 1)next
                            
                            combk3 <- combn(as.character(HYPO3), y, simplify = FALSE)
                            nb.comb <- nb.comb + length(combk3)
                            
                            for(u in 1:length(combk3)){
                              num <- num + 1
                              all_combinations[[num]] <- combk3[[u]]
                              names_comb[[num]] <- paste0(combk3[[u]], collapse  = "_and_")
                              new3 <- as.character(new.hyps[[num - 1]])
                              
                              for(d in 1:length(combk3[[u]])){
                                
                                new3[which(new3 == combk3[[u]][d])] <- names_comb[[num]]
                                
                              }
                              
                              new.hyps[[num]] <- factor(new3)
                              names(new.hyps)[num] <- paste(name.prior, names_comb[[num]], sep = "_with_")
                              
                              if(length(levels(new.hyps[[num]])) > 2){
                                
                                nn <- as.character(levels(new.hyps[[num]]))
                                del <- c()
                                
                                for(h in 1:length(nn)){
                                  
                                  str_h <- str_split_1(nn[h], pattern = "_and_")
                                  
                                  if(length(str_h) > 1){
                                    
                                    del <- c(del,h)
                                    
                                  }
                                  
                                }
                                
                                if(length(del) > 0) nn <- nn[-del]
                                if(length(nn) > 2){
                                  
                                  HYPO4 <- factor(nn)
                                  name.prior <- names(new.hyps)[num]
                                  #print("lev4")
                                  for(w in 1:length(HYPO4)){
                                    if(w == 1)next
                                    
                                    combk4 <- combn(as.character(HYPO4), w, simplify = FALSE)
                                    nb.comb <- nb.comb + length(combk4)
                                    
                                    for(s in 1:length(combk4)){
                                      num <- num + 1
                                      all_combinations[[num]] <- combk4[[s]]
                                      names_comb[[num]] <- paste0(combk4[[s]], collapse  = "_and_")
                                      new4 <- as.character(new.hyps[[num - 1]])
                                      
                                      for(e in 1:length(combk4[[s]])){
                                        
                                        new4[which(new4 == combk4[[s]][e])] <- names_comb[[num]]
                                        
                                      }
                                      
                                      new.hyps[[num]] <- factor(new4)
                                      names(new.hyps)[num] <- paste(name.prior, names_comb[[num]], sep = "_with_")
                                      
                                      if(length(levels(new.hyps[[num]])) > 2){
                                        
                                        nn <- as.character(levels(new.hyps[[num]]))
                                        del <- c()
                                        
                                        for(h in 1:length(nn)){
                                          
                                          str_h <- str_split_1(nn[h], pattern = "_and_")
                                          
                                          if(length(str_h) > 1){
                                            
                                            del <- c(del,h)
                                            
                                          }
                                          
                                        }
                                        
                                        if(length(del) > 0) nn <- nn[-del]
                                        if(length(nn) > 2){
                                          
                                          HYPO5 <- factor(nn)
                                          name.prior <- names(new.hyps)[num]
                                          #print("lev5")
                                          for(x in 1:length(HYPO5)){
                                            if(x == 1)next
                                            
                                            combk5 <- combn(as.character(HYPO5), x, simplify = FALSE)
                                            nb.comb <- nb.comb + length(combk5)
                                            
                                            for(g in 1:length(combk5)){
                                              num <- num + 1
                                              all_combinations[[num]] <- combk5[[g]]
                                              names_comb[[num]] <- paste0(combk5[[g]], collapse  = "_and_")
                                              new5 <- as.character(new.hyps[[num - 1]])
                                              
                                              for(l in 1:length(combk5[[g]])){
                                                
                                                new5[which(new5 == combk5[[g]][l])] <- names_comb[[num]]
                                                
                                              }
                                              
                                              new.hyps[[num]] <- factor(new5)
                                              names(new.hyps)[num] <- paste(name.prior, names_comb[[num]], sep = "_with_")
                                              
                                              if(length(levels(new.hyps[[num]])) > 2){
                                                
                                                nn <- as.character(levels(new.hyps[[num]]))
                                                del <- c()
                                                
                                                for(h in 1:length(nn)){
                                                  
                                                  str_h <- str_split_1(nn[h], pattern = "_and_")
                                                  
                                                  if(length(str_h) > 1){
                                                    
                                                    del <- c(del,h)
                                                    
                                                  }
                                                  
                                                }
                                                
                                                if(length(del) > 0) nn <- nn[-del]
                                                if(length(nn) > 2){
                                                  
                                                  HYPO6 <- factor(nn)
                                                  name.prior <- names(new.hyps)[num]
                                                  #print("lev6")
                                                  for(w in 1:length(HYPO6)){
                                                    if(w == 1)next
                                                    
                                                    combk6 <- combn(as.character(HYPO6), w, simplify = FALSE)
                                                    nb.comb <- nb.comb + length(combk6)
                                                    
                                                    for(z in 1:length(combk6)){
                                                      num <- num + 1
                                                      all_combinations[[num]] <- combk6[[z]]
                                                      names_comb[[num]] <- paste0(combk6[[z]], collapse  = "_and_")
                                                      new6 <- as.character(new.hyps[[num - 1]])
                                                      
                                                      for(x in 1:length(combk6[[z]])){
                                                        
                                                        new6[which(new6 == combk6[[z]][x])] <- names_comb[[num]]
                                                        
                                                      }
                                                      
                                                      new.hyps[[num]] <- factor(new6)
                                                      names(new.hyps)[num] <- paste(name.prior, names_comb[[num]], sep = "_with_")
                                                      
                                                      
                                                    }
                                                  }
                                                }
                                                
                                              }else{
                                                next
                                              }
                                            }
                                            
                                          }
                                        }
                                        
                                      }else{
                                        next
                                      }
                                    }
                                    
                                  }
                                }
                                
                              }else{
                                next
                              }
                            }
                            
                          }
                        }
                        
                      }else{
                        next
                      }
                    }
                    
                  }
                }
                
              }else{
                next
              }
            }
            
            }
          }
          
        }else{
          next
        }
      }
      setTxtProgressBar(pb,k)
    }
      close(pb)
    
    
    cat("","\n")
    cat("Number of combination excluding m = 1 and m = max:",nb.comb,"\n")
    cat("","\n")
  }
  
  

  mod.hyp <- new.hyps
  print(length(new.hyps))
  if(extract == T){
    
    out.gen <- data.frame(HYP = NULL, CR = NULL, Z = NULL, P.Z = NULL, pls = NULL, P.pls = NULL, ES = NULL)
    
    
  }
  
  CR.list<- list()
  
  I.list <- list()
  
  pb = txtProgressBar(min = 0, max = length(mod.hyp), initial = 0)
  
  # Brain
  ind = 0
  nb.null.int <- c()
  for(a in 1:length(mod.hyp)){
    if(length(levels(mod.hyp[[a]])) < 2){
      
      setTxtProgressBar(pb,a)
      next
    }
    
    ind <- ind + 1
    new <- phylo.modularity(A = shape, partition.gp = mod.hyp[[a]],CI = CI, phy = phy, print.progress = F)
    
    if(integration == T){
      
      indi <- "good"
      
      trying <- try(new.i <- phylo.integration(A = shape, partition.gp = mod.hyp[[a]], phy = phy, print.progress = F), silent = T)
      
      if(class(trying) == "try-error")indi <- "bad" else new.i <- phylo.integration(A = shape, partition.gp = mod.hyp[[a]], phy = phy, print.progress = F)
      
    }
    
    setTxtProgressBar(pb,a)
    
    CR.list[[ind]] <- new
    names(CR.list)[ind] <- names(mod.hyp)[a]
    
    if(integration == T & extract == T & indi == "good"){
      
      n.out <- data.frame(HYP = names(mod.hyp)[a], CR = new$CR, Z = new$Z, P.Z = new$P.value[1], pls = new.i$r.pls, P.pls = new.i$P.value[1], ES = new.i$Z)
      
    }else if(integration == T & extract == T & indi == "bad"){
      n.out <- data.frame(HYP = names(mod.hyp)[a], CR = new$CR, Z = new$Z, P.Z = new$P.value[1], pls = NA, P.pls = NA, ES = NA)
      
    }else if(integration == F & extract == T){
      
      n.out <- data.frame(HYP = names(mod.hyp)[a], CR = new$CR, Z = new$Z, P.Z = new$P.value[1], pls = NA, P.pls = NA, ES = NA)
    }
    
    if(extract == T)out.gen <- rbind(out.gen, n.out)
    
    if(integration == T & indi == "good")I.list[[ind]] <- new.i
    if(integration == T & indi == "bad"){
      I.list[[ind]] <- NA
      nb.null.int <- c(nb.null.int, ind)
    }
    if(integration == T)names(I.list)[ind] <- names(mod.hyp)[a]
  }
  
  if(comp == T){
    
    if(length(CR.list) > 1){
      
      comp.f <- compare.CR(CR.list, two.tailed = F)
      if(integration == T){
        
        if(length(nb.null.int) > 0){
          
          I.list <- I.list[-nb.null.int]
          
        }
        
        comp.i <- compare.pls(I.list, two.tailed = F)
      }
      
      CR.list$comp.CR <- comp.f
      if(integration == T)I.list$comp.pls <- comp.i
      
      
    }
    
  }
  close(pb)
  
  
  
  
  if(integration == T){
    
    out.list <- list(Modularity = CR.list, Integration = I.list)
  }else{
    
    out.list <- list(Modularity = CR.list)
  }
  
  if(extract == T)out.list$df <- out.gen
  out.list$new.hyps <- new.mod
  return(out.list)
}

# generate.comb.mod.hyp() ####

generate.comb.mod.hyp <- function(mod.hyp){
  
  library(knitr)
  library(geomorph)
  library(progress)
  library(stringr)
  
  nb.hyps <- length(mod.hyp)
  
  all_combinations <- list()
  names_comb <- list()
  new.hyps <- list()
  num <- 0
  
  cat("","\n")
  cat("Number of hypotheses:",nb.hyps,"\n")
  
  for(i in 1:nb.hyps){
    
    cat("","\n")
    cat("Hypothesis:",names(mod.hyp)[i],"\n")
    cat("","\n")
    cat("Number of modules:",length(levels(mod.hyp[[i]])),"\n")
    if(length(levels(mod.hyp[[i]])) == 1)next
    nb.comb <- 0
    
    num <- num + 1
    new.hyps[[num]] <- mod.hyp[[i]]
    names(new.hyps)[num] <- names(mod.hyp)[i]
    nb.comb <- nb.comb + 1
    
    pb = txtProgressBar(min = 0, max = length(levels(mod.hyp[[i]])), initial = 0)
    
    for(k in 1:length(levels(mod.hyp[[i]]))){
      
      if(k == 1 | k == length(levels(mod.hyp[[i]])))next
      
      combk <- combn(as.character(levels(mod.hyp[[i]])), k, simplify = FALSE)
      nb.comb <- nb.comb + length(combk)
      
      for(j in 1:length(combk)){
        num <- num + 1
        all_combinations[[num]] <- combk[[j]]
        names_comb[[num]] <- paste0(combk[[j]], collapse  = "_and_")
        new <- as.character(mod.hyp[[i]])
        
        for(a in 1:length(combk[[j]])){
          
          new[which(new == combk[[j]][a])] <- names_comb[[num]]
          
        }
        
        new.hyps[[num]] <- factor(new)
        names(new.hyps)[num] <- paste(names(mod.hyp)[i],names_comb[[num]],sep = "_with_")
        
        if(length(levels(new.hyps[[num]])) > 2){
          
          nn <- as.character(levels(new.hyps[[num]]))
          del <- c()
          
          for(h in 1:length(nn)){
            
            str_h <- str_split_1(nn[h], pattern = "_and_")
            
            if(length(str_h) > 1){
              
              del <- c(del,h)
              
            }
            
          }
          
          if(length(del) > 0) nn <- nn[-del]
          if(length(nn) > 2){
            
            HYPO1 <- factor(nn)
            name.prior <- names(new.hyps)[num]
            
            #print("lev1")
            for(r in 1:length(HYPO1)){
              if(r == 1)next
              
              combk1 <- combn(as.character(HYPO1), r, simplify = FALSE)
              combk1 <- combk1
              
              nb.comb <- nb.comb + length(combk1)
              
              for(q in 1:length(combk1)){
                num <- num + 1
                all_combinations[[num]] <- combk1[[q]]
                names_comb[[num]] <- paste0(combk1[[q]], collapse  = "_and_")
                new1 <- as.character(new.hyps[[num - 1]])
                
                for(b in 1:length(combk1[[q]])){
                  
                  new1[which(new1 == combk1[[q]][b])] <- names_comb[[num]]
                  
                }
                
                new.hyps[[num]] <- factor(new1)
                names(new.hyps)[num] <- paste(name.prior, names_comb[[num]], sep = "_with_")
                
                
              }
              
            }
          }
          
        }else{
          next
        }
      }
      setTxtProgressBar(pb,k)
    }
    close(pb)
    
    
    cat("","\n")
    cat("Number of combination excluding m = 1 and m = max:",nb.comb,"\n")
    cat("","\n")
  }
  cat("","\n")
  cat("Number of hypotheses to test", length(new.hyps), "\n")
  cat("","\n")
  cat("End of line","\n")
  cat("","\n")
  return(new.hyps)
}

# generate.comb.mod.hyp.comp() ####

generate.mod.hyp.comp <- function(mod.hyp.1, mod.hyp.2){
  
  library(knitr)
  library(geomorph)
  library(progress)
  library(stringr)
  
  nb.hyps.1 <- length(mod.hyp.1)
  nb.hyps.2 <- length(mod.hyp.2)
  
  new.hyps <- list()
  num <- 0
  
  cat("","\n")
  cat("Number of hypotheses in group 1:",nb.hyps.1,"\n")
  cat("Number of hypotheses in group 2:",nb.hyps.2,"\n")
  
  for(i in 1:nb.hyps.1){
    
    hyp1 <- mod.hyp.1[[i]]
    nb.mod.1 <- length(levels(hyp1))
    
    for(j in 1:nb.hyps.2){
      
      hyp2 <- mod.hyp.2[[j]]
      nb.mod.2 <- length(levels(hyp2))
      name.hypsij <- paste(names(mod.hyp.1)[i], names(mod.hyp.2)[j])
      
      num <- num + 1
      new.hyp <- c(mod.hyp.1[[i]], mod.hyp.2[[j]])
      new.hyp <- factor(new.hyp)
      new.hyps[[num]] <- new.hyp
      names(new.hyps)[num] <- name.hypsij
      
      for(k in 1:nb.mod.1){
        
        modk1 <- levels(hyp1)[k]
        
        for(l in 1:nb.mod.2){
          
          num <- num + 1
          
          modl2 <- levels(hyp2)[l]
          
          new.mod <- paste(modk1, modl2, sep = "_AND_")
          
          new.hyp <- as.character(c(mod.hyp.1[[i]], mod.hyp.2[[j]]))
          new.hyp[which(new.hyp == modk1)] <- new.mod
          new.hyp[which(new.hyp == modl2)] <- new.mod
          new.hyp <- factor(new.hyp)
          
          new.hyps[[num]] <- new.hyp
          names(new.hyps)[num] <- paste(name.hypsij, new.mod, sep = "_WITH_")
        }
        
      }
      
    }
    
  }
  
  cat("","\n")
  cat("Number of hypotheses to test", length(new.hyps), "\n")
  cat("","\n")
  cat("End of line","\n")
  cat("","\n")
  return(new.hyps)
}

# generate.comb.mod.hyp.compV2() ####

generate.mod.hyp.compV2 <- function(mod.hyp.1, mod.hyp.2, depth = 6){
  
  library(knitr)
  library(geomorph)
  library(progress)
  library(stringr)
  library(rapportools)
  
  nb.hyps.1 <- length(mod.hyp.1)
  nb.hyps.2 <- length(mod.hyp.2)
  
  new.hyps <- list()
  num <- 0
  
  cat("","\n")
  cat("Number of hypotheses in group 1:",nb.hyps.1,"\n")
  cat("Number of hypotheses in group 2:",nb.hyps.2,"\n")
  
  for(i in 1:nb.hyps.1){
    
    hyp1 <- mod.hyp.1[[i]]
    nb.mod.1 <- length(levels(hyp1))
    
    for(j in 1:nb.hyps.2){
      
      hyp2 <- mod.hyp.2[[j]]

      nb.mod.2 <- length(levels(hyp2))
      name.hypsij <- paste(names(mod.hyp.1)[i], names(mod.hyp.2)[j])

      num <- num + 1
      new.hyp <- c(mod.hyp.1[[i]], mod.hyp.2[[j]])
      new.hyp <- factor(new.hyp)
      new.hyps[[num]] <- new.hyp
      names(new.hyps)[num] <- name.hypsij
      
      for(k in 1:nb.mod.1){
        
        modk1 <- levels(hyp1)[k]

        for(l in 1:nb.mod.2){
          
          num <- num + 1
          
          modl2 <- levels(hyp2)[l]

          new.mod <- paste(modk1, modl2, sep = "_AND_")

          new.hyp <- as.character(c(mod.hyp.1[[i]], mod.hyp.2[[j]]))
          new.hyp[which(new.hyp == modk1)] <- new.mod
          new.hyp[which(new.hyp == modl2)] <- new.mod
          new.hyp <- factor(new.hyp)

          
          new.hyps[[num]] <- new.hyp
          names(new.hyps)[num] <- paste(name.hypsij, new.mod, sep = "_WITH_")
          
          rest.1 <- levels(mod.hyp.1[[i]])[-which(levels(mod.hyp.1[[i]]) == modk1)]
          rest.2 <- levels(mod.hyp.2[[j]])[-which(levels(mod.hyp.2[[j]]) == modl2)]
          
          if(length(rest.1) >= 1 & length(rest.2) >= 1 & depth >= 2){
            
            for(g in 1:length(rest.1)){
              
              modg1 <- rest.1[g]
              if(is.empty(modg1) == T)next
              if(is.na(modg1) == T)next
              for(h in 1:length(rest.2)){
                
                modh2 <- rest.2[h]
                if(is.empty(modh2) == T)next
                if(is.na(modh2) == T)next
                num <- num + 1
                
                new.mod <- paste(modg1, modh2, sep = "_AND_")
                
                new.hyp.2 <- as.character(new.hyp)
                new.hyp.2[which(new.hyp.2 == modg1)] <- new.mod
                new.hyp.2[which(new.hyp.2 == modh2)] <- new.mod
                new.hyp.2 <- factor(new.hyp.2)
                
                new.hyps[[num]] <- new.hyp.2
                names(new.hyps)[num] <- paste(names(new.hyps)[num - 1], new.mod, sep = "_WITH_")
                
                rest.1 <- rest.1[-which(rest.1 == modg1)]
                rest.2 <- rest.2[-which(rest.2 == modh2)]
                
                if(length(rest.1) >= 1 & length(rest.2) >= 1 & depth >= 3){
                  
                  for(q in 1:length(rest.1)){
                    
                    modq1 <- rest.1[q]
                    if(is.empty(modq1) == T)next
                    if(is.na(modq1) == T)next
                    for(w in 1:length(rest.2)){
                      
                      modw2 <- rest.2[w]
                      if(is.empty(modw2) == T)next
                      if(is.na(modw2) == T)next
                      num <- num + 1
                      
                      new.mod <- paste(modq1, modw2, sep = "_AND_")
                      
                      new.hyp.3 <- as.character(new.hyp)
                      new.hyp.3[which(new.hyp.3 == modq1)] <- new.mod
                      new.hyp.3[which(new.hyp.3 == modw2)] <- new.mod
                      new.hyp.3 <- factor(new.hyp.3)
                      
                      new.hyps[[num]] <- new.hyp.3
                      names(new.hyps)[num] <- paste(names(new.hyps)[num - 1], new.mod, sep = "_WITH_")
                      
                      rest.1 <- rest.1[-which(rest.1 == modq1)]
                      rest.2 <- rest.2[-which(rest.2 == modw2)]
                      
                      if(length(rest.1) >= 1 & length(rest.2) >= 1 & depth >= 4){
                        
                        for(e in 1:length(rest.1)){
                          
                          mode1 <- rest.1[e]
                          if(is.empty(mode1) == T)next
                          if(is.na(mode1) == T)next
                          for(r in 1:length(rest.2)){
                            
                            modr2 <- rest.2[r]
                            if(is.empty(modr2) == T)next
                            if(is.na(modr2) == T)next
                            num <- num + 1
                            
                            new.mod <- paste(mode1, modr2, sep = "_AND_")
                            
                            new.hyp.4 <- as.character(new.hyp)
                            new.hyp.4[which(new.hyp.4 == mode1)] <- new.mod
                            new.hyp.4[which(new.hyp.4 == modr2)] <- new.mod
                            new.hyp.4 <- factor(new.hyp.4)
                            
                            new.hyps[[num]] <- new.hyp.4
                            names(new.hyps)[num] <- paste(names(new.hyps)[num - 1], new.mod, sep = "_WITH_")
                            
                            rest.1 <- rest.1[-which(rest.1 == mode1)]
                            rest.2 <- rest.2[-which(rest.2 == modr2)]
                            
                            if(length(rest.1) >= 1 & length(rest.2) >= 1 & depth >= 5){
                              
                              for(t in 1:length(rest.1)){
                                
                                modt1 <- rest.1[t]
                                if(is.empty(modt1) == T)next
                                if(is.na(modt1) == T)next
                                for(y in 1:length(rest.2)){
                                  
                                  mody2 <- rest.2[y]
                                  if(is.empty(mody2) == T)next
                                  if(is.na(mody2) == T)next
                                  num <- num + 1
                                  
                                  new.mod <- paste(modt1, mody2, sep = "_AND_")
                                  
                                  new.hyp.5 <- as.character(new.hyp)
                                  new.hyp.5[which(new.hyp.5 == modt1)] <- new.mod
                                  new.hyp.5[which(new.hyp.5 == mody2)] <- new.mod
                                  new.hyp.5 <- factor(new.hyp.5)
                                  
                                  new.hyps[[num]] <- new.hyp.5
                                  names(new.hyps)[num] <- paste(names(new.hyps)[num - 1], new.mod, sep = "_WITH_")
                                  
                                  rest.1 <- rest.1[-which(rest.1 == modt1)]
                                  rest.2 <- rest.2[-which(rest.2 == mody2)]
                                  
                                  if(length(rest.1) >= 1 & length(rest.2) >= 1 & depth >= 6){
                                    
                                    for(o in 1:length(rest.1)){
                                      
                                      modo1 <- rest.1[o]
                                      if(is.empty(modo1) == T)next
                                      if(is.na(modo1) == T)next
                                      for(p in 1:length(rest.2)){
                                        
                                        modp2 <- rest.2[p]
                                        if(is.empty(modp2) == T)next
                                        if(is.na(modp2) == T)next
                                        num <- num + 1
                                        
                                        new.mod <- paste(modo1, modp2, sep = "_AND_")
                                        
                                        new.hyp.6 <- as.character(new.hyp)
                                        new.hyp.6[which(new.hyp.6 == modo1)] <- new.mod
                                        new.hyp.6[which(new.hyp.6 == modp2)] <- new.mod
                                        new.hyp.6 <- factor(new.hyp.6)
                                        
                                        new.hyps[[num]] <- new.hyp.6
                                        names(new.hyps)[num] <- paste(names(new.hyps)[num - 1], new.mod, sep = "_WITH_")
                                        
                                      }
                                      
                                    }
                                    
                                  }
                                  
                                }
                                
                              }
                              
                            }
                            
                          }
                          
                        }
                        
                      }
                      
                    }
                    
                  }
                  
                }
                
              }
              
            }
            
          }
        }
        
      }
      
    }
    
  }
  
  cat("","\n")
  cat("Number of hypotheses to test", length(new.hyps), "\n")
  cat("","\n")
  cat("End of line","\n")
  cat("","\n")
  return(new.hyps)
}



# pairwise.mod.analysesV3() ####

pairwise.mod.analysesV3 <- function(shape, phy, mod.hyp, CI = T, comp = T, integration = T, extract = T){
  
  library(knitr)
  library(geomorph)
  library(progress)
  library(stringr)
  
  if(extract == T){
    
    out.gen <- data.frame(HYP = NULL, CR = NULL, Z = NULL, P.Z = NULL, pls = NULL, P.pls = NULL, ES = NULL)
    
    
  }
  
  CR.list<- list()
  
  I.list <- list()
  
  cat("","\n")
  cat("Initiating modularity and morphological integration analyses","\n")
  cat("","\n")
  pb = txtProgressBar(min = 0, max = length(mod.hyp), initial = 0)
  
  # Brain
  ind = 0
  nb.null.int <- c()
  for(a in 1:length(mod.hyp)){
    if(length(levels(mod.hyp[[a]])) < 2){
      
      setTxtProgressBar(pb,a)
      next
    }
    
    ind <- ind + 1
    new <- phylo.modularity(A = shape, partition.gp = mod.hyp[[a]],CI = CI, phy = phy, print.progress = F)
    
    if(integration == T){
      
      indi <- "good"
      
      trying <- try(new.i <- phylo.integration(A = shape, partition.gp = mod.hyp[[a]], phy = phy, print.progress = F), silent = T)
      
      if(class(trying) == "try-error")indi <- "bad" else new.i <- phylo.integration(A = shape, partition.gp = mod.hyp[[a]], phy = phy, print.progress = F)
      
    }
    
    setTxtProgressBar(pb,a)
    
    CR.list[[ind]] <- new
    names(CR.list)[ind] <- names(mod.hyp)[a]
    CR.list <-CR.list
    
    if(integration == T & extract == T & indi == "good"){
      
      n.out <- data.frame(HYP = names(mod.hyp)[a], CR = new$CR, Z = new$Z, P.Z = new$P.value[1], pls = new.i$r.pls, P.pls = new.i$P.value[1], ES = new.i$Z)
      
    }else if(integration == T & extract == T & indi == "bad"){
      n.out <- data.frame(HYP = names(mod.hyp)[a], CR = new$CR, Z = new$Z, P.Z = new$P.value[1], pls = NA, P.pls = NA, ES = NA)
      
    }else if(integration == F & extract == T){
      
      n.out <- data.frame(HYP = names(mod.hyp)[a], CR = new$CR, Z = new$Z, P.Z = new$P.value[1], pls = NA, P.pls = NA, ES = NA)
    }
    
    if(extract == T)out.gen <- rbind(out.gen, n.out)
    
    if(integration == T & indi == "good"){I.list[[ind]] <- new.i}
    if(integration == T & indi == "bad"){
      I.list[[ind]] <- NA
      nb.null.int <- c(nb.null.int, ind)
    }
    if(integration == T){names(I.list)[ind] <- names(mod.hyp)[a]}
  }
  
  if(comp == T){
    
    if(length(CR.list) > 1){
      
      comp.f <- compare.CR(CR.list, two.tailed = F)
      if(integration == T){
        
        if(length(nb.null.int) > 0){
          
          I.list <- I.list[-nb.null.int]
          
        }
        
        comp.i <- compare.pls(I.list, two.tailed = F)
      }
      
      CR.list$comp.CR <- comp.f
      if(integration == T){
        I.list$comp.pls <- comp.i
        I.list <- I.list
        }
      
      
    }
    
  }
  close(pb)
  
  cat("","\n")
  cat("Creating output object","\n")
  cat("","\n")
  
  
  if(integration == T){
    
    out.list <- list(Modularity = CR.list, Integration = I.list)
  }else{
    
    out.list <- list(Modularity = CR.list)
  }
  
  if(extract == T)out.list$df <- out.gen
  cat("","\n")
  cat("End of line","\n")
  cat("","\n")
  return(out.list)
}

# read_all_trees() ####

read_all_trees <- function(filenames, keep.multi = T, name = F, keep = T, tips = NULL){
  
  all <- list()
  
  cat("","\n")
  cat("- Importing")
  if(keep == T)cat(" and pruning -", "\n") else cat(" -","\n")
  
  pb = txtProgressBar(min = 0, max = length(filenames), initial = 0)
  
  for(i in 1:length(filenames)){
    
    new <- read.tree(file = filenames[i], keep.multi = keep.multi)
    
    if(keep == T)new <- keep.tip(new, tip = tips)
    
    setTxtProgressBar(pb,i)
    
    all[[i]] <- new
    if(name == T)names(all)[i] <- filenames[i]
    
  }
  
  close(pb)
  
  return(all)
  
}

add.consensus.allom <- function(resids, consensus, nb.land, nb.dims){
  
  nb.ind <- length(dimnames(resids)[[3]])
  
  vec <- c()
  for(i in 1:nb.ind){
    
    mat <- resids[,,i] + consensus
    
    vec <- c(vec, mat)
    
  }
  
  out <- array(vec, dim = c(nb.land,nb.dims, nb.ind))
  
  return(out)
}

# pairwise_congruence() ####

pairwise_congruence <- function(aligned_landmarks) {
  n_landmarks <- dim(aligned_landmarks)[1]
  n_specimens <- dim(aligned_landmarks)[3]
  
  # Initialize a matrix to store the congruence correlation values
  congruence_matrix <- matrix(NA, n_landmarks, n_landmarks)
  
  # Loop over all pairs of landmarks
  for (i in 1:(n_landmarks - 1)) {
    for (j in (i + 1):n_landmarks) {
      # Extract the x, y, and z coordinates across all specimens for landmarks i and j
      vec_i <- c(aligned_landmarks[i, 1, ], aligned_landmarks[i, 2, ], aligned_landmarks[i, 3, ])
      vec_j <- c(aligned_landmarks[j, 1, ], aligned_landmarks[j, 2, ], aligned_landmarks[j, 3, ])
      
      # Calculate the correlation between the vectors (all x, y, z coordinates across specimens)
      corr <- cor(vec_i, vec_j)
      
      # Fill in the congruence matrix symmetrically
      congruence_matrix[i, j] <- corr
      congruence_matrix[j, i] <- corr
    }
  }
  
  # Set diagonal to 1 (self-correlation)
  diag(congruence_matrix) <- 1
  return(congruence_matrix)
}


# unique.names.array() ####

unique.names.array <- function(shape, ID){
  
  for(i in 1:length(shape[1,1,])){
    
    dni <- dimnames(shape)[[3]][i]
    
    if(i < length(shape[1,1,])){
      
      for(j in (i+1):length(shape[1,1,])){
        
        dnj <- dimnames(shape)[[3]][j]
        
        if(dni == dnj){
          
          dimnames(shape)[[3]][i] <- paste0(dni,ID[i])
          dimnames(shape)[[3]][j] <- paste0(dnj,ID[j])
          
        }
        
      }
    }
    
    
  }
  return(shape)
}
 
# qgraph.mod.int() ####
qgraph.mod.int <- function(mat, levels, export = "pdf", prefix = getwd(), hyp.name, width = 5, height = 5, plot = T){
  library(qgraph)
  
  if(export == "pdf"){
    
    filename <- paste0(prefix,hyp.name,".pdf")
    filename.leg <- paste0(prefix,hyp.name,"_leg.pdf")
    pdf(file = filename, width = width, height = height)
  }
  (res <- qgraph(mat, 
                 layout = "spring", 
                 groups = levels, 
                 labels = F,
                 label.cex = 3,
                 edge.color = "grey50", 
                 palette = "colorblind",
                 vsize = 5,
                 esize = 15,
                 legend = F))
  if(export == "pdf")dev.off()
  if(export == "pdf")pdf(file = filename.leg, width = width, height = height)
  (res.leg <- qgraph(mat, 
                     layout = "spring", 
                     label.cex = 3,
                     edge.color = "grey50", 
                     vsize = 5,
                     esize = 15,
                     legend = F))
  if(export == "pdf")dev.off()
  
  if(plot == T){
    par(mfrow = c(1, 2))
    plot(res)
    plot(res.leg)
  }
}

# qgraph.mod.mult() ####

qgraph.mod.mult <- function(res.list, all.levels, hyp.selected, export = "pdf", prefix = getwd(), width = 5, height = 5, plot = T, EMMLi = F, valeur.vsize = 3, vsize.scale = 2, vsize = NULL, treshold = 0.5, color = NULL){
  library(qgraph)
  if(EMMLi == F){
    for(i in 1:length(hyp.selected)){
      print(hyp.selected[i])
      tmat <- try(res.list[[hyp.selected[i]]][["CR.mat"]],silent = T)
      
      tlev <- try(levels(all.levels[[hyp.selected[i]]]), silent = T)
      if(class(tmat) == "NULL"){
        print("ERROR in hyp.selected, no CR.mat found")
        next
      }  
      if(class(tlev) == "NULL"){
        print("ERROR in hyp.selected, no CR.mat found")
        next
      }
      mat <- res.list[[hyp.selected[i]]][["CR.mat"]]
      max.val <- max(mat)
      min.val <- min(mat)
      ntreshold <- ((max.val - min.val) * treshold) + min.val
      mat[which(mat <= ntreshold)] <- 0.01
      mat <- as.matrix(mat)[levels(all.levels[[hyp.selected[i]]]), levels(all.levels[[hyp.selected[i]]])]
      mat <- as.dist(mat)

      
      levels.i <- levels(all.levels[[hyp.selected[i]]])
      
      if(class(color) == "list"){
        color.i <- color[[hyp.selected[i]]]
        
        color.i <- color.i[match(levels.i,names(color.i))]
        
      }else{color.i <- NULL}
      
      if(export == "pdf"){
        
        filename <- paste0(prefix,hyp.selected[i],".pdf")
        filename.leg <- paste0(prefix,hyp.selected[i],"_leg.pdf")
        pdf(file = filename, width = width, height = height)
      }
      (res <- qgraph(mat, 
                     layout = "circle", 
                     groups = levels.i, 
                     labels = F,
                     label.cex = 3,
                     edge.color = "grey50", 
                     color = color.i,
                     vsize = 5,
                     esize = 15,
                     legend = F))
      if(export == "pdf")dev.off()
      if(export == "pdf")pdf(file = filename.leg, width = width, height = height)
      (res.leg <- qgraph(mat, 
                         layout = "circle", 
                         label.cex = 3,
                         edge.color = "grey50", 
                         vsize = 5,
                         esize = 15,
                         legend = F))
      if(export == "pdf")dev.off()
      
      if(plot == T){
        par(mfrow = c(1, 2))
        plot(res)
        plot(res.leg)
      }
      
    }
  }
  if(EMMLi == T){
    
    for(i in 1:length(hyp.selected)){
      
      hyp.n <- names(hyp.selected)[i]
      
      tab.i <- res.list[[hyp.selected[i]]][["corr.modules.df"]]
      
      matrix <- corr.mod(tab.i)
      
      if(class(matrix)[1] == "character") next
      colnames(matrix) <- levels(all.levels[[hyp.selected[i]]])
      rownames(matrix) <- levels(all.levels[[hyp.selected[i]]])
      mat <- as.dist(matrix)
      max.val <- max(mat)
      min.val <- min(mat)

      ntreshold <- ((max.val - min.val) * treshold) + min.val

      mat[which(mat <= ntreshold)] <- 0.001
      
      mat <- as.matrix(mat)[levels(all.levels[[hyp.selected[i]]]), levels(all.levels[[hyp.selected[i]]])]
      mat <- as.dist(mat)


      
      levels.i <- levels(all.levels[[hyp.selected[i]]])
      
      vsize.mod.n <- str_split_i(colnames(vsize[[i]][["corr.modules.df"]]),pattern = "Module ", i = 2)
      val <- c()
      val.corr <- c()
      for(j in 1:length(vsize.mod.n)){
        
        if(is.na(vsize.mod.n[[j]]) == F){
          
          val <- c(val, as.numeric(vsize.mod.n[[j]]))
          val.corr <- c(val.corr, as.numeric(vsize[[i]][["corr.modules.df"]][2,j]))
          
        }
        
      }
      
      df.vsize <- data.frame(mod = val, corr = val.corr)
      df.vsize <- df.vsize[order(df.vsize$mod),]
      minvsize <- min(df.vsize$corr)
      maxvsize <- max(df.vsize$corr)
      df.vsize$corr <- (((df.vsize$corr - min(df.vsize$corr))/(max(df.vsize$corr) - min(df.vsize$corr)))*vsize.scale)+valeur.vsize
      
      if(class(color) == "list"){
        color.i <- color[[hyp.selected[i]]]
        
        color.i <- color.i[match(levels.i,names(color.i))]
        
      }else{color.i <- NULL}
      
      if(export == "pdf"){
        
        filename <- paste0(prefix,hyp.selected[i],".pdf")
        filename.leg <- paste0(prefix,hyp.selected[i],"_leg.pdf")
        pdf(file = filename, width = width, height = height)
      }
      (res <- qgraph(mat, 
                     layout = "circle", 
                     groups = levels.i, 
                     labels = F,
                     label.cex = 3,
                     edge.color = "grey50", 
                     color = color.i,
                     vsize = df.vsize$corr,
                     esize = 15,
                     legend = F))
      if(export == "pdf")dev.off()
      if(export == "pdf")pdf(file = filename.leg, width = width, height = height)
      (res.leg <- qgraph(mat, 
                         layout = "circle", 
                         label.cex = 3,
                         edge.color = "grey50", 
                         vsize = df.vsize$corr,
                         esize = 15,
                         legend = F))
      if(export == "pdf")dev.off()
      
      if(plot == T){
        par(mfrow = c(1, 2))
        plot(res)
        plot(res.leg)
      }
    }
    
  }
  
}

# qgraph.int.mult() ####

qgraph.int.mult <- function(res.list, all.levels, hyp.selected, export = "pdf", prefix = getwd(), width = 5, height = 5, plot = T, treshold = 0.5, color = NULL){
  library(qgraph)
  
  for(i in 1:length(hyp.selected)){
    print(hyp.selected[i])
    tmat <- try(res.list[[hyp.selected[i]]][["r.pls.mat"]],silent = T)
    
    tlev <- try(levels(all.levels[[hyp.selected[i]]]), silent = T)
    if(class(tmat) == "NULL"){
      print("ERROR in hyp.selected, no r.pls.mat found")
      next
    }  
    if(class(tlev) == "NULL"){
      print("ERROR in hyp.selected, no r.pls.mat found")
      next
    }  
    mat <- res.list[[hyp.selected[i]]][["r.pls.mat"]]
    max.val <- max(mat)
    min.val <- min(mat)
    ntreshold <- ((max.val - min.val) * treshold) + min.val
    mat[which(mat <= ntreshold)] <- 0.01
    mat <- as.matrix(mat)[levels(all.levels[[hyp.selected[i]]]), levels(all.levels[[hyp.selected[i]]])]
    mat <- as.dist(mat)

    levels.i <- levels(all.levels[[hyp.selected[i]]])
    
    if(class(color) == "list"){
      color.i <- color[[hyp.selected[i]]]
      
      color.i <- color.i[match(levels.i,names(color.i))]
      
    }else{color.i <- NULL}
    
    if(export == "pdf"){
      
      filename <- paste0(prefix,hyp.selected[i],".pdf")
      filename.leg <- paste0(prefix,hyp.selected[i],"_leg.pdf")
      pdf(file = filename, width = width, height = height)
    }
    (res <- qgraph(mat, 
                   layout = "circle", 
                   groups = levels.i, 
                   labels = F,
                   label.cex = 3,
                   edge.color = "grey50", 
                   color = color.i,
                   vsize = 5,
                   esize = 15,
                   legend = F))
    if(export == "pdf")dev.off()
    if(export == "pdf")pdf(file = filename.leg, width = width, height = height)
    (res.leg <- qgraph(mat, 
                       layout = "circle", 
                       label.cex = 3,
                       edge.color = "grey50", 
                       vsize = 5,
                       esize = 15,
                       legend = F))
    if(export == "pdf")dev.off()
    
    if(plot == T){
      par(mfrow = c(1, 2))
      plot(res)
      plot(res.leg)
    }
    
  }
  
}

# EMMli.mod() ####

EMMLi.mod <- function(shape, N_sample, landmark.names, hypotheses, show.comparison = F, mult = F, EMMLiv2 = F, phylo = NULL){
  library(stringr)
  
  corr <- as.data.frame(pairwise_congruence(shape))
  
  
  if(mult == F){
    hyps <- data.frame(name = landmark.names)
    k = 0
    nnam <- c()
    hyp.nam <- c()
    for(i in 1:length(hypotheses)){
      if(length(levels(hypotheses[[i]])) > 1){
        k <- k+1
        new <- paste0("mod",k)
        nn <- names(hypotheses)[i]
        hyp.nam <- c(hyp.nam, nn)
        nnam <- c(nnam, new)
        hyps <- cbind(hyps, hypotheses[[i]])
      }
      
      
    }
    
    colnames(hyps) <- c("names", nnam)
    cat("","\n")
    cat("Performing EMMLi analysis","\n")
    cat("","\n")
    if(EMMLiv2 == T){
      mods <- phyloEmmli(shape, phylo, N_sample = N_sample, mod = hyps, EMMLi = T)
    }else{
      mods <- EMMLi(corr, N_sample = N_sample, mod = hyps)
    }
    
    model.comp <- as.data.frame(mods$EMMLi$results)
    model.comp <- model.comp[order(model.comp$dAICc),]
    winner <- str_split_1(names(mods$EMMLi$rho)[1], pattern = "mod")[2]
    winner <- as.numeric(str_split_1(winner, pattern = ".sep")[1])
    name.winner <- hyp.nam[winner]
    cat("","\n")
    cat("Best supported hypothesis:",name.winner,"\n")
    cat("","\n")
    
    if(length(mods$EMMLi$rho) == 1)mods.df <- as.data.frame(mods$EMMLi$rho[[1]])
    if(length(mods$EMMLi$rho) > 1)mods.df <- mods$rho
    output <- list(model.comparison = model.comp,
                   EMMLi.analysis = mods,
                   corr.modules.df = mods.df,
                   hyps = hyps)
    if(show.comparison == T)print(model.comp)
    cat("","\n")
    cat("End of line","\n")
    cat("","\n")
    return(output)
    
  }
  if(mult == T){
    
    out.list <- list()

    for(i in 1:length(hypotheses)){
      hyps <- data.frame(name = landmark.names)
      
      nnam <- c()
      hyp.nam <- c()


      hyp.nam <- names(hypotheses)[i]
      nnam <- paste0("mod",i)
      hyps <- cbind(hyps,hypotheses[[i]])
 
        
      
      colnames(hyps) <- c("names", nnam)
      
      if(EMMLiv2 == T){
        
        mods <- phyloEmmli(shape, phylo, N_sample = N_sample, mod = hyps, EMMLi = T)
      }else{
        mods <- EMMLi(corr, N_sample = N_sample, mod = hyps)
      }
      
      model.comp <- as.data.frame(mods$EMMLi$results)
      model.comp <- model.comp[order(model.comp$dAICc),]
      
      winner <- str_split_1(names(mods$EMMLi$rho)[1], pattern = "mod")[2]
      winner <- as.numeric(str_split_1(winner, pattern = ".sep")[1])
      name.winner <- hyp.nam[winner]
      cat("","\n")
      cat("Hypothesis:",hyp.nam,"\n")
      cat("","\n")
      
      if(length(mods$EMMLi$rho) == 1)mods.df <- as.data.frame(mods$EMMLi$rho[[1]])
      if(length(mods$EMMLi$rho) > 1)mods.df <- as.data.frame(mods$EMMLi$rho[[1]])
      output <- list(model.comparison = model.comp,
                     EMMLi.analysis = mods,
                     corr.modules.df = mods.df,
                     hyps = hyps)

      out.list[[i]] <- output
      names(out.list)[i] <- hyp.nam
    }
    cat("","\n")
    cat("End of line","\n")
    cat("","\n")
    return(out.list)
  }
  
  
}

# move_values_below_diagonal() ####

move_values_below_diagonal <- function(mat) {
  # Check if the matrix is square
  if (nrow(mat) != ncol(mat)) {
    stop("Matrix must be square!")
  }
  
  n <- nrow(mat)
  
  # Create an empty matrix for the result
  result <- matrix(0, nrow = n, ncol = n)
  
  # Iterate through all elements above the diagonal
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      # If a value exists above the diagonal, move it below
      if (mat[i, j] != 0) {
        result[j, i] <- mat[i, j]
      } else if (mat[j, i] != 0) {
        result[j, i] <- mat[j, i]
      }
    }
  }
  
  # Copy the diagonal values as they are
  diag(result) <- diag(mat)
  
  return(result)
}


# corr.mod() ####

corr.mod <- function(data, pattern = " ", pos.row = 1, pos.col = 3){
  
  require(stringr)
  names.mod <- str_split(colnames(data), pattern = pattern)
  
  nb.modules <- 0
  for(q in 1:length(names.mod)){
    
    if(names.mod[[q]][1] == "Module")nb.modules <- nb.modules + 1
    
  }
  if(nb.modules <= 2)return("FAILED")
  nb.comb <- (nb.modules * (nb.modules - 1))/2

  
  
  data <- data[,(nb.modules + 1):(length(data[1,])-1)]

  str <- str_split(colnames(data), pattern = pattern)
  
  nb.comp <- length(str)
  
  rr <- c()
  cc <- c()
  
  for(i in 1:nb.comp){
    
    rr <- c(rr, str[[i]][pos.row])
    cc <- c(cc, str[[i]][pos.col])
  }
  
  mat <- matrix(nrow = length(unique(rr)), ncol = length(unique(cc)), dimnames = list(unique(rr), unique(cc)))

  k <- 1
  
  for(j in 1:length(mat[,1])){
    
    lgt.line <- length(mat[j,])
    nb.NA <- rep(NA,j-1)
    
    if(length(nb.NA) > 0){
      
      new <- nb.NA
      
    }else{
      new <- c()
    }
    
    diff.line <- lgt.line - length(new)
    new <- c(new, as.numeric(data[2,k:(k+(diff.line-1))]))
    k <- k + diff.line
    mat[j,] <- new
  }
  
  
  num <- as.numeric(factor(unique(c(colnames(mat), rownames(mat)))))
  fnum <- factor(unique(c(colnames(mat), rownames(mat))))
  names(num) <- fnum
  for(i in 1:length(colnames(mat))){
    
    for(j in 1:length(num)){
      
      if(names(num)[j] == colnames(mat)[i]){
        colnames(mat)[i] <- num[j]
      }
      
    }
  }
  
  for(i in 1:length(rownames(mat))){
    
    for(j in 1:length(num)){
      
      if(names(num)[j] == rownames(mat)[i]){
        rownames(mat)[i] <- num[j]
      }
      
    }
  }

  mat <- cbind(mat, as.numeric(rownames(mat)))
  colnames(mat)[length(mat[1,])] <- "ORDER"
  mat <- mat[order(mat[,length(mat[1,])]),]
  mat <- mat[,-length(mat[1,])]
  mat <- t(mat)
  mat <- cbind(mat, as.numeric(rownames(mat)))
  colnames(mat)[length(mat[1,])] <- "ORDER"
  mat <- mat[order(mat[,length(mat[1,])]),]
  mat <- mat[,-length(mat[1,])]
  
  
  mx.num <- max(num)
  suite <- 1:mx.num

  nums.col <- as.numeric(colnames(mat))
  for(k in suite){

    if(is.na(as.numeric(colnames(mat)[k])) == T){
      
      if(max(nums.col) == k)next
      mat <- cbind(mat[,1:(k-1)], rep(NA, length(mat[,1])))
      next
    }
    if(as.numeric(colnames(mat)[k]) != k){
      
      if(k == 1){
        mat <- cbind(rep(NA, length(mat[,1])), mat)
      }else{
        
        rest <- mat[,k:length(mat[1,])]
        mat <- cbind(mat[,1:(k-1)], rep(NA, length(mat[,1])), rest)
        
      }
        
      
    }
    
  }

  #print(colnames(mat))
  colnames(mat) <- suite

  for(k in suite){
    if(is.na(as.numeric(rownames(mat)[k])) == T){
      
      mat <- rbind(mat[1:(k-1),], rep(NA, length(mat[1,])))
      next
    }
    if(as.numeric(rownames(mat)[k]) != k){
      
      if(k == 1){
        mat <- rbind(rep(NA, length(mat[1,])), mat)
      }else{
        rest <- mat[k:length(mat[,1]),]
        mat <- rbind(mat[1:(k-1),], rep(NA, length(mat[1,])), rest)
        
      }
      
      
    }
    
  }
  
  #print(rownames(mat))
  rownames(mat) <- suite
  #print(mat)
  #print(rownames(mat))
  mat[is.na(mat)] <- 0
  mat <- as.matrix(mat)
  #print(mat)
  organized_matrix <- move_values_below_diagonal(mat)
  #print(organized_matrix)
  diag(organized_matrix) <- 1
  #print(organized_matrix)
  return(organized_matrix)
}

# EVRnRS.plot() ####

EVRnRS.plot <- function(color_palette = "viridis", tree, breaks = 6, node_cex, show.axis = T, edge_variable, edge.width = 3, cex.label = 0.6, label.offset = 3, export = F, export.params = list(file = paste0(getwd(),"plot.pdf"), width = 5, height = 5),  results){
  
  if(export == T){
    width.pl <- as.numeric(export.params$width[[1]])
    height.pl <- as.numeric(export.params$height[[1]])
    pdf(file = export.params$file[[1]], width = width.pl, height = height.pl, compress = F)
  }
  
  if(color_palette == "viridis")colors <- viridis(breaks)
  
  edge_colors <- colors[cut(edge_variable, breaks = breaks, labels = FALSE)]
  
  plot(tree, edge.color = edge_colors, edge.width = edge.width, show.tip.label = TRUE, type = "fan", cex = cex.label, label.offset = label.offset)
  nodelabels(pch = 24, cex = node_cex, col = "black", bg = "grey")
  
  if(show.axis == T)axisPhylo()
  
  if(export == T){
    dev.off()
  }
  
}

# CHEV.plot() ####

CHEV.plot <- function(simmap, groups, palette = "viridis", type = "fan",show.tip = T, cex = 0.2,fsize = 0.5, offset = 7, ftype = "i", export = F, export.params = list(file = paste0(getwd(), "/plot.pdf"), width = 5, height = 5)){
  
  require(viridis)
  
  if(export == T){
    width.pl <- as.numeric(export.params$width[[1]])
    height.pl <- as.numeric(export.params$height[[1]])
    pdf(file = export.params$file[[1]], width = width.pl, height = height.pl)
  }
  
  if(palette[1] == "viridis")palette <- viridis(length(levels(factor(groups))))
  
  
  plot(summary(simmap),colors = palette,ftype = ftype, type = type, cex = cex, fsize = fsize, offset = offset)
  if(export == T)dev.off()
}

# density.rate.plot() ####

density.rate.plot <- function(values, ylim1 = c(0,120),ylim2 = c(0,80), xlim = c(-22,-10), export = F, export.params = list(filename = "plot.pdf",width = 5, height = 5, path = getwd())){
  
  data <- data.frame(values = values, y = "NAME")
  data$y <- as.factor(data$y)
  
  g1 <- ggplot(data, aes(x = values)) +
    geom_histogram(aes(y = ..count.., fill = ..x..), bins = 6, color = "black") +
    coord_cartesian(ylim = ylim1, xlim = xlim)+
    scale_fill_gradientn(colors = viridis(6), 
                         guide = guide_colorbar(title = "X-axis Values")) +
    labs(x = "Values", 
         y = "count") +
    theme_bw() + theme(legend.position = "none")
  
  g2 <- ggplot(data, aes(x = values, y = y,fill=..x..)) + geom_density_ridges_gradient()+ scale_fill_gradientn(colors = viridis(6))+
    theme_bw()+theme(legend.position = "none")
  
  gar <- ggarrange(g1,g2)
  
  if(export == T){
    width.pl <- as.numeric(export.params$width[[1]])
    height.pl <- as.numeric(export.params$height[[1]])
    ggsave(export.params$filename[[1]],plot = gar, path = export.params$path[[1]], width = width.pl, height = height.pl, units = "cm")
  }else{
    print(gar)
  }
  
  
  
  return(list(plot.colors = g1, plot.density = g2))
}

# PP.node.plot() ####

PP.node.plot <- function(tree, PP.node, type = "fan",show.tip = T, bg = "grey", frame = "r", cex = 0.5, round.PP = 2, export = F, export.params = list(file = paste0(getwd(), "/plot.pdf"), width = 5, height = 5)){
  
  if(export == T){
    width.pl <- as.numeric(export.params$width[[1]])
    height.pl <- as.numeric(export.params$height[[1]])
    pdf(file = export.params$file[[1]], width = width.pl, height = height.pl)
  }
  
  plot(tree, type = type, show.tip.label = show.tip, label.offset = 3)
  nodelabels(round(c(0,PP.node), round.PP),bg = bg, frame = frame, cex = cex)
  if(export == T)dev.off()
}

# merge.results.CR.EMMLi() ####

merge.results.CR.EMMLi <- function(CR, EMMLi, names.hyp){
  
  library(tidyverse)
  
  tab.cr <- CR$df
  tab.em <- EMMLi$model.comparison
  tab.out <- cbind(tab.cr[NULL,], tab.em[NULL,])

  for(i in 1:length(names(names.hyp))){
    
    new.hyp <- names(names.hyp)[i]
    new.cr <- filter(tab.cr, HYP == new.hyp)
    new.em <- dplyr::filter(tab.em, Model.name2 == new.hyp)
    tab.na <- new.cr
    tab.na$CR <- NA
    tab.na$Z <- NA
    tab.na$P.Z <- NA
    tab.na$pls <- NA
    tab.na$P.pls <- NA
    tab.na$ES <- NA
    new.cr <- rbind(new.cr, tab.na,tab.na,tab.na)

    
    rownames(new.cr) <- NULL
    new.both <- cbind(new.cr, new.em)
    rownames(new.both) <- NULL
    tab.out <- rbind(tab.out, new.both)

  }
  rownames(tab.out) <- NULL
  return(tab.out)
}

# dtt.mult() ####

dtt.mult <- function(coords, tree, nsim = 1000, plot = T, calculateMDIp = T, warnings = F){
  
  if(class(coords) != "list")stop("ERROR - coords must be a list")
  
  out <- list()
  for(i in 1:length(coords)){
    
    cat(names(coords)[i],"\n")
    if(warnings == F){
      suppressWarnings(dtt.i <- dtt(phy=tree, data=two.d.array(coords[[i]]), nsim=nsim, plot=plot, calculateMDIp = calculateMDIp))
    }else{
      dtt.i <- dtt(phy=tree, data=two.d.array(coords[[i]]), nsim=nsim, plot=plot, calculateMDIp = calculateMDIp)
    }
    
    
    
    # plot
    
    sim_data <<- as.data.frame(dtt.i$sim)
    times <- dtt.i$times

    
    sim_df <- cbind(time = times, sim_data)  # Ensure time is the first column
    
    # Convert to long format
    sim_long <- sim_df %>%
      pivot_longer(-time, names_to = "simulation", values_to = "disparity")
    sim_long <- sim_long
    # Compute summary statistics
    summary_df <- sim_long %>%
      group_by(time) %>%
      dplyr::summarise(mean_disp = median(disparity),lower_disp = quantile(disparity, 0.025),  # 2.5% quantile
        upper_disp = quantile(disparity, 0.975))
    
    # Prepare data for geom_polygon (create lower and upper paths)
    polygon_df <- bind_rows(
      summary_df %>% select(time, disparity = lower_disp),
      summary_df %>% select(time, disparity = upper_disp) %>% arrange(desc(time))
    )

    polygon_df$time <- 200 - (polygon_df$time*200)
    dat.obs <- data.frame(time = (200-(dtt.i$times*200)), DTT = dtt.i$dtt)
    summary_df$time <- 200-(summary_df$time*200)


    GR <- ggplot() +
      geom_polygon(data = polygon_df, aes(x = time, y = disparity, group = 1), fill = "grey70", alpha = 0.5) +
      geom_line(data = summary_df, aes(x = time, y = mean_disp), color = "blue",linetype = 2, linewidth = 1) +
      geom_line(data = dat.obs, aes(x = time, y = DTT), linewidth = 1)+
      scale_x_reverse()+
      labs(x = "Time", y = "DTT Value", title = names(coords)[i]) +
      theme_bw() + 
      theme(plot.title = element_text(hjust = 0.5))

    out[[i]] <- list(DTT = dtt.i, plot = GR)
    names(out)[i] <- names(coords)[i]
    
  }
  return(out)
}

# comp.evol.r.multiple() ####
comp.evol.r.multiple <- function(list.A, gp, Subset = T, phy){
  
  out <- list()
  
  for(i in 1:length(list.A)){
    
    print(names(list.A)[i])
    
    phy.pruned <- keep.tip(phy, dimnames(list.A[[i]])[[3]])
    
    EMR <- compare.multi.evol.rates(A = list.A[[i]], gp = gp,Subset = T, phy = phy.pruned)
    out[[i]] <- EMR
    names(out)[i] <- names(list.A)[i]
  }
  return(out)
}

# morphol.disp.multiple() ####

morphol.disp.multiple <- function(gdf, phy, shape.name, variable.name = NULL, module,size.name = NULL, size = F){
  
  module <- factor(module)
  shape <- gdf[[shape.name]]
  if(is.null(variable.name) == F)variable <- gdf[[variable.name]]
  if(size == T)size.v <- gdf[[size.name]]
  out <- list()
  
  for(i in 1:length(levels(module))){
    
    print(levels(module)[i])
    
    shp.i <- shape[which(module == levels(module)[i]),,]
    
    phy.i <- keep.tip(phy,dimnames(shape)[[3]])
    
    if(size == T & is.null(variable.name) == F)gdf.i <- geomorph.data.frame(shape = shp.i, variable = variable ,phy = phy.i, size.v = size.v)
    if(size == T & is.null(variable.name) == T)gdf.i <- geomorph.data.frame(shape = shp.i,phy = phy.i, size.v = size.v)
    if(size == F)gdf.i <- geomorph.data.frame(shape = shp.i, variable = variable ,phy = phy.i)
    
    if(size == T & is.null(variable.name) == F)PGLS.res <- procD.pgls(shape ~ variable * size.v, data = gdf.i, phy = phy)
    if(size == T & is.null(variable.name) == T)PGLS.res <- procD.pgls(shape ~ size.v, data = gdf.i, phy = phy)
    if(size == F)PGLS.res <- procD.pgls(shape ~ variable, data = gdf.i, phy = phy)
    
    if(is.null(variable.name) == F)md.res <- morphol.disparity(PGLS.res, groups = ~ variable)
    if(is.null(variable.name) == T)md.res <- morphol.disparity(PGLS.res)
    
    out[[i]] <- md.res
    names(out)[i] <- levels(module)[i]
  }
  return(out)
}

# physignal.mult() ####
physignal.mult <- function(modules, coords, tree, lambda = "burn", PAC.no = NULL, threshold = 0.95){
  
  modules <- factor(modules)
  out <- list()
  out.df <- data.frame(module = NULL, K = NULL, lambda = NULL, Z = NULL, p = NULL)
  
  relevance_mod <- data.frame(mod = levels(modules))
  relevance_mod$rel <- NA
  
  for(i in 1:length(levels(modules))){
    
    mod.i <- levels(modules)[i]

    
    coords.i <- coords[which(modules == mod.i),,]
    if(is.null(PAC.no) == T)pz <- physignal.z(A = coords.i, phy = tree, lambda = lambda)
    if(is.null(PAC.no) == F){
      
      if(PAC.no == "find"){
        pzt <- physignal.z(A = coords.i, phy = tree, lambda = lambda)
        
        summ <- sum(pzt$PACA$d)
        eigen.perc <- pzt$PACA$d/summ
        nb.dims <- NA
        max.perc <- 0
        for(j in 1:length(eigen.perc)){
          max.perc <- max.perc + eigen.perc[j]
          if(max.perc > threshold){
            nb.dims <- j
            break
          }
        }
        
        pz <- physignal.z(A = coords.i, phy = tree, lambda = lambda, PAC.no = nb.dims)
        
      }else if(class(PAC.no) == "numeric" | class(PAC.no) == "integer"){
        pz <- physignal.z(A = coords.i, phy = tree, lambda = lambda, PAC.no = PAC.no)
      }else{
        cat("ERROR - PAC.no must be numeric or integer or = find","\n")
      }
      
      
    }
      

    if(is.na(pz$Z) == T){
      out.df <- rbind(out.df, data.frame(module = mod.i,K = 0, lambda = 0, Z = pz$Z, p = NA))
    }else{
      out.df <- rbind(out.df, data.frame(module = mod.i,K = pz$K, lambda = pz$lambda, Z = pz$Z, p = pz$pvalue))
      relevance_mod$rel[i] <- "YES"
    }
      
    
    
    out[[i]] <- pz
    names(out)[i] <- mod.i
    
  }

  out.rel <- out[which(relevance_mod$rel == "YES")]
  CPZ <- compare.physignal.z(out.rel)
  
  out[[length(out)+1]] <- out.df
  names(out)[length(out)] <- "Results.data.frame"
  out[[length(out)+1]] <- CPZ
  names(out)[length(out)] <- "Comparison"
  
  return(out)
}

comp.group.physignal <- function(shape, groups, tree, lambda = "burn", PAC.no = NULL, threshold = 0.95, pairwise = TRUE){
  
  out <- list()
  out.df <- data.frame(group = NULL, K = NULL, lambda = NULL, Z = NULL, p = NULL)
  
  relevance_gr <- data.frame(gr = unique(groups))
  relevance_gr$rel <- NA
  
  for(i in 1:length(unique(groups))){
    
    gr.i <- unique(groups)[i]
    tree.i <- keep.tip(tree, tip = dimnames(shape)[[3]][which(groups == unique(groups)[i])])
    shape.i <- shape[,,which(groups == unique(groups)[i])]
    
    if(is.null(PAC.no) == T)pz <- physignal.z(A = shape.i, phy = tree.i, lambda = lambda)
    if(is.null(PAC.no) == F){
      
      if(PAC.no == "find"){
        pzt <- physignal.z(A = shape.i, phy = tree.i, lambda = lambda)
        
        summ <- sum(pzt$PACA$d)
        eigen.perc <- pzt$PACA$d/summ
        nb.dims <- NA
        max.perc <- 0
        for(j in 1:length(eigen.perc)){
          max.perc <- max.perc + eigen.perc[j]
          if(max.perc > threshold){
            nb.dims <- j
            break
          }
        }
        
        pz <- physignal.z(A = shape.i, phy = tree.i, lambda = lambda, PAC.no = nb.dims)
        
      }else if(class(PAC.no) == "numeric" | class(PAC.no) == "integer"){
        pz <- physignal.z(A = shape.i, phy = tree.i, lambda = lambda, PAC.no = PAC.no)
      }else{
        cat("ERROR - PAC.no must be numeric or integer or = find","\n")
      }
      
      if(is.na(pz$Z) == T | length(tree.i$tip.label) <= 3){
        out.df <- rbind(out.df, data.frame(group = gr.i,K = 0, lambda = 0, Z = pz$Z, p = NA))
      }else{
        out.df <- rbind(out.df, data.frame(group = gr.i,K = pz$K, lambda = pz$lambda, Z = pz$Z, p = pz$pvalue))
        relevance_gr$rel[i] <- "YES"
      }
      
      
      
      out[[i]] <- pz
      names(out)[i] <- gr.i
      
    }
  }
  
  if(pairwise == T){
    out.rel <- out[which(relevance_gr$rel == "YES")]
    CPZ <- compare.physignal.z(out.rel)
  }
  
  out[[length(out)+1]] <- out.df
  names(out)[length(out)] <- "Results.data.frame"
  if(pairwise == T){
    out[[length(out)+1]] <- CPZ
    names(out)[length(out)] <- "Comparison"
    
  }
    
  
  return(out)
}
