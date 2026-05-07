# EMMLiv2 package

###########################################################################################
#' checkInput
#'
#' Checks the input of EMMLi to make sure it's all what is expected. Throws errors with
#' a descriptive error when it is not.
#' @name checkInput
#' @param corr Correlation matrix
#' @param mod Model file (landmarks classified into modules)
#' @param abs Logical - whether to use abosulte values of correlations
#' @param pprob posterior probability cut off
#' @param saveAs filename to save output too
#' @keywords internal

checkInput <- function(corr, mod, N_sample, abs, pprob, saveAs) {
  # Check inputs
  if(!is.numeric(corr) & !is.data.frame(corr)){
    stop('corr should be a (square) matrix or data frame.')
  }
  
  if(!(is.factor(mod[, 1]) | is.character(mod[, 1]))){
    stop('The first column of mod should be landmark names (factor or character).')
  }
  
  # Make elements in mod, after the first column, integers or NAs
  modClasses <- sapply(mod[, -1], class)
  modClasses[modClasses == 'integer'] <- 'numeric'
  if(!all(modClasses == 'factor') & !all(modClasses == 'numeric') & !all(modClasses == 'character')){
    stop('mod should contain landmark names in the first column and integers, factors or character vectors in subsequent columns')
  }
  
  # check that numerics are integers
  if(all(modClasses == 'numeric')){
    if(!all(sapply(mod[, -1], function(x) all(abs(stats::na.omit(x) - round(stats::na.omit(x))) < .Machine$double.eps^0.5)))){
      stop('mod should contain landmark names in the first column and integers, factors or character vectors in subsequent columns')
    }
  }
  
  if(all(modClasses == 'factor')){
    mod[, -1] <- sapply(mod[, -1], function(x) as.numeric(x))
  }
  
  if(all(modClasses == 'character')){
    mod[, -1] <- sapply(mod[, -1], function(x) as.numeric(factor(x)))
  }
  
  # Check corr
  if(!dim(corr)[1] == dim(corr)[2]) stop('corr should be a square matrix')
  
  # Check other parameters
  stopifnot(is.numeric(N_sample), N_sample > 0, is.logical(abs), pprob > 0, pprob < 1)
  
  if(!is.null(saveAs)){
    stopifnot(is.character(saveAs))
    if(!grepl('\\.csv$', saveAs)){
      warning("Output will be saved as a csv but saveAs does not end in '.csv'")
    }
  }
}

###########################################################################################
#' getModules
#'
#' Breaks correlation matrix up into discrete modules (within and between) based on
#' a model table (landmarks classified into modules)
#' @name getModules
#' @param varlist internal variable from EMMLi
#' @param symmet internal variable from EMMLi
#' @param corr correlation matrix
#' @param mod Landmark classifications
#' @param corr_list internal variable from EMMLi
#' @keywords internal

getModules <- function(corr, mod) {
  # Create varlist variable
  varlist = paste0('mod$', names(mod)[-1])
  
  # make the upper triangle of corr NA, we only use the lower triangle.
  corr[upper.tri(corr, diag = T)] = NA
  
  # array of coefficient matrix, NAs are removed.
  corr_list = (as.array(corr[!is.na(corr)]))
  symmet = corr
  
  # a symmetric matrix formed from the coefficient matrix, only used to find intermodular coefficients.
  symmet[upper.tri(symmet)] = t(symmet)[upper.tri(symmet)]
  
  all_modules = list()
  for(colnam in varlist){
    col = array(eval(parse(text = colnam)))
    
    # na.omit() used to remove unclassified landmarks.
    modNF = stats::na.omit(cbind(1:nrow(mod), col))
    w = unique(modNF[, 2])
    
    modules = list()
    for(i in seq(length(w))){
      # identify landmarks within a class
      fg = modNF[modNF[, 2] == w[i], ]
      
      # coefficients between identified landmarks.
      l <- corr[as.numeric(fg[, 1]), as.numeric(fg[, 1])]
      modules[[i]] = (as.array(l[!is.na(l)]))
    }
    names(modules) = paste("Module", w)
    
    between_mod = list()
    betweenModules = list()
    withinModules = list()
    unintegrated = list()
    betweenFloat = list()
    
    if (length(w) > 1){ #check that the num. of modules is greater than 1
      
      # all possible combination of modules
      cb = utils::combn(w, 2)
      for (i in seq(dim(cb)[2])){
        fg1 = modNF[modNF[, 2] == cb[1, i], ]
        fg2 = modNF[modNF[, 2] == cb[2, i], ]
        
        # setdiff(A,B) - present in A but not in B
        between = symmet[as.integer(setdiff(fg2[, 1], fg1[, 1])), as.integer(setdiff(fg1[, 1], fg2[, 1]))]
        between_mod[[i]] = between[!is.na(between)]
      }
      names(between_mod) = paste(cb[1, ], "to", cb[2, ])
      
      betweenModules['betweenModules'] = list(as.vector(rle(unlist(between_mod))$values))
      withinModules['withinModules'] = list(as.vector(rle(unlist(modules))$values))
    }
    
    unintegrated_list = setdiff(corr_list,unlist(c(modules,between_mod)))
    unintegrated['unintegrated'] = list(as.vector(rle(unlist(unintegrated_list))$values))
    # unintegrated['unintegrated+between'] = list(as.vector(rle(unlist(c(unintegrated_list,betweenModules)))$values))
    
    if (length(unintegrated_list) != 0){
      betweenFloat['betweenFloat'] = list(as.vector(rle(unlist(c(betweenModules['betweenModules'], unintegrated['unintegrated'])))$values))
      all_modules[[colnam]] = c(modules, between_mod, withinModules, betweenModules, betweenFloat, unintegrated)
    } else {
      all_modules[[colnam]] = c(modules, between_mod, withinModules, betweenModules, unintegrated)
    }
    
  }
  return(all_modules)
}

###########################################################################################
#' moduleLikelihoods
#'
#' Calculates the maximum likelihood value of rho for each module (between/within) from
#' getModules
#' @name moduleLikelihoods
#' @param all_modules output from getModules - internal to EMMLi
#' @param abs use absolute values or not
#' @param N_sample sample size
#' @param correction Experimental - whether to calculate number of parameters normally,
#' or number of parameters + number of modules - 1
#' @keywords internal

moduleLikelihoods <- function(all_modules, abs, N_sample, correction = "normal") {
  
  LogL = function(z_r, z_p) {
    -0.5 * log(var) - ((z_r - z_p)^2) / (2 * var)
  }
  
  output = list()
  maxlogL = list()
  logp = list()
  for (m in seq(length(all_modules))){
    # table of max. likelihood and rho
    maxres = matrix(, nrow = 2, ncol = (length(all_modules[[m]])))
    dimnames(maxres) = list(c("MaxL", "MaxL_p"), c(names(all_modules[[m]][seq((length(all_modules[[m]])))])))
    
    ###################################
    # Coulds for loop be a seperate function? Seems a good target.
    
    for(g in seq((length(all_modules[[m]])))){
      r = unlist(unname(all_modules[[m]][g]))
      n_value = length(r)
      
      # Calculate z_r differently depending on abs argument in original call.
      if(abs){
        z_r = 0.5 * log((1 + abs(r)) / (1 - abs(r)))
        #rho
        p = seq(0, 0.99, 0.01)
      } else if(!abs) {
        z_r = 0.5 * log((1 + r) / (1 - r))
        #rho
        p = seq(-0.99, 0.99, 0.01)
      }
      n = N_sample
      var = 1 / (n - 3)
      
      #zeta
      z_p = 0.5 * log((1 + p)/(1 - p))
      
      LogL_table = outer(z_r, z_p, LogL)
      Likelihoods = colSums(LogL_table)
      MaxIndex = which.max(Likelihoods)
      
      MaxL = Likelihoods[MaxIndex]
      MaxL_p = p[MaxIndex]
      maxres[1, g] = MaxL
      maxres[2, g] = MaxL_p
    }
    
    # list of maximum likelihood for all modules.
    output[[names(all_modules)[m]]] = maxres
    
    #calculate sum of modular likelihood and k
    if (dim(output[[m]])[2] == 2){
      
      maxlogL[[names(all_modules)[m]]][['default']] = c(output[[m]][1], 2)
      logp[[names(all_modules)[m]]][['default']] = output[[m]][, 1, drop = FALSE]
      
    } else if (dim(output[[m]])[2] == 6){
      
      #if(output[[m]][1,'unintegrated'] == 0){output[[m]] = output[[m]][,colnames(output[[m]]) != 'unintegrated']}
      
      mod_between = output[[m]]['MaxL', grep('Module |betweenModules|unintegrated', names(output[[m]][1, ]))]
      mod_between_p = output[[m]][,grep('Module |betweenModules|unintegrated', names(output[[m]][1, ]))]
      if (mod_between['unintegrated'] == 0){
        K = length(mod_between) - 1
      } else {
        K = length(mod_between)
      }
      maxlogL[[names(all_modules)[m]]][['sep.Mod + same.between']] = c(sum(mod_between), K + 1)
      logp[[names(all_modules)[m]]][['sep.Mod + same.between']] = mod_between_p
      
      
      within_between = output[[m]]['MaxL', c('withinModules','betweenModules','unintegrated')]
      within_between_p = output[[m]][, c('withinModules','betweenModules','unintegrated')]
      if (within_between['unintegrated'] == 0){
        K = length(within_between) - 1
      } else {
        K = length(within_between)
      }
      maxlogL[[names(all_modules)[m]]][['same.Mod + same.between']] = c(sum(within_between), K+1)
      logp[[names(all_modules)[m]]][['same.Mod + same.between']] = within_between_p
    } else {
      
      #if(output[[m]][1,'unintegrated'] == 0){output[[m]] = output[[m]][,colnames(output[[m]]) != 'unintegrated']}
      
      mod_to = output[[m]]['MaxL', grep('Module |to |unintegrated', names(output[[m]][1, ]))]
      mod_to_p = output[[m]][, grep('Module |to |unintegrated', names(output[[m]][1,]))]
      if (mod_to['unintegrated'] == 0){
        K = length(mod_to)-1
      } else {
        K = length(mod_to)
      }
      maxlogL[[names(all_modules)[m]]][['sep.Mod + sep.between']] = c(sum(mod_to), K + 1)
      logp[[names(all_modules)[m]]][['sep.Mod + sep.between']] = mod_to_p
      
      within_between = output[[m]]['MaxL', c('withinModules', 'betweenModules', 'unintegrated')]
      within_between_p = output[[m]][, c('withinModules', 'betweenModules', 'unintegrated')]
      if (within_between['unintegrated'] == 0){
        K = length(within_between) - 1
      } else {
        K = length(within_between)
      }
      maxlogL[[names(all_modules)[m]]][['same.Mod + same.between']] = c(sum(within_between), K + 1)
      logp[[names(all_modules)[m]]][['same.Mod + same.between']] = within_between_p
      
      mod_between = output[[m]]['MaxL', grep('Module |betweenModules|unintegrated', names(output[[m]][1, ]))]
      mod_between_p = output[[m]][, grep('Module |betweenModules|unintegrated', names(output[[m]][1, ]))]
      if (mod_between['unintegrated'] == 0){
        K = length(mod_between) - 1
      } else {
        K = length(mod_between)
      }
      maxlogL[[names(all_modules)[m]]][['sep.Mod + same.between']] = c(sum(mod_between), K + 1)
      logp[[names(all_modules)[m]]][['sep.Mod + same.between']] = mod_between_p
      
      to_within = output[[m]]['MaxL', grep('to |withinModules|unintegrated', names(output[[m]][1, ]))]
      to_within_p = output[[m]][, grep('to |withinModules|unintegrated', names(output[[m]][1, ]))]
      if (to_within['unintegrated'] == 0){
        K = length(to_within) - 1
      } else {
        K = length(to_within)
      }
      maxlogL[[names(all_modules)[m]]][['same.Mod + sep.between']] = c(sum(to_within), K + 1)
      logp[[names(all_modules)[m]]][['same.Mod + sep.between']] = to_within_p
      
      if (output[[m]]['MaxL', grep('unintegrated', names(output[[m]][1, ]))] != 0){
        
        sepmod_samebetweenunintegrated = output[[m]][, grep('Module |betweenFloat',names(output[[m]][1, ]))]
        K = length(sepmod_samebetweenunintegrated['MaxL', ])
        maxlogL[[names(all_modules)[m]]][['sep.mod + same.between.unintegrated']] = c(sum(sepmod_samebetweenunintegrated['MaxL', ]), K + 1)
        logp[[names(all_modules)[m]]][['sep.mod + same.between.unintegrated']] = sepmod_samebetweenunintegrated
        
        samemod_samebetweenunintegrated = output[[m]][, grep('withinModules|betweenFloat', names(output[[m]][1, ]))]
        K = length(samemod_samebetweenunintegrated['MaxL', ])
        maxlogL[[names(all_modules)[m]]][['same.mod + same.between.unintegrated']] = c(sum(samemod_samebetweenunintegrated['MaxL', ]), K + 1)
        logp[[names(all_modules)[m]]][['same.mod + same.between.unintegrated']] = samemod_samebetweenunintegrated
        
      }
    }
  }
  
  # Experimental correction to parameter number. Adds the number of modules - 1 to K, the
  # parameter count. Extra penalisation, to see if EMMLi then has less of a tendency to
  # overparamaterise.
  if (correction == "new") {
    for (i in seq_along(maxlogL)) {
      n_modules <- length(grep("^Module", names(all_modules[[i]])))
      for (g in seq_along(maxlogL[[i]])) {
        maxlogL[[i]][[g]][2] <- maxlogL[[i]][[g]][2] + (n_modules - 1)
      }
    }
  }
  
  rets <- list(maxlogL = maxlogL, logp = logp)
  return(rets)
}

##### Title          #####
#' Evaluating modularity with maximum likelihood
#'
####  Description    #####
#' Calculates the AICc values, model likelihoods, and posterior probabilities of different models of modularity, as described in Goswami and Finarelli (2016).
#'
####  Details        #####
#'  The publication describing this analysis is A. Goswami and J. Finarelli
#'    (2016) EMMLi: A maximum likelihood approach to the analysis of modularity.
#'    Evolution \url{http://onlinelibrary.wiley.com/doi/10.1111/evo.12956/abstract}.
#'
#'@param corr Lower triangle or full correlation matrix. n x n square matrix for n landmarks.
#'@param N_sample The number of specimens
#'@param mod A data frame defining the models. The first column should contain the landmark names.
#'Subsequent columns should define which landmarks are contained within each module with integers,
#'factors or characters. If a landmark should be ignored for a specific model (i.e., it is
#'unintegrated in any module), the element should be NA.
#'@param saveAs A character string defining the filename and path for where to save output. If NULL,
#'the output is not saved to file
#'@param abs Logical denoting whether absolute values should be used. Default is TRUE, as in Goswami
#'and Finarelli (2016)
#'@param pprob posterior probability cutoff for reporting of models. Default is 0.05, as suggested in
#'Goswami and Finarelli (2016)
#'@param correction If "normal" then AIC is calculated normally, if "new" then the number of modules - 1
#' is added to the n parameter penalisation during AIC calculation. This is experimental and is not
#' recommended!
#'
#'@export
#'@return A list containing two elements. The first (results) gives the AIC results for each model.
#'  The second (rho) gives the within and between module correlations.
#'  Optionally, the output is saved to the file defined by the saveAs argument with only models with a
#'  posterior probability > 0.01 being saved.
#'
#'
#'@examples
#'  set.seed(1)
#'
#'  # Chose a filename and directory for output
#'  dir <- tempdir()
#'  file <- paste0(dir, 'EMMLiTest.csv')
#'
#'  # Examine a correlation matrix and model dataframe
#'  dim(macacaCorrel)
#'  head(macacaModels)
#'
#'  # run EMMLi
#'  output <- EMMLi(macacaCorrel, 20, macacaModels, file)
#'
#'  unlink(file)
#'
#'  # run EMMLi without writing output
#'  output <- EMMLi(macacaCorrel, 20, macacaModels)
#'
#'  # Raw data example to illustrate pitfalls
#'  corrPath <- system.file("extdata", "M1lmcorrel.csv", package = "EMMLi")
#'  corr <- read.csv(corrPath, header = FALSE)
#'
#'  modelPath <- system.file("extdata", "macaca_landmarklist.csv", package = "EMMLi")
#'  mod <- read.csv(modelPath, header = TRUE, row.names = 1)
#'
#'  # First column should be character or factor. Subsequent columns integer
#'  sapply(mod, class)
#'
#'  out <- EMMLi(corr, 42, mod)
#'

EMMLi <- function(corr, N_sample, mod, saveAs = NULL, abs = TRUE, pprob = 0.05,
                  all_rhos = FALSE, correction = "normal"){
  
  checkInput(corr, mod, N_sample, abs, pprob, saveAs)
  
  # Create null model
  mod$No.modules = 1
  
  # get correlation matrices for each modules.
  all_modules <- getModules(corr, mod)
  
  # maxlogL will have the log likelihood for each module.
  liks <- moduleLikelihoods(all_modules, abs, N_sample, correction = correction)
  maxlogL <- liks$maxlogL
  
  logp <- liks$logp
  
  results = matrix(unlist(maxlogL), ncol = 2, byrow = TRUE)
  a = names(unlist(maxlogL))[seq(1, dim(results)[1] * dim(results)[2], dim(results)[2])]
  dimnames(results) = list(a, c('MaxL', 'K'))
  
  n = length(which(!is.na(corr) == TRUE))
  AICc = apply(results, 1, function(x) -2 * x['MaxL'] + 2 * x['K'] + (2 * x['K'] * (x['K'] + 1)) / (n - x['K'] - 1))
  dAICc = AICc - min(AICc)
  Model_L = exp(-0.5 * dAICc)
  Post_Pob = Model_L / sum(Model_L)
  
  results = cbind(results, n, AICc, dAICc, Model_L, Post_Pob)
  
  n = names(all_modules)
  nm = unlist(strsplit(n, split = 'mod\\$'))[seq(2, 2 * length(n), 2)]
  
  a = names(unlist(maxlogL))[seq(1, dim(results)[1] * 2, 2)]
  b = unlist(strsplit(a, split = 'mod\\$'))
  b = unlist(strsplit(b, split = '1$'))
  rownames(results) = b
  
  s = 1
  i = length(nm)
  o = order(results[grep(nm[i], rownames(results)), 2])
  t = results[grep(nm[i], rownames(results)), , drop = FALSE][o, , drop = FALSE]
  s = s + length(o)
  
  for(i in 1:(length(nm) - 1)){
    o = order(results[grep(paste(nm[i], '.s', sep = ""), rownames(results)), 2])
    t = rbind(t, results[grep(paste(nm[i], '.s', sep = ""), rownames(results)), , drop = FALSE][o, , drop = FALSE])
    s = s + length(o)
  }
  
  results = t
  
  rholist = list()
  h = 1
  
  for (i in 1:(length(logp))){
    for (j in 1:length(logp[[i]])){
      
      rholist[h] = logp[[i]][j]
      
      h = h + 1
    }
  }
  
  if (all_rhos) {
    names(rholist) <- names(Post_Pob)
  }
  
  # build output for the csv.
  rho_output = rholist[which(Post_Pob > pprob)]
  
  rholist_name = names(which(Post_Pob > pprob))
  rholist_name = unlist(strsplit(rholist_name, split = 'mod\\$'))
  rholist_name = unlist(strsplit(rholist_name, split = '1$'))
  
  return_rho <- rho_output
  names(return_rho) <- rholist_name
  
  if(!is.null(saveAs)){
    utils::write.table(results, file = saveAs, row.names = TRUE, col.names = NA, sep = ",")
    cat("\n\n", file = saveAs, append = TRUE)
    for(q in 1:length(rho_output)){
      cat(rholist_name[q], "\n", file = saveAs, append = TRUE)
      write(paste(c('', colnames(rho_output[[q]])), collapse = ','), saveAs, append = TRUE)
      utils::write.table(rho_output[q], saveAs, row.names = TRUE, col.names = FALSE, sep = ",", append = TRUE)
      cat("\n", file = saveAs, append = TRUE)
    }
  }
  
  if (all_rhos) {
    res <- list(results = results, best_rho = return_rho, all_rhos = rholist)
  } else {
    res <- list(results = results, rho = return_rho)
  }
  
  return(res)
  
}

################################################################################
#' IRSAL
#'
#' @param atlas object of class "atlas" created by \link[Morpho]{createAtlas}
#' @param landmarks k x 3 x n array containing reference landmarks of the sample
#' or a matrix in case of only one target specimen.
#' @param initial_fixed The fixed points for the initial
#' \link[Morpho]{placePatch} analysis.
#' @param n The number of patched points to use in each iteration of IRSAL.
#' @param reps The number of times to repeat the IRSAL procedure.
#' @param write.out If TRUE then the patched points from each iteration, for
#' each species, will be written to the working directory. This will NOT include
#' the temporary anchor points.
#' @param sp If write.out is true, then this determines which species is
#' written out. Defaults to 1 (the first species in the dataset). 2 would be
#' the second species in the dataset etc.
#' @param maxit The maximum number of attempts to find improved bending energy.
#' @param ... Additonal arguments passed to \link[Morpho]{placePatch}
#' @importFrom magrittr %>%
#' @name IRSAL
#' @export

IRSAL <- function(atlas, landmarks, initial_fixed, n, reps, write.out = FALSE,
                  sp = NULL, maxit = 50, ...) {
  
  # Place patches to get a starting point for each specimen.
  print("Initial patching...")
  start <- original <- Morpho::placePatch(atlas = atlas,
                                          dat.array = landmarks,
                                          keep.fix = initial_fixed,
                                          silent = TRUE,
                                          ...)
  
  # Make a series of percentage increases to go through.
  pcs <- round(seq.int(from = 0.10, to = 0.90, length.out = reps), 2)
  
  res <- array(0, dim = dim(start))
  dimnames(res) <- dimnames(start)
  
  # Calculate the bending energy for the template to reference against.
  template <- rbind(atlas[["landmarks"]], atlas[["patch"]])
  L_int <- CreateL(template)
  
  for (j in seq_len(dim(landmarks)[3])) {
    prefix <- dimnames(landmarks)[[3]][j]
    print(paste0("Specimen ", j, ":", prefix))
    t <- start[,,j]
    
    # calculate initial bending energy
    be_p <- t(t) %*% L_int$Lsubk %*% t %>%
      as.matrix %>%
      diag %>%
      sum
    failure <- 0
    i <- 1
    repeat {
      
      patch_idx <- (dim(landmarks)[1] + 1):dim(original)[1]
      pps <- t[patch_idx, ]
      
      # Sample some of the patches.
      samples <- sort(sample(1:dim(pps)[1], round(dim(pps)[1] * pcs[i], 0)))
      
      # Take the relevant parts out of the original atlas...
      a_patch <- atlas[["patch"]]
      a_lms <- atlas[["landmarks"]]
      a_fixed <- atlas[["keep.fix"]]
      
      # Set up the new info. Pathces stays the same.
      new_patch <- a_patch
      # Landmarks have the sampled patches added to them.
      # Landmarks should now be their original length + n.
      new_lms <- rbind(a_lms, a_patch[samples, ])
      # And the fixed points are defined
      new_fixed <- c(a_fixed, (nrow(a_lms) + 1):(nrow(a_lms) + length(samples)))
      
      # Now set up the new atlas. It's different because it has some new
      # fixed points (landmarks is longer - patches moved there).
      new_atlas <- atlas
      new_atlas[["patch"]] <- new_patch
      new_atlas[["landmarks"]] <- new_lms
      new_atlas[["keep.fix"]] <- new_fixed
      
      # Now move the patches into the landmarks data to make new landmark
      # data to correspond to the definition in the new atlas.
      new_landmarks <- abind::abind(landmarks[,,j], pps[samples, ], along = 1)
      # Now place the new patch
      t_new <- Morpho::placePatch(atlas = new_atlas,
                                  dat.array = new_landmarks,
                                  keep.fix = new_fixed,
                                  prefix = prefix,
                                  silent = TRUE,
                                  ...)
      
      # Compare bending energy to previous iteration.
      # Remove the "IRSALed" landmarks...
      f_lms <- (nrow(landmarks[,,j]) + 1):(nrow(landmarks[,,j]) + length(samples))
      t_test <- t_new[-f_lms, ]
      
      # Then compare the bending energy. Then compare this to the bending
      # energy of the last iteration (be_p)
      be_n <- t(t_test) %*% L_int$Lsubk %*% t_test %>%
        as.matrix %>%
        diag %>%
        sum
      
      # Now, if be_n (bending energy new) is smaller than be_p (bending
      # energy previous), t becomes t_test (the new patches, with the false
      # anchor points removed), and be_p becomes be_n.
      if (be_n < be_p) {
        print(paste0("Bending energy improved at rep ", i))
        t <- t_test
        be_p <- be_n
        # Increment i
        i <- i + 1
        failure <- 0
        # Then check if it equals reps + 1 (which means it has fully gone
        # through all the sample percentages).
        if (i == (reps + 1)) {
          break
        }
        
      } else {
        failure <- failure + 1
        
        if (failure == maxit) {
          print(paste0("No resolution at rep ", i, " after ", maxit,
                       " attempts."))
          break
        }
        
      }
    }
    # Remove anchor points (patched points that are now landmarks) and store
    # result. At this point t will be the patch that has the lowest bending
    # energy during IRSAL - and it will also have had the landmarks removed
    # alread (it will either be the initial patch that hasn't ever changed,
    # if bending energy was never reduced, OR it will be a t_test (that has
    # had the IRSAL landmarks removed in order to compare bending energy) that
    # had a lower bending energy. Pop that patch into the overall results...
    print("Storing best patch.")
    res[, , j] <- t
  }
  return(res)
}

###########################################################################################
#' compareModules
#'
#' This takes a correlation matrix or 3D landmark array, a model definition, and then two module
#' numbers or names to compare. It plots (if plot = TRUE) a figure of three boxplots - the first
#' two are the correltions within each of the two modules, and the third is the between-module
#' correlations. These boxes are coloured such that matching colours are not significantly
#' different according to a Tukey HSD test. The results of the anova and tukey HSD test are
#' also returned.
#' @name compareModules
#' @param corr A correlation matrix or a 3D array of landmarks. If 3D then a correlation matrix
#' is calculated with \link[paleomorph]{dotcorr}
#' @param model Either a vector of numbers describing a model of modules, or a 2 column dataframe
#' with the first bein landmark names and the second being the module definitions.
#' @param test_modules A vector of two module numbers to compare, or if the modules are named, the names
#' of those two modules.
#' @param plot Logical - if TRUE the plot is drawn.
#' @return A list with two elements - the first is the result of an ANOVA compaing the mean
#' correlations within- and between-modules, and the second is the results of a TukeyHSD test on
#' that ANOVA. If plot = TRUE a plot showing these results is called.
#' @export
#' @examples
#' data(macacaCorrel)
#' data(macacaModels)
#' # Pick a model to draw modules from - as a vector.
#' model <- macacaModels$Goswami
#' compareModules(corr = macacaCorrel, model = model, test_modules = c(2, 5))
#'
#' # Or as a 2 column dataframe...
#' model <- macacaModels[ , c(1, 4)]
#' compareModules(corr = macacaCorrel, model = model, test_modules = c(2, 5))

compareModules <- function(corr, model, test_modules, plot = TRUE) {
  
  # If landmarks supplied
  if (length(dim(corr)) == 3) {
    require(paleomorph)
    corr <- paleomorph::dotcorr(corr)
  }
  
  # If model has landmark names
  if (is.vector(model)) {
    lms <- model
  } else if (ncol(model) == 2) {
    lms <- array(model[ , 2])
    # If model is just a vector
  } else if (ncol(model) > 2) {
    stop("Model must either be a vector of model definitions or a data frame with the first column containing landmark names and the second of module definitions.")
  }
  
  symmet = corr
  symmet[upper.tri(symmet)] = t(symmet)[upper.tri(symmet)]
  
  if (is.vector(model)) {
    modNF <- stats::na.omit(cbind(1:length(model), lms))
  } else {
    modNF <- stats::na.omit(cbind(1:nrow(model), lms))
  }
  
  w <- unique(modNF[, 2])
  w <- w[w %in% test_modules]
  all_modules <- list()
  modules <- list()
  btw_mod = list()
  betweenModules = list()
  withinModules = list()
  unintegrated = list()
  betweenFloat = list()
  for(i in seq(length(w))){
    fg <- modNF[modNF[, 2] == w[i], ]
    l <- corr[as.numeric(fg[, 1]), as.numeric(fg[, 1])]
    modules[[i]] <- (as.array(l[!is.na(l)]))
  }
  
  names(modules) <- paste("Module", w)
  if (length(w) > 1) {
    cb <- utils::combn(w, 2)
    for (i in seq(dim(cb)[2])){
      fg1 <- modNF[modNF[, 2] == cb[1, i], ]
      fg2 <- modNF[modNF[, 2] == cb[2, i], ]
      btw <- symmet[
        as.integer(setdiff(fg2[, 1], fg1[, 1])),
        as.integer(setdiff(fg1[, 1], fg2[, 1]))]
      btw_mod[[i]] <- btw[!is.na(btw)]
    }
    names(btw_mod) <- paste(cb[1, ], "to", cb[2, ])
    betweenModules['betweenModules'] = list(as.vector(rle(unlist(btw_mod))$values))
    withinModules['withinModules'] = list(as.vector(rle(unlist(modules))$values))
  }
  all_modules = c(modules, btw_mod, withinModules, betweenModules)
  
  # Prepare data.
  td <- c(all_modules[[1]], all_modules[[2]], all_modules[[3]])
  group <- as.factor(c(
    rep(names(all_modules)[1], length(all_modules[[1]])),
    rep(names(all_modules)[2], length(all_modules[[2]])),
    rep(names(all_modules)[3], length(all_modules[[3]]))
  ))
  td <- data.frame(corrs = td, group = group)
  a <- aov(corrs ~ group, data = td)
  t <- TukeyHSD(a)
  
  res <- list(anova = a, tukeyhsd = t)
  
  if (plot) {
    # function to group variables that are not different.
    t_labs <- function(t, v){
      levs <- t[[v]][,4]
      labs <- data.frame(multcompView::multcompLetters(levs)['Letters'])
      labs$treatment <- rownames(labs)
      labs <- labs[order(labs$treatment) , ]
      return(labs)
    }
    labels <- t_labs(t, "group")
    
    # add labels to td for colouring.
    td$color <- NA
    for (i in seq_len(nrow(labels))) {
      td$color[td$group == labels$treatment[i]] <- as.character(labels$Letters[i])
    }
    
    pall <- c("#E69F00", "#56B4E9", "#7BB31A")
    p <- ggplot2::ggplot(td, aes(x = group, y = corrs, fill = color)) +
      geom_boxplot() +
      scale_fill_manual(values = pall) +
      xlab("") +
      ylab("Correlation") +
      guides(fill = guide_legend(title = "Significance")) +
      theme(
        legend.position = "none"
      )+
      stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                   width = .75, linetype = "dashed")
    print(p)
  }
  return(res)
}

###########################################################################################
#' phyloEmmli
#'
#' Takes landmarks and a phylogeny and then corrects the landmarks for the phylogeny according
#' to one of two methods, and then either returns the corrected landmarks, or calculates the
#' correlation matrix (using dotcorr) and fits EMMLi. Species missing from data or tree are
#' automatically dropped.
#' @name phyloEmmli
#' @param landmarks The landmarks. These can be in a 2D format with species as rownames,
#' and x, y, z as columns, or a 3D array with species names in the 3rd ([,,x]) dimension.
#' @param phylo A phylogeny describing the relationship between species in landmarks.
#' @param method Either "pgls" or "ic". If PGLS then corrected data are calculated as the
#' residuals of a phylogenetic least squares regression against 1 (\link[caper]{pgls}),
#' if IC then independent contrasts (\link[ape]{pic}).
#' @param EMMLi Logical - if TRUE EMMLi is fit and the results returned as well as the
#' phylogenetically correceted landmarks
#' @param ... Extra arguments required for EMMLi (at minimum models, and N_sample)
#' @export
#' @return A 2D or 3D array (depending on input) containing phylogenetically corrected
#' landmarks. If EMMLi = TRUE then the results of the EMMLi model are also returned along
#' with the phylogenetically corrected data in a two element list (data first element,
#' EMMLi output second element).

phyloEmmli <- function(landmarks, phylo, method = "pgls", EMMLi = FALSE, ...) {
  if (class(phylo) != "phylo") {
    stop("Tree must be an object of class 'phylo'.")
  }
  
  if(length(dim(landmarks)) == 3) {
    if (is.null(dimnames(landmarks)[[3]])) {
      stop("Landmarks must have species names in the 3rd dimension.")
    } else {
      dims <- 3
    }
  } else if (length(dim(landmarks)) == 2) {
    if (is.null(rownames(landmarks))) {
      stop("Landmarks must have species names as rownames.")
    } else {
      dims <- 2
    }
  } else {
    stop("Landmarks must be either a 2D or 3D array.")
  }
  
  if (EMMLi) {
    x <- list(...)
    if (!"mod" %in% names(x)) {
      stop("mod must be provided to fit EMMLi to phylo landmarks.")
    }
    if (!"N_sample" %in% names(x)) {
      stop("N_sample must be provided to fit EMMLi to phylo landmarks.")
    }
  }
  
  if (dims == 3) {
    sp_lm <- dimnames(landmarks)[[3]]
  } else if (dims == 2) {
    sp_lm <- rownames(landmarks)
  }
  
  sp_tr <- phylo$tip.label
  
  if (sum(sp_lm %in% sp_tr) != length(sp_lm)) {
    missing <- sum(!sp_lm %in% sp_tr)
    if (missing == length(sp_lm)) {
      stop("No species in dataset found on tree.")
    }
    print("Dropping", missing, "species from dataset - not in tree.")
    landmarks <- landmarks[sp_lm %in% sp_tr, ]
  }
  
  if (sum(sp_tr %in% sp_lm) != length(phylo$tip.label)) {
    missing <- sum(!sp_tr %in% sp_lm)
    if (missing == length(phylo$tip.label)) {
      stop("No species on tree found in dataset.")
    }
    print("Dropping", missing, "species from tree - not in dataset.")
    xtips <- phylo$tip.label[!sp_tr %in% sp_lm]
    phylo <- ape::drop.tip(phylo, xtips)
  }
  
  if (method == "pgls") {
    if (dims == 2) {
      landmarks <- as.data.frame(landmarks)
      landmarks$names <- rownames(landmarks)
      comp_data <- caper::comparative.data(phylo, landmarks, names = names)
      lms <- head(colnames(landmarks), -1)
      cl <- parallel::makeCluster(parallel::detectCores() - 2)
      parallel::clusterExport(cl, varlist = c("comp_data"), envir = environment())
      x <- parallel::parLapply(cl, lms, function(x) caper::pgls(formula(paste(x, "~", 1)), comp_data)$phyres)
      parallel::stopCluster(cl)
      phy_landmarks <- do.call(cbind, x)
      rownames(phy_landmarks) <- rownames(landmarks)
    } else if (dims == 3) {
      phy_landmarks <- landmarks
      cl <- parallel::makeCluster(parallel::detectCores() - 2)
      for (i in seq_len(dim(landmarks)[2])) {
        m_lms <- as.data.frame(t(landmarks[,i,]))
        m_lms$names <- rownames(m_lms)
        comp_data <- caper::comparative.data(phylo, m_lms, names = names)
        lms <- head(colnames(m_lms), -1)
        parallel::clusterExport(cl, varlist = c("comp_data"), envir = environment())
        x <- parallel::parLapply(cl, lms, function(x) caper::pgls(formula(paste(x, "~", 1)), comp_data)$phyres)
        tmp_landmarks <- do.call(cbind, x)
        rownames(tmp_landmarks) <- rownames(m_lms)
        phy_landmarks[,i,] <- t(tmp_landmarks)
      }
      parallel::stopCluster(cl)
    }
  } else if (method == "ic") {
    if (dims == 2) {
      landmarks <- landmarks[match(phylo$tip.label, rownames(landmarks)), ]
      phy_landmarks <- apply(landmarks, 2, function(x) ape::pic(x, phylo))
    } else if (dims == 3) {
      x <- dim(landmarks)
      phy_landmarks <- array(dim = c(x[1], x[2], x[3] - 1))
      pic_landmarks <- vector(mode = "list", length = dim(landmarks)[2])
      for (i in seq_len(dim(landmarks)[2])) {
        m_lms <- landmarks[,i,]
        m_lms <- m_lms[ , match(phylo$tip.label, colnames(m_lms))]
        phy_landmarks[,i,] <- t(apply(m_lms, 1, function(x) ape::pic(x, phylo)))
      }
    }
  }
  
  # Fit or don't fit EMMLi.
  if (!EMMLi) {
    res <- phy_landmarks
  } else if (EMMLi) {
    if (dims == 2) {
      arr <- geomorph::arrayspecs(phy_landmarks, ncol(phy_landmarks) / 3, 3)
    } else if (dims == 3) {
      arr <- landmarks
    }
    print("Computing correlation matrix...")
    corr <- paleomorph::dotcorr(arr)
    N_sample <- dim(landmarks)[3]
    emm <- EMMLi(corr = corr, ...)
    res <- list(EMMLi = emm, phy_landmarks = phy_landmarks)
  }
  
  return(res)
}

###########################################################################################
#' plotNetwork
#'
#' Plots a network diagram of the output of an EMMLi analysis. Nodes are proportionally
#' sized to within-module rho, and lines are larger and darker according to between-module
#' rho.
#' @name plotNetwork
#' @param rhos The rhos that come out of an \link[EMMLiv2]{EMMLi} analysis
#' @param module_names The names of the modules - if absent generic numbers are used.
#' If "rhos" then the nodes are named with the within-module rho.
#' @param linecolour The colour of the joining lines. If a single colour, then lines will
#' be that colour with width and transparency adjusted to the strength of the correlation.
#' If "viridis" then lines are coloured according to the viridis colour palette, with darker
#' colours corresponding to stronger correlations.
#' @param title Title for the plot.
#' @param layout A matrix describing the positions of each module on the canvas.
#' See \link[qgraph]{qgraph} for details.
#' @return A plotted network of the relationships between and within modules.
#' @examples
#' data("macacaCorrel")
#' data("macacaModels")
#' emm <- EMMLi(corr = macacaCorrel, mod = macacaModels, N_sample = 20)
#' plotNetwork(emm$rho[[1]], linecolour = "viridis")
#' @export

plotNetwork <- function(rhos, module_names = NULL, linecolour = "#56B4E9",
                        title = NULL, layout = NULL) {
  
  withins <- grep("Module*", colnames(rhos))
  nmodule <- length((grep("Module*", colnames(rhos))))
  rholist <- t(rhos)
  
  words <- strsplit(rownames(rholist), " ")
  
  plotcorr <- matrix(data = NA, nrow = nmodule, ncol = nmodule)
  modnums <- unlist(lapply(withins, function(x) strsplit(colnames(rhos)[x], " ")[[1]][2]))
  
  mods <- sapply(words[withins], function(x) paste(x, collapse = " "))
  mods <- gsub("Module", "M", mods)
  mods <- mods[order(sapply(mods, function(x)
    as.numeric(strsplit(x, "M ")[[1]][[2]])))]
  colnames(plotcorr) <- rownames(plotcorr) <- mods
  
  for (i in 1:(length(words) - 1)) {
    if (length(words[[i]]) == 2) {
      module <- paste("M", words[[i]][2])
      plotcorr[module, module] <- rholist[i, "MaxL_p"]
    }
    
    if (length(words[[i]]) == 3) {
      from_module <- paste("M", (words[[i]][1]))
      to_module <- paste("M", words[[i]][3])
      plotcorr[from_module, to_module] <- rholist[i, "MaxL_p"]
      plotcorr[to_module, from_module] <- rholist[i, "MaxL_p"]
    }
  }
  
  within <- diag(plotcorr)
  between <- plotcorr
  
  if (is.null(module_names)) {
    mod.names <- mods
  } else if (module_names == "rhos") {
    mod.names <- within
  }
  
  if (linecolour == "viridis") {
    cls <- viridis::viridis(100, direction = -1)
    vcols <- apply(between * 100, 1, function(x) cls[x])
    linecolour <- NULL
  } else {
    vcols <- NULL
  }
  
  qgraph::qgraph(between,
                 shape = "circle",
                 posCol = linecolour,
                 edge.color = vcols,
                 labels = mod.names,
                 vsize = within * 10,
                 diag = FALSE,
                 title = title,
                 layout = layout)
}

###########################################################################################
#' plotRandomSubsamples
#'
#' Generates network plots for a given number of randomly selected subsampled EMMLi analyses,
#' as returned by subSampleEmmli.
#' @name plotRandomSubsamples
#' @param subsasmples The output of the execution of \link[EMMLiv2]{subSampleEMMLi}
#' @param n The number of random subsamples to plot.
#' @param ... Additional arguments for plotNetwork (see \link[EMMLiv2]{plotNetwork})
#' @return A gridded plot with a network plot for each of n random subsamples.
#' @examples
#' ssemm <- subsampleEmmli(landmarks = landmarks, models = models, fractions = 0.4, min_landmark = 5, nsim = 25)
#' plotRandomSubsamples(subsamples = ssemm, n = 9, linecolour = "viridis")
#' @export

plotRandomSubsamples <- function(subsamples, n, ...) {
  samples <- sample(1:length(subsamples$results), n)
  par(mfrow = n2mfrow(n))
  for (i in samples) {
    rh <- subsamples$results[[i]]$rhos_best[[1]]
    plotNetwork(rh, ...)
  }
}

###########################################################################################
#' plotMeanNetwork
#'
#' Plots the mean network inferred from the rhos returned from a set of subsampled analyses.
#' @name plotMeanNetwork
#' @param summary The output of \link[EMMLiv2]{summariseResults}
#' @param ... Additional arguments for plotNetwork (see \link[EMMLiv2]{plotNetwork})
#' @export

plotMeanNetwork <- function(summary, ...) {
  bestRhos <- summary$bestRho
  n <- length(bestRhos)
  par(mfrow = n2mfrow(n))
  mod_names <- names(bestRhos)
  for (i in seq_along(bestRhos)) {
    yy <- colMeans(bestRhos[[i]])
    yy <- yy[1:(length(yy) - 3)]
    yy <- t(data.frame(MaxL_p = yy))
    plotNetwork(yy, title = mod_names[i], ...)
  }
}

###########################################################################################
#' identifyCandidates
#'
#' Identifies candidates for merging in simplifyEMMLi. Candidates are selected if their
#' between-module rho is larger than either of their within-module rhos by a factor of
#' 2*sd of the combined within-module correlations of the landmarks in the modules when
#' pooled.
#' @name simplifyEMMLi
#' @param fitted_emmli Fitted EMMLi model output.
#' @param corr_matrix The correlation matrix that EMMLi was fitted to.
#' @param models The original models that EMMLi tested.
#' @keywords internal

identifyCandidates <- function(fitted_emmli, corr_matrix, models) {
  fit_mods <- fitted_emmli$results
  best_mod <- names(fitted_emmli$rho)
  start_mod <- data.frame(data.point = models$Data.point,
                          original_best = models[ , grep(strsplit(best_mod, "\\.")[[1]][1], colnames(models))])
  rhos <- fitted_emmli$rho[[1]]
  
  # Seperate between from within.
  within_rhos <- rhos[ , grep("^Module", colnames(rhos))]
  between_rhos <- rhos[ , grep("to", colnames(rhos))]
  
  all_corrs <- getCorrs(fitted_emmli, models, corr_matrix)[[1]]
  between_rhos <- between_rhos[,order(between_rhos[2,], decreasing = TRUE)]
  
  # Select pairs that have a between rho greater than either of the withins by
  # a factor of 2*SD of the combined within-correlations of both modules pooled.
  
  pairs <- colnames(between_rhos)
  candidates <- list()
  for (i in seq_along(pairs)) {
    modules <- strsplit(pairs[1], " to ")[[1]]
    mod_1_rho <- rhos["MaxL_p" , grep(paste("Module", modules[1]), colnames(rhos))]
    mod_2_rho <- rhos["MaxL_p" , grep(paste("Module", modules[2]), colnames(rhos))]
    between_rho <- rhos["MaxL_p", grep(paste(modules[1], "to", modules[2]), colnames(rhos))]
    combined_mods <- c(
      all_corrs[[paste("Module", modules[1])]],
      all_corrs[[paste("Module", modules[2])]]
    )
    sd_comb <- sd(combined_mods)
    if (between_rho - mod_1_rho > 2 * sd_comb) {
      candidates[[length(candidates) + 1]] <- pairs[i]
    } else if (between_rho - mod_2_rho > 2 * sd_comb) {
      candidates[[length(candidates) + 1]] <- pairs[i]
    }
  }
  
  candidates <- lapply(candidates, function(x) as.numeric(strsplit(x, " to ")[[1]]))
  # return candidates, and the starting model (the best model from the fitted EMMLi)
  return(list(candidates = candidates, start_mod = start_mod))
}

###########################################################################################
#' simplifyEMMLi
#'
#' Experimental function - use with caution. Attempts to simplify a fitted EMMLi model,
#' looking for a simpler model that fits the data better by merging modules. Pairs of
#' modules are identified as candidates for merging if their between-module rho is higher
#' than either of the within-module rhos by more than 2 * SD of the modules combined. This
#' repeats until the model does not change or improve.
#'
#' Alternatively, pairs of modules can be offered as candidates for merging - in this
#' instance there is not exploration and just the suggested pairs are tested. Pairs are
#' offered as a list, where each element is a vector of two module numbers.
#' @name simplifyEMMLi
#' @param fitted_emmli Fitted EMMLi model output.
#' @param corr_matrix The correlation matrix that EMMLi was fitted to.
#' @param models The original models that EMMLi tested.
#' @param candidates A list of pairs of modules to test merging (each element is a vector of
#' length 2 with module numbers).
#' @param N_sample The sample size for the original EMMLi fit.
#' @param correction if "normal" uses the normal EMMLi calculation for K, if "new" then uses
#' the experimental adjustment to AICc (adding nmodules - 1 to K). Not reccommended.
#' @export

simplifyEMMLi <- function(fitted_emmli, corr_matrix, models, candidates = NULL,
                          correction = "normal", N_sample) {
  
  if (is.null(candidates)) {
    x <- identifyCandidates(fitted_emmli, corr_matrix, models)
    candidates <- x$candidates
    start_mod <- start_mod
    if (length(candidates) == 0) {
      stop("No candidate pairs found.")
    }
  } else {
    # Use user-supplied candidates and identify starting model.
    candidates <- candidates
    best_mod <- names(fitted_emmli$rho)
    start_mod <- data.frame(data.point = models$Data.point,
                            original_best = models[ , grep(strsplit(best_mod, "\\.")[[1]][1], colnames(models))])
  }
  
  repeat {
    if (length(candidates) == 0) {
      break
    }
    
    new_emms <- list()
    # test all candidate pairs.
    for (i in seq_along(candidates)) {
      test_mods <- start_mod
      # Combine the two modules in a new_model
      test_mods$new_model <- test_mods$original_best
      test_mods$new_model[test_mods$new_model == candidates[[i]][1]] <- candidates[[i]][2]
      # fit EMMLi
      new_emms[[i]] <- EMMLi(corr = corr_matrix, mod = test_mods, N_sample = N_sample, correction = correction)
      
      # Now if the best model is "original" go to the next i, otherwise break this loop.
      # if (!grepl("original_best", rownames(emm$results)[emm$results[ , "dAICc"] == 0])) {
      #   print("Better model found.")
      #   models <- test_mods
      #   new_emm <- emm
      #   break
      # }
      
    }
    best_liks <- rep(NA, length(candidates))
    for (i in seq_along(new_emms)) {
      # If the original model isn't the best record the likelihood of it.
      if (!grepl("original_best",
                 rownames(new_emms[[1]]$results)[new_emms[[1]]$results[ , "dAICc"] == 0])) {
        best_liks[i] <- new_emms[[i]]$results[new_emms[[i]]$results[ , "dAICc"] == 0, "MaxL"]
      }
    }
    
    # If best_liks is all NA then none of the merges helped, and break the loop.
    # Else take the one with the best likelihood and calculate new canidate pairs.
    # If there are no candidates the loop breaks, if there are, it repeats.
    if (all(is.na(best_liks))) {
      break
    } else {
      new_emm <- new_emms[which.max(best_liks)]
      x <- identifyCandidates(fitted_emmli = new_emm, models = test_mods,
                              corr = corr_matrix)
      candidates <- x$candidates
      start_mod <- x$start_mod
    }
  }
  # If new_emm exists it will be the simplest version of the model, so return it
  # and the test_mods. This will be the most recent test mods.
  # If it doesn't say so, and break.
  if (exists("new_emm")) {
    return (list(emmli = new_emm, model = test_mods))
  } else {
    stop("No simplification possible.")
  }
}

###########################################################################################
#' getCorrs
#'
#' This function takes a fitted EMMLi model and the models used to fit it, and then returns
#' all the correaltions for the best fitting model. This is recycled from the main EMMLi
#' function - Possible that this will change with the refactoring of EMMLi.
#' @name getCorrs
#' @param emm A fitted EMMLi model (output of EMMLi).
#' @param models A data frame defining the models. The first column should contain the landmark names
#' as factor or character. Subsequent columns should define which landmarks are contained within each
#' module with integers, factors or characters. If a landmark should be ignored for a specific model
#' (i.e., it is unintegrated in any module), the element should be NA.
#' @param corr The original correlation matrix used in the EMMLi analysis that generated
#' the input.
#' @return The correlations within and between modules of the best model.
#' @export

getCorrs <- function(emm, models, corr) {
  best_mod <- rownames(emm$results)[which(emm$res[ , "dAICc"] == 0)]
  best_mod <- strsplit(best_mod, " ")[[1]][1]
  best_mod <- paste(head(strsplit(best_mod, "\\.")[[1]], n = -2), collapse = ".")
  
  symmet = corr
  symmet[upper.tri(symmet)] = t(symmet)[upper.tri(symmet)]
  
  model <- paste0("models$", best_mod)
  lms <- array(eval(parse(text = model)))
  modNF <- stats::na.omit(cbind(1:nrow(models), lms))
  w <- unique(modNF[, 2])
  
  all_modules <- list()
  modules <- list()
  btw_mod = list()
  betweenModules = list()
  withinModules = list()
  
  for(i in seq(length(w))){
    # identify landmarks within a class
    fg <- modNF[modNF[, 2] == w[i], ]
    
    # coefficients between identified landmarks.
    l <- corr[as.numeric(fg[, 1]), as.numeric(fg[, 1])]
    modules[[i]] <- (as.array(l[!is.na(l)]))
  }
  names(modules) <- paste("Module", w)
  
  if (length(w) > 1) {
    # make combinations of modules for between module.
    cb <- utils::combn(w, 2)
    for (i in seq(dim(cb)[2])){
      fg1 <- modNF[modNF[, 2] == cb[1, i], ]
      fg2 <- modNF[modNF[, 2] == cb[2, i], ]
      
      # setdiff(A,B) - present in A but not in B
      btw <- symmet[as.integer(setdiff(fg2[, 1], fg1[, 1])), as.integer(setdiff(fg1[, 1], fg2[, 1]))]
      btw_mod[[i]] <- btw[!is.na(btw)]
    }
    names(btw_mod) <- paste(cb[1, ], "to", cb[2, ])
    betweenModules['betweenModules'] = list(as.vector(rle(unlist(btw_mod))$values))
    withinModules['withinModules'] = list(as.vector(rle(unlist(modules))$values))
  }
  
  all_modules[[model]] = c(modules, btw_mod, withinModules, betweenModules)
  
  return(all_modules)
}


###########################################################################################
#' subsampleLandmarks
#'
#' This function randomly removes landmarks from a dataset, and adjusts the corresponding
#' models object to match the newly subsampled landmark data. Internal, called by
#' subsampleEMMLi.
#' @name subsampleLandmarks
#' @param landmarks A 3D array of xyz landmarks to subsample from. Will be turned into a correlation
#' matrix using \link[paleopmorph]{dotcorr} for EMMLi analysis after subsampling.
#' @param fraction The decimal fraction to subsample down to. e.g. 0.2 will return 20% of the
#' original landmarks.
#' @param models A data frame defining the models. The first column should contain the landmark names
#' as factor or character. Subsequent columns should define which landmarks are contained within each
#' module with integers, factors or characters. If a landmark should be ignored for a specific model
#' (i.e., it is unintegrated in any module), the element should be NA.
#' @param min_landmark The minimum number of landmarks to subsample to. This ensures that a module isn't
#' totally removed during random subsampling. When subsampling causes a landmark to be removed or to be
#' subsampled below this threshold landmarks are drawn from the original module at random and added back in.
#' This means that sometimes (especially with low subsampling fractions and/or the presence of small
#' modules in a model) the actual subsampling level is higher than the requested subsampling. In these
#' cases a warning is printed to the screen.
#' @keywords internal

subsampleLandmarks <- function(landmarks, fraction, models, min_landmark) {
  sampleSpecimen <- function(lms, kps) {
    return(lms[kps, ])
  }
  
  x <- lapply(models[2:length(models)], table)
  x <- sapply(x, length)
  kps <- sort(sample(nrow(landmarks[,,1]), round(nrow(landmarks[,,1]) * fraction, 0)))
  m <- models[ , 2:ncol(models)]
  for (i in 1:ncol(m)) {
    mod_original <- unique(m[ , i])
    ss_mod_count <- table(m[kps, i])
    
    if (length(ss_mod_count) != length(mod_original)) {
      missing <- mod_original[!mod_original %in% names(ss_mod_count)]
      for (j in missing) {
        kps <- c(kps, sample(which(m[ , i] == j), min_landmark))
      }
      kps <- sort(unique(kps))
    }
    
    if (any(ss_mod_count < min_landmark)) {
      short <- names(ss_mod_count)[ss_mod_count < min_landmark]
      for (j in short) {
        kps <- c(kps, sample(which(m[ , i] == j), min_landmark))
      }
      kps <- sort(unique(kps))
    }
  }
  
  res <- array(dim = c(length(kps), 3, dim(landmarks)[3]))
  for (i in 1:dim(landmarks)[3]) {
    res[,,i] <- sampleSpecimen(landmarks[,,i], kps)
  }
  
  true_subsample = round(dim(res)[1] / dim(landmarks)[1], 2)
  if (true_subsample > fraction) {
    print(paste("+++ WARNING: Actual subsample is ", true_subsample, "+++"))
  }
  
  new_models <- models[kps, ]
  return(list(landmarks = res, models = new_models, true_subsample = true_subsample))
}


###########################################################################################
#' subsampleEMMLi
#'
#' Analyse random subsamples of a dataset repeatedly using EMMLi. This function subsamples
#' a 2D array of landmarks to a given fraction of it's original size and then fits EMMLi
#' to the smaller dataset, returning the results. This can be either done repeatedly for
#' a single subsampling fraction, or for a range of subsampling fractions (e.g. to investigate
#' the effects of increasing subsampling). In all cases EMMLi is fit to a correlation matrix
#' calculated from the subsampled landmarks using \link[paleomorph]{dotcorr}.
#'
#' If a single fraction is provided then an nrep argument is also required, and the function
#' will subsample the data to the given fraction nrep times. Alternatively, if a range of
#' fractions is given (in a vector) the nrep argument is not required, and the function will
#' subsample the data and fit EMMLi once for each fraction in the given fractions vector.
#' @name subsampleEMMLi
#' @param landmarks A 2D array of xyz landmarks to subsample from. Will be turned into a correlation
#' matrix using \link[paleopmorph]{dotcorr} for EMMLi analysis after subsampling.
#' @param fractions Either a single subsampling fraction (in which case nrep is required) or a
#' vector of fractions. Specified in decimal format, i.e. a fraction of 0.2 will subsample down
#' to 20\% of the original number of landmarks.
#' @param models A data frame defining the models. The first column should contain the landmark names
#' as factor or character. Subsequent columns should define which landmarks are contained within each
#' module with integers, factors or characters. If a landmark should be ignored for a specific model
#' (i.e., it is unintegrated in any module), the element should be NA.
#' @param min_landmark The minimum number of landmarks to subsample to. This ensures that a module isn't
#' totally removed during random subsampling. When subsampling causes a landmark to be removed or to be
#' subsampled below this threshold landmarks are drawn from the original module at random and added back in.
#' This means that sometimes (especially with low subsampling fractions and/or the presence of small
#' modules in a model) the actual subsampling level is higher than the requested subsampling. In these
#' cases a warning is printed to the screen.
#' @param aic_cut This is the threshold of dAICc below which two models are considered to be not different.
#' When this occurs multiple models are be returned as the best. Defaults to 0 (i.e., only the best
#' fitting model is returned, regardless of how close other models are.)
#' @param return_corr Logical - if TRUE then the full correlations of within and between modules are
#' returned for the single best model in addition to the EMMLi results. Defaults to FALSE.
#' @param nrep If a single subsampling fraction, this is the number of times that subsampling fraction
#' is repeated.
#' @return A list of n elements, where n is either the number of subsampling fractions, or nrep. Each
#' element contains the output of an EMMLi analysis on the subasampled data, consisting of four or five
#' elements:
#'
#'   - the best model(s) description
#'
#'   - the rho list(s) for the best model(s)
#'
#'   - the data used in the EMMLi analysis (two elements - the subsampled data, and the corresponding models)
#'
#'   - the true level of subsampling
#'
#'   - (optional) the correlation matrix of the single best model.
#'
#' @export

subsampleEMMLi <- function(landmarks, fractions, models, min_landmark, aic_cut = 0,
                           return_corr = FALSE, nrep = NULL) {
  
  if (!is.null(nrep)) {
    if (length(fractions) > 1) {
      stop("Only one subsampling fraction is allowed when multiple subsampling simulations are run.")
    }
    fractions <- rep(fractions, nrep)
  }
  
  fitEmmli <- function(landmarks, models, fraction, min_landmark, aic_cut) {
    dat <- subsampleLandmarks(landmarks = landmarks, fraction = fraction, models = models,
                              min_landmark = min_landmark)
    c <- paleomorph::dotcorr(dat$landmarks)
    
    emm <- EMMLi(corr = c, mod = dat$models, N_sample = dim(landmarks)[3], all_rhos = TRUE)
    emmli <- as.data.frame(emm$results)
    best_models <- emmli[emmli$dAICc <= aic_cut, ]
    best_names <- rownames(best_models)
    best_models <- do.call(rbind, best_models)
    colnames(best_models) <- best_names
    names(emm$all_rhos) <- sapply(names(emm$all_rhos), function(x) strsplit(x, "\\$")[[1]][2])
    names(emm$all_rhos) <- unlist(strsplit(names(emm$all_rhos), "1$"))
    rhos <- emm$all_rhos[names(emm$all_rhos) %in% best_names]
    all_data <- list(landmarks = dat$landmarks, models = dat$models)
    
    if (return_corr) {
      best_corr <- getCorrs(emm, dat$models, c)
      res <- list(best_models = best_models, rhos_best = rhos, best_corr = best_corr,
                  all_data = all_data, true_subsample = dat$true_subsample)
    } else {
      res <- list(best_models = best_models, rhos_best = rhos, all_data = all_data,
                  true_subsample = dat$true_subsample)
    }
    # print(ncol(best_models))
    # if (ncol(best_models) > 1) {
    #   print("multiple best")
    # }
    return(res)
  }
  
  res <- pbapply::pblapply(fractions, function(x) fitEmmli(fraction = x, landmarks = landmarks,
                                                           models = models, min_landmark = min_landmark, aic_cut = aic_cut))
  
  names(res) <- fractions
  
  return(res)
}

###########################################################################################
#' getBestMods
#'
#' Collates the best fitting models from each of the analyses in a subsampledEMMLi analysis
#' @name getBestMods
#' @param rs The output of subsampleEMMLi
#' @keywords internal

getBestMods <- function(rs) {
  bestMods <- vector(mode = "list", length = length(rs))
  for (i in seq_along(rs)) {
    bestMods[[i]] <- cbind(t(rs[[i]][[1]]),
                           subsample = as.numeric(names(rs)[i]),
                           true_subsample = as.numeric(rs[[i]]$true_subsample))
  }
  bestMods <- do.call(rbind, bestMods)
  return(bestMods)
}

###########################################################################################
#' getRhos
#'
#' Collates the rhos for each of the analyses in a subsampleEMMLi output. Collects them
#' together based on the model that derived them, like with like.
#' @name getRos
#' @param rs The output of subsampleEMMLi
#' @keywords internal

getRhos <- function(rs) {
  rho_names <- unlist(lapply(rs, function(x) names(x[[2]])))
  allRhos <- vector(mode = "list", length = length(unique(rho_names)))
  names(allRhos) <- unique(rho_names)
  for (i in seq_along(unique(rho_names))) {
    allRhos[[i]] <- vector(mode = "list", length = sum(rho_names == names(allRhos)[[i]]))
  }
  
  for (i in seq_along(rs)) {
    c_rs <- rs[[i]]$rhos_best
    for (j in seq_along(c_rs)) {
      slt <- which(names(allRhos) == names(c_rs)[j])
      first_null <- which(sapply(allRhos[[slt]], is.null))[[1]]
      xx <- cbind(c_rs[[j]],
                  subsample = as.numeric(names(rs)[i]),
                  true_subsample = as.numeric(rs[[i]]$true_subsample))
      bets <- grep("to", colnames(xx))
      ordered <- sapply(bets, function(x)
        paste(gtools::mixedsort(strsplit(colnames(xx)[x], " to ")[[1]]), collapse = " to ")
      )
      colnames(xx)[bets] <- ordered
      xx <- xx[ , order(colnames(xx), decreasing = TRUE)]
      xx <- xx[ , c(4:ncol(xx), 1:3)]
      allRhos[[slt]][[first_null]] <- xx
    }
  }
  
  for (i in seq_along(allRhos)) {
    allRhos[[i]] <- do.call(rbind, allRhos[[i]])
    allRhos[[i]] <- allRhos[[i]][rownames(allRhos[[i]]) != "MaxL", ]
  }
  return(allRhos)
}

###########################################################################################
#' subsampleSummary
#'
#' Summarises the output of subsampleEMMLi to allow easier comparison between multiple
#' subsampled analyses
#' @name subsampleSummary
#' @param subsamples The output of subsampleEMMLi
#' @return A list of two elements. The first of these is a matrix detailing the best-fitting
#' model(s) from each subsample, containing the model name, maximum likelihood, number of
#' parameters, number of landmarks, AICc, dAICc, posterior probabiliry, requested subsample
#' level and true subsample level. The second element is a list of n elements, where n is
#' the number of different models that emerged as the best fitting. Each of these elements
#' contains a matrix detailing the maximum likihood value of rho for all within- and
#' between-module correlations (columns) for each EMMLi analysis (rows).
#' @export

subsampleSummary <- function(subsamples) {
  bestMods <- getBestMods(subsamples)
  bestRhos <- getRhos(subsamples)
  return(list(bestModels = bestMods, bestRho = bestRhos))
}