######################################
# Quality control and transformation #
######################################
library(phangorn)
library(phyloseq)
library(zCompositions)
library(plotly)
library(dplyr)
library(forestplot)
library(quantreg)
library(fossil)
library(picante)
library(entropart)
library(gridGraphics)
library(gridExtra)
library(ggplot2)
library(patchwork)
library(ggthemes)
library(erer)
library(DiagrammeR)
library(stringr)
library(devtools)
library(betareg)

library(caret)
library(randomForest)
library(glmnet)

# surv.random.forest <- function(taxa.out, taxa.names.out, surv.dat, n.tree = 1000, n.fold = 10, tra.ratio = 0.8, ranks.upto = 5) { # taxa.out <- taxa.out$clr
#   rf.fit <- list()
#   mse <- numeric()
#   if (is.null(ncol(surv.dat$survtime))) {
#     survtime <- surv.dat$survtime
#   } else {
#     survtime <- matrix(ncol = 2)
#     survtime <- surv.dat$survtime
#   }
#   
#   censor <- surv.dat$status
#   
#   for(i in 1:ranks.upto){
#     p <- ncol(taxa.out[[i]])
#     dat <- taxa.out[[i]][rownames(surv.dat),]
#     imp <- numeric()
#     
#     tra.ind <- sample(1:nrow(dat), ceiling(nrow(dat)*tra.ratio))
#     
#     tra.dat <- dat[tra.ind,]
#     tes.dat <- dat[-tra.ind,]
#     
#     cens.tra.dat <- data.frame( "censor" = censor[tra.ind])
#     cens.tes.dat <- data.frame( "censor" = censor[-tra.ind])
#     
#     if (is.null(ncol(surv.dat$survtime))) {
#       time.tra.dat <- data.frame( "survtime" = survtime[tra.ind])
#       time.tes.dat <- data.frame( "survtime" = survtime[-tra.ind])
#       tra.dat <- cbind(time.tra.dat, cens.tra.dat, tra.dat)
#       tes.dat <- cbind(time.tes.dat, cens.tes.dat, tes.dat)
#     } else {
#       time.tra.dat1 <- data.frame( "survtime1" = survtime[tra.ind,1])
#       time.tes.dat1 <- data.frame( "survtime1" = survtime[-tra.ind,1])
#       time.tra.dat2 <- data.frame( "survtime2" = survtime[tra.ind,2])
#       time.tes.dat2 <- data.frame( "survtime2" = survtime[-tra.ind,2])
#       tra.dat <- cbind(time.tra.dat1, time.tra.dat2, cens.tra.dat, tra.dat)
#       tes.dat <- cbind(time.tes.dat1, time.tes.dat2, cens.tes.dat, tes.dat)
#     }
#     if(i == 1){
#       print(cbind(time.tra.dat, cens.tra.dat))
#       print(cbind(time.tes.dat, cens.tra.dat))
#     }
#     
#     
#     if (is.null(ncol(survtime))) {
#       model <- tune(Surv(survtime, censor) ~., data = tra.dat, sampsize = function(x) {x * .632})
#       # print(c(tune(Surv(survtime, censor) ~., data = tra.dat, folds = 10), tune(Surv(survtime, censor) ~., data = tra.dat, folds = 5)))
#       # print(tune(Surv(survtime, censor) ~., data = tra.dat, folds = 10))
#       mtry <- model$optimal[[2]]
#       rf.fit[[i]] <- rfsrc(Surv(survtime, censor) ~ ., data = tra.dat, mtry = mtry, ntree = n.tree, importance = TRUE)
#     } else {
#       model <- tune(Surv(survtime2, censor) ~., data = tra.dat, sampsize = "????")
#       # print(c(tune(Surv(survtime, censor) ~., data = tra.dat, folds = 10), tune(Surv(survtime, censor) ~., data = tra.dat, folds = 5)))
#       # print(tune(Surv(survtime, censor) ~., data = tra.dat, folds = 10))
#       mtry <- model$optimal[[2]]
#       rf.fit[[i]] <- rfsrc(Surv(survtime1, survtime2, censor) ~ ., data = tra.dat, mtry = mtry, ntree = n.tree, importance = TRUE)
#     }
#     survtime.hat <- predict(rf.fit[[i]], newdata = tes.dat)$survival
#     mse[i] <- mean((survtime.hat - survtime[-tra.ind])^2)
#     names(rf.fit[[i]]$importance) <- taxa.names.out$names[[i]]
#   }
#   if(ranks.upto == 6) {
#     names(rf.fit) <- c("phylum", "class", "order", "family", "genus", "species")
#   } else if(ranks.upto == 5) {
#     names(rf.fit) <- c("phylum", "class", "order", "family", "genus")
#   }
#   return(list(rf.fit = rf.fit, MSE = mse))
# }

surv.random.forest <- function(taxa.out, taxa.names.out, surv.dat, n.tree = 1000, bag.size = function(x){min(x * .632, max(150, x ^ (3/4)))}, ranks.upto = 5) { # taxa.out <- taxa.out$clr
  rf.fit <- list()
  mse <- numeric()
  if (is.null(ncol(surv.dat$survtime))) {
    survtime <- surv.dat$survtime
  } else {
    survtime <- matrix(ncol = 2)
    survtime <- surv.dat$survtime
  }
  
  colnames(surv.dat)[2] <- "censor"
  
  for(i in 1:ranks.upto){
    p <- ncol(taxa.out[[i]])
    dat <- taxa.out[[i]][rownames(surv.dat),]
    imp <- numeric()
    
    if (is.null(ncol(survtime))) {
      model <- tune(Surv(survtime, censor) ~., data = cbind(surv.dat, dat), sampsize = bag.size)
      mtry <- model$optimal[[2]]
      rf.fit[[i]] <- rfsrc(Surv(survtime, censor) ~ ., data = cbind(surv.dat, dat), mtry = mtry, sampsize = bag.size, ntree = n.tree, importance = TRUE)
    } else {
      model <- tune(Surv(survtime2, censor) ~., data = dat, sampsize = bag.size)
      mtry <- model$optimal[[2]]
      rf.fit[[i]] <- rfsrc(Surv(survtime1, survtime2, censor) ~ ., data = dat, mtry = mtry, ntree = n.tree, importance = TRUE)
    }
    names(rf.fit[[i]]$importance) <- taxa.names.out$names[[i]]
  }
  if(ranks.upto == 6) {
    names(rf.fit) <- c("phylum", "class", "order", "family", "genus", "species")
  } else if(ranks.upto == 5) {
    names(rf.fit) <- c("phylum", "class", "order", "family", "genus")
  }
  
  return(list(rf.fit = rf.fit))
}

surv.en.reg.fit <- function(taxa.out, surv.dat, taxa.names.rank.out, rank, n.folds = 10, loss.func = "MSE", seed.num = 1785) {
  
  set.seed(seed.num)
  
  if (is.null(ncol(surv.dat$survtime))) {
    y <- Surv(surv.dat$survtime, surv.dat$status)
  } else {
    y <- Surv(surv.dat$survtime[,1], surv.dat$survtime[,2], surv.dat$status)
  }
  
  if (n.folds == "LOOCV") {
    n.folds = length(y)
  }
  if (loss.func == "MSE") {
    loss.func = "mse"
  }
  
  X <- taxa.out[[rank]][rownames(surv.dat),]
  print(dim(X))
  alpha <- seq(0.05, 0.95, 0.05)
  
  lambda <- numeric()
  mse <- numeric()
  cv.kcv <- list()
  
  for (i in 1:length(alpha)) {
    cv.kcv[[i]] <- cv.glmnet(x = as.matrix(X), y = y, family = "cox", alpha = alpha[i], nfolds = n.folds, type.measure = loss.func)
    lambda[i] <- cv.kcv[[i]]$lambda.min
    mse[i] <- min(cv.kcv[[i]]$cvm)
  }
  ind.min.mse <- which.min(mse)
  elnet.kcv.err <- min(mse)
  
  lasso.kcv <- glmnet(as.matrix(X), y, family = "cox", alpha = alpha[ind.min.mse], standardize = TRUE, nfolds = n.folds, type.measure = loss.func)
  
  elnet.kcv <- glmnet(as.matrix(X), y, family = "cox", alpha = alpha[ind.min.mse], standardize = TRUE, lambda = lambda[ind.min.mse], nfolds = n.folds, type.measure = loss.func)
  
  coef.est <- as.numeric(elnet.kcv$beta)
  short.taxa <- taxa.names.rank.out$names[[rank]]
  full.taxa <- names(taxa.out[[rank]])
  cv.err <- elnet.kcv.err
  
  ind <- order(coef.est)
  ord.coef.est <- coef.est[ind]
  ord.short.taxa <- short.taxa[ind]
  ord.full.taxa <- full.taxa[ind]
  sig.ind <- which(ord.coef.est != 0)
  sig.ord.coef.est <- ord.coef.est[sig.ind]
  sig.ord.short.taxa <- ord.short.taxa[sig.ind]
  sig.ord.full.taxa <- ord.full.taxa[sig.ind]
  
  sig.coef.out <- as.data.frame(cbind(sig.ord.short.taxa, sig.ord.full.taxa, sig.ord.coef.est))
  colnames(sig.coef.out) <- c("Taxon", "Full Path", "Est.")
  
  cv.tab <- as.data.frame(rbind(format(alpha[ind.min.mse], digit = 3), format(lambda[ind.min.mse], digit = 3), format(cv.err, digit = 3)))
  colnames(cv.tab) <- ""
  rownames(cv.tab) <- c("Alpha", "Lambda", "CV Error")
  
  coef <- as.matrix(elnet.kcv$beta)
  colnames(coef) <- "Coefficient"
  return(list(sig.coef.out = sig.coef.out, coefficient = coef, cv.err = cv.err, alpha = alpha[ind.min.mse], lambda = lambda[ind.min.mse], 
              cv.tab = cv.tab, cv.kcv = cv.kcv[[ind.min.mse]], lasso.kcv = lasso.kcv))
}


lasso.en.barplot <- function(out) {

    sig.ord.coef.est <- as.numeric(out$sig.coef.out$Est.)
    #sig.ord.coef.est <- rep(sig.ord.coef.est, rep)
    sig.ord.taxa <- out$sig.coef.out$Taxon
    #sig.ord.taxa <- rep(sig.ord.taxa, rep)
    ind.pos <- which(sig.ord.coef.est > 0)
    ind.neg <- which(sig.ord.coef.est < 0)
    
    if (length(sig.ord.coef.est) == 0) {
      plot.new()
      text(x = 0.5, y = 0.5, "No significant taxa are found.", 
           cex = 1.2, col = "black")
    } else {
      
      # ind.long <- which(nchar(sig.ord.taxa) > 43)
      # sig.ord.taxa <- substr(sig.ord.taxa, 1, 43)
      # sig.ord.taxa[ind.long] <- paste("*", sig.ord.taxa[ind.long])
      
      if (length(sig.ord.coef.est) > 100) {
        sig.ord.coef.est <- c(sig.ord.coef.est[1:50], sig.ord.coef.est[(length(sig.ord.coef.est) - 49):length(sig.ord.coef.est)])
        sig.ord.taxa <- c(sig.ord.taxa[1:50], sig.ord.taxa[(length(sig.ord.coef.est) - 49):length(sig.ord.coef.est)])
        ind.pos <- which(sig.ord.coef.est > 0)
        ind.neg <- which(sig.ord.coef.est < 0)
      }
      
      if (length(sig.ord.coef.est) <= 50) {
        si <- (1.5-0.1*length(sig.ord.coef.est)^(0.5))*length(sig.ord.coef.est)^(0.8)
        par(mar=c(21 - si, 8, 19 - si, 5))
        barplot(height = sig.ord.coef.est, width = 2, cex.names = 0.7, cex.axis = 0.7, cex.lab = 0.7, space = 0.5, col = c(rep("lightblue", length(ind.neg)), rep("pink2", length(ind.pos))), horiz = TRUE, xlab = "Coefficient", names = sig.ord.taxa, las = 2)
      }
      
      if (length(sig.ord.coef.est) > 50 && length(sig.ord.coef.est) <= 100) {
        par(mfrow = c(1,2))
        #par(mar=c(3,17,1,2))
        par(mar=c(3,8,3,3))
        ind.1 <- round(length(sig.ord.coef.est)/2)
        sig.ord.coef.est.2 <- sig.ord.coef.est[(ind.1+1):length(sig.ord.coef.est)]
        sig.ord.taxa.2 <- sig.ord.taxa[(ind.1+1):length(sig.ord.coef.est)]
        ind.pos.2 <- which(sig.ord.coef.est.2 > 0)
        ind.neg.2 <- which(sig.ord.coef.est.2 < 0)
        barplot(height = sig.ord.coef.est.2, space = 0.5, cex.names = 0.6, cex.axis = 0.6, cex.lab = 0.7, col = c(rep("lightblue", length(ind.neg.2)), rep("pink2", length(ind.pos.2))), horiz = TRUE, xlab = "Coefficient", names = sig.ord.taxa.2, las = 2)
        sig.ord.coef.est.1 <- sig.ord.coef.est[1:ind.1]
        sig.ord.taxa.1 <- sig.ord.taxa[1:ind.1]
        ind.pos.1 <- which(sig.ord.coef.est.1 > 0)
        ind.neg.1 <- which(sig.ord.coef.est.1 < 0)
        barplot(height = sig.ord.coef.est.1, space = 0.5, cex.names = 0.6, cex.axis = 0.6, cex.lab = 0.7, col = c(rep("lightblue", length(ind.neg.1)), rep("pink2", length(ind.pos.1))), horiz = TRUE, xlab = "Coefficient", names = sig.ord.taxa.1, las = 2)
      }
    }
}


plot.importance <- function(out, num.display) {
  ind <- order(abs(out$importance), decreasing = TRUE)[1:num.display]
  out$importance <- out$importance[ind][!is.na(out$importance[ind])]
  plot.rfsrc(out, verbose = TRUE, plots.one.page = T)
  
}

surv.en.duplicate.list <- function(duplicate, duplicate.full.list){
  if(is.null(duplicate)) {
    text(x=0.5, y=0.5, "")
  } else {
    duplicate.taxa <- duplicate.full.list[duplicate]
    par(mar=c(0, 0.5, 0, 0.5))
    text(x=0, y=0.5, paste(duplicate.taxa, collapse = "\n"), cex = 0.75, adj = c(0, NA))
  }
}

surv.elastic.net <- function(taxa.out, survtime) {
  out <- list()
  
  for (i in 1:6) {
    dat <- taxa.out[[i]]
    alpha <- seq(0, 1, 0.1)
    lambda <- numeric()
    mse <- numeric()
    cv.loocv <- list()
    
    for (k in 1:length(alpha)) {
      cv.loocv[[k]] <- cv.glmnet(x = as.matrix(dat), y = survtime, alpha = alpha[k], nfolds = nrow(dat), family = "cox",
                                 type.measure = "mse")
                                 #type.measure = "mse")
      lambda[k] <- cv.loocv[[k]]$lambda.min
      mse[k] <- min(cv.loocv[[k]]$cvm)
    }
    
    ind.min.mse <- which.min(mse)
    #print(ind.min.mse)
    elnet.loocv <- glmnet(dat, survtime, alpha = alpha[ind.min.mse], family = "cox",
                          lambda = lambda[ind.min.mse], standardize = TRUE)
    
    tra.ind <- sample(1:nrow(dat), ceiling(nrow(dat)/2))
    
    tra.dat <- dat[tra.ind,]
    tes.dat <- dat[-tra.ind,]
    
    elnet.loocv.pr <- predict(elnet.loocv, s = lambda, newx = as.matrix(tes.dat))
    cv.error <- mean((elnet.loocv.pr - survtime[-tra.ind])^2)
    
    out[[i]] <- list(coef = as.matrix(coef(elnet.loocv)), CV_Error = cv.error, cv.kcv = cv.loocv[[ind.min.mse]]$cv.kcv, lasso.kcv = elnet.loocv)
  }
  names(out) <- c("phylum", "class", "order", "family", "genus", "species")
  return(out)
}

surv.lasso.reg.fit <- function(taxa.out, survtime, censoring, taxa.names.rank.out, rank, n.folds = c(5, 10, "LOOCV"), loss.func = "MSE", seed.num = 1785){
  
  set.seed(seed.num)
  
  y = Surv(survtime, censoring)
  
  if (n.folds == "LOOCV") {
    n.folds = length(y)
  }
  if (loss.func == "MSE") {
    loss.func = "mse"
  }
  
  out <- list()
  
  #for (i in 1:6) {
  
  X = taxa.out[[rank]]
  
  cv.kcv <- cv.glmnet(x = as.matrix(X), y = y, family = "cox", alpha = 1, nfolds = n.folds, type.measure = loss.func)
  lambda.kcv <- cv.kcv$lambda.min
  lasso.kcv.err <- min(cv.kcv$cvm)
  
  lasso.kcv.data <- glmnet(x = as.matrix(X), y = y, family = "cox", alpha = 1, type.measure = loss.func)
  
  # par(mar=c(6,5,6,2))
  # par(mfrow=c(2,1))
  # 
  # plot(cv.kcv, ylab = "MSE", xlab = expression(paste("log(", lambda, ")", sep="")))
  # plot(lasso.kcv.data, "lambda", ylab = "Coefficient", xlab = expression(paste("log(", lambda, ")", sep="")))
  # abline(v = log(lambda.kcv), col = "black")
  
  lasso.kcv <- glmnet(x = X, y = y, family = "cox", alpha = 1, type.measure = loss.func, lambda = lambda.kcv)
  coef.est <- as.numeric(lasso.kcv$beta)
  short.taxa <- taxa.names.rank.out$names[[rank]]
  full.taxa <- names(taxa.out[[rank]])
  cv.err <- lasso.kcv.err
  
  ind <- order(coef.est)
  ord.coef.est <- coef.est[ind]
  ord.short.taxa <- short.taxa[ind]
  ord.full.taxa <- full.taxa[ind]
  sig.ind <- which(ord.coef.est != 0)
  sig.ord.coef.est <- ord.coef.est[sig.ind]
  sig.ord.short.taxa <- ord.short.taxa[sig.ind]
  sig.ord.full.taxa <- ord.full.taxa[sig.ind]
  
  sig.coef.out <- as.data.frame(cbind(sig.ord.short.taxa, sig.ord.full.taxa, sig.ord.coef.est))
  colnames(sig.coef.out) <- c("Taxon", "Full Path", "Est.")
  #out[[i]] <- list(sig.coef.out = sig.coef.out, cv.err = cv.err, plot = plot)
  #}
  
  #names(out) <- c("phylum", "class", "order", "family", "genus", "species")
  return(list(sig.coef.out = sig.coef.out, cv.err = cv.err, cv.kcv = cv.kcv, lasso.kcv = lasso.kcv.data, lambda.kcv = lambda.kcv))
}

surv.alasso.reg.fit <- function(taxa.out, survtime, censoring, taxa.names.rank.out, rank, n.folds = c(5, 10, "LOOCV"), loss.func = "MSE", seed.num = 1785) {
  
  set.seed(seed.num)
  
  y = Surv(survtime, censoring)
  
  if (n.folds == "LOOCV") {
    n.folds = length(y)
  }
  if (loss.func == "MSE") {
    loss.func = "mse"
  }
  
  X = taxa.out[[rank]]
  
  ridge.cv.kcv <- cv.glmnet(x = as.matrix(X), y = y, family = "cox", alpha = 0, nfolds = n.folds, type.measure = loss.func)
  ridge.lambda.kcv <- ridge.cv.kcv$lambda.min
  
  ridge.kcv <- glmnet(as.matrix(X), y, family = "cox", alpha = 0, lambda = ridge.lambda.kcv, type.measure = loss.func)
  ridge.weight <- 1/abs(coef(ridge.kcv))
  
  cv.kcv <- cv.glmnet(x = as.matrix(X), y = y, family = "cox", alpha = 1, nfolds = n.folds, type.measure = loss.func, penalty.factor = ridge.weight)
  lambda.kcv <- cv.kcv$lambda.min
  lasso.kcv.err <- min(cv.kcv$cvm)
  
  lasso.kcv.data <- glmnet(x = as.matrix(X), y = y, family = "cox", alpha = 1, type.measure = loss.func, penalty.factor = ridge.weight)
  
  # par(mar=c(6,5,6,2))
  # par(mfrow=c(2,1))
  # 
  # plot(cv.kcv, ylab = "MSE", xlab = expression(paste("log(", lambda, ")", sep="")))
  # plot(lasso.kcv.data, "lambda", ylab = "Coefficient", xlab = expression(paste("log(", lambda, ")", sep="")))
  # abline(v = log(lambda.kcv), col = "black")
  
  lasso.kcv <- glmnet(x = as.matrix(X), y = y, family = "cox", alpha = 1, type.measure = loss.func, lambda = lambda.kcv, penalty.factor = ridge.weight)
  coef.est <- as.numeric(lasso.kcv$beta)
  short.taxa <- taxa.names.rank.out$names[[rank]]
  full.taxa <- names(taxa.out[[rank]])
  cv.err <- lasso.kcv.err
  
  ind <- order(coef.est)
  ord.coef.est <- coef.est[ind]
  ord.short.taxa <- short.taxa[ind]
  ord.full.taxa <- full.taxa[ind]
  sig.ind <- which(ord.coef.est != 0)
  sig.ord.coef.est <- ord.coef.est[sig.ind]
  sig.ord.short.taxa <- ord.short.taxa[sig.ind]
  sig.ord.full.taxa <- ord.full.taxa[sig.ind]
  
  sig.coef.out <- as.data.frame(cbind(sig.ord.short.taxa, sig.ord.full.taxa, sig.ord.coef.est))
  colnames(sig.coef.out) <- c("Taxon", "Full Path", "Est.")
  
  return(list(sig.coef.out = sig.coef.out, cv.err = cv.err, cv.kcv = cv.kcv, lasso.kcv = lasso.kcv.data, lambda.kcv = lambda.kcv))
}

plot.rfsrc <- function (x, m.target = NULL, plots.one.page = TRUE, sorted = TRUE, 
          verbose = TRUE, ...) 
{
  sf.flag <- FALSE
  if (sum(inherits(x, c("rfsrc", "synthetic"), TRUE) == c(1, 
                                                          2)) == 2) {
    if (sum(inherits(x, c("rfsrc", "synthetic", "oob"), 
                     TRUE) == c(1, 2, 3)) != 3) {
      sf.flag <- TRUE
      sf.message <- "OOB was not used for synthetic forests, error rates/VIMP will be unreliable"
    }
    x <- x$rfSyn
  }
  if (sum(inherits(x, c("rfsrc", "grow"), TRUE) == c(1, 2)) != 
      2 & sum(inherits(x, c("rfsrc", "predict"), TRUE) == 
              c(1, 2)) != 2) {
    stop("this function only works for objects of class `(rfsrc, grow)' or '(rfsrc, predict)'")
  }
  dots <- list(...)
  cex <- par("cex")
  if (!is.null(dots$cex)) {
    cex <- dots$cex
  }
  if (sum(inherits(x, c("rfsrc", "subsample"), TRUE) == c(1, 
                                                          4)) == 2 || sum(inherits(x, c("rfsrc", "bootsample"), 
                                                                                   TRUE) == c(1, 4)) == 2) {
    return(plot.subsample(x, m.target = m.target, cex = cex))
  }
  m.target <- get.univariate.target(x, m.target)
  x <- coerce.multivariate(x, m.target)
  if (is.null(x$err.rate)) {
    stop("object is devoid of performance values")
  }
  if (is.null(x$importance)) {
    x$importance <- NA
  }
  if (all(is.na(x$err.rate)) & all(is.na(x$importance))) {
    stop("performance values are all NA")
  }
  if (x$block.size != 1) {
    x$err.rate <- cbind(x$err.rate)
    fill.err.row <- x$block.size
    nullO <- lapply(1:x$ntree, function(i) {
      x$err.rate[i, ] <<- x$err.rate[fill.err.row, ]
      if (i == fill.err.row) {
        fill.err.row <<- min(fill.err.row + x$block.size, 
                             x$ntree)
      }
      NULL
    })
    rm(nullO)
  }
  if (!is.null(x$forest$perf.type) && (x$forest$perf.type == 
                                       "g.mean" || x$forest$perf.type == "g.mean.rfq")) {
    x$err.rate <- x$err.rate[, 1, drop = FALSE]
  }
  if (x$family == "surv-CR" | x$family == "surv") {
    plot.yvar.names <- ""
  }
  else {
    plot.yvar.names <- paste("(", x$yvar.names, ")", sep = "")
  }
  if (all(is.na(x$importance))) {
    if (x$ntree > 1 && !all(is.na(x$err.rate))) {
      err <- cbind(x$err.rate)
      par(mfrow = c(1, 1))
      plot.err(err, plot.yvar.names)
    }
  }
  else {
    err <- cbind(x$err.rate)
    imp <- cbind(x$importance)
    if (!is.null(x$forest$perf.type) && (x$forest$perf.type == 
                                         "g.mean" || x$forest$perf.type == "g.mean.rfq")) {
      imp <- imp[, 1, drop = FALSE]
    }
    x.var.names <- rownames(imp)
    n.pred <- nrow(imp)
    if (sorted) 
      pred.order <- order(imp[, 1])
    else pred.order <- n.pred:1
    if (ncol(imp) == 1) 
      max.pred <- 100
    else max.pred <- 80/ncol(imp)
    if (n.pred > max.pred) {
      dotchart.labels <- rep("", n.pred)
      pretty.pt <- pretty(1:n.pred, n = max.pred)
      dotchart.labels[pretty.pt] <- x.var.names[pred.order][pretty.pt]
    }
    else {
      dotchart.labels <- x.var.names[pred.order]
    }
    if (x$ntree > 1 & !all(is.na(x$err.rate)) & plots.one.page) {
      #par(mfrow = c(2, 1))
      par(mfrow = c(1, 1))
    }
    else {
      par(mfrow = c(1, 1))
    }
    # if (x$ntree > 1 & !all(is.na(x$err.rate))) {
    #   plot.err(err, plot.yvar.names)
    # }
    if (ncol(imp) > 1) {
      imp.out <- imp[rev(pred.order), , drop = FALSE]
      dotChart(imp[pred.order, , drop = FALSE], plot.yvar.names, 
               dotchart.labels, cex = cex)
    }
    if (ncol(imp) == 1) {
      dotChart(imp[pred.order, ], plot.yvar.names, dotchart.labels, 
               cex = cex)
      if (!is.null(x$xvar.wt) & length(unique(x$xvar.wt)) > 
          1) {
        if (length(unique(x$xvar.wt)) == 1) 
          x$xvar.wt <- 1
        imp.out <- as.data.frame(cbind(imp, imp/max(abs(imp), 
                                                    na.rm = TRUE), x$xvar.wt), row.names = x.var.names)[rev(pred.order), 
                                                    ]
        if (nrow(imp.out) == 1) 
          imp.out[1, 2] <- 1
        colnames(imp.out) <- c("Importance", "Relative Imp", 
                               "xvar weight")
      }
      else {
        imp.out = as.data.frame(cbind(imp, imp/max(abs(imp), 
                                                   na.rm = TRUE)), row.names = x.var.names)[rev(pred.order), 
                                                   ]
        if (nrow(imp.out) == 1) 
          imp.out[1, 2] <- 1
        colnames(imp.out) <- c("Importance", "Relative Imp")
      }
    }
    cat("\n")
    if (verbose) {
      print(round(imp.out[1:min(n.pred, max.pred), , drop = FALSE], 
                  4), justify = "right", print.gap = 3)
    }
    # if (x$ntree > 1 & !all(is.na(x$err.rate))) {
    #   plot.err(err, plot.yvar.names)
    # }
  }
  if (sf.flag) {
    message(sf.message)
  }
}

get.univariate.target <- function(x, outcome.target = NULL) {
  ## This function takes a grow, grow-equivalent, or predict object and returns a single coherent target.
  ## That is, if no target has been specified, the first regression outcome with statistics is chosen.
  ## If no regression outcome exists, the first classification outcome with statistics is chosen.
  ## If the target is specified, the object is verified to contain the target outcome statistics
  ## for that y-var.  If none exist, the function will error.
  if (x$family == "regr+" | x$family == "class+" | x$family == "mix+") {
    if (is.null(outcome.target)) {
      ## Check the y-vars against regression and then classification.
      ## We choose the "first" variable, favoring regression, then
      ## classification.
      target <- match(c("regrOutput", "classOutput"), names(x))
      target <- target[!is.na(target)]
      if(length(target) > 0) {
        do.break <- FALSE
        for (i in target) {
          for (j in 1:length(x[[i]])) {
            if (length(x[[i]][[j]]) > 0) {
              ## This is a non-null output.
              outcome.target <- names(x[[i]][j])
              ## Exit the loop.
              do.break <- TRUE
              break
            }
          }
          if (do.break == TRUE) {
            break
          }
        }
      }
      else {
        ## Something would have to be seriously wrong for this to happen.
        stop("No outcomes found in object.  Please contact technical support.")
      }
    }
    else {
      ## Check that one and only one target has been specified.
      if (sum(is.element(outcome.target, x$yvar.names)) != 1) {
        stop("Specified target was not found or too many target outcomes were supplied (only one is allowed).")
      }
      ## A target outcome has been specified.  Verify that it contains outcome statistics.
      target <- match(c("regrOutput", "classOutput"), names(x))
      target <- target[!is.na(target)]
      found = FALSE
      if(length(target) > 0) {
        do.break <- FALSE
        for (i in target) {
          for (j in 1:length(x[[i]])) {
            if (length(x[[i]][[j]]) > 0) {
              ## This is a non-null output.
              if (outcome.target == names(x[[i]][j])) {
                found = TRUE
                ## Exit the loop.
                do.break <- TRUE
                break
              }
            }
          }
          if (do.break == TRUE) {
            break
          }
        }     
      }
      if (!found) {
        stop("Target outcome has been correctly specified but it did not contain outcome statistics.")
      }
    }
  }
  ## This function will return NULL if the function is not
  ## multivariate.  Otherwise, the outcome and its associated statistics is
  ## guaranteed to exist in the object.
  outcome.target
}

coerce.multivariate <- function(x, outcome.target) {
  ## Warning:  This functon assumes that get.univariate.target has been called first, to
  ## verify the coherency of the target.  This means that the target exists in the forest object, and that
  ## it contains outcome statistics.
  ## If this is a multivarate family, we coerce the object, based on outcome.target
  ## into a univaritate regression or classification object.
  x$univariate <- TRUE
  if (x$family == "regr+" | x$family == "class+" | x$family == "mix+") {
    ## coerce the mulitvariate object into a univariate object
    x.coerced <- unlist(list(x$classOutput, x$regrOutput), recursive = FALSE)[[outcome.target]]
    x$univariate <- FALSE
    x$yvar <- x$yvar[, outcome.target]
    ## test for factors - ordered factors are treated as factors!
    if (is.factor(x$yvar) || is.ordered(x$yvar)) {
      x$family <- "class"
    }
    else {
      x$family <- "regr"
    }
    ## make various assignments to the coerced object.
    x$predicted <- x.coerced$predicted
    x$predicted.oob <- x.coerced$predicted.oob
    x$class <- x.coerced$class
    x$class.oob <- x.coerced$class.oob
    x$err.rate <- x.coerced$err.rate
    x$err.block.rate <- x.coerced$err.block.rate
    x$importance <- x.coerced$importance
    x$yvar.names <- outcome.target
    x$cse.num <- x.coerced$cse.num
    x$cse.den <- x.coerced$cse.den
    x$csv.num <- x.coerced$csv.num
    x$csv.den <- x.coerced$csv.den
  }
  x$outcome.target <- outcome.target
  x
}

plot.err <- function(err, yname = NULL) {
  matplot(1:nrow(err), err,
          xlab = "Number of Trees",
          ylab = paste("Error rate", yname),
          type = c("p", "l")[1 + 1 * (nrow(err) > 1)], pch = 16, lty = 1, lwd = 3)
  if (ncol(err) > 1) {
    legend("topright",
           legend = colnames(err), col = 1:ncol(err), lty = 1, lwd = 3)
  }
}

dotChart <- function(x, yname = NULL, labels = NULL, cex = cex) {
  if (!is.null(dim(x))) {
    ncol  <- ncol(x)
    x.dot <- NULL
    for (k in ncol(x):1) {x.dot <- c(x.dot, x[, k])}
    gcolor <- 1:ncol
  }
  else {
    x.dot <- x
    gcolor <- par("fg")
  }
  
  y.dot <- dot.chart.main(x, labels = labels, xlab = paste("Variable Importance", yname),
                          cex = cex, pch="", lwd = 2, lcolor = "white", gcolor = gcolor)
  segments(rep(max(0, min(x.dot, na.rm = TRUE)) - 1e-6, length(y.dot)),
           y.dot, x.dot, y.dot, col=c(4,2)[1 + 1 * (x.dot > 0)], lwd = 4)
  if (min(x.dot, na.rm = TRUE) < 0) abline(v=0, lwd = 2, lty = 2, col = 1)
}
## workhorse for dotchart
dot.chart.main <- function (x, labels = NULL, groups = NULL, gdata = NULL, cex = NULL,
                            pch = 21, gpch = 21, bg = par("bg"), color = par("fg"), gcolor = par("fg"), 
                            lcolor = "gray", xlim = range(x[is.finite(x)]), main = NULL, 
                            xlab = NULL, ylab = NULL, ...) 
{
  opar <- par("mai", "mar", "cex", "yaxs")
  on.exit(par(opar))
  par(yaxs = "i")
  if (!is.numeric(x)) 
    stop("'x' must be a numeric vector or matrix")
  n <- length(x)
  if (is.matrix(x)) {
    if (is.null(labels)) 
      labels <- rownames(x)
    if (is.null(labels)) 
      labels <- as.character(1:nrow(x))
    labels <- rep(labels, length.out = n)
    if (is.null(groups)) 
      groups <- col(x, as.factor = TRUE)
    glabels <- levels(groups)
  }
  else {
    if (is.null(labels)) 
      labels <- names(x)
    glabels <- if (!is.null(groups)) 
      levels(groups)
  }
  plot.new()
  linch <- if (!is.null(labels)) 
    max(strwidth(labels, "inch"), na.rm = TRUE)
  else 0
  if (is.null(glabels)) {
    ginch <- 0
    goffset <- 0
  }
  else {
    ginch <- max(strwidth(glabels, "inch"), na.rm = TRUE)
    goffset <- 0.4
  }
  if (!(is.null(labels) && is.null(glabels))) {
    nmai <- par("mai")
    nmai[2] <- nmai[4] + max(linch + goffset, ginch) + 0.1
    par(mai = nmai)
  }
  if (is.null(groups)) {
    o <- 1:n
    y <- o
    ylim <- c(0, n + 1)
  }
  else {
    o <- sort.list(as.numeric(groups), decreasing = TRUE)
    x <- x[o]
    groups <- groups[o]
    color <- rep(color, length.out = length(groups))[o]
    lcolor <- rep(lcolor, length.out = length(groups))[o]
    offset <- cumsum(c(0, diff(as.numeric(groups)) != 0))
    y <- 1:n + 2 * offset
    ylim <- range(0, y + 2)
  }
  plot.window(xlim = xlim, ylim = ylim, log = "")
  lheight <- par("csi")
  if (!is.null(labels)) {
    linch <- max(strwidth(labels, "inch"), na.rm = TRUE)
    loffset <- (linch + 0.1)/lheight
    labs <- labels[o]
    mtext(labs, side = 2, line = loffset, at = y, adj = 0, 
          col = color, las = 2, cex = cex, ...)
  }
  abline(h = y, lty = "dotted", col = lcolor)
  points(x, y, pch = pch, col = color, bg = bg)
  if (!is.null(groups)) {
    gpos <- rev(cumsum(rev(tapply(groups, groups, length)) + 
                         2) - 1)
    ginch <- max(strwidth(glabels, "inch"), na.rm = TRUE)
    goffset <- (max(linch + 0.2, ginch, na.rm = TRUE) + 0.1)/lheight
    mtext(glabels, side = 2, line = goffset, at = gpos, adj = 0, 
          col = gcolor, las = 2, cex = cex, ...)
    if (!is.null(gdata)) {
      abline(h = gpos, lty = "dotted")
      points(gdata, gpos, pch = gpch, col = gcolor, bg = bg, 
             ...)
    }
  }
  axis(1)
  box()
  title(main = main, xlab = xlab, ylab = ylab, ...)
  invisible(y)
}
