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
library(lme4)
library(lmerTest)
library(dirmult) 
library(robustbase)
library(robCompositions) 
library(BiasedUrn)
library(CompQuadForm)
library(GUniFrac) 
library(ecodist) 
library(MiRKAT)
#library(GLMMMiRKAT)
library(proxy)
library(vegan3d)
library(pca3d)
#####################
# Data manipulation #
#####################

beta.bin.cat.func <- function(sam.dat, sel.bin.var) {
  bin.cat <- levels(as.factor(unlist(sam.dat[,sel.bin.var])))
  return(bin.cat)
}

beta.bin.cat.ref.ori.func <- function(sam.dat, sel.bin.var = "ecig_status") {
  
  return(levels(as.factor(as.data.frame(as.matrix(sam.dat))[,sel.bin.var])))
}


beta.bin.cat.recode.func <- function(sam.dat, sel.bin.var = "ecig_status", ori.cat, rename.ref, rename.com) {
  ind.ref <- which(sam.dat[,sel.bin.var] == ori.cat[1])
  ind.com <- which(sam.dat[,sel.bin.var] == ori.cat[2])
  sam.dat[ind.ref,sel.bin.var] <- rename.ref
  sam.dat[ind.com,sel.bin.var] <- rename.com
  
  return(sam.dat)
}

betaS.bin.cat.ref.func <- function(sel.bin.var, sel.con.var, sel.ref, sel.com, sam.dat, Ds.Ks) {
  bin.var <- unlist(sam.dat[,sel.bin.var])
  ind.ref <- which(bin.var[1:nrow(Ds.Ks$Ds[[1]])] == sel.ref)
  ind.com <- which(bin.var[1:nrow(Ds.Ks$Ds[[1]])] == sel.com)
  bin.var <- factor(c(bin.var[ind.ref], bin.var[ind.com]), levels = c(sel.ref, sel.com))

  
  #X.SampleID <- rownames(Ds.Ks$Ds$W.UniFrac)
  X.SampleID <- rownames(data.frame(Ds.Ks$Ds[1]))
  con.var <- sam.dat[X.SampleID,][[sel.con.var]]
  
  
  Ds <- Ds.Ks$Ds
  Ks <- Ds.Ks$Ks
  for (i in 1:length(Ds)) {
    Ds[[i]] <- Ds[[i]][c(ind.ref, ind.com), c(ind.ref, ind.com)]
  }
  for (i in 1:length(Ks)) {
    Ks[[i]] <- Ks[[i]][c(ind.ref, ind.com), c(ind.ref, ind.com)]
  }
  
  return(list(bin.var = bin.var, con.var = con.var, Ds = Ds, Ks = Ks))
}

betaS.bin.cov.cat.ref.func <- function(sel.bin.var, sel.con.var, sel.ref, sel.com, sel.cov.var, sam.dat, Ds.Ks) {  
  bin.var <- unlist(sam.dat[,sel.bin.var])
  ind.ref <- which(bin.var[1:nrow(Ds.Ks$Ds[[1]])] == sel.ref)
  ind.com <- which(bin.var[1:nrow(Ds.Ks$Ds[[1]])] == sel.com)
  bin.var <- factor(c(bin.var[ind.ref], bin.var[ind.com]), levels = c(sel.ref, sel.com))
  
  cov.var <- sam.dat[,sel.cov.var]
  cov.var <- as.data.frame(rbind(cov.var[ind.ref,], cov.var[ind.com,]))
  
  
  #X.SampleID <- rownames(Ds.Ks$Ds$W.UniFrac)
  X.SampleID <- rownames(data.frame(Ds.Ks$Ds[1]))
  con.var <- sam.dat[X.SampleID,][[sel.con.var]]
  
  
  Ds <- Ds.Ks$Ds
  Ks <- Ds.Ks$Ks
  for (i in 1:length(Ds)) {
    Ds[[i]] <- Ds[[i]][c(ind.ref, ind.com), c(ind.ref, ind.com)]
  }
  for (i in 1:length(Ks)) {
    Ks[[i]] <- Ks[[i]][c(ind.ref, ind.com), c(ind.ref, ind.com)]
  }
  
  return(list(bin.var = bin.var, con.var = con.var, cov.var = cov.var, Ds = Ds, Ks = Ks))
}



##################
# Beta diversity #
##################

Ds.Ks.func <- function(rare.biom, biom.after.qc) {
  rare.otu.tab <- otu_table(rare.biom)
  no.rare.otu.tab <- otu_table(biom.after.qc)
  no.rare.tree <- phy_tree(biom.after.qc)
  
  jac <- as.matrix(proxy::dist(t(rare.otu.tab), method = "Jaccard"))
  bc <- as.matrix(bcdist(t(rare.otu.tab)))
  unifs <- GUniFrac(t(no.rare.otu.tab ), no.rare.tree, alpha = c(0.5, 1))$unifracs
  u.unif <- unifs[, , "d_UW"]
  g.unif <- unifs[, , "d_0.5"]
  w.unif <- unifs[, , "d_1"]
  
  jac.k <- D2K(jac)
  bc.k <- D2K(bc)
  u.unif.k <- D2K(u.unif)
  g.unif.k <- D2K(g.unif)
  w.unif.k <- D2K(w.unif)
  
  rownames(jac.k) <- colnames(jac.k) <- colnames(rare.otu.tab)
  rownames(bc.k) <- colnames(bc.k) <- colnames(rare.otu.tab)
  rownames(u.unif.k) <- colnames(u.unif.k) <- colnames(rare.otu.tab)
  rownames(g.unif.k) <- colnames(g.unif.k) <- colnames(rare.otu.tab)
  rownames(w.unif.k) <- colnames(w.unif.k) <- colnames(rare.otu.tab)
  
  return(
    list(Ds = list(Jaccard = jac, Bray.Curtis = bc, U.UniFrac = u.unif, G.UniFrac = g.unif, W.UniFrac = w.unif),
         Ks = list(Jaccard = jac.k, Bray.Curtis = bc.k, U.UniFrac = u.unif.k, G.UniFrac = g.unif.k, W.UniFrac = w.unif.k))
  )
}


#################
# Data analysis #
#################

### MiRKAT

mirkatS.bin <- function(beta.bin.out) {
  set.seed(487)
  out <- MiRKATS(obstime = as.numeric(beta.bin.out$con.var)-1,#beta.bin.out$con.var,
                 delta   = as.numeric(beta.bin.out$bin.var)-1,#beta.bin.out$bin.var,
                 X = NULL, 
                 Ks = beta.bin.out$Ks,
                 nperm = 3000)
  
  return(out)
}

# mirkatS.bin.plot <- function(out, beta.bin.out) {
#   par(mfrow = c(3, 2))
#   for (i in 1:length(beta.bin.out$Ds)) {
#     if (out$p_values[i] < 0.05) {
#       sub.tit <- paste("*p:", p.value.0.1(out$p_values[i]), sep="")
#     }
#     if (out$p_values[i] >= 0.05) {
#       sub.tit <- paste("p:", p.value.0.1(out$p_values[i]), sep="")
#     }
#     mod <- betadisper(as.dist(beta.bin.out$Ds[[i]]), beta.bin.out$bin.var)
#     plot(mod, ellipse = TRUE, hull = FALSE, main = names(beta.bin.out$Ds)[i], xlab="PC 1", ylab="PC 2",
#         sub=NA, col = c("blue2", "red2"), mgp=c(2.5,1,0), cex=1.7, label.cex=1.3, cex.lab=1.2, cex.main=1.7)
#     mtext(sub.tit, side=1, line=3.8, cex=1.0)
#   }
#   
#   plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
#   legend("center", title = NULL, legend = levels(beta.bin.out$bin.var), fil = c("blue2", "red2", cex=2.5, box.lty=0), bty = "n", cex=1.6)
#   legend("bottom", paste("Omnibus MiRKAT: ", p.value.0.1(out$omnibus_p), sep=""), bty = "n", cex=1.6)
# }

mirkatS.bin.plot1 <- function(out, beta.bin.out) {
  mod <- list()
  sub.tit <- list()
  
  bin.var <- c() 
  
  bin.var[beta.bin.out$bin.var == "0" | beta.bin.out$bin.var == 0] <- "Censored"
  bin.var[beta.bin.out$bin.var == "1" | beta.bin.out$bin.var == 1] <- "Event"
  
  for (i in 1:length(beta.bin.out$Ds)) {
    if (out$p_values[i] < 0.05) {
      sub.tit[[i]] <- paste("*p:", p.value.0.1(out$p_values[i]), sep="")
    }
    if (out$p_values[i] >= 0.05) {
      sub.tit[[i]] <- paste("p:", p.value.0.1(out$p_values[i]), sep="")
    }
    
    
    mod[[i]] <- betadisper(as.dist(beta.bin.out$Ds[[i]]), bin.var)
  }
  
  return(list(mod = mod, sub.tit = sub.tit))
}

mirkatS.bin.plot2 <- function(out, beta.bin.out, mod, sub.tit) {
  par(mfrow = c(3,2))
  for (i in 1:length(beta.bin.out$Ds)) {
    
    plot(mod[[i]], ellipse = TRUE, hull = FALSE, main = names(beta.bin.out$Ds)[i], xlab="PC 1", ylab="PC 2",
         sub=NA, col = c("blue2", "red2"), mgp=c(2.5,1,0), cex=1.7, label.cex=1.3, cex.lab=1.2, cex.main=1.7)
    mtext(sub.tit[[i]], side=1, line=3.8, cex=1.0)
  }
  
  legends <- levels(beta.bin.out$bin.var)
  legends[legends == "0"] <- "Censored"
  legends[legends == "1"] <- "Event"
  
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  
  legend("center", title = NULL, legend = legends, fil = c("blue2", "red2", cex=2.5, box.lty=0), bty = "n", cex=1.6)
  legend("bottom", paste("OMiRKAT-S: ", p.value.0.1(out$omnibus_p), sep=""), bty = "n", cex=1.6)
}

mirkatS.bin.plot3d <- function(out, beta.bin.out, mod, sub.tit) {
  
  par(mfrow = c(2,2))
  
  plot(mod[[1]], ellipse = TRUE, hull = FALSE, 
       main = paste(names(beta.bin.out$Ds)[1], "2D" ), 
       xlab="PC 1", ylab="PC 2",
       sub=NA, col = c("blue2", "red2") )#, mgp=c(2.5,1,0), cex=1.7, label.cex=1.3, cex.lab=1.2, cex.main=1.7)
  mtext(sub.tit[[1]], side=1, line=3.8, cex=1.0)
  
  
  p2 <- ordiplot3d(mod[[1]],
                   scaling = "symmetric", angle=30,
                   main = paste(names(beta.bin.out$Ds)[1], "3D" ),
                   xlab="PC 1", ylab="PC 2", zlab="PC 3",
                   col = c("blue2", "red2"))
  mtext(sub.tit[[1]], side=1, line=3.8, cex=1.0)
  
  
  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  legend("center", title = NULL, legend = levels(beta.bin.out$bin.var), fil = c("blue2", "red2", cex=2.5, box.lty=0), bty = "n", cex=1.6)
  legend("bottom", paste("Omnibus MiRKAT: ", p.value.0.1(out$omnibus_p), sep=""), bty = "n", cex=1.6)
  
  barplot(mod[[1]]$eig[1:5], main = "Significant PC Values", col = "green2",
          names.arg=c("PC1","PC2","PC3","PC4","PC5") )
  
}

mirkatS.bin.cov <- function(beta.bin.cov.out) {
  set.seed(487)
  #if ( is.numeric( beta.bin.cov.out$  )  )
  
  cov.ls <- names(beta.bin.cov.out$cov.var)
  
  library(dplyr)
  
  for ( i in 1:length(cov.ls)){
    if ( !is.numeric(beta.bin.cov.out$cov.var[[cov.ls[i]]] ) ){
      beta.bin.cov.out$cov.var[[cov.ls[i]]] <- ifelse(beta.bin.cov.out$cov.var[[cov.ls[i]]] == unique(beta.bin.cov.out$cov.var[[cov.ls[i]]])[1], 0, 1)
    }
  }
  
  out <- MiRKATS(#y = as.numeric(beta.bin.cov.out$bin.var)-1, 
                 obstime = as.numeric(beta.bin.cov.out$con.var)-1,#beta.bin.out$con.var,
                 delta   = as.numeric(beta.bin.cov.out$bin.var)-1,#beta.bin.out$bin.var,
                 X = as.vector(beta.bin.cov.out$cov.var),
                 beta = NULL,
                 Ks = beta.bin.cov.out$Ks,
                 nperm = 3000)
  
  return(out)
}

mirkat.bin.cov.plot <- function(out, beta.bin.cov.out) {
  par(mfrow = c(3, 2))
  for (i in 1:length(beta.bin.cov.out$Ds)) {
    if (out$p_values[i] < 0.05) {
      sub.tit <- paste("*p:", p.value.0.1(out$p_values[i]), sep="")
    }
    if (out$p_values[i] >= 0.05) {
      sub.tit <- paste("p:", p.value.0.1(out$p_values[i]), sep="")
    }
    mod <- betadisper(as.dist(beta.bin.cov.out$Ds[[i]]), beta.bin.cov.out$bin.var)
    plot(mod, ellipse = TRUE, hull = FALSE, main = names(beta.bin.cov.out$Ds)[i], xlab="PC 1", ylab="PC 2",
         sub=NA, col = c("blue2", "red2"), mgp=c(2.5,1,0), cex=1.7, label.cex=1.3, cex.lab=1.2, cex.main=1.7)
    mtext(sub.tit, side=1, line=3.8, cex=1.0)
  }

  plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
  legend("center", title = NULL, legend = levels(beta.bin.cov.out$bin.var), fil = c("blue2", "red2", cex=2.5, box.lty=0), bty = "n", cex=1.6)
  legend("bottom", paste("Omnibus MiRKAT: ", p.value.0.1(out$omnibus_p), sep=""), bty = "n", cex=1.6)
}


###################
# Other functions #
###################

q.func <- function(out, method = c("BH", "BY")) {
  Q.value <- p.adjust(out$P.value, method = method)
  
  return(cbind(out, Q.value))
}

p.value.0.1 <- function(x, round.x = 3) {
  x <- format(round(x, 3), digits = 3)
  ind.0 <- which(x == "0.000" | x == 0)
  x[ind.0] <- "<.001"
  ind.1 <- which(x == "1.000" | x == 1)
  x[ind.1] <- ">.999"
  
  return(x)
}

