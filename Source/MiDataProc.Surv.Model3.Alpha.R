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

###################
# Alpha diversity #
###################

alpha.pe.pqe.func <- function(x, tree, norm = TRUE) {
  ind <- which(x != 0)
  s.tree <- prune_taxa(names(x[ind]), tree)
  pe <- AllenH(x[ind], 1, s.tree, Normalize = norm, CheckArguments = FALSE)
  pqe <- AllenH(x[ind], 2, s.tree, Normalize = norm, CheckArguments = FALSE)
  return(c(pe, pqe))
}

alpha.v1.func <- function(biom) {
  
  otu.tab <- otu_table(biom)
  tree <- phy_tree(biom)
  biom <- merge_phyloseq(round(otu.tab), tree)
  prop.otu.tab <- apply(otu.tab, 2, function(x) x/sum(x))
  t.prop.otu.tab <- t(prop.otu.tab)
  t.otu.tab <- t(otu.tab)
  
  alpha.abu <- estimate_richness(biom, split = TRUE)
  alpha.ice <- apply(otu.tab, 2, ICE)
  alpha.pd <- pd(t.otu.tab, tree)$PD
  
  alpha.div <- as.data.frame(cbind(alpha.abu[,c("Observed", "Shannon", "Simpson", "InvSimpson", "Fisher", "Chao1", "ACE")], alpha.ice, alpha.pd))
  colnames(alpha.div) <- c("Observed", "Shannon", "Simpson", "InvSimpson", "Fisher", "Chao1", "ACE", "ICE", "PD")
  
  rownames(alpha.div) <- colnames(otu.tab)
  
  return(alpha.div)
}

alpha.v2.func <- function(biom) {
  
  otu.tab <- otu_table(biom)
  tree <- phy_tree(biom)
  biom <- merge_phyloseq(round(otu.tab), tree)
  prop.otu.tab <- apply(otu.tab, 2, function(x) x/sum(x))
  t.prop.otu.tab <- t(prop.otu.tab)
  t.otu.tab <- t(otu.tab)
  
  alpha.abu <- estimate_richness(biom, split = TRUE)
  alpha.ice <- apply(otu.tab, 2, ICE)
  alpha.pd <- pd(t.otu.tab, tree)$PD
  alpha.pe.pqe <- t(apply(t.prop.otu.tab, 1, function(x) alpha.pe.pqe.func(x, tree)))
  
  alpha.div <- as.data.frame(cbind(alpha.abu[,c("Observed", "Shannon", "InvSimpson", "Chao1", "ACE")], alpha.ice, alpha.pd, alpha.pe.pqe))
  colnames(alpha.div) <- c("Observed", "Shannon", "InvSimpson", "Chao1", "ACE", "ICE", "PD", "PE", "PQE")
  
  rownames(alpha.div) <- colnames(otu.tab)
  return(alpha.div)
}

############################
# Alpha diversity Survival analysis #
############################

alpha.surv.events.dat.func <- function(sam.dat, surv.con, surv.second.bin) {
  
  survtime    <- sam.dat[[surv.con]]
  status     <- sam.dat[[surv.second.bin]]
  
  stat.var.ordered <- unique(status)[order(unique(status))]
  
  if( !is.numeric(stat.var.ordered) ){
    status[ status == stat.var.ordered[1] ] <- 0
    status[ status == stat.var.ordered[2] ] <- 1
    status <- as.numeric(status)
  }
  else if( !(stat.var.ordered[1] == 0 || stat.var.ordered[1] == 1)  ) {
    if( !(stat.var.ordered[2] == 0 || stat.var.ordered[2] == 1)  ){
      status[ status == stat.var.ordered[1] ] <- 0
      status[ status == stat.var.ordered[2] ] <- 1
      status <- as.numeric(status)
    }
  }
  
  events    <- data.frame(survtime = survtime, status = status)
  
  return( events )
}
  
# Data Union
data.union.func <- function(sam_dat, alpha.div){
  X.SampleID <- rownames(alpha.div)
  new_alpha.div <- cbind(X.SampleID, alpha.div)
  new_dat <- merge( data.frame(cbind(X.SampleID, sam_dat)), new_alpha.div, by = "X.SampleID", all.y = TRUE)
  print(head(new_dat))

  new_sam.dat   <- new_dat[ , ! names(new_dat) %in% colnames(alpha.div)  ]
  new_alpha.div <- subset(new_dat, select = colnames(alpha.div)) #new_alpha.div <- new_dat[ ,   names(new_dat) %in% colnames(alpha.div) ]
  rownames(new_alpha.div) <- new_dat$X.SampleID
  
  return(list(new_sam.dat = new_sam.dat, new_alpha.div = new_alpha.div) )
}


# Table for Forest Plot
alpha.cox.test <- function(new_surv.dat, new_alpha.div, cov.dat = NULL) {
  n.alpha <- ncol(new_alpha.div)
  out <- matrix(NA, n.alpha, 6)
  
  if ( is.null(cov.dat) ){
    
    for (i in 1:n.alpha) {
      fit <- coxph(Surv(survtime, status) ~  scale(new_alpha.div[,i]), data = new_surv.dat)
      fit.info <- summary(fit)
      out[i,] <- c( round(fit.info$coefficients[1,1], 3), 
                    round(fit.info$coefficients[1,3], 3), 
                    round(fit.info$coefficients[1,5], 3),
                    round(fit.info$conf.int[1,1]    , 3),            #round(exp(fit$coefficients), 3)[1],
                    round(fit.info$conf.int[1,3]    , 3),
                    round(fit.info$conf.int[1,4]    , 3)
      )
    }
  }
  else{
    
    cov.list <- paste( names(cov.dat), collapse = ' + ') #paste(' + ', paste(input, collapse = ' + ') )

    new_dat <- data.frame( cbind(new_surv.dat, cov.dat) )
    
    form <- as.formula( paste("Surv(survtime, status) ~ scale(new_alpha.div[,i]) + ", cov.list ) )
    
    for (i in 1:n.alpha) {
      fit <- coxph( form, data = new_dat)
      fit.info <- summary(fit)
      out[i,] <- c( round(fit.info$coefficients[1,1], 3), 
                    round(fit.info$coefficients[1,3], 3), 
                    round(fit.info$coefficients[1,5], 3),
                    round(fit.info$conf.int[1,1]    , 3),            #round(exp(fit$coefficients), 3)[1],
                    round(fit.info$conf.int[1,3]    , 3),
                    round(fit.info$conf.int[1,4]    , 3)
      )
    }
  }
  
  out <- as.data.frame(out)
  rownames(out) <- colnames(new_alpha.div)
  colnames(out) <- c("Coef", "SE", "P.value", "HR", "Lower", "Upper")
  
  return(out)
}


# Forest Plot
alpha.cox.forest.plot <- function(out, multi.test = FALSE) {
  
  if(multi.test) {
    text.tab.all <- as.matrix(rbind(c("Alpha Diversity", "Coef", "SE", "HR", "P.value", "Q.value"), 
                                    cbind(rownames(out),
                                          format(round(out[, c(1, 2)], digits = 3), nsmall = 3), 
                                          p.value.0.1(out[,4]),
                                          format(round(out[, c(3,7)], digits = 3), nsmall = 3) )))
  } else {
    text.tab.all <- as.matrix(rbind(c("Alpha Diversity", "Coef", "SE", "HR", "P.value"), 
                                    cbind(rownames(out),
                                          format(round(out[, c(1, 2)], digits = 3), nsmall = 3), 
                                          p.value.0.1(out[,4]),
                                          format(round(out[, 3], digits = 3), nsmall = 3) )))
  }
  
    
  ci.tab.all <- as.matrix(rbind(c(NA, NA, NA), cbind(out[,4], out[,c(5,6)])))
    
    
  forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3], 
             
             hrzl_lines=TRUE, new_page=TRUE, boxsize=0.25, line.margin=0.3, grid=1, lineheight = unit(1, "cm"), colgap = unit(0.5, "cm"), graphwidth = unit(6, "cm"),
               
             col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval for HR",
               
             txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.7), gpar(fontfamily="", cex=0.7)),
                            ticks=gpar(fontfamily="", cex=0.7),   
                            xlab=gpar(fontfamily="", cex=0.7)))
    
}




