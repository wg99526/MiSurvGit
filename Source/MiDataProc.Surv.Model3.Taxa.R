
taxa.cox.test <- function(survtime, censoring, taxa, cov.dat = NULL, multi.method = "BH") {
  
  censoring <- censoring
  cens.var.ord <- unique(censoring[[1]])[order(unique(censoring[[1]]))]
  print(cens.var.ord)
  
  if( !is.numeric(cens.var.ord) ){
    censoring[ censoring == cens.var.ord[1] ] <- 0
    censoring[ censoring == cens.var.ord[2] ] <- 1
    censoring <- as.numeric(censoring)
  }
  else if( !(cens.var.ord[1] == 0 || cens.var.ord[1] == 1)  ) {
    if( !(cens.var.ord[2] == 0 || cens.var.ord[2] == 1)  ){
      censoring[ censoring == cens.var.ord[1] ] <- 0
      censoring[ censoring == cens.var.ord[2] ] <- 1
      censoring <- as.numeric(censoring)
    }
  }
  
  cox.test <- list()
  for(i in 1:6) {
    
    taxa[[i]] <- taxa[[i]][!is.na(match(rownames(taxa[[i]]), rownames(survtime))),]
    survtime1 <- survtime[match(rownames(survtime), rownames(taxa[[i]]))]
    censoring1 <- censoring[match(rownames(survtime), rownames(taxa[[i]]))]
    
    n.tax <- ncol(taxa[[i]])
    cox.out <- matrix(NA, n.tax, 6)
    for (j in 1:n.tax) {
      taxon <- taxa[[i]][,j]
      if (is.null(cov.dat)){
        fit <- coxph(Surv(survtime1[[1]], censoring1[[1]]) ~  scale(taxon))#, silent = TRUE)
      } else {
        dat <- as.data.frame(cbind(survtime, censoring[[1]], taxon, cov.dat))
        f <- formula(paste("Surv(survtime1[[1]], censoring1[[1]]) ~  scale(taxon) +",  paste(colnames(dat)[4:ncol(dat)], collapse = "+")))
        fit <- try(coxph(f, data = dat), silent = TRUE)
      }
      
      
      if (class(fit) == "try-error"){
        cox.out[j,] <- c(rep(NA, 6))
      } else {
        fit.info <- summary(fit)
        
        cox.out[j,] <- c( round(fit.info$coefficients[1,1], 3), 
                          round(fit.info$coefficients[1,3], 3), 
                          round(fit.info$conf.int[1,1]    , 3),            #round(exp(fit$coefficients), 3)[1],
                          round(fit.info$conf.int[1,3]    , 3),
                          round(fit.info$conf.int[1,4]    , 3),
                          round(fit.info$coefficients[1,5], 3))
      }
    }
    
    cox.out <- as.data.frame(cox.out)
    rownames(cox.out) <- colnames(taxa[[i]])
    colnames(cox.out) <- c("Coef", "SE", "HR", "Lower", "Upper", "P.value")
    cox.out <- bin.q.func(cox.out, method = multi.method)
    cox.test[[i]] <- cox.out
  }
  names(cox.test) <- names(taxa)
  return(cox.test)
}

surv.taxa.forest.plot.pages <- function(all.taxa.q.out, species.include, mult.test.cor = "TRUE") {
  sig.by.rank <- list()
  
  if(species.include){
    range = 6
  }else{
    range = 5
  }
  
  for(i in 1:range) {
    out <- all.taxa.q.out[[i]]
    
    if (mult.test.cor) {
      ind.sig <- which(out$Q.value < 0.05)  
    } else {
      ind.sig <- which(out$P.value < 0.05)
    }
    sig.by.rank[[i]] <- ind.sig
  }
  
  total <- length(unlist(sig.by.rank))
  
  if(total>80) {
    plot.per.page <- total
    num.pages <- 0
    while(plot.per.page > 80){
      num.pages <- num.pages +1
      plot.per.page <- total%/%num.pages
    }
  } else if(total<=80 & total>=0) {
    num.pages = 1
    plot.per.page <-total
  }
  
  return(num.pages)
}

surv.taxa.forest.plot.pages1 <- function(all.taxa.q.out, taxa.names.out, species.include, report.type = "Coef", mult.test.cor = "TRUE") {
  sig.by.rank <- list()
  
  if(species.include){
    range = 6
  }else{
    range = 5
  }
  
  if(report.type == "Est"){
    report.txt <- "Est."
  } else if(report.type == "OR") {
    report.txt <- "OR"
  } else if(report.type == "HR") {
    report.txt <- c("Coef", "HR")
  }
  
  for(i in 1:range) {
    out <- all.taxa.q.out[[i]]
    
    if (mult.test.cor) {
      ind.sig <- which(out$Q.value < 0.05)  
    } else {
      ind.sig <- which(out$P.value < 0.05)
    }
    
    sig.by.rank[[i]] <- ind.sig
  }
  
  total <- length(unlist(sig.by.rank))
  
  if(total>80) {
    capacity <- c(40:80)
    plot.per.page <- total
    num.pages <- 0
    num.mod <- 0
    while(plot.per.page > 80){
      num.pages <- num.pages +1
      plot.per.page <- total%/%num.pages
      num.mod <- mod(total, num.pages)
    }
  } else if(total<=80 & total>=1) {
    num.pages = 1
    plot.per.page <-total
    num.mod <- 0
  }
  
  all.text.tab <- list()
  all.ci.tab <- list()
  if(total >0) {
    plot.taxa <- list()
    rank <- 1
    tax.rank <- c("Phylum","Class","Order","Family","Genus","Species")
    current.page <- 1
    text.tab.all <- 0
    text.tab.save <- numeric()
    ci.tab.save <- numeric()
    info.taxa <- numeric()
    info.data <- numeric()
    info.rank <- numeric()
    info.ci <- numeric()
    sig.out <- numeric()
    actual.plot.per.page <- plot.per.page
    
    if (text.tab.all == 0) {
      if (mult.test.cor){
        text.tab.all <- matrix(c("ID", "Rank", "Taxon", report.txt, "P-value", "Q-value"), 1, (5+length(report.txt)))  
      } else {
        text.tab.all <- matrix(c("ID", "Rank", "Taxon", report.txt, "P-value"), 1, (4+length(report.txt)))
      }
      ci.tab.all <- matrix(c(NA, NA, NA), 1, 3)
    }
    
    for (rank in 1:range) {
      if(!is.null(all.taxa.q.out[[rank]])){
        ind.sig <- sig.by.rank[[rank]]
        sig.out <- all.taxa.q.out[[rank]][sig.by.rank[[rank]],]
        
        info.taxa <- c(info.taxa, taxa.names.out$names[[rank]][sig.by.rank[[rank]]])
        info.data <- rbind(info.data, sig.out)
        info.ci <- rbind(info.ci, cbind(sig.out[,report.type], sig.out[,"Lower"], sig.out[,"Upper"]))
        info.rank <- c(info.rank, rep(tax.rank[rank], length(sig.by.rank[[rank]])))
      }
    }
    
    if(mult.test.cor){
      info <- as.matrix(cbind(c(1:nrow(info.data)), info.rank, info.taxa, format(round(info.data[,c(report.txt, "P.value", "Q.value")],3), nsmall = 3)))
    } else {
      info <- as.matrix(cbind(c(1:nrow(info.data)), info.rank, info.taxa, format(round(info.data[,c(report.txt, "P.value")],3), nsmall = 3)))
    }
    
    
    for(p in 1:num.pages) {
      if(p <= num.mod){
        actual.plot.per.page <- plot.per.page+1
        initial <- 1+actual.plot.per.page*(p-1)
      } else {
        initial <- 1+actual.plot.per.page*(p-1)
        actual.plot.per.page <- plot.per.page
      }
      names(info) <- names(text.tab.all)
      all.text.tab[[p]] <- rbind(as.matrix(text.tab.all), info[c(initial:(initial+actual.plot.per.page-1)),])
      all.ci.tab[[p]] <- rbind(ci.tab.all, info.ci[c(initial:(initial+actual.plot.per.page-1)),])
    }
  } else {
    all.text.tab <- NULL
    all.ci.tab <- NULL
  }
  return(list(all.text.tab = all.text.tab, all.ci.tab = all.ci.tab))
}


taxa.surv.forest.plot.pages2 <- function(page.taxa.q.out, page) {
  
  text.tab.all <- page.taxa.q.out$all.text.tab[[page]]
  ci.tab.all <- page.taxa.q.out$all.ci.tab[[page]]
  
  if(is.null(text.tab.all) & is.null(ci.tab.all)){
    plot.new()
    text(x = 0.5, y = 0.5, "No significant taxa are found.", 
         cex = 1.2, col = "black")
  }else{
    str.max <- lapply(page.taxa.q.out$all.text.tab, nchar)
    for(i in 1:length(page.taxa.q.out$all.text.tab)){
      str.max[[i]] <- str.max[[i]][,3]
    }
    maxStr <- max(unlist(str.max))
    if(!is.numeric(maxStr)){
      maxStr <- 0
    }
    
    text.tab.all[,3] <- substr(text.tab.all[,3], 1, 55)
    text.tab.all[,6] <- p.value.0.1_char(text.tab.all[,6]) 
    text.tab.all[,7] <- p.value.0.1_char(text.tab.all[,7])
    
    #par(mar=c(0, 0.2, 0, 0.2))
    if(text.tab.all[1,5] == "Coef"){
      if(nrow(ci.tab.all) <= 5) {
        forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3],
                   hrzl_lines=TRUE, new_page=TRUE, boxsize=0.12, grid=0, colgap = unit(1.2, "cm"), graphwidth = unit(12, "cm"), lineheight = unit(1, "cm"), #line.margin = unit(0.2, "cm"),
                   col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval for HR", mar = unit(c(0.5,0,0.5,0), "cm"), #mar = unit(c(blank.space,0,0,0), "cm"),
                   txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
                                  ticks=gpar(fontfamily="", cex=0.75),
                                  xlab=gpar(fontfamily="", cex=0.75)))
      }
      else if(nrow(ci.tab.all) <= 45){
        forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3],
                   hrzl_lines=TRUE, new_page=TRUE, boxsize=0.12, grid=0, colgap = unit(1.2, "cm"), graphwidth = unit(12, "cm"), lineheight = "lines", line.margin = unit(0.12, "cm"),
                   col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval for HR", mar = unit(c(0.5,0,0.5,0), "cm"), #mar = unit(c(blank.space,0,0,0), "cm"),
                   txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
                                  ticks=gpar(fontfamily="", cex=0.75),
                                  xlab=gpar(fontfamily="", cex=0.75)))
      } else {
        forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3],
                   hrzl_lines=TRUE, new_page=TRUE, boxsize=0.12, grid=0, colgap = unit(1.2, "cm"), graphwidth = unit(12, "cm"), 
                   col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval for HR", mar = unit(c(0.5,0,0.5,0), "cm"), #mar = unit(c(blank.space,0,0,0), "cm"),
                   txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
                                  ticks=gpar(fontfamily="", cex=0.75),
                                  xlab=gpar(fontfamily="", cex=0.75)))
      }
      
    } else {
      if(nrow(ci.tab.all) <= 5){
        forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3],
                   zero = 1, hrzl_lines=TRUE, new_page=TRUE, boxsize=0.12, grid=0, colgap = unit(1.2, "cm"), graphwidth = unit(12, "cm"), lineheight = unit(1, "cm"), #line.margin = unit(0.08, "cm"),
                   col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval for HR", mar = unit(c(0.5,0,0.5,0), "cm"), #mar = unit(c(blank.space,0,0,0), "cm"),
                   txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
                                  ticks=gpar(fontfamily="", cex=0.75),
                                  xlab=gpar(fontfamily="", cex=0.75)))
      } else if(nrow(ci.tab.all) <= 45){
        forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3],
                   zero = 1, hrzl_lines=TRUE, new_page=TRUE, boxsize=0.12, grid=0, colgap = unit(1.2, "cm"), graphwidth = unit(12, "cm"), lineheight = "lines", line.margin = unit(0.12, "cm"),
                   col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval for HR", mar = unit(c(0.5,0,0.5,0), "cm"), #mar = unit(c(blank.space,0,0,0), "cm"),
                   txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
                                  ticks=gpar(fontfamily="", cex=0.75),
                                  xlab=gpar(fontfamily="", cex=0.75)))
      } else {
        forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3],
                   zero = 1, hrzl_lines=TRUE, new_page=TRUE, boxsize=0.12, grid=0, colgap = unit(1.2, "cm"), graphwidth = unit(12, "cm"), 
                   col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval for HR", mar = unit(c(0.5,0,0.5,0), "cm"), #mar = unit(c(blank.space,0,0,0), "cm"),
                   txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
                                  ticks=gpar(fontfamily="", cex=0.75),
                                  xlab=gpar(fontfamily="", cex=0.75)))
      }
    } 
  }
}


# 
# 
# taxa.surv.forest.plot.pages2 <- function(page.taxa.q.out, page) {
#   
#   text.tab.all <- page.taxa.q.out$all.text.tab[[page]]
#   ci.tab.all <- page.taxa.q.out$all.ci.tab[[page]]
#   
#   if(is.null(text.tab.all) & is.null(ci.tab.all)){
#     plot.new()
#     text(x = 0.5, y = 0.5, "No significant taxa are found.", 
#          cex = 1.2, col = "black")
#   }else{
#     str.max <- lapply(page.taxa.q.out$all.text.tab, nchar)
#     for(i in 1:length(page.taxa.q.out$all.text.tab)){
#       str.max[[i]] <- str.max[[i]][,3]
#     }
#     maxStr <- max(unlist(str.max))
#     if(!is.numeric(maxStr)){
#       maxStr <- 0
#     }
#     
#     text.tab.all[,3] <- substr(text.tab.all[,3], 1, 55)
# 
#     if(text.tab.all[1,5] == "Coef"){
#       if(nrow(ci.tab.all) <= 5) {
#         forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3],
#                    hrzl_lines=TRUE, new_page=TRUE, boxsize=0.25, grid=0, colgap = unit(1, "cm"), graphwidth = unit(6, "cm"), lineheight = unit(1, "cm"), #line.margin = unit(0.2, "cm"),
#                    col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval for HR", mar = unit(c(0.5,0,0.5,0), "cm"), #mar = unit(c(blank.space,0,0,0), "cm"),
#                    txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
#                                   ticks=gpar(fontfamily="", cex=0.75),
#                                   xlab=gpar(fontfamily="", cex=0.75)))
#       }
#       else if(nrow(ci.tab.all) <= 45){
#         forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3],
#                    hrzl_lines=TRUE, new_page=TRUE, boxsize=0.25, grid=0, colgap = unit(1, "cm"), graphwidth = unit(6, "cm"), lineheight = "lines", line.margin = unit(0.12, "cm"),
#                    col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval for HR", mar = unit(c(0.5,0,0.5,0), "cm"), #mar = unit(c(blank.space,0,0,0), "cm"),
#                    txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
#                                   ticks=gpar(fontfamily="", cex=0.75),
#                                   xlab=gpar(fontfamily="", cex=0.75)))
#       } else {
#         forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3],
#                    hrzl_lines=TRUE, new_page=TRUE, boxsize=0.25, grid=0, colgap = unit(1, "cm"), graphwidth = unit(6, "cm"), 
#                    col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval for HR", mar = unit(c(0.5,0,0.5,0), "cm"), #mar = unit(c(blank.space,0,0,0), "cm"),
#                    txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
#                                   ticks=gpar(fontfamily="", cex=0.75),
#                                   xlab=gpar(fontfamily="", cex=0.75)))
#       }
#       
#     } else {
#       if(nrow(ci.tab.all) <= 5){
#         forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3],
#                    zero = 1, hrzl_lines=TRUE, new_page=TRUE, boxsize=0.25, grid=0, colgap = unit(1, "cm"), graphwidth = unit(6, "cm"), lineheight = unit(1, "cm"), #line.margin = unit(0.08, "cm"),
#                    col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval for HR", mar = unit(c(0.5,0,0.5,0), "cm"), #mar = unit(c(blank.space,0,0,0), "cm"),
#                    txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
#                                   ticks=gpar(fontfamily="", cex=0.75),
#                                   xlab=gpar(fontfamily="", cex=0.75)))
#       } else if(nrow(ci.tab.all) <= 45){
#         forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3],
#                    zero = 1, hrzl_lines=TRUE, new_page=TRUE, boxsize=0.25, grid=0, colgap = unit(1, "cm"), graphwidth = unit(6, "cm"), lineheight = "lines", line.margin = unit(0.12, "cm"),
#                    col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval for HR", mar = unit(c(0.5,0,0.5,0), "cm"), #mar = unit(c(blank.space,0,0,0), "cm"),
#                    txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
#                                   ticks=gpar(fontfamily="", cex=0.75),
#                                   xlab=gpar(fontfamily="", cex=0.75)))
#       } else {
#         forestplot(labeltext=text.tab.all, mean=ci.tab.all[,1], lower=ci.tab.all[,2], upper=ci.tab.all[,3],
#                    zero = 1, hrzl_lines=TRUE, new_page=TRUE, boxsize=0.25, grid=0, colgap = unit(1, "cm"), graphwidth = unit(6, "cm"), 
#                    col=fpColors(box=rgb(1,0,0,0.5), line="black", summary="red3"), xlab="95% Confidence Interval for HR", mar = unit(c(0.5,0,0.5,0), "cm"), #mar = unit(c(blank.space,0,0,0), "cm"),
#                    txt_gp=fpTxtGp(label=list(gpar(fontfamily="", cex=0.75), gpar(fontfamily="", cex=0.75)),
#                                   ticks=gpar(fontfamily="", cex=0.75),
#                                   xlab=gpar(fontfamily="", cex=0.75)))
#       }
#     } 
#   }
# }
# 
# 
