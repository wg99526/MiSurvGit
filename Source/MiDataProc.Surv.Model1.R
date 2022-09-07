##############################################
##########   Required Libraries     ##########
##############################################

library(survival)
library(survminer)
library(coin)
library(dplyr)

##############################################
##########   Which Variable?        ##########
##############################################

surv.pri.func <- function(sam.dat, mon.sin.rev.bin.con) {
  is.bin <- colnames(sam.dat)[(mon.sin.rev.bin.con$is.bin) & !mon.sin.rev.bin.con$is.mon & !mon.sin.rev.bin.con$is.sin]
  is.con <- colnames(sam.dat)[(mon.sin.rev.bin.con$is.con) & !mon.sin.rev.bin.con$is.mon & !mon.sin.rev.bin.con$is.sin]
  return(list(is.bin = is.bin, is.con = is.con))
}

##############################################
##########   Numerical Values?      ##########
##############################################

surv.num.treat1 <- function(sam.dat, select.var) {
  #count( sam.dat[ sam.dat[[first.binary]] == "Tylosin" ] )$n
  select.var.list <- sam.dat[[select.var]]
  lvl <- unique( select.var.list  )
  return( length( select.var.list[  select.var.list == lvl[1] ]  ) )
}

surv.num.treat2 <- function(sam.dat, select.var) {
  select.var.list <- sam.dat[[select.var]]
  lvl <- unique( select.var.list  )
  return( length( select.var.list[  select.var.list == lvl[2] ]  ) )
}

surv.num.censr <- function(sam.dat) {
  dplyr::count( sam.dat[ sam.dat$status == 0, ])$n
}

surv.num.uncen <- function(sam.dat) {
  dplyr::count( sam.dat[ sam.dat$status == 1, ])$n
}


##############################################
##########   Necessary Data Sets    ##########
##############################################

# For KM and no covariate Cox PH model
surv.events.dat.func <- function(sam.dat, surv.first.bin = NA, surv.con = "T1Dweek", surv.second.bin = "T1D") {
  
  survtime    <- sam.dat[[surv.con]]
  status     <- sam.dat[[surv.second.bin]]
  
  if( !is.numeric(unique(status)) ){
    status[ status == unique(status)[1] ] <- 0
    status[ status == unique(status)[2] ] <- 1
    status <- as.numeric(status)
  }
  else if( !(unique(status)[1] == 0 || unique(status)[1] == 1)  ) {
    if( !(unique(status)[2] == 0 || unique(status)[2] == 1)  ){
      status[ status == unique(status)[1] ] <- 0
      status[ status == unique(status)[2] ] <- 1
      status <- as.numeric(status)
    }
  }
  
  if (is.na(surv.first.bin)) {
    events    <- data.frame(survtime = survtime, status = status)
  } else {
    treat  <- sam.dat[[surv.first.bin]]
    events    <- data.frame(survtime = survtime, status = status, treat = treat)
  }
  rownames(events) <- rownames(sam.dat)
  
  return( events )
}

follow.time.func <- function(surv.dat, follow.time) {
  
  surv.dat$survtime <- surv.dat$survtime - follow.time[1]
  ind1 <- surv.dat$survtime <= 0
  surv.dat <- surv.dat[!ind1,]
  
  ind2 <- surv.dat$survtime >= (follow.time[2]-follow.time[1])
  surv.dat$status[ind2] <- 0
  
  surv.dat$survtime[ind2] <- follow.time[2]-follow.time[1]
  
  return(surv.dat)
}

sam.follow.time.func <- function(sam.dat, time.var, censor.var, follow.time) {
  
  sam.dat[[time.var]] <- sam.dat[[time.var]] - follow.time[1]
  
  ind1 <- sam.dat[[time.var]] <= 0
  sam.dat <- sam.dat[!ind1,]
  
  ind2 <- sam.dat[[time.var]] >= (follow.time[2]-follow.time[1])
  
  if(!is.na(censor.var)) {
    sam.dat[[censor.var]][ind2] <- 0
  }
  
  sam.dat[[time.var]][ind2] <- follow.time[2]-follow.time[1]
  
  return( sam.dat )
}

# For covariate Cox PH model
surv.events.dat.func2 <- function(sam.dat,
                                  surv.first.bin = "Antibiotics", 
                                  surv.con = "T1Dweek", 
                                  surv.second.bin = "T1D",
                                  input) {
  
  treat  <- sam.dat[[surv.first.bin]]
  survtime    <- sam.dat[[surv.con]]
  status     <- sam.dat[[surv.second.bin]]
  
  if( !is.numeric(unique(status)) ){
    status[ status == unique(status)[1] ] <- 0
    status[ status == unique(status)[2] ] <- 1
    status <- as.numeric(status)
  }
  else if( !(unique(status)[1] == 0 || unique(status)[1] == 1)  ) {
    if( !(unique(status)[2] == 0 || unique(status)[2] == 1)  ){
      status[ status == unique(status)[1] ] <- 0
      status[ status == unique(status)[2] ] <- 1
      status <- as.numeric(status)
    }
  }
  
  events    <- data.frame(survtime = survtime, status = status, treat = treat, subset( sam.dat, select = input))
  
  ind <- complete.cases(events)
  events <- events[ind,]
  
  return( events )
}


##############################################
##############   Methods     #################
##############################################

# Log-rank Test
surv.logrank.test.func <- function( surv.dat ) {
  surv.logrank.test <- survdiff( Surv( survtime, status) ~ treat, data = surv.dat )
  logrank.p         <- 1 - pchisq( surv.logrank.test$chisq, 1)
  
  return(logrank.p)
}

# Wilcoxon Test
surv.wilcoxon.test.func <- function( surv.dat ) {
  surv.wilcoxon.test <- survdiff( Surv( survtime, status) ~ treat, rho = 1, data = surv.dat )
  wilcoxon.p <- 1 - pchisq( surv.wilcoxon.test$chisq, 1)
  
  return(wilcoxon.p)
}

# Other Log-rank Family Tests
surv.Other.test.func <- function( surv.dat, input) {
  
  treat_new <- factor(surv.dat$treat)
  treat_new <- unclass(treat_new)
  surv.dat$treat_new <- treat_new
  surv.dat$treat_new <- as.factor(surv.dat$treat_new)
  
  test_res <- logrank_test(Surv(survtime, status) ~ treat_new,
                           data = surv.dat,
                           type = input)
  
  return ( test_res@distribution@pvalue(test_res@statistic@teststatistic)  )
}


##############################################
##########   Kaplan-Meier Model     ########## & Plot
##############################################

# Kaplan-Meier Model
surv.KM.fit.func <- function( surv.dat ) {
  survfit( Surv( survtime, status )~treat, data = surv.dat, conf.type = c("plain"), type = 'kaplan-meier')
}

surv.KM.plot.func <- function( surv.KM.fit, surv.dat, input = "Log-rank test", follow.start, follow.end ) {
  ls <- unique(surv.dat$treat)
  if ( input != "Cox PH Model"){
    ggsurv <- ggsurvplot( surv.KM.fit,
                          data = surv.dat,
                          pval = FALSE, conf.int = TRUE,
                          linetype = c(1,2), 
                          #title = "Kaplan Meier curve",
                          risk.table = FALSE,
                          risk.table.col = "strata",
                          surv.median.line = "hv",
                          ggtheme = theme_bw(),
                          palette = c("#E7B800", "#2E9FDF"),
                          legend.title = "Treatment",
                          legend = "right",
                          legend.labs = c(ls[2], ls[1]) )
    
    if (input == "Log-rank test" || input == "Wilcoxon test"){
      output <- ifelse( input == "Log-rank test",  surv.logrank.test.func(surv.dat), surv.wilcoxon.test.func(surv.dat) )
    }
    else{
      output <- surv.Other.test.func( surv.dat, input)
    }
    
    output <- ifelse( input == "Log-rank test",  surv.logrank.test.func(surv.dat), surv.wilcoxon.test.func(surv.dat) )
    #output <- paste( input, "p-value =", format(round(output,3), nsmall = 3))
    p.val <- format(round(output,3), nsmall = 3)
    
    if (p.val < 0.001){
      output <- "P-value: <.001"
    }
    else{ output <- paste("P-value: ", p.val) }
    
    ggsurv$plot <- ggsurv$plot +
      ggplot2::annotate("text",
                        x = 8, y = 0.1, 
                        label = output, size = 4)
  }
  else {
    ggsurv <- ggsurvplot( survfit(surv.KM.fit), data = surv.dat,
                          color = "#2E9FDF",
                          ggtheme = theme_minimal() )
  }

  
  
    
  return(ggsurv)
}

##############################################
#####   Cox Proportional Hazards Model   ##### & Plot
##############################################

# Univariate Cox Model
surv.Cox.fit.func1 <- function( surv.dat) {
  coxph( Surv(survtime, status) ~ treat, data = surv.dat)
}

# Multivariate Cox Model
surv.Cox.fit.func2 <- function( surv.dat, input) {
  
  f <- paste(input, collapse = ' + ')
  cox.model <- coxph( as.formula( paste("Surv(survtime, status) ~ treat + ", f ) ), data = surv.dat )
  
  return( cox.model )
}

# Cox Model Plot
surv.Cox.plot.func <- function( cox.fit, surv.dat, input = NULL) {
  ls <- unique(surv.dat$treat)
  if ( is.null(input) ){
    new_dat <- with(surv.dat, data.frame(treat = c(ls[2],ls[1])))
  }
  else{
    new_dat <- with(surv.dat, data.frame(treat = c(ls[2],ls[1])))
    surv.pri <- surv.pri.func(surv.dat, is.mon.sin.rev.bin.con(surv.dat))
    for ( i in 1:length(input)){
      if ( input[i] %in% surv.pri$is.bin) { # var is binary
        ls_ <- unique(surv.dat[[input[i]]])
        new <- c(ls_[1], ls_[1])}#ls[2]) }
      else {
        new <- rep(mean(surv.dat[[input[i]]], na.rm = TRUE), 2)
      }
      new_dat[ , input[i] ] <- new 
    }
  }
  
  fit <- try(survfit(cox.fit, newdata = new_dat),silent = T)
  
  ggsurv <- ggsurvplot( fit,
                          data = surv.dat,
                          pval = FALSE, conf.int = TRUE,
                          linetype = c(1,2),
                          surv.median.line = "hv",
                          ggtheme = theme_bw(),
                          palette = c("#E7B800", "#2E9FDF"),
                          legend.title = "Treatment",
                          legend = "right",
                          legend.labs = c(ls[2], ls[1]) )
    
    
    HR <- round(exp(cox.fit$coefficients), 3)[1]# round( summary(cox.fit.1)$coefficient[2], 3)
    
    output <- paste("HR: ", HR)
    p.val <- round( summary(cox.fit)$coefficients[1,5], 3)# 
    
    if (!is.na(p.val)) {
      if (p.val < 0.001){
        output1 <- "P-value: <.001"
      }
      else{
        output1<- paste("P-value: ", p.val ) 
      }
      
      ggsurv$plot <- ggsurv$plot +
        ggplot2::annotate("text",
                          x = 8, y = 0.1, 
                          label = output, size = 4)
      
      ggsurv$plot <- ggsurv$plot +
        ggplot2::annotate("text",
                          x = 8, y = 0.06, 
                          label = output1, size = 4)
    }
    return( ggsurv  )
}



##############################################
##############################################




