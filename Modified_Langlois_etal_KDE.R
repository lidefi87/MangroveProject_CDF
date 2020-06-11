#Code from Langlois, T.J., Fitzpatrick, B., Wakefield, C.B., Fairclough, D.V., Hesp, A., McLean, D., 
#Meeuwig, J.J., Harvey, E.S., in prep. Similarities between fish length-frequency distributions estimated 
#from line fishing and baited stereo-video: novel application of kernel density estimates.

#Modified by Denisse Fierro Arcos (DFA): Using ggplot to obtain graphs as outputs for further manipulation 
#in addition to estimates and p values originally given by the sm.density.compare() function
#Date: 2020-03-26

kde.compare <- function(length = Pagrus_auratus$Length,
                        group = Pagrus_auratus$Method,
                        align = c('no','by.median')[1],
                        nboot = 500,
                        xlab = 'Fish Length (mm)',
                        ylab = 'Probability Density',
                        main = ''){
  
  #check for required packages ('sm' and 'KernSmooth') and load if they are already installed 
  #or install and load them if they are not already installed
  check <- require(sm)
  if(check == FALSE) {install.packages('sm')}
  ks.check <- require(KernSmooth, quietly = TRUE)
  if(ks.check == FALSE) {install.packages('KernSmooth')}
  #Adding tidyverse to be able to use pipes and ggplot - DFA
  tidy.check <- require(tidyverse, quietly = TRUE)
  if(tidy.check == FALSE) {install.packages('tidyverse')}
  print(paste(c('Dependencies Installed:', require(sm,quietly = TRUE), 
                require(KernSmooth, quietly = TRUE),
                require(tidyverse, quietly = T))))
  
  library(sm) #should load without error
  library(KernSmooth) ##should load without error
  library(tidyverse)
  
  Data <- data.frame(length, group)
  
  #distribution alignment
  method <- levels(Data$group)
  length.aligner <- function(Data){
    for(i in 1:length(method)){
      yi = Data[Data$group == method[i],]$length
      med = median(yi)
      sc = diff(quantile(yi,c(0.25,0.75)))/1.349
      Data[Data$group==method[i],]$length <- (yi - med)/sc
    }
    return(Data)}    
  
  if(align == c('by.median')){Data <- length.aligner(Data)}
  
  #bandwidth calculation
  
  h.mean <- function(Data){
    hl = dpik(Data$length[Data$group == method[1]], kernel='normal')
    hv = dpik(Data$length[Data$group == method[2]], kernel='normal')
    h.m = (hl+hv)/2 #arithmetic mean of the KDEs
    return(h.m)}
  
  #####################
  
  sm.density.compare3 <- function (x, group, h, model = "none", ...) 
  {if (!is.vector(x)) 
      stop("sm.density.compare can handle only 1-d data")
    opt <- sm.options(list(...))
    sm:::replace.na(opt, ngrid, 50)
    sm:::replace.na(opt, display, "line")
    sm:::replace.na(opt, xlab, deparse(substitute(x)))
    sm:::replace.na(opt, ylab, "Density")
    sm:::replace.na(opt, xlim, c(min(x) - diff(range(x))/4, max(x) + 
                                   diff(range(x))/4))
    sm:::replace.na(opt, eval.points, seq(opt$xlim[1], opt$xlim[2], 
                                          length = opt$ngrid))
    if (is.na(opt$band)) {
      if (model == "none") 
        opt$band <- FALSE
      else opt$band <- TRUE}
    
    if ((model == "none") && opt$band){
      opt$band <- FALSE}
    band <- opt$band
    ngrid <- opt$ngrid
    xlim <- opt$xlim
    nboot <- opt$nboot
    y <- x
    
    if (is.na(opt$test)) {
      if (model == "none") 
        opt$test <- FALSE
      else opt$test <- TRUE}
    
    if ((model == "none") && opt$test) 
      opt$test <- FALSE
    test <- opt$test
    
    if (opt$display %in% "none") 
      band <- FALSE
    fact <- factor(group)
    fact.levels <- levels(fact)
    nlev <- length(fact.levels)
    ni <- table(fact)
    
    if (band & (nlev > 2)) {
      cat("Reference band available to compare two groups only.", 
          "\n")
      band <- FALSE}
    
    if (length(opt$lty) < nlev) 
      opt$lty <- 1:nlev
    
    if (length(opt$col) < nlev) 
      opt$col <- 2:(nlev + 1)
    
    if (missing(h)) 
      h <- h.select(x, y = NA, group = group, ...)
    opt$band <- band
    opt$test <- test
    #Estimate1 added as data frame to ease in the creation of graphs in ggplot - DFA
    estimate1 <- data.frame(est = rep(0, opt$ngrid*nlev), 
                            fact = rep(fact.levels, each = opt$ngrid))
    se <- matrix(0, ncol = opt$ngrid, nrow = nlev)
    
    for (i in 1:nlev) {
      sm <- sm.density(y[fact == fact.levels[i]], h = h, display = "none", 
                       eval.points = opt$eval.points)
      estimate1$est[estimate1$fact == fact.levels[i]] <- sm$estimate
      se[i, ] <- sm$se}
    
    #Estimate1 converted into matrix as set in original code by Bowan & Azzalini - DFA
    estimate <- estimate1 %>% mutate(id = rep(1:opt$ngrid, nlev)) %>% 
      pivot_wider(names_from = fact, values_from = est) %>% 
      select(-id) %>% t()
    eval.points <- sm$eval.points
    
    if (!(opt$display %in% "none" | band)) {
      #Swapping ggplot for base graphic commands to create graphs - DFA
      fig <- ggplot(mapping = aes(x = rep(eval.points, nlev), y = estimate1$est, 
                                  linetype = estimate1$fact))+
        xlim(xlim)+ylim(0, 1.1*max(as.vector(estimate1$est)))+
        labs(x = opt$xlab, y = opt$ylab, title = main)+
        geom_line()+theme_bw()+
        theme(legend.title = element_blank(),
              legend.position = c(0.85,0.85),
              axis.title = element_text(family = "sans", size = 12), 
              axis.text = element_text(family = "sans", size = 12),
              legend.text = element_text(family = "sans", size = 12))} 
    
    est <- NULL
    p <- NULL
    
    if (model == "equal" & test) {
      if (nlev == 2) {
        ts <- sum((estimate[1, ] - estimate[2, ])^2)
      }else {
        sm.mean <- sm.density(y, h = h, xlim = opt$xlim, 
                              ngrid = opt$ngrid, display = "none")$estimate
        ts <- 0
        for (i in 1:nlev) ts <- ts + ni[i] * sum((estimate[i,] - sm.mean)^2)}
      
      p <- 0
      est.star <- matrix(0, ncol = opt$ngrid, nrow = nlev)
      
      for (iboot in 1:nboot) {
        ind <- (1:length(y))
        for (i in 1:nlev) {
          indi <- sample((1:length(ind)), ni[i])
          est.star[i, ] <- sm.density(y[ind[indi]], h = h, 
                                      ngrid = opt$ngrid, xlim = opt$xlim, display = "none")$estimate
          ind <- ind[-indi]}
        
        if (nlev == 2) {
          ts.star <- sum((est.star[1, ] - est.star[2, ])^2)
        }else {
          sm.mean <- sm.density(y, h = h, xlim = opt$xlim, 
                                ngrid = opt$ngrid, display = "none")$estimate
          ts.star <- 0
          
          for (i in 1:nlev) {
            ts.star <- ts.star + ni[i] * sum((est.star[i,] - sm.mean)^2)}}
        
        if (ts.star > ts) 
          p <- p + 1
        
        if (opt$verbose > 1) {
          cat(iboot)
          cat(" ")}}
      
      p <- p/nboot
      cat("\nTest of equal densities:  p-value = ", round(p,3), "\n")
      est <- list(p = p, h = h)}
    
    if (model == "equal" & band) {
      av <- (sqrt(estimate[1, ]) + sqrt(estimate[2, ]))/2
      se <- sqrt(se[1, ]^2 + se[2, ]^2)
      upper <- (av + se)^2
      lower <- pmax(av - se, 0)^2
      #Swapping ggplot for base graphic commands to create graphs - DFA
      fig <- ggplot(mapping = aes(x = rep(eval.points, nlev), y = estimate1$est, 
                                  linetype = estimate1$fact))+
        xlim(xlim)+
        ylim(0, 1.1*max(as.vector(estimate1$est), upper))+
        labs(x = opt$xlab, y = opt$ylab, title = main)+
        geom_ribbon(aes(x = rep(eval.points, 2),
                        ymin = rep(lower,2), ymax = rep(upper,2)),
                    fill = opt$col.band, show.legend = F)+
        geom_line()+theme_bw()+
        theme(legend.title = element_blank(),
              legend.position = c(0.85,0.85),
              axis.title = element_text(family = "sans", size = 12), 
              axis.text = element_text(family = "sans", size = 12),
              legend.text = element_text(family = "sans", size = 12))
      est <- list(p = p, upper = upper, lower = lower, h = h)}
    
    invisible(est)
    #Creating output list containing two elements: ggplot and estimates - DFA
    return(list(graph = fig, est = est))}
  
  with(Data,sm.density.compare3(x = length, group = group, h = h.mean(Data),
                                nboot = nboot, model = "equal", ngrid=500,
                                col = c('black','black'),
                                xlab = paste(xlab),
                                ylab = paste(ylab), col.band = "grey"))} 
# End of Function: kde.compare

#Examples of how to use kde.compare() 
#Simply choose species, align method and number of bootstraps.

# kde.compare(length = Choerodon_rubescens$Length,
#             group = Choerodon_rubescens$Method,
#             align = 'no',nboot = 500,
#             xlab = 'Fish Length (mm)',
#             ylab = 'Probability Density',
#             main = 'Choerodon rubescens')
# 
# kde.compare(length = Choerodon_rubescens$Length,
#             group = Choerodon_rubescens$Method,
#             align = 'by.median',
#             nboot = 500,
#             xlab = 'Standardised Fish Lengths',
#             ylab = 'Probability Density',
#             main = 'Choerodon rubescens')

#######################################
