
source('source/vc_treat.R')

if(!require('cowplot')) install.packages('cowplot'); library(cowplot)

# library(extrafont)
# font_import()
library(showtext)
library(tidyverse)

font_add("Myriad Pro", "MyriadPro-Regular.otf")
showtext_auto()





# Function to set up data for Figure 4 Vc x R0 plot -----------------------

Vc_curve_plot_data <- function(delay=1, alpha=.9){

  if(alpha==0) {
    
    # LOAD STAN DATA ----------------------------------------------------------
    
    load("results/R0_res_knownGT_2000iters.RData")
    #load("Vc_with_75pct_treatment_3days.RData")
    
    ## Extract data from stan
    knGT_samps <- rstan::extract(fit_knGT_fixCoV)
    ## Make a clean dataframe  
    data_raw <- data.frame(R0 = knGT_samps$R0_preds, Vc = knGT_samps$Vc*100)
    rm(knGT_samps)
    
  } else {
    # Get Treatment-Adjusted Vc estimates by R0
    vc_treat_data <- vc_treat(delay = delay, alpha = alpha, restrict_quantiles=c(0,1))
    data_raw <- data.frame(R0 = vc_treat_data$R0_preds, Vc = vc_treat_data$Vc.treat*100)
    data_raw <- data_raw[order(data_raw$R0),]
    rm(vc_treat_data)
  }  
    
  # Make smooth data for polygon where Vc<100% ------------------------------
  
  ## Calculate the 95% CI and mean and run a regression to smooth them
  summaryQ <- data_raw %>% mutate(R0=round(R0, 1)) %>%
    group_by(R0) %>% 
    summarize(pLo = quantile(Vc, probs = 0.025, na.rm = TRUE),
              pMid = quantile(Vc, probs = 0.5, na.rm = TRUE),
              pHi = quantile(Vc, probs = 0.975, na.rm = TRUE)) 

  ## Store the mean in a data frame
  mid.smooth <- data.frame(R0= seq(1,5,0.001), 
                           Vc=predict(loess(summaryQ$pMid ~ summaryQ$R0), newdata = seq(1,5,0.001)),
                           type='mid') 
  low.smooth <- data.frame(R0= seq(1,5,0.001), 
                           Vc=predict(loess(summaryQ$pLo ~ summaryQ$R0), newdata = seq(1,5,0.001)),
                           type='low')
  high.smooth <- data.frame(R0= seq(1,5,0.001), 
                            Vc=predict(loess(summaryQ$pHi ~ summaryQ$R0), newdata = seq(1,5,0.001)),
                            type='high')
  smooth.data <- rbind(mid.smooth, low.smooth, high.smooth)
  smooth.data <- smooth.data[which(!is.na(smooth.data$Vc)),]
  smooth.data<- smooth.data %>% filter(Vc >= 0) ## Filter to be greater than Vc=0
  

  ## Filter to Vc < 100%
  smooth.data.lessthan100 <- smooth.data %>% filter(Vc <= 100) ## Filter to be lower than Vc=100
  
  smooth.data.lessthan100.wide <- reshape(smooth.data.lessthan100, idvar = 'R0', timevar = 'type', direction = 'wide')
  # Fill max Vc as 100
  smooth.data.lessthan100.wide$Vc.high[is.na(smooth.data.lessthan100.wide$Vc.high)] <- 100
  # Fill min Vc as 0
  smooth.data.lessthan100.wide$Vc.low[is.na(smooth.data.lessthan100.wide$Vc.low)] <- 0
  
  
  # Make the data for the Vc > 100 % (red) scatter --------------------------
  data.morethan100 <- data_raw %>% filter(Vc>=100)
  
  
  # Return list of 2 data.frames for plotting
  return(list(data.morethan100=data.morethan100, 
              smooth.data.lessthan100.wide=smooth.data.lessthan100.wide))
}



# Function to plot Vc by R0 curve -----------------------------------------

plot_Vc_curve_plot_data_untreated <- function(restrict_quantiles=c(0.025,0.975)){
  
  
  # Untreated data
  data_untreated <- Vc_curve_plot_data(delay=1, alpha=0)
  r0_iqr <- get_R0_IQR(delay=1, alpha=0, restrict_quantiles=c(0,1))
    
  ggplot() + 
    annotate("rect", xmin=0.9, xmax=4.2, ymin=100, ymax=Inf, alpha=0.12, fill="maroon") +
    #annotate("rect", xmin=r0_iqr[1], xmax=r0_iqr[2], ymin=0, ymax=Inf, alpha=0.2, fill="grey") +
    coord_cartesian(xlim = c(1.09, 4)) +
    geom_point(data=data_untreated$data.morethan100, aes(x=R0, y=Vc), alpha=0.3, shape=1, color="indianred1") +
    geom_line(data=data_untreated$smooth.data.lessthan100.wide, aes(x=R0, y=Vc.mid)) +
    geom_ribbon(data=data_untreated$smooth.data.lessthan100.wide, aes(x=R0, ymin=Vc.low, ymax=Vc.high), fill="Blue", alpha=0.5) +
    geom_hline(yintercept=100, color="red") + 
    ylim(0,200) +
    ylab(expression("Critical vaccination threshold (%)")) +
    #ylab(expression(paste(italic(V[c]), "  (%)"))) +
    background_grid(major = "xy", minor = "none") +
    annotate(geom="text", x=1.1, y=200, label='A. No treatment', color="black", hjust=0) +
    theme(text=element_text(size = 10, family="MyriadPro"),
          axis.title.x=element_blank(), 
          axis.title.y=element_text(face="italic")) 
}



# Function to plot Vc by R0 curve, with Vc adjusted for antibiotic --------

plot_Vc_curve_plot_data_treated <- function(delay=1, alpha=0, panel='A', restrict_quantiles=c(0.025,0.975)){
  
  # Untreated data
  data_untreated <- Vc_curve_plot_data(delay=1, alpha=0)
  # Treated data
  data_treated <- Vc_curve_plot_data(delay=delay, alpha=alpha)

  r0_iqr <- get_R0_IQR(delay=1, alpha=0, restrict_quantiles=restrict_quantiles)
  treat_label <- paste0(panel, '. ', round(alpha*100,0),'% treated, ', delay,'-day delay')
  
  ggplot() + 
    annotate("rect", xmin=.9, xmax=4.2, ymin=100, ymax=Inf, alpha=0.12, fill="maroon") +
    #annotate("rect", xmin=r0_iqr[1], xmax=r0_iqr[2], ymin=0, ymax=Inf, alpha=0.2, fill="grey") +
    coord_cartesian(xlim = c(1.09,4)) +
    geom_point(data=data_untreated$data.morethan100, aes(x=R0, y=Vc), alpha=0.15, shape=4, color="grey") +
    geom_ribbon(data=data_untreated$smooth.data.lessthan100.wide, aes(x=R0, ymin=Vc.low, ymax=Vc.high), fill="steelblue3", alpha=0.3) +
    geom_line(data=data_untreated$smooth.data.lessthan100.wide, aes(x=R0, y=Vc.mid), colour="grey")+
    geom_point(data=data_treated$data.morethan100, aes(x=R0, y=Vc), alpha=0.3, shape=1, color="indianred1") +
    geom_line(data=data_treated$smooth.data.lessthan100.wide, aes(x=R0, y=Vc.mid)) +
    geom_ribbon(data=data_treated$smooth.data.lessthan100.wide, aes(x=R0, ymin=Vc.low, ymax=Vc.high), fill="Blue", alpha=0.5) +
    geom_hline(yintercept=100, color="red") + 
    ylim(0,200) +
    ylab(expression("Critical vaccination threshold (%)")) +
    #ylab(expression(paste(italic(V[c]), "  (%)"))) +
    background_grid(major = "xy", minor = "none") +
    annotate(geom="text", x=1.1, y=200, label=treat_label, color="black", hjust=0) +
    theme(text=element_text(size = 10, family="MyriadPro"),
          axis.title.x=element_blank(), 
          axis.title.y=element_text(face="italic"))
}




# Function to generate and plot R0 density curve, shaded by Vc < 100% -------

plot_R0_density <- function(delay=1, alpha=0.9, alpha_density=0.5, blank=FALSE, font_fam="MyriadPro"){
  # Source code to adjust Vc for treatment
  source('source/vc_treat.R')
  
  # Calculate Treated Vc
  data_<- vc_treat(delay = delay, alpha = alpha, restrict_quantiles = c(0,1))
  
  KnGT_R0 <- data.frame(R0 = data_$R0_preds, Vc = data_$Vc.treat*100)
  KnGT_R0 <- KnGT_R0[order(KnGT_R0$R0),]
  
  ## Make a dataframe of Just R0 Values
  KnGT_R0_df <- as.data.frame(KnGT_R0)
  KnGT_R0_df <- KnGT_R0_df %>% filter(R0<=4)
  KnGT_R0_df <- KnGT_R0_df %>% filter(R0>=1)
  
    ## Plot a blank geom density to extract informationf nrom
  p <- ggplot(KnGT_R0_df, aes(x=R0)) + 
    geom_density() + 
    scale_x_continuous(breaks=seq(0,5,by=1)) +
    scale_y_continuous(labels = c("0.0"="200","0.2"=" ","0.4"=" ","0.6"=" ","0.8"=" ")) +
    #xlim(1,4) +
    xlab(expression(paste(italic(R[0])))) +
    ylab(expression("Freq")) +
    #ylab(expression("Critical vaccination threshold (%)")) +
    #ylab(expression(paste(italic(V[c]), "  (%)"))) +
    theme(axis.text.y=element_text(colour = "white"),
          axis.title.y=element_text(size = 10, colour = "white"),
          text=element_text(family=font_fam),
          axis.title.x=element_text(face="italic"))
  
  ## Extract the height etc to fill
  q <- ggplot_build(p)$data[[1]] %>% dplyr::select(x, y)

  ## Calculate the color based on the proportion of R0 above 100%
  shade<-NULL
  color<-NULL
  index<-seq(from=1, to=3.5, by=0.05)
  for(i in 1:length(index)){
    low<-index[i]
    high<-index[i]+0.08
    shade[i]<-mean(KnGT_R0_df$Vc[KnGT_R0_df$R0 >= low & KnGT_R0_df$R0 < high]<100)
  }
  shade[which(is.na(shade))]<-1
  for(i in 1:length(index)){
    color[i]<-rgb(1-shade[i], 0, shade[i])
  }

  ## Add the different filled in areas to the plot
  for(i in 1:length(index)){
    p <- p + geom_area(data=subset(q, x>index[i] & x<=(index[i]+0.05572)),
                       aes(x=x, y=y), fill=color[i], alpha=alpha_density)
  }
  p <- p + geom_area(data = subset(q, x >= 3.55 & x <= 5), aes(x=x, y=y), fill=color[51], alpha=alpha_density)

  ## Add it all to on plot
  p + 
    geom_density() +
    #geom_vline(xintercept = 3.21, linetype=2) + 
    background_grid(major = "xy", minor = "none")
}


plot_combined_VcR0 <- function(delay=1, alpha=0, alpha_density_R0=0.5, panel='A'){
  plot_grid(
    # Top Plot
    if(alpha==0){
      plot_Vc_curve_plot_data_untreated()
    }else{
      plot_Vc_curve_plot_data_treated(delay=delay, alpha=alpha, panel=panel)
    },
    # Bottom Plot
    plot_R0_density(delay=delay, alpha=alpha, alpha_density=alpha_density_R0),
    
    nrow=2, ncol = 1, rel_heights = c(1, .5))
}


blankplot1_fn <- function(restrict_quantiles=c(0.025,.975), alpha_plot=0.35){
  # Source code to adjust Vc for treatment
  source('source/vc_treat.R')
  
  # Calculate Treated Vc
  data_<- vc_treat(delay = 1, alpha = 0, restrict_quantiles = restrict_quantiles)
  
  KnGT_R0 <- data.frame(R0 = data_$R0_preds, Vc = data_$Vc.treat*100)
  KnGT_R0 <- KnGT_R0[order(KnGT_R0$R0),]
  
  ## Make a dataframe of Just R0 Values
  KnGT_R0_df <- as.data.frame(KnGT_R0)
  KnGT_R0_df <- KnGT_R0_df %>% filter(R0<=4)
  KnGT_R0_df <- KnGT_R0_df %>% filter(R0>=1)
  
  r0_iqr <- get_R0_IQR(delay=1, alpha=0, restrict_quantiles=restrict_quantiles)
  
    ggplot(KnGT_R0_df) +
        annotate("rect", xmin=r0_iqr[1], xmax=r0_iqr[2], ymin=-Inf, ymax=Inf, alpha=alpha_plot, fill="gray90") +
        scale_x_continuous(breaks=seq(0,5,by=1), limits=c(1.09,4)) +
        ylim(0,200) +
        xlab(expression(paste(italic(R[0])))) +
        ylab(expression("Critical vaccination threshold (%)")) +
        #ylab(expression(paste(italic(V[c]), "  (%)"))) +
        theme(text=element_text(family="Myriad Pro"),
              axis.text.y=element_text(colour = "white"), 
              axis.title.y=element_text(size = 10, colour = "white"),
              axis.ticks.y=element_blank(),  
              axis.line.y= element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.x=element_text(colour = "white"), 
              axis.line.x= element_blank(),
              axis.text.x=element_text(colour = "white"))
}

blankplot2_fn <- function(restrict_quantiles = c(0.025,.975)){
  
  # Source code to adjust Vc for treatment
  source('source/vc_treat.R')
  
  # Calculate Treated Vc
  data_<- vc_treat(delay = 1, alpha = 0, restrict_quantiles = restrict_quantiles)
  
  KnGT_R0 <- data.frame(R0 = data_$R0_preds, Vc = data_$Vc.treat*100)
  KnGT_R0 <- KnGT_R0[order(KnGT_R0$R0),]
  
  ## Make a dataframe of Just R0 Values
  KnGT_R0_df <- as.data.frame(KnGT_R0)
  KnGT_R0_df <- KnGT_R0_df %>% filter(R0<=4)
  KnGT_R0_df <- KnGT_R0_df %>% filter(R0>=1)
  
    ggplot(KnGT_R0_df) +
        scale_x_continuous(breaks=seq(0,5,by=1), limits=c(1.09,4)) +
        xlab(expression(paste(italic(R[0])))) +
        ylab(expression("Critical vaccination threshold (%)")) +
        #ylab(expression(paste(italic(V[c]), "  (%)"))) +
        theme(text=element_text(family="Myriad Pro"),
              axis.text.y=element_text(colour = "white"), 
              axis.title.y=element_text(size = 10, colour = "white"),
              axis.ticks.y=element_blank(),  
              axis.line.y= element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.x=element_text(colour = "white"), 
              axis.line.x= element_blank(),
              axis.text.x=element_text(colour = "white"))
}