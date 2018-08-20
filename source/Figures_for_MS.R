if(!require('reshape2')) install.packages('reshape2'); library(reshape2)
if(!require('ggplot2')) install.packages('ggplot2'); library(ggplot2)
if(!require('dplyr')) install.packages('dplyr'); library(dplyr)
if(!require('gridExtra')) install.packages('gridExtra'); library(gridExtra)
if(!require('knitr')) install.packages('knitr'); library(knitr)
if(!require('ggbeeswarm')) install.packages('ggbeeswarm'); library(ggbeeswarm)
if(!require('cowplot')) install.packages('cowplot'); library(cowplot)
if(!require('ggpubr')) install.packages('ggpubr'); library(ggpubr)
if(!require('grid')) install.packages('grid'); library(grid)
if(!require('latex2exp')) install.packages('latex2exp'); library(latex2exp)
if(!require("devtools")) install.packages("devtools");library(devtools)
if(!require('quantreg')) install.packages('quantreg'); library(quantreg)
if(!require('grid')) install.packages('grid'); library(grid)

#devtools::install_github("kassambara/ggpubr")

# FIGURES FOR PAPER -------------------------------------------------------

# FIGURES FOR PAPER -------------------------------------------------------
# Figure 2: Epidemiology Figure -------------------------------------------
serial_interval_curve <- read.csv(file='results/serial_interval_curve.csv', header=TRUE) ## Load Generation Time 
incubation_curve <- read.csv(file='results/incubation_curve.csv', header=TRUE)
clearance_curve_shared <- read.csv(file='results/clearance_curve_shared.csv', header=TRUE)


MyColour <- c("#756bb1", "#31a354") 
## Fig. 2a: Incubation Period
Incub <- ggplot(data=incubation_curve, aes(x=q)) +
  geom_ribbon(aes(ymin=plow, ymax=phigh, fill=type),  alpha=.35) +
  geom_line(aes(y=pmid, colour=type), size=0.5) + 
  scale_fill_manual(values="springgreen4") +
  scale_colour_manual(values="springgreen4") +
  scale_x_continuous(limits=c(0, 18), expand = c(0, 0)) + 
  theme(legend.position="none", aspect.ratio=0.9) +
  ylab("Proportion of Cases") + 
  xlab("Incubation Period (days)") +
  background_grid(major = "xy", minor = "none") +
  theme(text=element_text(family="Myriad Pro"))  

## Fig. 2b: Clearance
Clear <- ggplot(clearance_curve_shared, aes(x=q), group=type) +
  geom_ribbon(aes(ymin=vcplow*100, ymax=vcphigh*100, fill=type),  alpha=.35) +
  geom_line(aes(y=vcpmid*100, colour=type), size=0.5) + 
  scale_fill_manual(values=c("steelblue3"),name = "Infection Type") +
  scale_colour_manual(values=c("steelblue3"),name = "Infection Type") +
  scale_x_continuous(limits=c(0, 70), expand = c(0, 0)) + 
  ylab("Percent Colonized (%)") + xlab("Duration (days)") + 
  theme(legend.position="none", aspect.ratio=0.9) +
  background_grid(major = "xy", minor = "none") +
  theme(text=element_text(family="Myriad Pro")) 

## Fig. 2c: Serial Interval
GenTime <- ggplot(serial_interval_curve, aes(x=q), group=type) +
  geom_ribbon(aes(ymin=plow, ymax=phigh, fill=type),  alpha=.35) +
  geom_line(aes(y=pmid, colour=type), size=0.5) +
  scale_fill_manual(values="red") +
  scale_colour_manual(values="red") +
  scale_x_continuous(limits=c(0, 18), expand = c(0, 0)) + 
  theme(legend.position="none", aspect.ratio = 0.9) +
  ylab("Proportion of Cases") + 
  xlab("Serial Interval (days)") +
  background_grid(major = "xy", minor = "none") +
  theme(text=element_text(family="Myriad Pro")) 

(p <- plot_grid(Incub, Clear, GenTime, labels=c("A", "B", "C"), nrow=1, ncol = 3))


## Figure 3: CFR -----------------------------------------------------------
source('source/cfr_outbreak_size.R')
source('source/cfr_overall.R')
#source('source/cfr_age.R')
source('source/cfr_decade.R')
cfr.dat.delay <- read.csv(file='results/cfr_antitoxin_delay_boots_results_raw.csv', header=TRUE)
cfr.age <- read.csv(file='results/cfr_age_results.csv', header = TRUE)
cfr.age.long <- melt(cfr.age)
cfr.vax <- read.csv(file='results/cfr_vax.csv', header = TRUE)
cfr.vax.long <- melt(cfr.vax)
cfr.antitoxin <- read.csv(file='results/cfr_antitoxin.csv', header=TRUE)
cfr.antitoxin.long <- melt(cfr.antitoxin)

## Fig 3a: CFR x Date
q0<-ggplot() + 
    geom_point(data=CFR.data, aes(x=Mid.Date, y=CFR.adjusted, group=Region, color=Region, size=(Outbreak.Size))) +
    xlab("Date (year)") +
    ylab("Case Fatality Ratio (%)") + 
    geom_line(data=pred.frame, aes(x=Decade,y=CFR)) + 
    scale_x_continuous(breaks=seq(1875,2017,25)) + 
    labs(size="Outbreak Size") + 
    scale_colour_discrete(labels=c("NIS" = "FSU")) +
    background_grid(major = "xy", minor = "none") +
    theme(text=element_text(family="Myriad Pro")) 

## Fig 3b: CFR x Age
q1 <- ggplot(cfr.age.long, aes(x=variable, y=value)) +
  geom_boxplot()+ 
  xlab("Age (y)") +
  ylab("Case Fatality Ratio (%)") + 
  scale_x_discrete(labels=c("age1"= "0-4 Years", "age2"="5-19 Years", "age3"="≥20 Years")) +
  theme(legend.position = "none",plot.margin = unit(c(0,0,0,0), "cm"), axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 0.8) + 
  background_grid(major = "xy", minor = "none") +
  theme(text=element_text(family="Myriad Pro")) 

## Fig 3c: CFR x DAT
q2<-ggplot(cfr.antitoxin.long, aes(x=variable, y=value)) +geom_boxplot() + 
    xlab("Antitoxin")+ylab("Case Fatality Ratio (%)") +xlab("Antitoxin Treatment") +
    scale_shape_manual(values=c(16, 24))+scale_colour_gradient( limits=c(1730, 2010)) + 
    guides(shape = guide_legend(order=1),colour = guide_colourbar(order = 2))  +
    theme(legend.position = "none",axis.title.y=element_blank(),plot.margin = unit(c(0,0,0,0), "cm")) + 
    scale_x_discrete(labels=c("no_antitoxin" = "No\nAntitoxin", "antitoxin" = "Antitoxin"))+ theme(axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 0.8) +
    background_grid(major = "xy", minor = "none") +
    theme(text=element_text(family="Myriad Pro")) 

## Fig 3e: CFR x Vax
q4<-ggplot(cfr.vax.long, aes(x=variable, y=value)) +geom_boxplot()+ 
    xlab("Vaccination Status")+ylab("Case Fatality Ratio (%)") + 
    scale_x_discrete(labels=c("NoVacc" = "None", "Partial" = "Imperfect\nDTP1 or 2", "Full" = "Full\nDTP3")) +
    theme(legend.position = "none",axis.title.y=element_blank(),plot.margin = unit(c(0,0,0,0), "cm"), axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 0.8) + 
    background_grid(major = "xy", minor = "none") +
    theme(text=element_text(family="Myriad Pro")) 

## Fig 3d: CFR x DAT Delay
q3 <- ggplot(cfr.dat.delay, aes(x=delay, y=prob)) +
    geom_boxplot() +
    xlab("Antitoxin Treatment Delay") +
    ylab("Case Fatality Ratio (%)") +
    theme(legend.position = "none",axis.title.y=element_blank(),plot.margin = unit(c(0,0,0,0), "cm"),axis.text.x = element_text(angle = 45, hjust = 1), aspect.ratio = 0.8) +
    scale_x_discrete(labels=c("1"= "1 Day", "2"="2 Days", "3"="3 Days", "4"="4 Days", "5+"="≥5 Days")) + 
    background_grid(major = "xy", minor = "none") +
    theme(text=element_text(family="Myriad Pro")) 


# arrange the three plots in a single row
bottomrow <- plot_grid( q1 + theme(legend.position="none"),
                     q2 + theme(legend.position="none"),
                     q3 + theme(legend.position="none"),
                     q4 + theme(legend.position="none"),
                     align = 'vh',
                     labels = c("B", "C", "D", "E") ,
                     #hjust = -1,
                     rel_widths = c(1.2,1.1,1.1,1.1),
                     nrow = 1)

toprow <- plot_grid(q0, labels = c('A'))

(p <- plot_grid(toprow, bottomrow, nrow=2,ncol = 1, rel_heights = c(0.5, .5)))





# Figure 4: Critical Vaccination Figure --------------------------------------

# Source Functions to make Figure 4
source('source/figure_4_functions.R')

# PLOT ALL TOGETHER -------------------------------------------------------

A <- plot_combined_VcR0(delay=1, alpha=0, alpha_density_R0=0.5)
B <- plot_combined_VcR0(delay=1, alpha=0.25, alpha_density_R0=0.5)
C <- plot_combined_VcR0(delay=2, alpha=0.5, alpha_density_R0=0.5)
D <- plot_combined_VcR0(delay=5, alpha=.9, alpha_density_R0=0.5)

plot_grid(A, B, C, D, nrow=2, ncol=2, labels=c('A','B','C','D'))

 
vp1 <- viewport(width = 0.5, height = 0.5, x = 0.25, y = 0.75)
vp2 <- viewport(width = 0.5, height = 0.5, x = 0.75, y = 0.75)
vp3 <- viewport(width = 0.5, height = 0.5, x = 0.25, y = 0.25)
vp4 <- viewport(width = 0.5, height = 0.5, x = 0.75, y = 0.25)


print(blankplot2_fn())
print(blankplot1_fn(), vp = vp1)
print(A, vp = vp1)
print(blankplot1_fn(), vp = vp2)
print(B, vp = vp2)
print(blankplot1_fn(), vp = vp3)
print(C, vp = vp3)
print(blankplot1_fn(), vp = vp4)
print(D, vp = vp4)





