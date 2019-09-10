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


# FIGURES FOR SUPPLEMENT --------------------------------------------------

## Figure S1: Cases vs DTP3
source('source/Vaccination_Coverage_Incidence.R')

ggplot() + 
  geom_line(data=VaxCoverage, aes(x=Year,y = value, color=variable)) + 
  geom_line(data=VaxCoverage2, aes(x=Year,y = value*1000, color=variable)) +
  ylab("Cases") + 
  scale_y_continuous(sec.axis = sec_axis(~.*.001, name = "DPT3 Coverage (%)")) + 
  theme(legend.position = "top") + scale_colour_manual(values=c("darkblue", "tomato"),name = " ") + 
  background_grid(major = "xy", minor = "none") +
  theme(text=element_text(family="Myriad Pro")) 


## Figure S3: Age Distribution of Cases
S3 <- ggplot(data=MeanAge, aes(x=Year.Mid, y=meanAge))+
  geom_point(aes(color=Region, size=TotalCases)) + 
  xlab("Outbreak Year") + 
  ylab("Average Age") +
  geom_smooth(method = "lm", color="Black", size=0.75)+ 
  scale_color_discrete(labels=c("NIS"= "FSU"))+
  scale_x_continuous(breaks=c(1890,1930,1970,2010)) + 
  labs(size = "Outbreak\nSize")

## Figure S4: Case fatality ratio by outbreak size
S4 <- ggplot(data=CFR.Size.data, aes(x=Outbreak.Size, y=CFR.adjusted)) + 
  geom_point(aes(group=Region, color=Region))  + 
  stat_smooth(method = "lm", col = "black",size=0.5) + 
  scale_x_log10() +
  xlab("Log Outbreak Size") + 
  ylab("Case Fatality Ratio (%)") + 
  scale_color_discrete(labels=c("NIS"= "FSU")) + 
  background_grid(major = "xy", minor = "none") +
  theme(text=element_text(family="Myriad Pro")) 

## Fig. S5: Effectiveness of diphtheria vaccination 
load("~/Documents/diphtheria/trunk/results/VE_results.RData") ## Load Vaccine efficacy
load("~/Documents/diphtheria/trunk/results/VE_age_boots.RData") ## Load Vaccine efficacy
source(paste0(trunk.dir, '/source/Severity_By_Vaccination_Staus.R')) ## Load Vaccine effectiveness at preventing severe cases
source(paste0(trunk.dir, '/source/Duration_Of_Protection_and_Boosting.R')) ## Load Boosting and Coverage by age
#source(paste0(trunk.dir, '/source/Antitoxin.R')) ## Load Antitoxin
source(paste0(trunk.dir, '/source/schick_test_analysis.R')) ## Load Schick Test

VE_ests_boot_summary_part<- VE_ests_boot_summary %>% filter(Dose=="Partial")
VE_ests_boot_summary_full<- VE_ests_boot_summary %>% filter(Dose=="Full")
load(file='./trunk/results/VE.partial.stan.RData') # Loads VE.partial.stan.fit
VE.partial <- rstan::extract(VE.partial.stan.fit)$VE
load(file='./trunk/results/VE.full.stan.RData') # Loads VE.full.stan.fit
VE.full <- rstan::extract(VE.full.stan.fit)$VE
VE_part<-data.frame(VE.partial)
VE_part$Population<-"Total"
VE_Full<-data.frame(VE.full)
VE_Full$Population<-"Total"
VE_part <- VE_part %>% filter(VE.partial>=0 & VE.partial<=1) # get rid of any less than 0 or greater than 1
VE_Full <- VE_Full %>% filter(VE.full>=0 & VE.full<=1) # get rid of any less than 0 or greater than 1
VE_ests_boot_summary_part<- VE_ests_boot_summary %>% filter(Dose=="Partial")
VE_ests_boot_summary_full<- VE_ests_boot_summary %>% filter(Dose=="Full")

## Fig. S5a: VE by study
q1<-ggplot(VE_ests_boot_summary_part, aes(x=VE, y=(Population))) + 
  geom_point( aes(group=(Population)), shape=15) + xlim(0, 1) + 
  geom_errorbarh( aes(xmax=CI_high, xmin=CI_low), height=0.1, size=0.6) + 
  #geom_vline(xintercept = 0, linetype="dashed") + 
  theme( axis.ticks.x = element_blank(),axis.line.x=element_blank(), 
         axis.title.x=element_blank(), axis.text.x = element_blank()) + 
  ylab("Partially\nVaccinated") +
  background_grid(major = "y", minor = "none") +
  theme(text=element_text(family="Myriad Pro")) 

q2<- ggplot()  + 
  geom_violin(data=VE_part, aes(Population, VE.partial), fill="black",trim=TRUE)+ ylim(0, 1) +
  coord_flip()  + 
  scale_x_discrete(labels=c("Total" = "San Antonio, TX 1970")) + 
  xlab("Partially\nVaccinated") +
  theme(axis.text.y=element_text(colour = "white"), axis.title.y=element_text(colour = "white"),axis.line.x=element_blank(), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.text.x = element_blank()) +
  geom_vline(xintercept = 0, linetype="dashed") + 
  background_grid(major = "y", minor = "none") +
  theme(text=element_text(family="Myriad Pro")) 

q3<-ggplot(VE_ests_boot_summary_full, aes(x=VE, y=(Population))) + 
  geom_point( aes(group=(Population)), shape=15) + 
  xlim(0, 1) + 
  geom_errorbarh( aes(xmax=CI_high, xmin=CI_low), height=0.1, size=0.6) + 
  #geom_vline(xintercept = 0, linetype="dashed") + 
  ylab("Fully\nVaccinated") + 
  theme(axis.line.x=element_blank(), axis.ticks.x = element_blank(), 
        axis.title.x=element_blank(), axis.text.x = element_blank()) + 
  background_grid(major = "y", minor = "none") +
  theme(text=element_text(family="Myriad Pro")) 

q4<-ggplot()  + 
  geom_violin(data=VE_Full, aes(Population, VE.full), fill="black",trim=FALSE)+ ylim(0, 1) +
  coord_flip()  + 
  scale_x_discrete(labels=c("Total" = "San Antonio, TX 1970")) + 
  xlab("Partially\nVaccinated") +
  theme(axis.text.y=element_text(colour = "white"), axis.title.y=element_text(colour = "white")) +
  geom_vline(xintercept = 0, linetype="dashed") + 
  ylab("Vaccine Effectiveness") +
  background_grid(major = "y", minor = "none") +
  theme(text=element_text(family="Myriad Pro")) 

(p1<-plot_grid(q1,q2,q3,q4, nrow=4,ncol = 1, rel_heights = c(0.8, 0.2, 1,0.5)))

q10<-ggplot(VE_ests_boot_summary_part, aes(x=VE, y=(Population))) + 
  geom_point( aes(group=(Population)), shape=15) + xlim(0, 1) + 
  geom_errorbarh( aes(xmax=CI_high, xmin=CI_low), height=0.1, size=0.6) + 
  #geom_vline(xintercept = 0, linetype="dashed") + 
  theme( axis.ticks.x = element_blank(),axis.line.x=element_blank(), axis.title.x=element_blank(), axis.text.x = element_blank(), axis.title.y=element_text(colour = "white"),axis.text.y=element_text(colour = "white")) + 
  ylab("Partially\nVaccinated") +
  background_grid(major = "y", minor = "none") +
  theme(text=element_text(family="Myriad Pro")) 

q20<- ggplot()  + 
  geom_violin(data=VE_part, aes(Population, VE.partial), fill="black",trim=TRUE)+ ylim(0, 1) +
  coord_flip()  + 
  scale_x_discrete(labels=c("Total" = "San Antonio, TX 1970")) + 
  xlab("Partially\nVaccinated") +
  theme(axis.text.y=element_text(colour = "white"), axis.title.y=element_text(colour = "white"),
        axis.line.x=element_blank(), axis.ticks.x = element_blank(), axis.title.x=element_blank(), 
        axis.text.x = element_blank()) +
  geom_vline(xintercept = 0, linetype="dashed") + 
  background_grid(major = "y", minor = "none") +
  theme(text=element_text(family="Myriad Pro")) 

q30<-ggplot(VE_ests_boot_summary_full, aes(x=VE, y=(Population))) + 
  geom_point( aes(group=(Population)), shape=15) + 
  xlim(0, 1) + 
  geom_errorbarh( aes(xmax=CI_high, xmin=CI_low), height=0.1, size=0.6) + 
  #geom_vline(xintercept = 0, linetype="dashed") + 
  ylab("Fully\nVaccinated") + 
  theme(axis.line.x=element_blank(), axis.ticks.x = element_blank(), axis.title.x=element_blank(), 
        axis.text.x = element_blank(), axis.title.y=element_text(colour = "white"),
        axis.text.y=element_text(colour = "white")) + 
  background_grid(major = "y", minor = "none") +
  theme(text=element_text(family="Myriad Pro")) 

q40<-ggplot()  + 
  geom_violin(data=VE_Full, aes(Population, VE.full), fill="black",trim=FALSE)+ ylim(0, 1) +
  coord_flip()  + 
  scale_x_discrete(labels=c("Total" = "San Antonio, TX 1970")) + 
  xlab("Partially\nVaccinated") +
  theme(axis.text.y=element_text(colour = "white"), axis.title.y=element_text(colour = "white"), 
        axis.title.x=element_text(colour = "white"),axis.text.x=element_text(colour = "white")) +
  geom_vline(xintercept = 0, linetype="dashed") + 
  ylab("Vaccine Effectiveness") +
  background_grid(major = "y", minor = "none") +
  theme(text=element_text(family="Myriad Pro")) 

(p10<-plot_grid(q10,q20,q30,q40, nrow=4,ncol = 1, rel_heights = c(0.8, 0.2, 1,0.5)))

## Fig. S5b: VE by age
load(file='./trunk/Results/VE_age_boots.RData')  # load VE_age_boots_summary, VE_age_boots, VE_age_boots_long
#View(VE_age_boots_summary)
VE_age_boots_summary_part <-VE_age_boots_summary %>% filter(Dose=="Partial")
VE_age_boots_summary_full <-VE_age_boots_summary %>% filter(Dose=="Full")
VE_age_boots_summary_part$AC<-NULL
VE_age_boots_summary_part$AC<-c("3","2","1")
VE_age_boots_summary_full$AC<-NULL
VE_age_boots_summary_full$AC<-c("3","2","1")

q10<-ggplot(VE_age_boots_summary_part, aes(x=Mean, y=(AC))) + 
  geom_point( aes(group=(AC)), shape=15) + xlim(0, 1) + 
  geom_errorbarh(aes(xmax=CI.UB, xmin=CI.LB), height=0.1, size=0.6) + 
  #geom_vline(xintercept = 0, linetype="dashed") + 
  theme(axis.line.x=element_blank(), axis.ticks.x = element_blank(), axis.title.x=element_blank(), 
        axis.text.x = element_blank(), axis.title.y=element_text(colour = "white"),
        axis.text.y=element_text(colour = "white")) + 
  ylab("Partially\nVaccinated") +
  scale_y_discrete(labels=c("3"="0-4y","2"="5-19y","1"="≥20"))+ 
  background_grid(major = "y", minor = "none") +
  theme(text=element_text(family="Myriad Pro")) 

q30<-ggplot(VE_age_boots_summary_full, aes(x=Mean, y=(AC))) + 
  geom_point( aes(group=(AC)), shape=15) + 
  xlim(0, 1) + 
  geom_errorbarh(aes(xmax=CI.UB, xmin=CI.LB), height=0.1, size=0.6) + 
  #geom_vline(xintercept = 0, linetype="dashed") + 
  xlab("Vaccine Effectiveness")+
  ylab("Fully\nVaccinated") + 
  scale_y_discrete(labels=c("3"="0-4y","2"="5-19y","1"="≥20")) + 
  background_grid(major = "y", minor = "none") +
  theme(text=element_text(family="Myriad Pro"), axis.title.y=element_text(colour = "white"),
        axis.text.y=element_text(colour = "white"), axis.title.x=element_text(colour = "white"),
        axis.text.x=element_text(colour = "white")) 


(p20<-plot_grid(q10,q30, nrow=2,ncol = 1, rel_heights = c(0.8, 1)))


q1<-ggplot(VE_age_boots_summary_part, aes(x=Mean, y=(AC))) + 
  geom_point( aes(group=(AC)), shape=15) + xlim(0, 1) + 
  geom_errorbarh(aes(xmax=CI.UB, xmin=CI.LB), height=0.1, size=0.6) + 
  #geom_vline(xintercept = 0, linetype="dashed") + 
  theme(axis.line.x=element_blank(), axis.ticks.x = element_blank(), 
        axis.title.x=element_blank(), axis.text.x = element_blank()) + 
  ylab("Partially\nVaccinated") +
  scale_y_discrete(labels=c("3"="0-4y","2"="5-19y","1"="≥20"))+ 
  #background_grid(major = "xy", minor = "none") +
  theme(text=element_text(family="Myriad Pro")) 

q3<-ggplot(VE_age_boots_summary_full, aes(x=Mean, y=(AC))) + 
  geom_point( aes(group=(AC)), shape=15) + 
  xlim(0, 1) + 
  geom_errorbarh(aes(xmax=CI.UB, xmin=CI.LB), height=0.1, size=0.6) + 
  #geom_vline(xintercept = 0, linetype="dashed") + 
  xlab("Vaccine Effectiveness")+
  ylab("Fully\nVaccinated") + 
  scale_y_discrete(labels=c("3"="0-4y","2"="5-19y","1"="≥20")) + 
  #background_grid(major = "xy", minor = "none") +
  theme(text=element_text(family="Myriad Pro")) 

(p2<-plot_grid(q1,q3, nrow=2,ncol = 1, rel_heights = c(0.8, 1)))



top_row <- plot_grid(p10, p20, labels = c('A', 'B'),nrow = 1, ncol=2, align = 'h') 

vp1 <- viewport(width = 0.5, height = 1, x = 0.25, y = 0.5)
vp2 <- viewport(width = 0.5, height = 1, x = 0.75, y = 0.5)
vp <- viewport(width = 1, height = 1, x = 0, y = 0)

blankplot1<-ggplot(VE_ests_boot_summary_part, aes(x=VE, y=(Population))) +  
  geom_point(color="white", alpha=0) +
  ylab("Fully\nVax") +
  xlab(" ") +
  xlim(0, 1) + 
  theme(text=element_text(family="Myriad Pro")) + 
  scale_y_discrete(labels=c("Ukraine 1992"=" ","Lao 2012"=" ","Russia 1995"= "San Antonio, TX 1970",
                            "Viet Nam 1999"= " ", "San Antonio, TX 1970" = " ", "Elgin, TX 1970"= " ", "Russia 1993"= " ")) + 
  theme(axis.title.y=element_text(colour = "white"), axis.ticks.y=element_blank(),axis.ticks.x=element_blank(), 
        axis.title.x=element_text(colour = "white"), axis.line.x= element_blank(),axis.text.y=element_text(colour = "white")) +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 1, linetype="dashed") +
  background_grid(major = "x", minor = "none") 


blankplot2<-ggplot(VE_age_boots_summary_full,  aes(x=Mean, y=(AC))) +  
  geom_point(color="white", alpha=0) +
  scale_y_discrete(labels=c("3"=" ","2"="5-19y","1"=" ")) + 
  ylab("Fully\nVax") +
  xlab(" ") +
  xlim(0, 1) + 
  theme(text=element_text(family="Myriad Pro")) + 
  theme(axis.title.y=element_text(colour = "white"),axis.text.y=element_text(colour = "white"), 
        axis.ticks.y=element_blank(),axis.ticks.x=element_blank(), axis.title.x=element_text(colour = "white"), 
        axis.line.x= element_blank()) +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 1, linetype="dashed") +
  background_grid(major = "x", minor = "none") 

print(top_row)
print(blankplot1, vp=vp1)
print(p1, vp=vp1)
print(blankplot2, vp=vp2)
print(p2, vp=vp2)
grid.text("Total", 0.125, 0.645, gp=gpar(fontfamily="Myriad Pro"))
grid.text("Total", 0.125, 0.135, gp=gpar(fontfamily="Myriad Pro"))

source('source/Duration_Of_Protection_and_Boosting.R') ## Load Boosting and Coverage by age
source(paste0(trunk.dir, '/source/schick_test_analysis.R')) ## Load Schick Test
## Figure S6: Boosting and Coverage by age
S6 <- DurProt %>% 
  ggplot(aes(x, y, group=g, color=g, fill=g)) + 
  theme(legend.position="none")+
  annotate("text", x = 34, y = 15, label = "Full Protection")+ 
  annotate("text", x = 34, y = 60, label = "Basic Protection")+ 
  stat_smooth(
    geom = 'area', method = 'loess', position="stack", span = 1/3,
    alpha = 1/2) +
  labs(y="Percent Protected (%)", x="Age") + 
  background_grid(major = "xy", minor = "none") +
  theme(text=element_text(family="Myriad Pro")) 

# Figure S7: Schick
S7 <- ggplot(schick.dat1, aes(x=Age, y=Prct.Exposed)) + 
  geom_point(aes(color=Loc)) + 
  geom_smooth(data=schick.dat1, aes(x=Age, y=Prct.Exposed),color="black", size=0.5, alpha=0.2) + 
  xlab("Age (year)") + 
  ylab("Percent Previously Exposed (%)") + 
  geom_hline(aes(yintercept=100), linetype="dashed") + 
  labs(color="Location") + 
  theme(legend.position = c(0.65, 0.28)) + 
  background_grid(major = "xy", minor = "none") +
  theme(text=element_text(family="Myriad Pro")) #+ 


##  Vaccine effectiveness at preventing severe cases
dodge <- position_dodge(width=0.7)

p1<-ggplot() + 
  geom_point(data=VaxSever, aes(x=Immunization.Status, y=Percent.Severe, group=Location, 
                                size = Study.Size, color=Location), shape=15, position=dodge) + 
  geom_errorbar(data=VaxSever, mapping=aes(x=Immunization.Status, ymin=Percent.Lower.CI, 
                                           ymax=Percent.Upper.CI,color=Location), width=0.1, size=0.6, position=dodge) + 
  geom_violin(data=VaxTotalLong, aes(x=Type, y=value, group=interaction(Type)),position=position_nudge(x=0.41), 
              width=0.15, fill="gray30",trim = FALSE) + 
  xlab("Immunization Status") + 
  ylab("Percent Severe Cases (%)")+ 
  scale_x_discrete(labels=c("Full" = "Full\nDTP3", "Imperfect" = "Imperfect\nDTP 1 or 2","None" = "None")) +
  labs(size="Study\nSize") +
  ylim(0,100)+theme(aspect.ratio = 0.6) + 
  background_grid(major = "xy", minor = "none") +
  theme(text=element_text(family="Myriad Pro")) 


# Figure S8: Vaccine effectiveness 
load("~/Documents/diphtheria/trunk/results/VE_results.RData") ## Load Vaccine efficacy
load("~/Documents/diphtheria/trunk/results/VE_age_boots.RData") ## Load Vaccine efficacy
source(paste0(trunk.dir, '/source/Severity_By_Vaccination_Staus.R')) ## Load Vaccine effectiveness at preventing severe cases

VE_ests_boot_summary_part<- VE_ests_boot_summary %>% filter(Dose=="Partial")
VE_ests_boot_summary_full<- VE_ests_boot_summary %>% filter(Dose=="Full")
load(file='./trunk/results/VE.partial.stan.RData') # Loads VE.partial.stan.fit
VE.partial <- rstan::extract(VE.partial.stan.fit)$VE
load(file='./trunk/results/VE.full.stan.RData') # Loads VE.full.stan.fit
VE.full <- rstan::extract(VE.full.stan.fit)$VE
VE_part<-data.frame(VE.partial)
VE_part$Population<-"Total"
VE_Full<-data.frame(VE.full)
VE_Full$Population<-"Total"
VE_part <- VE_part %>% filter(VE.partial>=0 & VE.partial<=1) # get rid of any less than 0 or greater than 1
VE_Full <- VE_Full %>% filter(VE.full>=0 & VE.full<=1) # get rid of any less than 0 or greater than 1
VE_ests_boot_summary_part<- VE_ests_boot_summary %>% filter(Dose=="Partial")
VE_ests_boot_summary_full<- VE_ests_boot_summary %>% filter(Dose=="Full")

## Fig. 4a: VE by study
q1<-ggplot(VE_ests_boot_summary_part, aes(x=VE, y=(Population))) + 
  geom_point( aes(group=(Population)), shape=15) + xlim(0, 1) + 
  geom_errorbarh( aes(xmax=CI_high, xmin=CI_low), height=0.1, size=0.6) + 
  #geom_vline(xintercept = 0, linetype="dashed") + 
  theme( axis.ticks.x = element_blank(),axis.line.x=element_blank(), axis.title.x=element_blank(), 
         axis.text.x = element_blank()) + 
  ylab("Partially\nVaccinated") +
  background_grid(major = "y", minor = "none") +
  theme(text=element_text(family="Myriad Pro")) 

q2<- ggplot()  + 
  geom_violin(data=VE_part, aes(Population, VE.partial), fill="black",trim=TRUE)+ ylim(0, 1) +
  coord_flip()  + 
  scale_x_discrete(labels=c("Total" = "San Antonio, TX 1970")) + 
  xlab("Partially\nVaccinated") +
  theme(axis.text.y=element_text(colour = "white"), axis.title.y=element_text(colour = "white"),
        axis.line.x=element_blank(), axis.ticks.x = element_blank(), axis.title.x=element_blank(), 
        axis.text.x = element_blank()) +
  geom_vline(xintercept = 0, linetype="dashed") + 
  background_grid(major = "y", minor = "none") +
  theme(text=element_text(family="Myriad Pro")) 

q3<-ggplot(VE_ests_boot_summary_full, aes(x=VE, y=(Population))) + 
  geom_point( aes(group=(Population)), shape=15) + 
  xlim(0, 1) + 
  geom_errorbarh( aes(xmax=CI_high, xmin=CI_low), height=0.1, size=0.6) + 
  #geom_vline(xintercept = 0, linetype="dashed") + 
  ylab("Fully\nVaccinated") + 
  theme(axis.line.x=element_blank(), axis.ticks.x = element_blank(), axis.title.x=element_blank(), 
        axis.text.x = element_blank()) + 
  background_grid(major = "y", minor = "none") +
  theme(text=element_text(family="Myriad Pro")) 

q4<-ggplot()  + 
  geom_violin(data=VE_Full, aes(Population, VE.full), fill="black",trim=FALSE)+ ylim(0, 1) +
  coord_flip()  + 
  scale_x_discrete(labels=c("Total" = "San Antonio, TX 1970")) + 
  xlab("Partially\nVaccinated") +
  theme(axis.text.y=element_text(colour = "white"), axis.title.y=element_text(colour = "white")) +
  geom_vline(xintercept = 0, linetype="dashed") + 
  ylab("Vaccine Effectiveness") +
  background_grid(major = "y", minor = "none") +
  theme(text=element_text(family="Myriad Pro")) 

(p1<-plot_grid(q1,q2,q3,q4, nrow=4,ncol = 1, rel_heights = c(0.8, 0.2, 1,0.5)))

q10<-ggplot(VE_ests_boot_summary_part, aes(x=VE, y=(Population))) + 
  geom_point( aes(group=(Population)), shape=15) + xlim(0, 1) + 
  geom_errorbarh( aes(xmax=CI_high, xmin=CI_low), height=0.1, size=0.6) + 
  #geom_vline(xintercept = 0, linetype="dashed") + 
  theme( axis.ticks.x = element_blank(),axis.line.x=element_blank(), axis.title.x=element_blank(), 
         axis.text.x = element_blank(), axis.title.y=element_text(colour = "white"),
         axis.text.y=element_text(colour = "white")) + 
  ylab("Partially\nVaccinated") +
  background_grid(major = "y", minor = "none") +
  theme(text=element_text(family="Myriad Pro")) 

q20<- ggplot()  + 
  geom_violin(data=VE_part, aes(Population, VE.partial), fill="black",trim=TRUE)+ ylim(0, 1) +
  coord_flip()  + 
  scale_x_discrete(labels=c("Total" = "San Antonio, TX 1970")) + 
  xlab("Partially\nVaccinated") +
  theme(axis.text.y=element_text(colour = "white"), axis.title.y=element_text(colour = "white"),
        axis.line.x=element_blank(), axis.ticks.x = element_blank(), axis.title.x=element_blank(), 
        axis.text.x = element_blank()) +
  geom_vline(xintercept = 0, linetype="dashed") + 
  background_grid(major = "y", minor = "none") +
  theme(text=element_text(family="Myriad Pro")) 

q30<-ggplot(VE_ests_boot_summary_full, aes(x=VE, y=(Population))) + 
  geom_point( aes(group=(Population)), shape=15) + 
  xlim(0, 1) + 
  geom_errorbarh( aes(xmax=CI_high, xmin=CI_low), height=0.1, size=0.6) + 
  #geom_vline(xintercept = 0, linetype="dashed") + 
  ylab("Fully\nVaccinated") + 
  theme(axis.line.x=element_blank(), axis.ticks.x = element_blank(), axis.title.x=element_blank(), 
        axis.text.x = element_blank(), axis.title.y=element_text(colour = "white"),axis.text.y=element_text(colour = "white")) + 
  background_grid(major = "y", minor = "none") +
  theme(text=element_text(family="Myriad Pro")) 

q40<-ggplot()  + 
  geom_violin(data=VE_Full, aes(Population, VE.full), fill="black",trim=FALSE)+ ylim(0, 1) +
  coord_flip()  + 
  scale_x_discrete(labels=c("Total" = "San Antonio, TX 1970")) + 
  xlab("Partially\nVaccinated") +
  theme(axis.text.y=element_text(colour = "white"), axis.title.y=element_text(colour = "white"), 
        axis.title.x=element_text(colour = "white"),axis.text.x=element_text(colour = "white")) +
  geom_vline(xintercept = 0, linetype="dashed") + 
  ylab("Vaccine Effectiveness") +
  background_grid(major = "y", minor = "none") +
  theme(text=element_text(family="Myriad Pro")) 

(p10<-plot_grid(q10,q20,q30,q40, nrow=4,ncol = 1, rel_heights = c(0.8, 0.2, 1,0.5)))

## Fig. 4b: VE by age
load(file='./trunk/Results/VE_age_boots.RData')  # load VE_age_boots_summary, VE_age_boots, VE_age_boots_long
#View(VE_age_boots_summary)
VE_age_boots_summary_part <-VE_age_boots_summary %>% filter(Dose=="Partial")
VE_age_boots_summary_full <-VE_age_boots_summary %>% filter(Dose=="Full")
VE_age_boots_summary_part$AC<-NULL
VE_age_boots_summary_part$AC<-c("3","2","1")
VE_age_boots_summary_full$AC<-NULL
VE_age_boots_summary_full$AC<-c("3","2","1")

q10<-ggplot(VE_age_boots_summary_part, aes(x=Mean, y=(AC))) + 
  geom_point( aes(group=(AC)), shape=15) + xlim(0, 1) + 
  geom_errorbarh(aes(xmax=CI.UB, xmin=CI.LB), height=0.1, size=0.6) + 
  #geom_vline(xintercept = 0, linetype="dashed") + 
  theme(axis.line.x=element_blank(), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.text.x = element_blank(), axis.title.y=element_text(colour = "white"),axis.text.y=element_text(colour = "white")) + 
  ylab("Partially\nVaccinated") +
  scale_y_discrete(labels=c("3"="0-4y","2"="5-19y","1"="≥20"))+ 
  background_grid(major = "y", minor = "none") +
  theme(text=element_text(family="Myriad Pro")) 

q30<-ggplot(VE_age_boots_summary_full, aes(x=Mean, y=(AC))) + 
  geom_point( aes(group=(AC)), shape=15) + 
  xlim(0, 1) + 
  geom_errorbarh(aes(xmax=CI.UB, xmin=CI.LB), height=0.1, size=0.6) + 
  #geom_vline(xintercept = 0, linetype="dashed") + 
  xlab("Vaccine Effectiveness")+
  ylab("Fully\nVaccinated") + 
  scale_y_discrete(labels=c("3"="0-4y","2"="5-19y","1"="≥20")) + 
  background_grid(major = "y", minor = "none") +
  theme(text=element_text(family="Myriad Pro"), axis.title.y=element_text(colour = "white"),axis.text.y=element_text(colour = "white"), axis.title.x=element_text(colour = "white"),axis.text.x=element_text(colour = "white")) 

(p20<-plot_grid(q10,q30, nrow=2,ncol = 1, rel_heights = c(0.8, 1)))

q1<-ggplot(VE_age_boots_summary_part, aes(x=Mean, y=(AC))) + 
  geom_point( aes(group=(AC)), shape=15) + xlim(0, 1) + 
  geom_errorbarh(aes(xmax=CI.UB, xmin=CI.LB), height=0.1, size=0.6) + 
  #geom_vline(xintercept = 0, linetype="dashed") + 
  theme(axis.line.x=element_blank(), axis.ticks.x = element_blank(), axis.title.x=element_blank(), axis.text.x = element_blank()) + 
  ylab("Partially\nVaccinated") +
  scale_y_discrete(labels=c("3"="0-4y","2"="5-19y","1"="≥20"))+ 
  #background_grid(major = "xy", minor = "none") +
  theme(text=element_text(family="Myriad Pro")) 

q3<-ggplot(VE_age_boots_summary_full, aes(x=Mean, y=(AC))) + 
  geom_point( aes(group=(AC)), shape=15) + 
  xlim(0, 1) + 
  geom_errorbarh(aes(xmax=CI.UB, xmin=CI.LB), height=0.1, size=0.6) + 
  #geom_vline(xintercept = 0, linetype="dashed") + 
  xlab("Vaccine Effectiveness")+
  ylab("Fully\nVaccinated") + 
  scale_y_discrete(labels=c("3"="0-4y","2"="5-19y","1"="≥20")) + 
  #background_grid(major = "xy", minor = "none") +
  theme(text=element_text(family="Myriad Pro")) 

(p2<-plot_grid(q1,q3, nrow=2,ncol = 1, rel_heights = c(0.8, 1)))



top_row <- plot_grid(p10, p20, labels = c('A', 'B'),nrow = 1, ncol=2, align = 'h') 
#bottom_row <- plot_grid(p3, p4, labels = c('C', 'D'), align = 'h')
#plot1<-plot_grid(top_row, bottom_row,  ncol = 1, rel_heights = c(1, 1))

vp1 <- viewport(width = 0.5, height = 1, x = 0.25, y = 0.5)
vp2 <- viewport(width = 0.5, height = 1, x = 0.75, y = 0.5)

blankplot1<-ggplot(VE_ests_boot_summary_part, aes(x=VE, y=(Population))) +  
  geom_point(color="white", alpha=0) +
  ylab("Fully\nVax") +
  xlab(" ") +
  xlim(0, 1) + 
  theme(text=element_text(family="Myriad Pro")) + 
  scale_y_discrete(labels=c("Ukraine 1992"=" ","Lao 2012"=" ","Russia 1995"= "San Antonio, TX 1970", "Viet Nam 1999"= " ", "San Antonio, TX 1970" = " ", "Elgin, TX 1970"= " ", "Russia 1993"= " ")) + 
  theme(axis.title.y=element_text(colour = "white"), axis.ticks.y=element_blank(),axis.ticks.x=element_blank(), 
        axis.title.x=element_text(colour = "white"), axis.line.x= element_blank(),axis.text.y=element_text(colour = "white")) +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 1, linetype="dashed") +
  background_grid(major = "x", minor = "none") 

blankplot2<-ggplot(VE_age_boots_summary_full,  aes(x=Mean, y=(AC))) +  
  geom_point(color="white", alpha=0) +
  scale_y_discrete(labels=c("3"=" ","2"="5-19y","1"=" ")) + 
  ylab("Fully\nVax") +
  xlab(" ") +
  xlim(0, 1) + 
  theme(text=element_text(family="Myriad Pro")) + 
  theme(axis.title.y=element_text(colour = "white"),axis.text.y=element_text(colour = "white"), axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(), axis.title.x=element_text(colour = "white"), axis.line.x= element_blank()) +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 1, linetype="dashed") +
  background_grid(major = "x", minor = "none") 


print(top_row)
print(blankplot1, vp=vp1)
print(p1, vp=vp1)
print(blankplot2, vp=vp2)
print(p2, vp=vp2)
grid.text("Total", 0.125, 0.645, gp=gpar(fontfamily="Myriad Pro"))
grid.text("Total", 0.125, 0.135, gp=gpar(fontfamily="Myriad Pro"))






