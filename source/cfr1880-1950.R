##########################################################################
##########################################################################
##                                                                      ##
## Calculate the Case Fatility Ratio in 1880 and in 1940-1950           ## 
##                                                                      ##
##########################################################################
##########################################################################

if(!require('dplyr')) install.packages('dplyr'); library(dplyr)


CFR.data <- read_csv("data/CFR.data.csv")

# CFR 1880 compared to 1940-1950 ------------------------------------------

CFR1880<-CFR.data %>% filter(Mid.Date<1890)
CFR1880<-CFR1880 %>% filter(Mid.Date>1880)
CFR1880.mean<-CFR1880 %>%  summarize(mean=mean(CFR.adjusted))

CFR194050<-CFR.data %>% filter(Mid.Date<1950)
CFR194050<-CFR194050 %>% filter(Mid.Date>1940)
CFR194050.mean<-CFR194050 %>%  summarize(mean=mean(CFR.adjusted))


CFR.compare.df<-data.frame(CFR1880.mean,CFR194050.mean)
colnames(CFR.compare.df)<-c("CFR 1880", "CFR 1940-1950")
