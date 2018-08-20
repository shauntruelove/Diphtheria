if(!require('readr')) install.packages('readr'); library(readr)
if(!require('dplyr')) install.packages('dplyr'); library(dplyr)

CFR.data <- read_csv("data/CFR.data.csv")

CFR.noTreat.noVacc<-CFR.data %>% filter(Mid.Date<1928, Antitoxin=="No")
CFR.noTreat.noVacc %>% summarise(mean=mean(CFR.adjusted), ql = quantile(CFR.adjusted,0.025),qu = quantile(CFR.adjusted,0.975))