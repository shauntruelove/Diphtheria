if(!require('readr')) install.packages('readr'); library(readr)
if(!require('dplyr')) install.packages('dplyr'); library(dplyr)

CFR.data <- read_csv("data/CFR.data.csv")

# CFR ---------------------------------------------------------------------
CFR.df<-NULL
ConfIntCFR<-CFR.data %>% filter(!is.na(Cases))
CFR.df$mean<-ConfIntCFR %>% summarise(mean=mean(CFR.adjusted),LowCI=100*(mean(CFR.adjusted)/100) - 1.96*sqrt(((mean(CFR.adjusted)/100)*(1-(mean(CFR.adjusted)/100)))/mean(Cases)), HighCI=100*(mean(CFR.adjusted)/100) + 1.96*sqrt(((mean(CFR.adjusted)/100)*(1-(mean(CFR.adjusted)/100)))/sum(Cases)))

CFR.df$median<-CFR.data %>% summarize(med=median(CFR.adjusted), ql = quantile(CFR.adjusted,.25),qu = quantile(CFR.adjusted,.75))
