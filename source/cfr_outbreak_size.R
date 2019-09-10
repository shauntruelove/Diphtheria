if(!require('readr')) install.packages('readr'); library(readr)
if(!require('dplyr')) install.packages('dplyr'); library(dplyr)
if(!require('ggplot2')) install.packages('ggplot2'); library(ggplot2)


CFR.data <- read_csv("data/CFR.data.csv")

# CFR by outbreak size ----------------------------------------------------

CFR.Size.data <- CFR.data[!is.na(CFR.data$Cases),]

if (!is.na(CFR.Size.data$Deaths)){ 
  CFR.Size.data$CFR <- round(CFR.Size.data$Deaths / CFR.Size.data$Cases, 3) * 100
}

lm_size <- lm(CFR ~ log(Cases/100) + Mid.Date, data=CFR.Size.data)
summary(lm_size)
exp(coef(lm_size)[2])

predict(lm_size, data.frame(Cases=c(10, 100, 1000, 10000, 100000), Mid.Date=1990))
predict(lm_size, data.frame(Cases=c(10, 100, 1000, 10000, 100000), Mid.Date=2000))
predict(lm_size, data.frame(Cases=c(10, 100, 1000, 2000, 3000, 10000, 20000, 30000, 100000), Mid.Date=2000))

tmp <- predict(lm_size, data.frame(Cases=c(10, 100, 1000, 10000, 100000), Mid.Date=2000))
(tmp[-5] - tmp[-1]) / tmp[-5]
tmp <- predict(lm_size, data.frame(Cases=c(10, 100, 101, 102, 103), Mid.Date=2000))
(tmp[-5] - tmp[-1]) / tmp[-5]
tmp <- predict(lm_size, data.frame(Cases=c(10, 100, 101, 110, 1100), Mid.Date=2000))
(tmp[-5] - tmp[-1]) / tmp[-5]
tmp 

confint(lm_size)


ggplot(data=CFR.Size.data, aes(x=Outbreak.Size, y=CFR.adjusted)) + 
  geom_point(aes(group=Region, color=Region),size=2)  + 
  stat_smooth(method = "lm", col = "black",size=0.5) + scale_x_log10() 
