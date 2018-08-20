
# Install packages
if(!require('reshape2')) install.packages('reshape2'); library(reshape2)
if(!require('devtools')) install.packages('devtools'); library(devtools)
if(!require('ggbeeswarm')) install.packages('ggbeeswarm'); library(ggbeeswarm)
if(!require('ggpubr')) install.packages('ggpubr'); library(rjags)
if(!require('cowplot')) install.packages('cowplot'); library(cowplot)
if(!require('dplyr')) install.packages('dplyr'); library(dplyr)
if(!require('boot')) install.packages('boot'); library(boot)
library(ggplot2)
library(gridExtra)

# Get rid of scientific notation numbers
options(scipen=999)


# CFR ---------------------------------------------------------------------

CFR.data <- read.csv(file="./data/cfr_all.csv", header = TRUE, stringsAsFactors = FALSE)

CFR.data <- CFR.data[-which(is.na(CFR.data$CFR.adjusted)),] # Get rid of NAs
CFR.data <- CFR.data[-which(CFR.data$Mid.Date<1800),]       # Remove stuff before 1800

CFR.data$Deaths[is.na(CFR.data$Deaths)] <- round(CFR.data$CFR.adjusted[is.na(CFR.data$Deaths)] * CFR.data$Cases[is.na(CFR.data$Deaths)], 0)
CFR.data$Outbreak.Size <- CFR.data$Cases
CFR.data$CFR.adjusted <- CFR.data$CFR.adjusted*100


# CFR 1880 compared to 1940-1950 ------------------------------------------

CFR1880 <- CFR.data %>% filter(Mid.Date<1890)
CFR1880 <- CFR1880 %>% filter(Mid.Date>1880)
CFR1880 %>%  summarize(mean=mean(CFR.adjusted))

CFR194050 <- CFR.data %>% filter(Mid.Date<1950)
CFR194050 <- CFR194050 %>% filter(Mid.Date>1940)
CFR194050 %>%  summarize(mean=mean(CFR.adjusted))

sum(CFR.data$Deaths, na.rm = T) / sum(CFR.data$Cases, na.rm = T)


# CFR by decade -----------------------------------------------------------

fit1 <- smooth.spline(CFR.data$Mid.Date,CFR.data$CFR.adjusted,df=5) #5 degrees of freedom

pred <- predict(fit1, x=seq(min(CFR.data$Mid.Date),max(CFR.data$Mid.Date), by=1))
pred.frame <- data.frame(pred)
(p1 <- ggplot() + geom_point(data=CFR.data, aes(x=Mid.Date, y=CFR.adjusted, group=Region, color=Region, size=Outbreak.Size)) +xlab("Date (year)")+ylab("Case Fatality Ratio") + geom_line(data=pred.frame, aes(x,y)))


# CFR by outbreak size ----------------------------------------------------

CFR.Size.data <- CFR.data[!is.na(CFR.data$Cases),]
ggplot(data=CFR.Size.data, aes(x=Outbreak.Size, y=CFR.adjusted)) + 
  geom_point(aes(group=Region, color=Region),size=2)  + 
  stat_smooth(method = "lm", col = "black",size=0.5) + scale_x_log10() 

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




