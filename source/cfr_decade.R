##########################################################################
##########################################################################
##                                                                      ##
## Calculate the Case Fatility Ratio by decade using a spline with 5    ## 
## knots                                                                ##
##                                                                      ##
##########################################################################
##########################################################################


if(!require('readr')) install.packages('readr'); library(readr)
if(!require('dplyr')) install.packages('dplyr'); library(dplyr)

CFR.data <- read_csv("data/CFR.data.csv")

# CFR by decade -----------------------------------------------------------

fit1<-smooth.spline(CFR.data$Mid.Date,CFR.data$CFR.adjusted,df=5) #5 degrees of freedom

pred<-predict(fit1, x=seq(min(CFR.data$Mid.Date),max(CFR.data$Mid.Date), by=1))
pred.frame<-data.frame(pred)
colnames(pred.frame)<-c("Decade", "CFR")

