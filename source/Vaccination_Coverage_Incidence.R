library(reshape2)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(knitr)


# Load Data ---------------------------------------------------------------

VaxCov <- read.csv(file="data/who_vacc_cov_and_incid.csv", header = TRUE, stringsAsFactors = FALSE)

VaxCoverage<-melt(VaxCov[,-3], id.vars="Year")
VaxCoverage2<-melt(VaxCov[,-2], id.vars="Year")

# Save multiple objects
save(VaxCoverage, VaxCoverage2, file = "data/VacCovIncid.RData")


