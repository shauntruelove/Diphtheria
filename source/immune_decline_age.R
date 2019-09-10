
#...........................................................................................
# Estimate decline in protection from vaccination with Age
#...........................................................................................


# Packages and Source Code ------------------------------------------------

if(!require('dplyr')) install.packages('dplyr'); library(dplyr)
if(!require('lme4')) install.packages('lme4'); library(lme4)
if(!require('effects')) install.packages('effects'); library(effects)


# Load and Clean Data -----------------------------------------------------

dat <- read.csv(file="./data/immune_decline_age.csv", header = TRUE, stringsAsFactors = FALSE)

dat.clean <- dat[,c("Ref_index","Study","Location","Year.Mid",'Age.min','Age.max',"prct.Full","prct.Basic","prct.Suscept")]
dat.clean <- dat.clean[!is.na(dat.clean$prct.Full),]

dat.clean$Age.mean <- (dat.clean$Age.max + dat.clean$Age.min) / 2
#View(dat.clean[is.na(dat.clean$Age.mean),])

#dat.clean$abs.decr.per.agegroup <- NA
dat.clean$abs.decr.per.year <- NA
dat.clean$prct.decr.per.year <- NA

pop.ind <- with(dat.clean, paste(Location, Study, Ref_index, Year.Mid, sep='-'))
#unique(pop.ind)
mean.study <- sd.study <- numeric(length(unique(pop.ind)))
for (i in 1:length(unique(pop.ind))) {
    tmp_inds <- which(pop.ind == unique(pop.ind)[i])
    dat.tmp <- dat.clean[tmp_inds,]
    diffs_ <- diff(dat.tmp$prct.Full) 
    dat.clean$abs.decr.per.year[tmp_inds[-length(tmp_inds)]] <- diffs_ / diff(dat.tmp$Age.mean)
    prct.diffs_ <- diffs_ / dat.tmp$prct.Full[-nrow(dat.tmp)]
    prct.diffs.per.year_ <- prct.diffs_ / diff(dat.tmp$Age.mean)
    dat.clean$prct.decr.per.year[tmp_inds[-length(tmp_inds)]] <- prct.diffs.per.year_
    prct.diffs.per.year_ <- prct.diffs.per.year_[!is.infinite(prct.diffs.per.year_) & !is.na(prct.diffs.per.year_) & !is.nan(prct.diffs.per.year_)]
    mean.study[i] <- mean(prct.diffs.per.year_, na.rm=T)
    sd.study[i] <- sd(prct.diffs.per.year_, na.rm=T)
}

prct.decr.per.year <- dat.clean$prct.decr.per.year
prct.decr.per.year <- prct.decr.per.year[!is.na(prct.decr.per.year) & !is.infinite(prct.decr.per.year)]
mean(prct.decr.per.year, na.rm=T)
hist(prct.decr.per.year, breaks=100)

abs.decr.per.year <- dat.clean$abs.decr.per.year
abs.decr.per.year <- abs.decr.per.year[!is.na(abs.decr.per.year) & !is.infinite(abs.decr.per.year)]
mean(abs.decr.per.year, na.rm=T)
hist(abs.decr.per.year, breaks=100)

mean.study <- mean.study[!is.infinite(mean.study) & !is.na(mean.study) & !is.nan(mean.study)]
sd.study <- sd.study[!is.infinite(sd.study) & !is.na(sd.study) & !is.nan(sd.study)]
mean(mean.study)
hist(mean.study, breaks=50)




# Linear model ------------------------------------------------------------
dat.clean$pop.ind <- with(dat.clean, paste(Location, Study, Ref_index, Year.Mid, sep='-'))
dat.clean$prct.Full[dat.clean$prct.Full==0] <- 0.00001
mod <- lmer(log(prct.Full) ~ Age.mean + (1 | pop.ind), data=dat.clean)

cbind("% increase per age yr" = fixef(mod), confint(mod, method="Wald", level = 0.95))



# How many and Which Studies
length(unique(dat.clean$Study))
length(unique(dat.clean$Ref_index))
unique(dat.clean$Ref_index)
  