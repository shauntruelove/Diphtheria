# Install & Load packages
if(!require('reshape2')) install.packages('reshape2'); library(reshape2)
if(!require('ggplot2')) install.packages('ggplot2'); library(ggplot2)
if(!require('dplyr')) install.packages('dplyr'); library(dplyr)
if(!require('tidyr')) install.packages('tidyr'); library(tidyr)
if(!require('lme4')) install.packages('lme4'); library(lme4)
if(!require('meta')) install.packages('meta'); library(meta)
if(!require('readr')) install.packages('readr'); library(readr)
if(!require('effects')) install.packages('effects'); library(effects)

# Load Age Data -----------------------------------------------------------
Age.mult.imput <- read_csv("data/Age.mult.imput.csv")

# Sum cases to split ages back to the new groups
Age.grouped <- Age.mult.imput %>% group_by(Ref_index, Location, Region, Decade, Group) %>%
    summarise(Cases.group=round(sum(Cases, na.rm=TRUE),1)) %>% as.data.frame
Age.grouped <- rename(Age.grouped, Cases=Cases.group)

Age.grouped2 <- Age.mult.imput %>% group_by(Ref_index, Location, Region, Decade, Group2) %>%
    summarise(Cases.group=round(sum(Cases, na.rm=TRUE),1)) %>% as.data.frame
Age.grouped2 <- rename(Age.grouped2, Cases=Cases.group)



# Get age proportions before and after 1980 -------------------------------
Age.grouped %>% group_by(Group, Decade>=1980) %>% summarize(sum(Cases))
Age.grouped %>% group_by(Decade>=1980) %>% summarize(sum(Cases))


# Get Proportions for Each population/outbreak/study
Age.grouped <- Age.grouped %>% 
    group_by(Ref_index, Location, Region, Decade) %>%
    mutate(N = sum(Cases), Group.prop=round(Cases / sum(Cases),3)) %>%
    ungroup %>% as.data.frame
# Get Proportions for Each population/outbreak/study
Age.grouped2 <- Age.grouped2 %>% 
    group_by(Ref_index, Location, Region, Decade) %>%
    mutate(N = sum(Cases), Group.prop=round(Cases / sum(Cases),3)) %>%
    ungroup %>% as.data.frame
Age.grouped$Decade.fact <- factor(Age.grouped$Decade)
Age.grouped2$Decade.fact <- factor(Age.grouped2$Decade)

# Get Proportions for Each Decade
Age.grouped.decade <- Age.grouped %>% 
    group_by(Decade, Group) %>%
    summarise(N = sum(Cases)) %>%
    mutate(Group.prop=round(N / sum(N),3))
# Get Proportions for Each Decade
Age.grouped.decade2 <- Age.grouped2 %>% 
    group_by(Decade, Group2) %>%
    summarise(N = sum(Cases)) %>%
    mutate(Group.prop=round(N / sum(N),3))



# Model with mixed effect for Decade  ---------------------------------------

mod1 <- lm((Cases/N) ~ Group*Decade.fact, data=Age.grouped)
summary(mod1)

fit1980 <- predict(mod1, data.frame(Group=unique(Age.grouped$Group), Decade.fact="1980"), se.fit=TRUE)
fit1990 <- predict(mod1, data.frame(Group=unique(Age.grouped$Group), Decade.fact="1990"), se.fit=TRUE)
fit2000 <- predict(mod1, data.frame(Group=unique(Age.grouped$Group), Decade.fact="2000"), se.fit=TRUE)
fit2010 <- predict(mod1, data.frame(Group=unique(Age.grouped$Group), Decade.fact="2010"), se.fit=TRUE)

fits <- rbind(fit1980$fit, fit1990$fit, fit2000$fit, fit2010$fit)
row.names(fits) <- c(1980, 1990, 2000, 2010)
colnames(fits) <- unique(Age.grouped$Group)
t(round(fits,2))



fits_ci <- rbind(c(fit1980$fit[1], fit1980$fit[1] - 1.96*fit1980$se.fit[1], fit1980$fit[1] + 1.96*fit1980$se.fit[1]),
              c(fit1990$fit[1], fit1990$fit[1] - 1.96*fit1990$se.fit[1], fit1990$fit[1] + 1.96*fit1990$se.fit[1]),
              c(fit2000$fit[1], fit2000$fit[1] - 1.96*fit2000$se.fit[1], fit2000$fit[1] + 1.96*fit2000$se.fit[1]),
              c(fit2010$fit[1], fit2010$fit[1] - 1.96*fit2010$se.fit[1], fit2010$fit[1] + 1.96*fit2010$se.fit[1]))
row.names(fits_ci) <- c(1980, 1990, 2000, 2010)
fits_ci


# USE THIS FOR PAPER ------------------------------------------------------



dat <- Age.grouped2[Age.grouped2$Group2=='>=20y' & Age.grouped2$Decade<1980,]
m1 <- metaprop(round(dat$Cases,0), round(dat$N,0), sm="PLOGIT"); summary(m1)
summary(dat$Group.prop)
range(dat$Group.prop)

dat <- Age.grouped2[Age.grouped2$Group2=='>=20y' & Age.grouped2$Decade>=1980,]
m1 <- metaprop(round(dat$Cases,0), round(dat$N,0), sm="PLOGIT"); summary(m1)
summary(dat$Group.prop)
range(dat$Group.prop)


length(unique(paste0(Age.grouped2$Location, Age.grouped2$Decade)))





# Model with mixed effect for Decade  ---------------------------------------
library(lme4)

Age.grouped$Decade.fact <- factor(Age.grouped$Decade)

mod2 <- lmer(Group.prop ~ Group + (1|Decade.fact), data=Age.grouped)
summary(mod2)
fixef(mod2)
exp(fixef(mod2))
exp(confint(mod2, method="Wald", level = 0.95)[-1,])

# Convert to Probabilities
library(effects)
(ef.1 <- data.frame(effect(c("Group"), mod2)))

predict(mod2, data.frame(GroupDecade.fact=c(1950,1950,1970,1970,1990,1990)), type="response")





mod2 <- glmer(died ~ TreatDelay + (1|Ref_index), data=CFR.DAT.long)
dat.indiv <- CFR.DAT.long

set.seed(20)
dat.indiv$ID <- 1:nrow(dat.indiv)
dat.indiv <- as.data.frame(dat.indiv)
tmp <- sampler(dat.indiv, "Ref_index", reps = 1000)
bigdata <- cbind(tmp, dat.indiv[tmp$RowID,])
save(dat.indiv, bigdata, mod2, file = './trunk/data/CFRxDATdelay_boots.RData')
f <- fixef(mod2)
r <- getME(mod2, "theta")

library(parallel)
cl <- makeCluster(6)
clusterExport(cl, c("bigdata", "f", "r"))
clusterEvalQ(cl, require(lme4))






# Adults and Children

Age.grouped2$Decade.fact <- factor(Age.grouped2$Decade)
Age.grouped2$LocDecade <- factor(paste0(Age.grouped2$Location, Age.grouped2$Decade))

Age.grouped2 <- Age.grouped2 %>% mutate(row = 1:nrow(Age.grouped2))
Age2.long <- Age.grouped2 <- Age.grouped2 %>% group_by(row) %>%
  do({ left_join(data_frame(row=.$row, case = rep(1,.$Cases)),.,by='row') })
Age2.long$adult <- ifelse(Age2.long$Group2=="0-19y", 0,1)

mod2 <- glmer(adult ~ I(Decade>=1980) + (1|LocDecade), data=Age2.long, family = binomial(link = "logit"))
mod2 <- glm(adult ~ I(Decade>=1980), data=Age2.long, family = binomial(link = "logit"))

summary(mod2)
fixef(mod2)
exp(cbind("Odds ratio" = fixef(mod2), confint(mod2, method="Wald", level = 0.95)[-1,]))

# Convert to Probabilities
library(effects)
(ef.1 <- data.frame(effect(c("Group2"), mod2)))
probs1

predict(mod2, data.frame(adult=c(0,1,0,1,0,1), Decade=c(1950,1950,1970,1970,1990,1990)), type="response")























#### Make a fake population using the weights given.
#### For a given size, distribute the percentages over that population
#### using the cumulative distribution function

## How large a population to simulate
population_size = 1000
## Convert percentages to number of members of our fake population

Age.Data.integer <- Age.mult.imput %>%
  mutate_all(funs({
    . = ./sum(.) * (population_size - .2)
    ceiling(cumsum(.))-ceiling(c(0,lag(cumsum(.))[-1]))
  }))

Age.Data.integer <- Age.mult.imput %>%
  mutate(funs({
    . = ./sum(.) * (population_size - .2)
    ceiling(cumsum(.))-ceiling(c(0,lag(cumsum(.))[-1]))
  }))

## Age doesn't need to sum to population, so we convert that back
Age.Data.integer$Age = Age.Data$Age


## For each population count, make an actual population of that count
all_names = names(Age.Data.integer)
##Pre-allocate a population data frame.  We now have to go back and populate it
##We subtract 1 here because we don't need age
Age.Data.population = as.data.frame(matrix(NA,nrow=population_size,ncol=length(all_names)-1))
names(Age.Data.population)  <-  all_names[-1]
for(name in all_names[-1]){
  ### this is essentially an unlist(rep(age,count)), but modified so the syntax works.
  Age.Data.population[[name]]  <-  unlist(mapply(count = Age.Data.integer[[name]],age = Age.Data.integer[['Age']], function(age,count){c(rep(age,count))}))
}
## Melt for ggplot
Age.Data.integer.melted = melt(Age.Data.population)

## Once we've melted, pull out the years separately
##    (This is the reason the variables have an @ symbol in them)
## We want to keep the locations 
Age.Data.melted <- Age.Data.integer.melted %>%
  arrange(as.character(variable)) %>%
  separate(col='variable',sep='@ *',into= c('Year','Location')) %>%
  separate(col='Year',sep='-',into=c('start','end')) %>%
  mutate(
    Range = paste(start,end,sep='-'),
    Decade = floor(ifelse(is.na(end),as.numeric(start),.5*(as.numeric(start)+as.numeric(end)))/10)*10
  ) %>%
  mutate(xaxis = paste(Range,Location)) %>%
  arrange(xaxis)



ggplot(Age.mult.imput) +
  geom_violin(aes(x=Location,y=Age), bw=10/3) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) + ylab("Age")+ scale_fill_brewer(palette = "Greens")






Age.Data <- data.frame(AgeData)
MeanAge.Data <- data.frame(MeanAgeData)
colnames(Age.Data) <- c("Age",
                        "1990-1999@ Thailand", 
                        "1993-1996@ Georgia", 
                        "1967-1969@ Austin, TX, USA", 
                        "1947@      Esher County, UK", 
                        "1961-1962@ Omaha, NE, USA", 
                        "1970@      San Antonio, TX, USA", 
                        "1939-1951@ Jamacia", 
                        "1940-1941@ Halifax, NS, Canada", 
                        "1993-1997@ Georgia", 
                        "1908-1917@ Maryland, USA",
                        "1892-1897@ London County, UK", 
                        "1981-1982@ Yemen", 
                        "1977-1978@ Jordan", 
                        "1982-1983@ Jordan", 
                        "1978@      Sudan", 
                        "1988@      Sudan", 
                        "1989@      Lesotho", 
                        "1940-1950@ Los Angeles, CA, USA", 
                        "2011@      Maharashtra, India", 
                        "1991-1997@ Ukraine", 
                        "2008@      Assam, India", 
                        "1947@      Utah", 
                        "1990-1996@ Belarus", 
                        "1994-1996@ Moldova", 
                        "1990-1996@ Armenia", 
                        "1990-1996@ Azerbaijan")

Alph.Age.Data <- Age.Data[,sort(colnames(Age.Data))]
AgeMelted <- melt(Alph.Age.Data, id.vars = "Age")
AgeMelted$Decade <- rep(c("1890", "1900", "1910", "1940","1940", "1940", "1940","1940","1960","1960", "1970","1970", "1970","1980", "1980","1980", "1990","1990","1990", "1990","1990", "1990","1990","1990", "2000", "2010"), each=9)
#AgeMelted1 <- melt(Alph.Age.Data1, id.vars = "Age")
colnames(AgeMelted) <- c("Age", "Location", "value", "Decade")
MeanAgeData  <-  MeanAge.Data[ -c(2:25) ]



ggplot(Age.Data.melted) +
  geom_violin(aes(x=xaxis,y=value, fill=as.factor(Decade)), bw=10/3) +
  theme(axis.text.x=element_text(angle=90,hjust=1)) + ylab("Age")+ scale_fill_brewer(palette = "Greens") + scale_x_discrete(labels=c("1892-1897 London County, UK" = "London County, UK", "1908-1917 Maryland, USA" = "Maryland, USA", "1939-1951 Jamacia" ="Jamacia","1940-1941 Halifax, NS, Canada" =" Halifax, NS, Canada",  "1940-1950 Los Angeles, CA, USA" ="Los Angeles, CA, USA", "1947-NA Esher County, UK" = "Esher County, UK", "1947-NA Utah" ="Utah, USA", "1961-1962 Omaha, NE, USA" ="Omaha, NE, USA", "1967-1969 Austin, TX, USA" ="Austin, TX, USA", "1970-NA San Antonio, TX, USA" ="San Antonio, TX, USA",   "1977-1978 Jordan"=" Jordan", "1978-NA Sudan" = "Sudan", "1981-1982 Yemen" ="Yemen",  "1982-1983 Jordan"="Jordan", "1988-NA Sudan" ="Sudan", "1989-NA Lesotho" ="Lesotho", "1990-1996 Armenia" ="Armenia", "1990-1996 Azerbaijan"= "Azerbaijan", "1990-1996 Belarus" ="Belarus", "1990-1999 Thailand"="Thailand", "1991-1997 Ukraine"="Ukraine", "1993-1996 Georgia"="Georgia", "1993-1997 Georgia"= "Georgia", "1994-1996 Moldova"="Moldova", "2008-NA Assam, India"="Assam, India", "2011-NA Maharashtra, India"="Maharashtra, India")) +labs(fill="Decade") +xlab("Location")+ylab("Age (Years)")






# Median by year  ---------------------------------------------------------

Age <- gs_read(ss=diph.sheet, ws="AgeMedian")
Age <- data.frame(Age)
Age <- Age[-c(1,3)]

ggplot(Age, aes(x=Decade, y=Median.Age.range, color=X.Region))+geom_point()



