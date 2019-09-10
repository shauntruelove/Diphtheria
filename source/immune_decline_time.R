
# Immunity Decline Over Time



# Packages and Source Code ------------------------------------------------

if(!require('dplyr')) install.packages('dplyr'); library(dplyr)
if(!require('lme4')) install.packages('lme4'); library(lme4)



# LOAD AND CLEAN DATA -----------------------------------------------------

dat.ab.time <- read.csv(file="./data/immune_decline_time.csv", header = TRUE, stringsAsFactors = FALSE)

dat.ab.time <- dat.ab.time[-which(dat.ab.time$Time_since_vacc_month=='Pre'),]
dat.ab.time$Time_since_vacc_month <- as.integer(dat.ab.time$Time_since_vacc_month)
dat.ab.time$Time_since_vacc_year[is.na(dat.ab.time$Time_since_vacc_year)] <- round(dat.ab.time$Time_since_vacc_month/12,2)[is.na(dat.ab.time$Time_since_vacc_year)]
dat.ab.time$Time_since_vacc_year <- as.integer(dat.ab.time$Time_since_vacc_year)
dat.ab.time <- dat.ab.time[!is.na(dat.ab.time$Time_since_vacc_year),]
dat.ab.time$pop_id <- paste(dat.ab.time$Study, dat.ab.time$Location, sep="-")
dat.ab.time <- dat.ab.time[!is.na(dat.ab.time$prop_protected),]


# Individual data ---------------------------------------------------------

dat.ab.time.long <-  dat.ab.time %>% mutate(row = 1:nrow(dat.ab.time)) %>%
    group_by(row) %>%
    do({
        left_join(data_frame(row=.$row, protected=c(rep(1,.$no.protected), rep(0, .$N-.$no.protected))),., by='row')
    })

# Aggregate data back to groups

dat.ab.time.long$time_year <- round(dat.ab.time.long$Time_since_vacc_year,0)
dat.ab.time.long <- dat.ab.time.long %>% mutate( year_group = cut(time_year, breaks = seq(0,100,10), labels = F))
dat.ab.time.long$year_group <- dat.ab.time.long$year_group*10 - 5
dat.ab.time.long$year_group[dat.ab.time.long$time_year<=10] <- dat.ab.time.long$time_year[dat.ab.time.long$time_year<=10]


dat.ab.time.grouped <- dat.ab.time.long %>% group_by(Ref_index, Location, Study, Year, year_group) %>%
    summarise(n.protected=round(sum(protected, na.rm=TRUE),0), n.suscept=round(sum(protected==0, na.rm=TRUE),0)) %>% as.data.frame
dat.ab.time.grouped$prop.protected <- dat.ab.time.grouped$n.protected / (dat.ab.time.grouped$n.protected + dat.ab.time.grouped$n.suscept)
dat.ab.time.grouped$prop.protected[dat.ab.time.grouped$prop.protected==0] <- 0.001

# Linear model to estimate decline over time.
mod <- lmer(log(prop.protected) ~ year_group + (1|Study), data=dat.ab.time.grouped)
cbind("% increase per age yr" = fixef(mod), confint(mod, method="Wald", level = 0.95))


# How many and Which Studies
length(unique(dat.ab.time.grouped$Study))
length(unique(dat.ab.time.grouped$Ref_index))
unique(dat.ab.time.grouped$Ref_index)
unique(dat.ab.time.grouped$Study)

