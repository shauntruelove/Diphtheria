# Clean Data for R0 STAN model to Source



# LOAD PACKAGES ----------------------------------------------------------
if(!require('ggplot2')) install.packages('ggplot2'); library(ggplot2)
if(!require('lubridate')) install.packages('lubridate'); library(lubridate)
if(!require('dplyr')) install.packages('dplyr'); library(dplyr)
if(!require('pspline')) install.packages('pspline'); library(pspline)



get_STAN_data <- function(return.fake.test.data=FALSE){

    # FUNCTIONS ---------------------------------------------------------------
    
    # Function Name:    incid_week_to_day
    # Purpose:          Translate weekly case series into daily case series through fitting CDF of cases.
    incid_week_to_day <- function(week.cases) {
        
        require(pspline)
        
        # Set up week start and week number variables
        weeknum <- 0:(length(week.cases)-1)
        week.start <- weeknum*7
        
        # Get Empirical CDF of Data ( in long format)
        ecdf_fit <- ecdf(rep(week.start, week.cases))
        
        # Fit Spline to CDF
        x <- c(week.start, max(week.start)+7)
        y <- c(0, ecdf_fit(week.start))
        
        x_pred <- 0:(max(x)-1)
        #pspl <- smooth.Pspline(x, y, method=4, spar=.1) # requires pspline package
        pspl <- smooth.spline(x, y, spar=.1) 
        
        f1 <- predict(pspl, x_pred, deriv=1)
        f1 <- f1$y
        f1 <- f1 - min(f1) # normalize to 0 minimum
        f1 <- f1/sum(f1)
        
        # Impute daily case, on a weekly basis
        daily.cases <- data.frame(day=x_pred, cases=NA)
        for (w in 1:length(week.start)) {
            
            week.cases_ <- week.cases[w]
            week.days_ <- seq(week.start[w]+1, week.start[w]+7, by=1) # add 1 to make easier for for loop
            
            p_ <- f1[week.days_]
            cases.by.day <- round((p_/sum(p_)) * week.cases_)
            
            # If too many imputed cases, remove some
            if (sum(cases.by.day) > week.cases_) { 
                print(paste0(w, ', Too Many Imputed Cases: ', sum(cases.by.day)-week.cases_))
                tmp <- sample((1:7)[which(cases.by.day>0)], sum(cases.by.day)-week.cases_, replace=FALSE, prob=1-p_[which(cases.by.day>0)])
                cases.by.day[tmp] <- cases.by.day[tmp]-1
            }
            # If too few imputed cases, add some
            if (sum(cases.by.day) < week.cases_) { 
                print(paste0(w, ', Too Few Imputed Cases: ', week.cases_-sum(cases.by.day)))
                tmp <- sample(1:7, week.cases_-sum(cases.by.day), replace=FALSE, prob=p_)
                cases.by.day[tmp] <- cases.by.day[tmp]+1
            }
            
            daily.cases$cases[week.days_] <- cases.by.day
        }
        
        # Trim initial 0's
        daily.cases <- daily.cases[min(which(daily.cases$cases>0)):nrow(daily.cases),]
        daily.cases$day <- daily.cases$day - daily.cases$day[1]
        
    
        # Check that they are equal
        if ( sum(daily.cases$cases) == sum(week.cases) ) {
            return(daily.cases)
        } else {
            return('Failed -- Did not produce equal case totals.')   
        }
    }
    
    
    
    # Function Name:    incid_week_to_day_stoch
    # Purpose:          Stochastically mpute daily cases from weekly case series.
    incid_week_to_day_stoch <- function(week.cases) {
        
        # Set up week start and week number variables
        weeknum <- 0:(length(week.cases)-1)
        week.start <- weeknum*7
        
        daily.cases <- data.frame(day=0:(max(week.start)+6), cases=0)
        
        d <- sapply(X=1:length(week.start), FUN=function(day=X) { as.character(sample(size=week.cases[day], x=0:6, replace=T) + week.start[day]) } )
        d <- sort(as.integer(unlist(d)))
        hist.d<- hist(d, breaks=0:max(d), plot=FALSE)$counts
        daily.cases$cases[1:length(hist.d)] <- hist.d 
        
        # Trim initial 0's
        daily.cases <- daily.cases[min(which(daily.cases$cases>0)):nrow(daily.cases),]
        daily.cases$day <- daily.cases$day - daily.cases$day[1]
        
        return(daily.cases)
    }
    
    
    
    # Match the length between outbreak case series so they can easily be merged into a matrix
    match.l <- function(cases, length.needed=max(outbreak_T, na.rm=TRUE)+1) {
        if (length(cases)>length.needed) { 
            return(cases[1:length.needed])
        } else {
            return(c(cases, rep(0, length.needed-length(cases))))
        }
    }
    
    # Match the length between outbreak case series so they can easily be merged into a matrix
    match.l.matrix <- function(cases, length.needed=max(outbreak_T, na.rm=TRUE)+1) {
        if (ncol(cases)==length.needed) { 
            return(cases)
        } else if (ncol(cases)>length.needed) { 
            return(cases[,1:length.needed])
        } else {
            return(cbind(cases, matrix(0, ncol=length.needed-ncol(cases), nrow=nrow(cases))))
        }
    }
    
    
    
    # SET UP OBJECTS TO HOLD DATA & VARIABLES ---------------------------------
    outbreak_T <- NULL           # Max length of outbreak, in days
    outbreak_T_stoch <- NULL           # Max length of outbreak, in days
    
    
    # SET UP OUTBREAKS --------------------------------------------------------
    
    #.............................................
    # Daily data
    
    dailydata <- read.csv(file="data/outbreak_data_daily.csv", header = TRUE, stringsAsFactors = FALSE)
    row.names(dailydata) <- dailydata[,1]
    dailydata <- dailydata[,-1]
    
    # Extract matrix of cases
    start.row_ <- min(grep('Day', row.names(dailydata)))
    dailycases <- as.matrix(dailydata[-(1:(start.row_-1)),])
    dailycases <- apply(dailycases, 2, as.integer)
    row.names(dailycases) <- 0:(nrow(dailycases)-1)
    
    # Remove outbreaks to be skipped (i.e., without peak weeks, less than 15 cases)
    to_remove_ <- ( is.na(dailydata["Peak_day",]) | 
                    is.na(dailydata["use_schick",]) |
                    colSums(dailycases, na.rm=TRUE)<15 )
    dailycases <- dailycases[,!to_remove_]
    dailydata <- dailydata[,!to_remove_]               
    
    # Outbreak variables needed for STAN Model
    daily.vars <- list(outbreak_id=dailydata["Location",], 
                       outbreak_id2=paste0(1:ncol(dailydata), "-", dailydata["Location",]),
                       year.start=as.integer(dailydata["Year.start",]), 
                       use_schick=as.integer(dailydata["use_schick",]),
                       est_Vs_from_cases=as.integer(dailydata["est_Vs_from_cases",]),
                       est_Vs=as.integer(dailydata["est_Vs",]),
                       NoATB=as.integer(dailydata["NoATB",]),
                       vacc_mat=(dailydata[c("V_0","V_1&2","V_3"),]),
                       vacc_case_mat=(dailydata[c("V0_cases","V2_cases","V3_cases"),]),
                       peak.day=as.integer(dailydata["Peak_day",]))
    
    # Trim daily cases to peak days
    dailycases <- dailycases[1:(max(daily.vars$peak.day)+1), ]
    max_row_ <- nrow(dailycases)
    
    for (s in 1:ncol(dailycases)) {
        if ((daily.vars$peak.day[s]+2) > max_row_) next
        if ((daily.vars$peak.day[s]+2) == max_row_) {
            dailycases[max_row_, s] <- NA
        } else {
            dailycases[((daily.vars$peak.day[s]+2):max_row_), s] <- NA
        }    
    }
    
    daily.vars$outbreak_T <- daily.vars$peak.day
    
    # Remove extra rows with all NAs & Save max.T
    #dailycases <- dailycases[1:max(unlist(apply(X=dailycases, 2, FUN=function(x=X) which(!is.na(x))))), ]
    daily.vars$max.T <- max(as.integer(row.names(dailycases)))

    
    #.............................................
    # Weekly Data

    weeklydata <- read.csv(file="data/outbreak_cases_weekly.csv", header = TRUE, stringsAsFactors = FALSE)
    row.names(weeklydata) <- weeklydata[,1]
    weeklydata <- weeklydata[,-1]
    
    # Extract matrix of cases
    start.row_ <- grep('Week 0', row.names(weeklydata))
    weeklycases <- as.matrix(weeklydata[-(1:(start.row_-1)),])
    weeklycases <- apply(weeklycases, 2, as.integer)
    row.names(weeklycases) <- 0:(nrow(weeklycases)-1)

    # Remove outbreaks to be skipped (i.e., without peak weeks, less than 12 cases)
    to_remove_ <- ( is.na(weeklydata["Peak_week",]) | 
                        is.na(weeklydata["use_schick",]) |
                        colSums(weeklycases, na.rm=TRUE)<12 )
    weeklycases <- weeklycases[,!to_remove_]
    weeklydata <- weeklydata[,!to_remove_]               
    
    # Get outbreak variables
    weekly.vars <- list(outbreak_ids=as.character(weeklydata["Outbreak",]),
                        outbreak_ids2=paste0(1:ncol(weeklydata),"-",weeklydata["Outbreak",]),
                        year.start=as.integer(weeklydata["Year.start",]),
                        use_schick=as.integer(weeklydata["use_schick",]),
                        est_Vs_from_cases=as.integer(weeklydata["est_Vs_from_cases",]),
                        est_Vs=as.integer(weeklydata["est_Vs",]),
                        NoATB=as.integer(weeklydata["NoATB",]),
                        vacc_mat=(weeklydata[c("V_0","V_1&2","V_3"),]),
                        vacc_case_mat=(weeklydata[c("V0_cases","V2_cases","V3_cases"),]),
                        peak.week=as.integer(weeklydata["Peak_week",]))
    colnames(weeklycases) <- weekly.vars$outbreak_ids

    weeklycases <- weeklycases[1:max(unlist(apply(X=weeklycases, 2, FUN=function(x=X) which(!is.na(x))))), ]
    weekly.vars$max.week.T <- as.integer(max(nrow(weeklycases))-1)
    weekly.vars$max.T <- (weekly.vars$max.week.T*7)+6

    # Make case series for each outbreak
    weeklydata <- weeklycases
    weeklycases <- weeklycases_stoch <- weeklycases_all <- NULL
    #weeklycases <- 0:max.T
    for (s in 1:ncol(weeklydata)) {
        
        weeklycases_ <- as.integer(weeklydata[,s])
        weeklycases_ <- weeklycases_[!is.na(weeklycases_)]
        # Restrict to exponential growth phase
        weeklycases_ <- weeklycases_[1:(weekly.vars$peak.week[s]+1)]
        
        # Convert to daily case series
        weeklydaily_ <- incid_week_to_day(week.cases=weeklycases_) # using CDF
        weeklydaily.stoch_ <- incid_week_to_day_stoch(week.cases=weeklycases_) # stochastically
     
        weekday_ <- 0:(length(weeklycases_)-1)*7

        # Plot exponential growth only
        plot(weekday_, weeklycases_, col='red', type='l', main=weekly.vars$outbreak_ids[s])
        # Plot daily cases
        lines(weeklydaily_$day, weeklydaily_$cases, col='purple', type='l', lwd=1.5)
        lines(weeklydaily.stoch_$day, weeklydaily.stoch_$cases, col='green', type='l', lwd=1.5)
        
        # # cut off front of data before cases (again, for daily cases)
        # weeklydaily_ <- weeklydaily_$cases[min(which(weeklydaily_$cases>0)):nrow(weeklydaily_)]
        # weeklydaily.stoch_ <- weeklydaily.stoch_$cases[-(1:(min(which(weeklydaily.stoch_$cases>0))-1))]
        # Save max T for each school outbreak
        weekly.vars$outbreak_T[s] <- max(weeklydaily_$day[weeklydaily_$cases > 0])
        weekly.vars$outbreak_T_stoch[s] <- max(weeklydaily.stoch_$day[weeklydaily.stoch_$cases > 0])
        
        weeklycases_all <- cbind(weeklycases_all, 
                             as.integer(match.l(weeklydaily_$cases, length.needed=weekly.vars$max.T+1)), 
                             as.integer(match.l(weeklydaily.stoch_$cases, length.needed=weekly.vars$max.T+1)))
        weeklycases <- cbind(weeklycases, as.integer(match.l(weeklydaily_$cases, length.needed=weekly.vars$max.T+1))) 
        weeklycases_stoch <- cbind(weeklycases_stoch, as.integer(match.l(weeklydaily.stoch_$cases, length.needed=weekly.vars$max.T+1)))
    }
    colnames(weeklycases_all) <- c(paste0(c('cases_','cases_stoch_'), sort(rep((1:ncol(weeklydata)),2))))
    colnames(weeklycases) <- c(paste0(c('cases_'), (1:ncol(weeklydata))))
    colnames(weeklycases_stoch) <- c(paste0(c('cases_'), (1:ncol(weeklydata))))
    

    # End of Exponential growth period
    outbreak_T <- c(daily.vars$outbreak_T[!is.na(daily.vars$outbreak_T)], 
                    weekly.vars$outbreak_T[!is.na(weekly.vars$outbreak_T)])
    outbreak_T_stoch <- c(daily.vars$outbreak_T[!is.na(daily.vars$outbreak_T)], 
                          weekly.vars$outbreak_T_stoch[!is.na(weekly.vars$outbreak_T_stoch)])
    
    
    # FULL STAN Data ---------------------------------------------------------------

    full.cases <- rbind(match.l.matrix(cases=t(dailycases), length.needed=max(outbreak_T, na.rm=TRUE)+1), 
                        match.l.matrix(cases=t(weeklycases), length.needed=max(outbreak_T, na.rm=TRUE)+1))
    
    full_data <- list(M=as.integer(nrow(full.cases)),
                      T=as.integer(outbreak_T),
                      T_max=as.integer(max(outbreak_T)),
                      N_mat=as.matrix(full.cases[,-1]),
                      N0=as.integer(full.cases[,1]),
                      use_schick=as.integer(c(daily.vars$use_schick, weekly.vars$use_schick)),
                      est_Vs_from_cases=as.integer(c(daily.vars$est_Vs_from_cases, weekly.vars$est_Vs_from_cases)),
                      est_Vs=as.integer(c(daily.vars$est_Vs, weekly.vars$est_Vs)),
                      NoATB=as.integer(c(daily.vars$NoATB, weekly.vars$NoATB)),
                      Vk_pop=t(data.matrix(bind_cols(daily.vars$vacc_mat, weekly.vars$vacc_mat))),
                      Vk_cases=t(data.matrix(bind_cols(daily.vars$vacc_case_mat, weekly.vars$vacc_case_mat))))

    full_data$Vk_pop[is.na(full_data$Vk_pop)] <- 0
    full_data$Vk_cases[is.na(full_data$Vk_cases)] <- 0
    
    save(full_data, file='data/R0_STAN_data.RData')
    
    
    
    # FULL STAN Data - STOCHASTIC ---------------------------------------------------------------
    
    full.cases.stoch <- rbind(match.l.matrix(cases=t(dailycases), length.needed=max(outbreak_T_stoch, na.rm=TRUE)+1), 
                        match.l.matrix(cases=t(weeklycases_stoch), length.needed=max(outbreak_T_stoch, na.rm=TRUE)+1))
    
    full_data_stoch <- list(M=as.integer(nrow(full.cases.stoch)),
                      T=as.integer(outbreak_T_stoch),
                      T_max=as.integer(max(outbreak_T_stoch)),
                      N_mat=as.matrix(full.cases.stoch[,-1]),
                      N0=as.integer(full.cases.stoch[,1]),
                      N_R0pred=100,
                      use_schick=as.integer(c(daily.vars$use_schick, weekly.vars$use_schick)),
                      est_Vs_from_cases=as.integer(c(daily.vars$est_Vs_from_cases, weekly.vars$est_Vs_from_cases)),
                      est_Vs=as.integer(c(daily.vars$est_Vs, weekly.vars$est_Vs)),
                      NoATB=as.integer(c(daily.vars$NoATB, weekly.vars$NoATB)),
                      Vk_pop=t(data.matrix(bind_cols(daily.vars$vacc_mat, weekly.vars$vacc_mat))),
                      Vk_cases=t(data.matrix(bind_cols(daily.vars$vacc_case_mat, weekly.vars$vacc_case_mat))))
    
    full_data_stoch$Vk_pop[is.na(full_data_stoch$Vk_pop)] <- 0
    full_data_stoch$Vk_cases[is.na(full_data_stoch$Vk_cases)] <- 0
    
    
    save(full_data_stoch, file='data/R0_STAN_data_stoch.RData')
    
    
    if (return.fake.test.data){
        return(list(full_data_test, full_data, full_data_stoch))
    } else {
        return(list(full_data, full_data_stoch))
    }
}


