

# Vaccine Effectiveness - Against Colonization ----------------------------


# Data from Miller et al. 1972
data <- data.frame(Group=c('Symptomatic Inf','Asymptomatic Inf','No Infection'), 
                   Full=c(2,71,132), Partial=c(10,18,59), None=c(3,0,7))

# RR by colonization - Full Vacc vs None
a <- 2 + 71
b <- 132
c <- 3 + 0
d <- 7

(RR <- (a / (a+b)) / (c / (c+d)))
1 - RR # VE

# CI of RR
RR.ci <- c(exp(log(RR) - 1.96 * sqrt( (b/a)/(a+b) + (c/d)/(c+d))), 
           exp(log(RR) + 1.96 * sqrt( (b/a)/(a+b) + (c/d)/(c+d))))
# CI of VE
1 - RR.ci[c(2,1)]



# RR by colonization - Partial Vacc vs None -------------------------------
a <- 10 + 18
b <- 59
c <- 3 + 0 
d <- 7

(RR <- (a / (a+b)) / (c / (c+d)))
1 - RR # VE

# CI of RR
RR.ci <- c(exp(log(RR) - 1.96 * sqrt( (b/a)/(a+b) + (c/d)/(c+d))), 
           exp(log(RR) + 1.96 * sqrt( (b/a)/(a+b) + (c/d)/(c+d))))
# CI of VE
1 - RR.ci[c(2,1)]

