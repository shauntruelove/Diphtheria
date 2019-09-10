
load_VE_long <- function(){
    
    VE.raw <- read.csv(file="./data/VE_data_all.csv", header = TRUE, stringsAsFactors = FALSE)
    
    # Case Control Method -----------------------------------------------------
    
    VE.cc <- VE.raw[!(is.na(VE.raw$Cases_total) | is.na(VE.raw$Controls)),]
    VE.cc <- VE.cc[,c('Study','Location', "Year", 'Dose',"Age.min","Age.max","Cases_total","Controls")]
    VE.cc$Population <- factor(paste(VE.cc$Location, VE.cc$Year,sep=" "))
    
    VE.cc <- VE.cc[!(VE.cc$Cases_total==0 & VE.cc$Controls==0),]
    VE.cc$Study <- factor(VE.cc$Study)
    # Group Doses into 3 categories: Full (3+), Partial (1-2), or None
    VE.cc$Dose013 <- NA
    VE.cc <- VE.cc[!is.na(VE.cc$Dose),] # remove NAs. These should be fixed and not exist
    VE.cc <- VE.cc %>% mutate(Dose013 = ifelse(Dose=="0", "None", 
                                               ifelse(Dose=="1" | Dose=="2" | Dose=="1-2", "Partial",    "Full")))
    VE.cc$Dose013 <- factor(VE.cc$Dose013, levels=c("None","Partial","Full"))

    
    # Group Age into 3 categories: 0-4y, 5-19y, 20+y
    VE.cc$AgeGroup3 <- NA
    VE.cc <- VE.cc %>% mutate(AgeGroup3 = ifelse((!is.na(Age.max) & Age.max<=4), "0-4y",
                                                 ifelse((!is.na(Age.max) & Age.min>=5 & Age.max<=19), "5-19y",
                                                        ifelse((!is.na(Age.max) & Age.min>=20), "20+", NA))))
    VE.cc$AgeGroup3 <- factor(VE.cc$AgeGroup3, levels = c("0-4y","5-19y","20+"))
    
    VE.cc <- VE.cc %>% mutate(row = 1:nrow(VE.cc))
    VE.cc.long <- VE.cc %>% group_by(row) %>%
        do({ left_join(data_frame(row=.$row, case = c(rep(1,.$Cases_total), rep(0,.$Controls))),.,by='row') })

    return(VE.cc.long)
}
    