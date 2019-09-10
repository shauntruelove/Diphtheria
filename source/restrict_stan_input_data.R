

# Function to Restrict data

restrict_stan_input_data <- function(data=full_data, outbreaks_remove=1){
  vectors_to_change <- which(lapply(data, length)==length(data$T))
  mats <- lapply(data, nrow)
  mats_to_change <- which(!is.na(as.integer(lapply(X=mats, FUN=function(x=X) which(x==length(data$T))))))

  for (i in vectors_to_change){
    data[[i]] <- data[[i]][-outbreaks_remove]
  }
  for (r in mats_to_change){
    data[[r]] <- data[[r]][-outbreaks_remove,]
  }
  # check that max T is not being changed
  data$T_max <- max(data$T)
  data$M <- length(data$T)
  data$N_mat <- data$N_mat[,1:data$T_max]
  
  return(data)
}



# EXAMPLE RUN:
#  - To exclude the Rohingya outbreak
# full_data_noRoh <- restrict_stan_input_data(data=full_data, outbreaks_remove=c(1))

