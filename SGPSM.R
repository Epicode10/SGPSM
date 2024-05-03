####################################################################################
# Description: generalized propensity score matching for multiple treatment groups
# Algorithm: similarity measures for generalized propensity score matching (SGPSM) 
# Version: 1.0
# Lab version: beta 1.0
# Created: 2024/3/27
####################################################################################



SGPSM <- function(
  df,
  treat_group,
  gpsm,
  algo,
  caliper){



################################################################################
# require packages 
################################################################################

library(proxyC)
library(dplyr)

################################################################################
# data preparation 
################################################################################

row_data <- df

row_data$PID <- 1:nrow(row_data)

df$PID <- 1:nrow(df)

df <- df[,append(append("PID", treat_group), gpsm)]

# number of treatment groups

ng_df <- data.frame(table(df[,treat_group]))
ng_df$Var1 <- as.integer(as.character(ng_df$Var1)) 

ng <- nrow(ng_df)

min_trt <- ng_df$Var1[ng_df$Freq == min(ng_df$Freq)]

n_trt <- min(ng_df$Freq)

################################################################################  
# similarity measures for GPSM (function)
################################################################################

GPSM_cos_simil <- function(x){

  set.seed(0)
  
  
  ### similarity measures
  
  vec <- df[,gpsm]
  vec <- data.matrix(vec)
  
  dist_mx <- proxyC::simil(vec, method = algo)
  dist_mx <- data.matrix(dist_mx)
  dist_mx <- data.frame(dist_mx)
  dist_mx$PID <- as.numeric(rownames(dist_mx))
  
  label_df <- df[,c("PID", treat_group)]
  names(label_df) <- c("PID", "group")
  
  dist_mx <- merge(label_df , dist_mx, by="PID")
  
  
  ### GPSM: base matched set
  
  matched_cases <- data.frame()
  
  base_mx <- dist_mx[dist_mx$group == min_trt,]
  base_mx_pid <- paste0("X", rownames(base_mx))
  
  compare_mx <- dist_mx[dist_mx$group != min_trt, append( c("PID", "group") ,base_mx_pid)]
  
  for(i in 1:n_trt){
    cat("\r","Process rate: ", round(i/n_trt*100,2), '%     ') # process rate
    
    base_arm_match <- base_mx[i, c("PID", "group")] %>%
      mutate(dist_measure = 1) %>%
      mutate(matched = i)
    
    
    compare_arm_match <- compare_mx %>%
      group_by(group) %>%
      dplyr::select(PID, group, base_mx_pid[i]) %>%
      rename(dist_measure = base_mx_pid[i]) %>%
      filter(dist_measure >= (1-caliper)) %>%
      filter(dist_measure == max(dist_measure)) %>%
      unique() %>%
      mutate(matched = i)
    
    
    ### GPSM: check other arms
    
    check_n_group <- nrow(compare_arm_match)
    
    if(check_n_group != (ng-1)){
      check_index <- 1
    }else{
      check_mx <- dist_mx[dist_mx$PID %in% compare_arm_match$PID, paste0("X", compare_arm_match$PID)] 
      
      check_vec <- c()
      for(j in 1:nrow(check_mx)){
        check_vec[j] <- sum(check_mx[j,] < (1-caliper))
      }
      check_index <- sum(check_vec)    
    }
    
    
    
    ### matched cases
    
    if( (check_n_group == (ng-1)) & (check_index == 0) ){
      matched_cases <- rbind(matched_cases, base_arm_match, compare_arm_match)
      compare_mx <- compare_mx[!(compare_mx$PID %in% compare_arm_match$PID),]
    }else{
      next
    }
  }
  
  matched_data <- merge(row_data, matched_cases[,c("PID", "dist_measure", "matched")], by="PID")
  
  
return(matched_data) 

}


###  


GPSM_euc_dist <- function(x){
  
  set.seed(0)
  
  
  ### similarity measures
  
  vec <- df[,gpsm]
  vec <- data.matrix(vec)
  
  dist_mx <- proxyC::dist(vec, method = algo)
  dist_mx <- data.matrix(dist_mx)
  dist_mx <- data.frame(dist_mx)
  dist_mx$PID <- as.numeric(rownames(dist_mx))
  
  label_df <- df[,c("PID", treat_group)]
  names(label_df) <- c("PID", "group")
  
  dist_mx <- merge(label_df , dist_mx, by="PID")
  
  
  ### GPSM: base matched set
  
  matched_cases <- data.frame()
  
  base_mx <- dist_mx[dist_mx$group == min_trt,]
  base_mx_pid <- paste0("X", rownames(base_mx))
  
  compare_mx <- dist_mx[dist_mx$group != min_trt, append( c("PID", "group") ,base_mx_pid)]
  
  for(i in 1:n_trt){
    cat("\r","Process rate: ", round(i/n_trt*100,2), '%     ') # process rate
    
    base_arm_match <- base_mx[i, c("PID", "group")] %>%
      mutate(dist_measure = 1) %>%
      mutate(matched = i)
    
    
    compare_arm_match <- compare_mx %>%
      group_by(group) %>%
      dplyr::select(PID, group, base_mx_pid[i]) %>%
      rename(dist_measure = base_mx_pid[i]) %>%
      filter(dist_measure <= (caliper)) %>%
      filter(dist_measure == min(dist_measure)) %>%
      unique() %>%
      mutate(matched = i)
    
    
    ### GPSM: check other arms
    
    check_n_group <- nrow(compare_arm_match)
    
    if(check_n_group != (ng-1)){
      check_index <- 1
    }else{
      check_mx <- dist_mx[dist_mx$PID %in% compare_arm_match$PID, paste0("X", compare_arm_match$PID)] 
      
      check_vec <- c()
      for(j in 1:nrow(check_mx)){
        check_vec[j] <- sum(check_mx[j,] > (caliper))
      }
      check_index <- sum(check_vec)    
    }
    
    
    
    ### matched cases
    
    if( (check_n_group == (ng-1)) & (check_index == 0) ){
      matched_cases <- rbind(matched_cases, base_arm_match, compare_arm_match)
      compare_mx <- compare_mx[!(compare_mx$PID %in% compare_arm_match$PID),]
    }else{
      next
    }
  }
  
  matched_data <- merge(row_data, matched_cases[,c("PID", "dist_measure", "matched")], by="PID")
  
  
  return(matched_data) 
  
}


 

###  


GPSM_man_dist <- function(x){
  
  set.seed(0)
  
  
  ### similarity measures
  
  vec <- df[,gpsm]
  vec <- data.matrix(vec)
  
  dist_mx <- proxyC::dist(vec, method = algo)
  dist_mx <- data.matrix(dist_mx)
  dist_mx <- data.frame(dist_mx)
  dist_mx$PID <- as.numeric(rownames(dist_mx))
  
  label_df <- df[,c("PID", treat_group)]
  names(label_df) <- c("PID", "group")
  
  dist_mx <- merge(label_df , dist_mx, by="PID")
  
  
  ### GPSM: base matched set
  
  matched_cases <- data.frame()
  
  base_mx <- dist_mx[dist_mx$group == min_trt,]
  base_mx_pid <- paste0("X", rownames(base_mx))
  
  compare_mx <- dist_mx[dist_mx$group != min_trt, append( c("PID", "group") ,base_mx_pid)]
  
  for(i in 1:n_trt){
    cat("\r","Process rate: ", round(i/n_trt*100,2), '%     ') # process rate
    
    base_arm_match <- base_mx[i, c("PID", "group")] %>%
      mutate(dist_measure = 1) %>%
      mutate(matched = i)
    
    
    compare_arm_match <- compare_mx %>%
      group_by(group) %>%
      dplyr::select(PID, group, base_mx_pid[i]) %>%
      rename(dist_measure = base_mx_pid[i]) %>%
      filter(dist_measure <= (caliper)) %>%
      filter(dist_measure == min(dist_measure)) %>%
      unique() %>%
      mutate(matched = i)
    
    
    ### GPSM: check other arms
    
    check_n_group <- nrow(compare_arm_match)
    
    if(check_n_group != (ng-1)){
      check_index <- 1
    }else{
      check_mx <- dist_mx[dist_mx$PID %in% compare_arm_match$PID, paste0("X", compare_arm_match$PID)] 
      
      check_vec <- c()
      for(j in 1:nrow(check_mx)){
        check_vec[j] <- sum(check_mx[j,] > (caliper))
      }
      check_index <- sum(check_vec)    
    }
    
    
    
    ### matched cases
    
    if( (check_n_group == (ng-1)) & (check_index == 0) ){
      matched_cases <- rbind(matched_cases, base_arm_match, compare_arm_match)
      compare_mx <- compare_mx[!(compare_mx$PID %in% compare_arm_match$PID),]
    }else{
      next
    }
  }
  
  matched_data <- merge(row_data, matched_cases[,c("PID", "dist_measure", "matched")], by="PID")
  
  
  return(matched_data) 
  
}




################################################################################
# similarity measures for GPSM (perform analysis) 
################################################################################

if(algo == "cosine"){
  matched_data <- GPSM_cos_simil(x)
}else if(algo == "euclidean"){
  matched_data <- GPSM_euc_dist(x)
}else{
  matched_data <- GPSM_man_dist(x)
}


return(matched_data)
}











