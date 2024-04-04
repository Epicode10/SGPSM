####################################################################################
# Description: generalized propensity score matching for multiple treatment groups
# Algorithm: similarity measures for generalized propensity score matching (SGPSM) 
# Version: 1.0
# Lab version: beta 1.0
# Created: 2024/3/27
####################################################################################


data <- read.csv("~/data.csv")

source("~/SGPSM.R")

matched_data <- SGPSM(
  df = data, # file name
  treat_group = "trt", # column name of treatment
  gpsm = c("X0","X1","X2"),  # input column name of propensity scores vector
  algo = "manhattan",  # cosine, euclidean, manhattan
  caliper = 0.1
)

