# SGPSM
Similarity measures for generalized propensity score matching (SGPSM) is a propensity score method for multiple treatments matching. It can reduce confounding effects and provide unbiased average treatment effects for multiple treated groups, similar to other propensity score methods. This R code implements SGPSM for propensity score matching with (k >= 3) treated groups. 


# Usage

### Example:
```diff

data <- read.csv("~/data.csv")

source("~/SGPSM.R")

matched_data <- SGPSM(
  df = data,   # file name
  treat_group = "trt",   # column name of treatment
  gpsm = c("X0","X1","X2"),    # input column name of propensity scores vector
  algo = "manhattan",    # cosine, euclidean, manhattan
  caliper = 0.1   # criteria for selecting subjects in treatment A whose propensity score is "close" to that of a subject in treatment B
)

```

* df: &emsp; File name.<br>
* treat_group: &emsp; Variable for treatment groups.<br>
* gpsm: &emsp; Variables for generalized propensity score (input column name of propensity scores vector).<br>
* algo: &emsp; Algorithm of cosine, euclidean, manhattan.<br>
* caliper: &emsp; Caliper.<br>


# Dependencies
R package: dplyr >= 1.0.8, proxyC >= 0.2.4

# Papers
Yi-Chia Su, Chih-Chien Wu, Yi-Hau Chen, Pei-Ting Lee, **Chien-Chou Su**. *Similarity measures for generalized propensity score matching with three target therapies in cancer study: An idea from distance metrics in unsupervised learning* (2024; submission)

*Correspondent: Chien-Chou Su*
