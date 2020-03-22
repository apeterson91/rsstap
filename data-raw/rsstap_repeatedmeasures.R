## code to prepare `bbnet_repeatedmeasures` dataset goes here

library(dplyr)
library(tidyr)

set.seed(1214)
num_subj <- 5.5E2
num_bef <- 55


alpha <- 26.5
Z <- rbinom(n = num_subj, prob = .45, size = 1)
subject_intercepts <- rnorm(num_subj, mean = 0 ,sd= 1)
subject_slopes <- rnorm(num_subj,mean = 0,sd=.5)
time_effect <- 1
delta <- -.8
sigma <- 1

true_direct_effect <- function(x){ ((-.5 + x -.5 *x^2))*(x<=1) }

sub_df <- tibble(x = runif(n = num_subj, min = 0, max = 5.0),
                 y = runif(n = num_subj, min = 0, max = 5.0),
                 class = "Subject")



## HFS's
bef_df <- tibble(x = runif(n = num_bef, min = 0, max = 5.0),
                 y = runif(n = num_bef, min = 0, max = 5.0),
                 class = "Healthy Food Stores")

direct_dists <- fields::rdist(as.matrix(sub_df[,1:2]),
                              as.matrix(bef_df[,1:2]))
subject_bef_df <- as_tibble(direct_dists) %>% 
  mutate(subj_ID = 1:num_subj) %>% 
  gather(contains("V"),key="BEF",value="Distance") %>% 
  mutate(BEF = "HFS",measure_ID = 0L)

X_1 <- rowSums(true_direct_effect(direct_dists))

y <- alpha + Z*delta + X_1 + subject_intercepts + rnorm(n = num_subj, mean = 0, sd = sigma) ## timepoint 0


## move ten people
movers <- sample(1:num_subj,10)

sub_df[movers,c("x","y")] <- cbind(runif(n = 10,min = 0, max = 5.0),
                                   runif(n = 10, min = 0, max = 5.0))

direct_dists <- fields::rdist(as.matrix(sub_df[,1:2]),
                              as.matrix(bef_df[,1:2]))

X_1 <- rowSums(true_direct_effect(direct_dists))

y_2 <- alpha + Z*delta + X_1 + subject_intercepts + time_effect + subject_slopes +  rnorm(n = num_subj, mean = 0, sd = sigma) ## timepoint 1


movers <- sample(1:num_subj,10)

sub_df[movers,c("x","y")] <- cbind(runif(n = 10,min = 0, max = 5.0),
                                   runif(n = 10, min = 0, max = 5.0))

direct_dists <- fields::rdist(as.matrix(sub_df[,1:2]),
                              as.matrix(bef_df[,1:2]))

X_1 <- rowSums(true_direct_effect(direct_dists))

subject_bef_df <- subject_bef_df %>% 
  rbind(.,as_tibble(direct_dists) %>% 
          mutate(subj_ID = 1:num_subj) %>% 
          gather(contains("V"),key="BEF",value="Distance") %>% 
          mutate(BEF = "HFS",measure_ID = 1L))

y_3 <- alpha + Z*delta + X_1 + subject_intercepts + time_effect*2 + subject_slopes*2 +  rnorm(n = num_subj, mean = 0, sd = sigma) ## timepoint 2


movers <- sample(1:num_subj,10)

sub_df[movers,c("x","y")] <- cbind(runif(n = 10,min = 0, max = 5.0),
                                   runif(n = 10, min = 0, max = 5.0))

direct_dists <- fields::rdist(as.matrix(sub_df[,1:2]),
                              as.matrix(bef_df[,1:2]))

X_1 <- rowSums(true_direct_effect(direct_dists))

subject_bef_df <- subject_bef_df %>% 
  rbind(.,as_tibble(direct_dists) %>% 
          mutate(subj_ID = 1:num_subj) %>% 
          gather(contains("V"),key="BEF",value="Distance") %>% 
          mutate(BEF = "HFS",measure_ID = 2L))

y_4 <- alpha + Z*delta + X_1 + subject_intercepts + time_effect*3 + subject_slopes*3 +  rnorm(n = num_subj, mean = 0, sd = sigma) ## timepoint 3



subject_df <- tibble(subj_ID = 1:num_subj,
                     BMI = y,
                     sex = Z,
                     time = 0,
                     measure_ID = 0L) %>% 
  rbind(.,tibble(subj_ID = 1:num_subj,
                 BMI = y_2,
                 sex = Z,
                 time = 1,
                 measure_ID = 1L)) %>% 
  rbind(.,tibble(subj_ID = 1:num_subj,
                 BMI = y_3,
                 sex = Z, 
                 time = 2, 
                 measure_ID = 2L)) %>% 
  rbind(.,tibble(subj_ID = 1:num_subj,
                 BMI = y_4,
                 sex = Z,
                 time = 3,
                 measure_ID = 3L))

subject_bef_df <- subject_bef_df %>% 
  rbind(.,as_tibble(direct_dists) %>% 
  mutate(subj_ID = 1:num_subj) %>% 
  gather(contains("V"),key="BEF",value="Distance") %>% 
  mutate(BEF = "HFS",measure_ID = 3L))

rsstap_repeatedmeasures <- list(distance_df = subject_bef_df %>% filter(Distance<=2),
                                subject_df = subject_df,
                                exposure_function = true_direct_effect)

usethis::use_data(rsstap_repeatedmeasures, overwrite = TRUE)
