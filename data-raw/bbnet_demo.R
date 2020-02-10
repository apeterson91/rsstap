## code to prepare `bnet_demo` dataset goes here


# Libraries ---------------------------------------------------------------


library(dplyr)
library(tidyr)

# Initialize RNG, base parameters -----------------------------------------


set.seed(1214)
num_subj <- 5.5E2
num_bef <- 55


alpha <- 26.5
Z <- rbinom(n = num_subj, prob = .45, size = 1)
delta <- -.8
sigma <- 1

true_direct_effect <- function(x){ ((-.5 + x -.5 *x^2)/3)*(x<=1) }

# d <- seq(from = 0, to = 1,by=0.01)
# plot(d,true_direct_effect(d),type='l')

# true_indirect_effect <- function(x){ (-.25 + x - x^2)/20}

# Generate Position Locations ---------------------------------------------


sub_df <- tibble(x = runif(n = num_subj, min = 0, max = 5.0),
                 y = runif(n = num_subj, min = 0, max = 5.0))
## HFS's
bef_df <- tibble(x = runif(n = num_bef, min = 0, max = 5.0),
                 y = runif(n = num_bef, min = 0, max = 5.0))



direct_dists <- fields::rdist(as.matrix(sub_df[,1:2]),
                       as.matrix(bef_df[,1:2]))

X_1 <- rowSums(true_direct_effect(direct_dists))



indirect_dists <- fields::rdist(as.matrix(bef_df[,1:2]),
                                as.matrix(bef_df[,1:2]))



y <- alpha + Z*delta + X_1 + rnorm(n = num_subj, mean = 0, sd = sigma)


subject_data <- tibble(subj_id = 1:num_subj,
                           y = y,
                           sex = factor(Z,labels=c("M","F")))

subject_distance_data <- fields::rdist(as.matrix(sub_df[,1:2]),
                                       as.matrix(bef_df[,1:2])) %>%
  as_tibble() %>%
  mutate(subj_id = 1:num_subj) %>% 
  gather(contains("V"),key = 'BEF',value = 'Distance') %>% 
  mutate(BEF = 'HFS')

bef_bef_distance_data <- fields::rdist(as.matrix(bef_df[,1:2]),
                                       as.matrix(bef_df[,1:2])) %>% as_tibble() %>%
  mutate(from_bef = 1:num_bef) %>%
  gather(contains("V"),key = "to_bef",value = "Distance") %>%
  mutate(to_bef = as.integer(stringr::str_replace(from_bef,"V","")) )





bbnet_demo <- list(subject_data = subject_data,
                  direct_distance_data = subject_distance_data,
                  indirect_distance_data = bef_bef_distance_data,
                  alpha = alpha,
                  beta=beta,
                  delta = delta,
                  sigma = sigma,
                  true_effect = true_direct_effect)



usethis::use_data(bbnet_demo, overwrite = TRUE)
