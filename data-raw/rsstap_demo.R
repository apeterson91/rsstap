## code to prepare `rsstap_demo` dataset goes here


# Libraries ---------------------------------------------------------------


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


sub_df <- dplyr::tibble(x = runif(n = num_subj, min = 0, max = 5.0),
                 y = runif(n = num_subj, min = 0, max = 5.0),
                 class = "Subject")
## HFS's
bef_df <- dplyr::tibble(x = runif(n = num_bef, min = 0, max = 5.0),
                 y = runif(n = num_bef, min = 0, max = 5.0),
                 class = "Healthy Food Stores")



direct_dists <- fields::rdist(as.matrix(sub_df[,1:2]),
                       as.matrix(bef_df[,1:2]))

X_1 <- rowSums(true_direct_effect(direct_dists))



indirect_dists <- fields::rdist(as.matrix(bef_df[,1:2]),
                                as.matrix(bef_df[,1:2]))



y <- alpha + Z*delta + X_1 + rnorm(n = num_subj, mean = 0, sd = sigma)


subject_data <- dplyr::tibble(subj_id = 1:num_subj,
                           y = y,
                           sex = factor(Z,labels=c("M","F")))

subject_distance_data <- fields::rdist(as.matrix(sub_df[,1:2]),
                                       as.matrix(bef_df[,1:2])) %>%
  dplyr::as_tibble() %>%
  dplyr::mutate(subj_id = 1:num_subj) %>% 
  tidyr::gather(dplyr::contains("V"),key = 'BEF',value = 'Distance') %>% 
  dplyr::mutate(BEF = 'HFS')




rsstap_demo <- list(subject_data = subject_data,
                  direct_distance_data = subject_distance_data,
                  alpha = alpha,
                  delta = delta,
                  sigma = sigma,
                  true_effect = true_direct_effect,
                  subjectbef_positions = rbind(sub_df,bef_df)
                  )



usethis::use_data(rsstap_demo, overwrite = TRUE)
