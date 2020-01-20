## code to prepare `bnet_demo` dataset goes here


# Libraries ---------------------------------------------------------------


library(dplyr)
library(tidyr)

# Initialize RNG, base parameters -----------------------------------------


set.seed(1214)
num_subj <- 2.5E2
num_bef <- 55


alpha <- 26.5
Z <- rbinom(n = num_subj, prob = .45, size = 1)
delta <- -.8
sigma <- 4.3

true_direct_effect <- function(x){ (-.5 + x -.5 *x^2)/20 }

true_indirect_effect <- function(x){ (-.25 + x - x^2)/20}

# Generate Position Locations ---------------------------------------------


sub_df <- data_frame(x = runif(n = num_subj, min = 1, max = 2.0),
                     y = runif(n = num_subj, min = 1, max = 2.0))
## HFS's
bef_df <- data_frame(x = runif(n = num_bef, min = 0, max = 3.0),
                     y = runif(n = num_bef, min = 0, max = 3.0))



direct_dists <- fields::rdist(as.matrix(sub_df[,1:2]),
                       as.matrix(bef_df[,1:2]))

## create spatial scale of 1
direct_dists <- Matrix::Matrix((direct_dists<=1)*direct_dists, sparse = TRUE)

X_1 <- true_direct_effect(direct_dists) %*% rep(1,55)

nonzero_ics <- apply(direct_dists,1,function(x) which(x>0))

indirect_dists <- fields::rdist(as.matrix(bef_df[,1:2]),
                                as.matrix(bef_df[,1:2]))

indirect_dists <- Matrix::Matrix((indirect_dists<=1)*indirect_dists,sparse=TRUE)

all(Matrix::t(indirect_dists) == indirect_dists) ## symmetric
indirect_dists[lower.tri(indirect_dists)] <- 0
all(Matrix::t(indirect_dists) == indirect_dists) ## not symmetric; only upper triangle is nonzero
X_2 <- purrr::map_dbl(nonzero_ics,function(x){
  sum(true_indirect_effect(indirect_dists[x,x]) %*% rep(1,length(x)))
})



y <- as.numeric(alpha + Z*delta + X_1 + X_2 + rnorm(n = num_subj, mean = 0, sd = sigma))


subject_data <- data_frame(subj_id = 1:num_subj,
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


# Create-and-combine distance basis matrices  -----------------------------

## set knots along decile of included distances distribution
knots <- quantile(subject_distance_data %>% filter(Distance<=2) %>% pull(Distance),probs = seq(from=0.1,to=.9,by=0.1))

## using b-splines here, but could use something else
B_list <- subject_distance_data %>% filter(Distance<=2) %>% split(.$subj_id) %>% 
  purrr::map(.,function(x) splines::bs(x %>% pull(Distance),knots = knots[2:10],Boundary.knots = c(0,2)) )

max_q <- subject_distance_data %>% filter(Distance<=2) %>% 
  group_by(subj_id) %>% count() %>%  ungroup() %>% summarise(max = max(n)) %>% pull(max)

B_mat <- purrr::map(B_list,function(x){
  attributes(x) <- attributes(x)["dim"]
  z <- purrr::map(1:ncol(x),function(y){
    c(x[,y],rep(0,max_q - nrow(x)))
  })
  matrix(unlist(z),nrow=1)
})

B <- Matrix::Matrix(do.call(rbind,B_mat),sparse=TRUE) ## dimension n x num_basis * max_num_distances
D <- cbind(1,B) ## add intercept


## eleven betas from knots, 1 for intercept = 12 beta's total
## create binary matrix to transform 12 x 1 beta vector--> 606 x 1 vector with appropriately placed betas
T_mat <- matrix(0,nrow=ncol(D),ncol=ncol(B_list[[1]])+1)
T_mat[1,1] <- 1
end <- 1
for(col_ix in 2:ncol(T)){
  start <- end+1
  end <- start + max_q -1
  m <- matrix(0,nrow=ncol(D),ncol=ncol(B_list[[1]])+1)
  m[start:end,col_ix] <- 1
  T_mat <- T_mat + m
}



bnet_demo <- list(subject_data = subject_data,
                  direct_distance_data = subject_distance_data,
                  indirect_distance_data = bef_bef_distance_data,
                  D = D,
                  T_mat = T_mat)



usethis::use_data(bnet_demo, overwrite = TRUE)
