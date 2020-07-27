## code to prepare `poisson_benvo` dataset goes here

num_subj <- 300
sex <- rbinom(n = num_subj,size=1,prob = .5)
centered_income <- rlnorm(n=num_subj,mean= 0,log(1.5))
hasFFR_exp <- rbinom(num_subj,size = 1,prob = .9)
FFRcnt <- rpois(num_subj,10)*hasFFR_exp
ldists <- lapply(FFRcnt,function(x) runif(x) )
f <- function(x) .1*pweibull(x,shape=5,scale=.6,lower.tail = F)
FFRexposure <- sapply(ldists,function(x) sum(f(x)))
hasHFS_exp <- rbinom(num_subj,size=1,prob=.8)
HFScnt <- (rpois(num_subj,5)+1)*hasHFS_exp
HFSdists <- lapply(HFScnt,function(x) runif(x) )
f <- function(x) -.3*pweibull(x,shape=3.5,scale=.7,lower.tail = F)
HFSexposure <- sapply(HFSdists,function(x) sum(f(x)))
eta <- 2 + sex*1 + centered_income*-.1 + HFSexposure + FFRexposure
y <- rpois(n = num_subj,lambda = poisson()$linkinv(eta))



subj_df <- dplyr::tibble(ID= 1:num_subj,
                         NumObese = y,
                         sex = sex,
                         Centered_Scaled_Income = centered_income)

FFR <- purrr::map2_dfr(1:length(ldists),ldists,function(x,y) dplyr::tibble(ID=x,Distance=y))
HFS <- purrr::map2_dfr(1:length(HFSdists),HFSdists,function(x,y) dplyr::tibble(ID=x,Distance=y))

poisson_benvo <- benvo(subject_data = subj_df,
                        bef_data = list(FFR,HFS),
                        bef_names = c("FFR","HFS"),
                        distance_col = c("Distance","Distance"),
                        joining_id = c("ID"),
                        exposed_time_col = rep(NA,2))

usethis::use_data(poisson_benvo, overwrite = TRUE)
