## code to prepare `network_benvo` dataset goes here

set.seed(234131)
num_subj <- 2E3
sex <- rbinom(n = num_subj,size=1,prob = .5)
centered_income <- rlnorm(n=num_subj,mean= 0,log(1.5))
subj_pos <- cbind(runif(num_subj),runif(num_subj))
FFRpos <- cbind(runif(30),runif(30))
subj_FFR <- fields::rdist(subj_pos,FFRpos)
FFR_FFR <- fields::rdist(FFRpos,FFRpos)
f_direct <- function(x) .3*pweibull(x,shape=5,scale=.6,lower.tail = F)
f_indirect <- function(x) .1*pweibull(x,shape=4,scale=.3,lower.tail=F)
FFRexposure <- apply(subj_FFR,1,function(x) {sum(f_direct(x))})
FFR_sq_exposure <- sapply(1:num_subj,function(x) {
  ics <- which(subj_FFR[x,]<=.5)
  mat <- FFR_FFR[ics,]
  sum(f_indirect(mat[lower.tri(mat)]))
  })

y <- 25 + sex*-2  + centered_income*-2  + FFRexposure + FFR_sq_exposure + rnorm(num_subj,sd = .5)


subj_df <- dplyr::tibble(ID = 1:num_subj,
                         BMI = y,
                         sex = sex,
                         FFR_exposure = FFRexposure,
                         FFR_sq_exposure = FFR_sq_exposure)

FFR_df <- purrr::map_dfr(1:num_subj,function(x) tibble(ID = x,
                                                       Distances = subj_FFR[x,]))

FFR_df <- FFR_df %>% dplyr::filter(Distances<=1)

FFR_FFR_df <- purrr::map_df(1:num_subj,function(x) {
  ics <- which(subj_FFR[x,]<=.5)
  mat <- FFR_FFR[ics,]
  out <- tibble(ID = x,
         Distances = mat[lower.tri(mat)])
  return(out)
  })

FFR_FFR_df <- FFR_FFR_df %>% dplyr::filter(Distances<=1)

network_benvo <- benvo(subject_data = subj_df,
                       bef_data = list(FFR_df,FFR_FFR_df),
                       bef_names = c("Direct FFR","Indirect FFR"),
                       joining_id = "ID",
                       distance_col = c("Distances","Distances"))


usethis::use_data(network_benvo, overwrite = TRUE)
