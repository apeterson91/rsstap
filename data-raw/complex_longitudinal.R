## code to prepare `complex_longitudinal` dataset goes here

set.seed(3431)
num_subj <- 6E2
num_visits <- rpois(num_subj,.9)+1
visit_num <- purrr::map_dfr(1:num_subj,function(x) dplyr::tibble(id = x,measurement = 1:num_visits[x]))
num_obs <- sum(num_visits)
Z <- rbinom(num_subj,1,.5)
sjdf <- dplyr::tibble(id=1:num_subj,
                      sex = Z,
                      subj_effect = rnorm(num_subj,mean = 0,sd = .3),
                      subj_slope = rnorm(num_subj,mean=0,sd=0.3))
has_exp <- rbinom(num_obs,size = 1,prob = .95)
cnt <- rpois(num_obs,10)*has_exp
ldists <- lapply(cnt,function(x) runif(x) )
ltime <- lapply(cnt,function(x) 5*runif(x) )
f <- function(x,y) pweibull(x,shape=5,scale=.6,lower.tail = F)*pweibull(y,shape=1,scale=1.3,lower.tail=T)

HFS_distances_times <- purrr::map_dfr(1:length(ldists),function(x) {dplyr::tibble(ix=x,Distance=ldists[[x]],Time=ltime[[x]])}) %>% 
  dplyr::right_join(visit_num %>% dplyr::mutate(ix=1:dplyr::n())) %>%
  dplyr::filter(!is.na(Distance)) %>%
  dplyr::select(id,measurement,Distance,Time)

HFS_distances_times %>% dplyr::group_by(id,measurement) %>% 
  dplyr::summarise(exposure = sum(f(Distance,Time))) %>% 
  dplyr::mutate(between = -.85*mean(exposure),
         within = -.25 *(exposure- mean(exposure))
         ) ->edf

sjdf <- dplyr::left_join(sjdf,visit_num) %>% dplyr::left_join(edf) %>% 
  dplyr::mutate(year = measurement -1,
         between = tidyr::replace_na(between,0),
         within  = tidyr::replace_na(within,0))


sjdf$BMI <- 33 +  sjdf$sex* -2.2 + .1*sjdf$year + sjdf$between + sjdf$within + sjdf$subj_effect + sjdf$subj_slope*sjdf$year + rnorm(num_obs)



complex_longitudinal <- rbenvo::benvo(subject_data = sjdf,
                          bef_data = list(HFS=HFS_distances_times),
                          by = c("id","measurement"))

usethis::use_data(complex_longitudinal, overwrite = TRUE)
