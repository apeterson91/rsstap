## code to prepare `binomial_benvo` dataset goes here


num_subj <- 300
sex <- rbinom(n = num_subj,size=1,prob = .5)
centered_income <- rlnorm(n=num_subj,mean= 0,log(1.5))
hasFFR_exp <- rbinom(num_subj,size = 1,prob = .9)
FFRcnt <- rpois(num_subj,10)*hasFFR_exp
ldists <- lapply(FFRcnt,function(x) runif(x) )
f <- function(x) .8*pweibull(x,shape=5,scale=.6,lower.tail = F)
FFRexposure <- sapply(ldists,function(x) sum(f(x)))
hasHFS_exp <- rbinom(num_subj,size=1,prob=.8)
HFScnt <- (rpois(num_subj,5)+1)*hasHFS_exp
HFSdists <- lapply(HFScnt,function(x) runif(x) )
f <- function(x) -.4*pweibull(x,shape=3.5,scale=.7,lower.tail = F)
HFSexposure <- sapply(HFSdists,function(x) sum(f(x)))
eta <- -0.5462739 + sex*-1 + centered_income*-2 + HFSexposure + FFRexposure
nn <- rpois(num_subj,25)
y <- rbinom(n = num_subj,size=nn,prob = binomial()$linkinv(eta))



subj_df <- dplyr::tibble(ID= 1:num_subj,
                         NumObese = y,
                         NumTotal = nn,
                         Num_NotObese = nn-y,
                         sex = sex,
                         Centered_Scaled_Income = centered_income)

FFR <- purrr::map2_dfr(1:length(ldists),ldists,function(x,y) dplyr::tibble(ID=x,Distance=y))
HFS <- purrr::map2_dfr(1:length(HFSdists),HFSdists,function(x,y) dplyr::tibble(ID=x,Distance=y))

binomial_benvo <- rbenvo::benvo(subject_data = subj_df,
                        bef_data = list(FFR=FFR,HFS=HFS),
                        by = c("ID"))


usethis::use_data(binomial_benvo, overwrite = TRUE)
