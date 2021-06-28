## code to prepare `rvc1` dataset goes here
library(tidyverse)
set.seed(34234)
n <- 1000
num_dists <- rpois(n,20)
dists <- map2_dfr(1:n,num_dists,function(x,y) tibble(subj_id = x,
                                                     Distance = runif(y),
                                                     Rank = rank(Distance))) %>% 
  arrange(subj_id,Distance)


sink(file=tempfile(fileext=".txt"))
ddf <- map2_dfr(dists$subj_id,dists$Rank,function(x,y){
  dists %>% 
    filter(subj_id == {{x}},
           Rank == {{y}}) %>% 
    select(subj_id,Rank,Distance) %>% 
    crossing(dists %>% 
               filter(subj_id == {{x}},
                      Rank != {{y}}) %>% select(Distance) %>% 
               rename(OtherDistances=Distance)) %>% 
    mutate(gd_ten = (abs(Distance-OtherDistances)<=.1)*1,
           gd_twenty = (abs(Distance-OtherDistances)<=.2)*1,
           gd_fifty = (abs(Distance-OtherDistances)<=.5)*1
    ) %>% group_by(subj_id,Rank,Distance) %>% 
    summarise(gd_ten = sum(gd_ten),
              gd_twenty = sum(gd_twenty),
              gd_fifty = sum(gd_fifty))
  
}) %>% ungroup()
sink()


fexp <- function(x) pweibull(x,shape = 5, 
                             scale=.5,
                             lower.tail = FALSE)

qexp <- function(x,y,z){
  exp( x*-.4 )
}

sdf <- ddf %>% 
  group_by(subj_id) %>% 
  summarise(Exposure = sum(pmap_dbl(list(gd_ten,gd_twenty,gd_fifty,Distance),function(x,y,z,d) qexp(x,y,z)*fexp(d) ))) %>% 
  mutate(Sex = rbinom(n(),size = 1,prob = .5),
         BMI = 26 + -2.2*Sex + Exposure + rnorm(n()))

rvc_bdf <- rbenvo::benvo(subject_data = sdf,sub_bef_data = list(FFR=ddf %>% select(-Rank)))

usethis::use_data(rvc_bdf, overwrite = TRUE)
