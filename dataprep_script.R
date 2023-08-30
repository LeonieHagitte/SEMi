---
  title: "Dataprep"
author: "Leonie H."
date: "2023-08-16"
output: html_document
---
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
setwd("/Users/admin/Desktop/MPI/Share/ShareSEM")

library(haven)
library(tidyverse)

## Wave 1: Demographics (demo), Cognitive Function (cf), Mental Health (mf)

demo1 <- read_dta("sharew1_rel8-0-0_ALL_datasets_stata/sharew1_rel8-0-0_dn.dta")  ## read data
cf1 <- read_dta("sharew1_rel8-0-0_ALL_datasets_stata/sharew1_rel8-0-0_cf.dta")
mh1 <- read_dta("sharew1_rel8-0-0_ALL_datasets_stata/sharew1_rel8-0-0_mh.dta")
casp1 <- read_dta("sharew1_rel8-0-0_ALL_datasets_stata/sharew1_rel8-0-0_dropoff.dta")

demo1items<- demo1 %>% select("mergeid", "dn003_", "dn042_") %>%  ## select items
  rename_at(vars(-mergeid),function(x) paste0(x,"_1"))            ## add wave number as suffix to variable name

cf1items<- cf1 %>% select(mergeid, contains ("cf0"))%>%
  rename_at(vars(-mergeid),function(x) paste0(x,"_1"))

mh1items<- mh1 %>% select(mergeid, contains ("mh0")) %>%
  select(- c( "mh018_", "mh019_", "mh020_", "mh021_")) %>%
  rename_at(vars(-mergeid),function(x) paste0(x,"_1"))

casp1items <- casp1 %>% dplyr::select("mergeid", "q2_a", "q2_b", "q2_c", "q2_d","q2_e","q2_f","q2_g","q2_h","q2_i",
                                      "q2_j","q2_k","q2_l") %>%
  rename_at(vars(-mergeid),function(x) paste0(x,"_1"))


## Wave 2 

demo2 <- read_dta("sharew2_rel8-0-0_ALL_datasets_stata/sharew2_rel8-0-0_dn.dta")  
cf2 <- read_dta("sharew2_rel8-0-0_ALL_datasets_stata/sharew2_rel8-0-0_cf.dta")
mh2 <- read_dta("sharew2_rel8-0-0_ALL_datasets_stata/sharew2_rel8-0-0_mh.dta")
casp2 <- read_dta("sharew2_rel8-0-0_ALL_datasets_stata/sharew2_rel8-0-0_ac.dta")

demo2items<- demo2 %>% select("mergeid", "dn003_", "dn042_") %>%
  rename_at(vars(-mergeid),function(x) paste0(x,"_2"))

cf2items<- cf2 %>% select(mergeid, contains ("cf0")) %>%
  rename_at(vars(-mergeid),function(x) paste0(x,"_2"))

mh2items<- mh2 %>% select(mergeid, contains ("mh0")) %>%
  select (- c ("mh018_", "mh019_", contains ("mh02"))) %>%
  rename_at(vars(-mergeid),function(x) paste0(x,"_2"))

casp2items<- casp2 %>% dplyr::select("mergeid", "ac014_", "ac015_", "ac016_", "ac017_","ac018_","ac019_","ac020_", "ac021_","ac022_","ac023_","ac024_","ac025_") %>%
  rename_at(vars(-mergeid),function(x) paste0(x,"_2"))

## Wave 3 

#demo3 <- read_dta("")  
#cf3 <- read_dta("")
#mh3 <- read_dta("")
casp3 <- read_dta("sharew3_rel8-0-0_ALL_datasets_stata/sharew3_rel8-0-0_ac.dta")

#demo3items<- demo3 %>% select("mergeid", "dn003_", "dn042_") %>%
#  rename_at(vars(-mergeid),function(x) paste0(x,"3"))

#cf3items<- cf3 %>% select(mergeid, contains ("cf0")) %>%
#  rename_at(vars(-mergeid),function(x) paste0(x,"_3"))

#mh3items<- mh3 %>% select(mergeid, contains ("mh0")) %>%
#  select (- c ("mh018_", "mh019_", contains ("mh02"))) %>%
#  rename_at(vars(-mergeid),function(x) paste0(x,"_3"))

#casp3items<- casp3 %>% dplyr::select("mergeid", "ac014_", "ac015_", "ac016_", "ac017_","ac018_","ac019_","ac020_","ac021_","ac022_","ac023_","ac024_","ac025_") %>%
#  rename_at(vars(-mergeid),function(x) paste0(x,"_3"))


## Wave 4

demo4 <- read_dta("sharew4_rel8-0-0_ALL_datasets_stata/sharew4_rel8-0-0_dn.dta")
cf4 <- read_dta("sharew4_rel8-0-0_ALL_datasets_stata/sharew4_rel8-0-0_cf.dta")
mh4 <- read_dta("sharew4_rel8-0-0_ALL_datasets_stata/sharew4_rel8-0-0_mh.dta")
casp4 <- read_dta("sharew4_rel8-0-0_ALL_datasets_stata/sharew4_rel8-0-0_ac.dta")

demo4items<- demo4 %>% select("mergeid", "dn003_", "dn042_") %>%
  rename_at(vars(-mergeid),function(x) paste0(x,"_4"))

cf4items<- cf4 %>% select(mergeid, contains ("cf0")) %>%
  rename_at(vars(-mergeid),function(x) paste0(x,"_4"))

mh4items<- mh4 %>% select(mergeid, contains ("mh0")) %>%
  select( -c( "mh018_", "mh019_", contains ("mh02"), contains("mh03"))) %>%
  rename_at(vars(-mergeid),function(x) paste0(x,"_4"))

casp4items<- casp4 %>% dplyr::select("mergeid", "ac014_", "ac015_", "ac016_", "ac017_","ac018_","ac019_","ac020_", "ac021_","ac022_","ac023_","ac024_","ac025_") %>%
  rename_at(vars(-mergeid),function(x) paste0(x,"_4"))

## Wave 5

demo5 <- read_dta("sharew5_rel8-0-0_ALL_datasets_stata/sharew5_rel8-0-0_dn.dta")
cf5 <- read_dta("sharew5_rel8-0-0_ALL_datasets_stata/sharew5_rel8-0-0_cf.dta")
mh5 <- read_dta("sharew5_rel8-0-0_ALL_datasets_stata/sharew5_rel8-0-0_mh.dta")
casp5 <- read_dta("sharew5_rel8-0-0_ALL_datasets_stata/sharew5_rel8-0-0_ac.dta")

demo5items<- demo5 %>% select("mergeid", "dn003_", "dn042_") %>%
  rename_at(vars(-mergeid),function(x) paste0(x,"_5"))

cf5items<- cf5 %>% select(mergeid, contains ("cf0")) %>%
  rename_at(vars(-mergeid),function(x) paste0(x,"_5"))

mh5items<- mh5 %>% select(mergeid, contains ("mh0")) %>%
  select( -c(contains ("mh02"), contains("mh03"))) %>%
  rename_at(vars(-mergeid),function(x) paste0(x,"_5"))

casp5items<- casp5 %>% dplyr::select("mergeid", "ac014_", "ac015_", "ac016_", "ac017_","ac018_","ac019_","ac020_", "ac021_","ac022_","ac023_","ac024_","ac025_") %>%
  rename_at(vars(-mergeid),function(x) paste0(x,"_5"))

## Wave 6

demo6 <- read_dta("sharew6_rel8-0-0_ALL_datasets_stata/sharew6_rel8-0-0_dn.dta")
cf6 <- read_dta("sharew6_rel8-0-0_ALL_datasets_stata/sharew6_rel8-0-0_cf.dta")
mh6 <- read_dta("sharew6_rel8-0-0_ALL_datasets_stata/sharew6_rel8-0-0_mh.dta")
casp6 <- read_dta("sharew6_rel8-0-0_ALL_datasets_stata/sharew6_rel8-0-0_ac.dta")

demo6items<- demo6 %>% select("mergeid", "dn003_", "dn042_") %>%
  rename_at(vars(-mergeid),function(x) paste0(x,"_6"))

cf6items<- cf6 %>% select(mergeid, contains ("cf0")) %>%
  select( -"cf019_") %>%
  rename_at(vars(-mergeid),function(x) paste0(x,"_6"))

mh6items<- mh6 %>% select(mergeid, contains ("mh0")) %>%
  select( -c(contains ("mh02"), contains("mh03"))) %>%
  rename_at(vars(-mergeid),function(x) paste0(x,"_6"))

casp6items<- casp6 %>% dplyr::select("mergeid","ac014_", "ac015_", "ac016_", "ac017_","ac018_","ac019_","ac020_",  "ac021_","ac022_","ac023_","ac024_","ac025_") %>%
  rename_at(vars(-mergeid),function(x) paste0(x,"_6"))

## Wave 7

demo7 <- read_dta("sharew7_rel8-0-0_ALL_datasets_stata/sharew7_rel8-0-0_dn.dta")
cf7 <- read_dta("sharew7_rel8-0-0_ALL_datasets_stata/sharew7_rel8-0-0_cf.dta")
mh7 <- read_dta("sharew7_rel8-0-0_ALL_datasets_stata/sharew7_rel8-0-0_mh.dta")
casp7 <- read_dta("sharew7_rel8-0-0_ALL_datasets_stata/sharew7_rel8-0-0_ac.dta")


demo7items<- demo7 %>% select("mergeid", "dn003_", "dn042_") %>%
  rename_at(vars(-mergeid),function(x) paste0(x,"_7"))

cf7items<- cf7 %>% select(mergeid, contains ("cf0")) %>%
  rename_at(vars(-mergeid),function(x) paste0(x,"_7"))

mh7items<- mh7 %>% select(mergeid, contains ("mh0")) %>%
  select( -c(contains ("mh02"), contains("mh03"))) %>%
  rename_at(vars(-mergeid),function(x) paste0(x,"_7"))

casp7items<- casp7 %>% select("mergeid","ac014_", "ac015_", "ac016_", "ac017_","ac018_","ac019_","ac020_", "ac021_","ac022_","ac023_","ac024_","ac025_") %>%
  rename_at(vars(-mergeid),function(x) paste0(x,"_7"))

## Wave 8

demo8 <- read_dta("sharew8_rel8-0-0_ALL_datasets_stata/sharew8_rel8-0-0_dn.dta")
cf8 <- read_dta("sharew8_rel8-0-0_ALL_datasets_stata/sharew8_rel8-0-0_cf.dta")
mh8 <- read_dta("sharew8_rel8-0-0_ALL_datasets_stata/sharew8_rel8-0-0_mh.dta")
casp8 <- read_dta("sharew8_rel8-0-0_ALL_datasets_stata/sharew8_rel8-0-0_ac.dta")

demo8items<- demo8 %>% select("mergeid", "dn003_", "dn042_") %>%
  rename_at(vars(-mergeid),function(x) paste0(x,"_8"))

cf8items<- cf8 %>% select(mergeid, contains ("cf0")) %>%
  rename_at(vars(-mergeid),function(x) paste0(x,"_8"))

mh8items<- mh8 %>% select(mergeid, contains ("mh0")) %>%
  select( -c(contains ("mh02"), contains("mh03"))) %>%
  rename_at(vars(-mergeid),function(x) paste0(x,"_8"))

casp8items<- casp8 %>% select("mergeid", "ac014_", "ac015_", "ac016_", "ac017_","ac018_","ac019_","ac020_", "ac021_","ac022_","ac023_","ac024_","ac025_") %>%
  rename_at(vars(-mergeid),function(x) paste0(x,"_8"))


## Merge all Data for Cognitive Function including Demographics

cf_all <- full_join(demo1items, cf1items, by = "mergeid") %>%
  #  full_join (demo2items, by = "mergeid") %>%
  full_join (cf2items, by = "mergeid") %>%
  #  full_join (demo4items, by = "mergeid") %>%
  full_join (cf4items, by = "mergeid") %>%
  #  full_join (demo5items, by = "mergeid") %>%
  full_join (cf5items, by = "mergeid") %>%
  #  full_join (demo6items, by = "mergeid") %>%
  full_join (cf6items, by = "mergeid") %>%
  #  full_join (demo7items, by = "mergeid") %>%
  full_join (cf7items, by = "mergeid") %>%
  #  full_join (demo8items, by = "mergeid") %>%
  full_join (cf8items, by = "mergeid")

dim(cf_all)  
names(cf_all)

library(tidyr)

head(cf_all)
cf_long <- pivot_longer(cf_all, cols = -c(1:3), 
                        names_to = c('.value', 'wave'), 
                        names_pattern = '(.*)(\\d+)')


## Merge all Data for Mental Health including Demographics

mh_all <- full_join(demo1items, mh1items, by = "mergeid") %>%
  #  full_join (demo2items, by = "mergeid") %>%
  full_join (mh2items, by = "mergeid") %>%
  #  full_join (demo4items, by = "mergeid") %>%
  full_join (mh4items, by = "mergeid") %>%
  #  full_join (demo5items, by = "mergeid") %>%
  full_join (mh5items, by = "mergeid") %>%
  #  full_join (demo6items, by = "mergeid") %>%
  full_join (mh6items, by = "mergeid") %>%
  #  full_join (demo7items, by = "mergeid") %>%
  full_join (mh7items, by = "mergeid") %>%
  #  full_join (demo8items, by = "mergeid") %>%
  full_join (mh8items, by = "mergeid")

dim(mh_all)
mh_long <- pivot_longer(mh_all, cols = -c(1:3), 
                        names_to = c('.value', 'wave'), 
                        names_pattern = '(.*)(\\d+)')
head(mh_long)

cfmh_all <- full_join(cf_all, mh_all, by = "mergeid") 
dim(cfmh_all)

cfmh_long <- pivot_longer(cfmh_all, cols = -c(1:3), 
                          names_to = c('.value', 'wave'), 
                          names_pattern = '(.*)(\\d+)')
head(cfmh_long)


## Merge all Data for Casp including Demographics

casp_all <- full_join(demo1items, casp1items, by = "mergeid") %>%
  full_join (casp2items, by = "mergeid") %>%
  #full_join (casp3items, by = "mergeid") %>%
  full_join (casp4items, by = "mergeid") %>%
  full_join (casp5items, by = "mergeid") %>%
  full_join (casp6items, by = "mergeid") %>%
  full_join (casp7items, by = "mergeid") %>%
  full_join (casp8items, by = "mergeid")

dim(casp_all)  
names(casp_all)

library(tidyr)

head(casp_all)
casp_long <- pivot_longer(casp_all, cols = -c(1:3), 
                          names_to = c('.value', 'wave'), 
                          names_pattern = '(.*)(\\d+)')



# Sanity check
write.csv2(cf_all, "cf_all.csv")
write.csv2(mh_all, "mh.csv")
write.csv2(casp_all, "casp.csv")
write.csv2(cfmh_all, "cfmh.csv")

write.csv2(cf_long, "cf_long.csv")
write.csv2(mh_long, "mh_long.csv")
write.csv2(casp_long, "casp_long.csv")
write.csv2(cfmh_long, "cfmh_long.csv")


```
