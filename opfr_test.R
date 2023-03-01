library(tidyverse)
library(readxl)
library(janitor)
#==========================================================================
#Original data
dat <-read_csv2("E://BBK17/2. BBK17_pending lab tests_July19 2022//LAB//MIREC_MICROALIQOUT RESULTS table code 001-01.001-03//QueryData_code001-03.csv") |>
  clean_names()|>
  rename(metric = micro_aliquot_tests_results_test_name,
         value =micro_aliquot_tests_results_raw_result)|>
  filter(metric == "BClEtPUR" | 
           metric == "BDCliPrPUR" |
           metric == "DPhPUR" |
           metric== "DBPUR")|>
  filter(!is.na(as.numeric(value)))|>
  select(c(5,6))

dat$value<- as.numeric(dat$value)
summary(factor(dat$metric))


mirec_all<- dat |> 
  group_by(metric) |>
  summarise(median = median((value), na.rm = T), 
            q25 = quantile(value, .25, na.rm=T),
            q75 = quantile(value, .75, na.rm=T),
            max= max(value, na.rm = T)) |>
  mutate_if(is.numeric, round, 3)

ggplot(dat, aes(x=value))+geom_density()+facet_wrap(~metric, scales = "free")

#==========================================================================
#WPPSI subset
dat_wide <- read_csv("E:/BBK17/pj/wppsi_exp_out_comb/wppsi_fin.csv") |>
  mutate_if(is.numeric, round, 4) |>
  filter(cohort==2) |>
  drop_na() |>
  clean_names() |>
  rename(BCEtP=bcetp_sg, BDCPP=bdcpp_sg, DBuP=dbup_sg, DPhP=dphp_sg) |>
  select(subject_id, BCEtP,BDCPP,DBuP,DPhP)|>
  pivot_longer(!subject_id, names_to = "measure", values_to = "value")


mirec_wppsi_set<- dat_wide |> 
  group_by(measure) |>
  summarise(median = median((value), na.rm = T), q25 = quantile(value, .25, na.rm=T),
            q75 = quantile(value, .75, na.rm=T),
            max= max(value, na.rm = T)) |>
  mutate_if(is.numeric, round, 3)






