library(dplyr)
library(truncnorm)
library(readr)
library(tidyr)
library(here)
library(stringr)
select <- dplyr::select
rename <- dplyr::rename


# The single imputation process is carried out with the "impute" function below in step 7. Note that the impute function relies on the "fLOD" function. These functions could have been combined into 1, but they were created at different times so that's just how it is. 

fLOD <- function(test.result, meanlog2, sdlog2){ ##Create a new function, fLOD, which needs 3 inputs: the LOD value (represented by labres$test.result in FLAGGED rows); the mean (meanlog2); the sd (sdlog2)
  impute <- rtruncnorm(1, b=log2(test.result), mean=meanlog2, sd=sdlog2) #the function will then run the rtruncnorm command given the inputs (eg. b= upper limit = log2(LOD), mean, and sd) 
  return(impute) 
}

impute <- function(data, dist_data, name, short_name) {
  name <- as.character(name)
  meanlog2 <- dist_data$meanlog2[dist_data$test.name == name]
  sdlog2 <- dist_data$sdlog2[dist_data$test.name == name]
  
  data1 <- data %>% # make an updated dataset with a new column coresponding to the imputed data. 
    mutate(log2.imputed = ifelse(data$lod.flag==1, # if row is flagged for being below the LOD
                                 apply(data, MARGIN=1, # apply the imputation
                                       FUN = function(x) 
                                         fLOD(data$test.result, # lower LOD for a particular row 
                                              meanlog2,
                                              sdlog2)),
                                 data$log2.test.result)) # if row is NOT flagged, keep it as is. 
  data2 <- data1 %>%
    mutate(untrans.res = 2^(data1$log2.imputed)) %>% 
    select(c(specimen_bar_code, lod.flag, untrans.res))
  
  # fix column names
  name_untrans <- paste0(short_name, ".res")
  name_lod <- paste0(short_name, ".flag")
  
  colnames(data2)[2] <- name_lod
  colnames(data2)[3] <- name_untrans
  
  return(data2)
}

#Start here
###############################################################################
###############################################################################

labres <- mutate(labres, lod.flag = ifelse((substr(labres$tests_results_raw_result,1,1) == "<"), 1, 0))
labres$test.result <- str_replace_all(labres$tests_results_raw_result, "[^[:alnum:]\\.\\s]","")

labres$test.result <- as.numeric(labres$test.result)
labres$test.result <- ifelse(labres$lod.flag==1, labres$test.result/sqrt(2), 
                             labres$test.result)

#########################################
#Check percent of observations below LOD
#########################################
summ<- labres |>
  group_by(tests_results_test_name) |>
  summarise(tested=n(),
            bel_lod= sum(lod.flag)) |>
  mutate(pct_lod= round((bel_lod/tested)*100,2))


labres <-  mutate(labres, log2.test.result = log2(labres$test.result)) 

LODmeansd.all <- labres %>% 
  group_by(tests_results_test_name) %>% 
  summarise(meanlog2 = mean(log2.test.result, na.rm= TRUE), 
            sdlog2 = sd(log2.test.result, na.rm=TRUE)) %>% 
  mutate(test.name = (tests_results_test_name))


labres$test.result <- ifelse(labres$lod.flag==1, labres$test.result*(2/sqrt(2)), 
                             labres$test.result) 
labres$log2.test.result <- log2(labres$test.result) 

## Step 7: Complete the imputation for each chemical, and store the 
#results in a dataset for each individual chemical/ timepoint.  

analyte<-"DPhPUR"

dat <- labres |>
  filter(tests_results_test_name == analyte)

set.seed(1010) 
dat <- impute(dat, LODmeansd.all, analyte, "result") 

dat<- dat |>
  rename(subject_id=specimen_bar_code) |>
  mutate(visit="12w") |>
  mutate(analyte="DPhP") |>
  drop_na()

write_csv(dat[c(1,3:5)], "E://BBK17//pj//imputed_data//DPhP.csv")
