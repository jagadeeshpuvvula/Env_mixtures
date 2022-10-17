
library(tidyverse)
library(haven)
library(janitor)
library(truncnorm)
select <- dplyr::select
rename <- dplyr::rename

####imputation function for HOME cohort
fLOD <- function(result, meanlog2, sdlog2){ 
  impute <- rtruncnorm(1, b=log2(result), mean=meanlog2, sd=sdlog2) 
  return(impute) 
}

impute <- function(data, dist_data, name, short_name) {
  name <- as.character(name)
  meanlog2 <- dist_data$meanlog2[dist_data$test.name == name]
  sdlog2 <- dist_data$sdlog2[dist_data$test.name == name]
  
  data1 <- data %>% # make an updated dataset with a new column coresponding to the imputed data. 
    mutate(log2.imputed = ifelse(data$flg_lod==1, # if row is flagged for being below the LOD
                                 apply(data, MARGIN=1, # apply the imputation
                                       FUN = function(x) 
                                         fLOD(data$result, # lower LOD for a particular row 
                                              meanlog2,
                                              sdlog2)),
                                 data$log2.result)) # if row is NOT flagged, keep it as is. 
  data2 <- data1 %>%
    mutate(untrans.res = 2^(data1$log2.imputed)) %>% # un-transform the results (for easier interpretation)
    select(c(subject_id, flg_lod, untrans.res))
  
  # fix column names
  name_untrans <- paste0(short_name, ".res")
  name_lod <- paste0(short_name, ".flag")
  
  colnames(data2)[2] <- name_lod
  colnames(data2)[3] <- name_untrans
  
  return(data2)
}
################ Function ends here #############

#Impute values from data

labres$result <- ifelse(labres$flg_lod==1, labres$analyte_lod/sqrt(2), 
                        labres$result)
labres <-  mutate(labres, log2.result = log2(labres$result)) 

LODmeansd.all <- labres |>
  group_by(analyte_code, visit) |> 
  summarise(meanlog2 = mean(log2.result, na.rm= TRUE), 
            sdlog2 = sd(log2.result, na.rm=TRUE)) |>
              mutate(test.name = (paste0(analyte_code, ".t", visit)))

labres$result <- ifelse(labres$flg_lod==1, labres$result*(2/sqrt(2)), 
                        labres$result) 
labres$log2.result <- log2(labres$result) 

###############################################################################
#repeat per analyte per visit
dat <- labres |> 
  filter(analyte_code == "BP-3") |> 
  filter(visit == "M12")

set.seed(1010)
dat <- impute(dat, LODmeansd.all, "BP-3.tM12", "result")
###############################################################################
dat <- dat|>
  mutate(visit="m12") |>
  mutate(analyte="bp3")

write_csv(dat[c(1,3:5)], "/Users/jpuvvula/Documents/data/imputed/bp3_m12.csv")
###############################################################################
