
  ## Step 0: load libraries and functions

#load packages
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
    mutate(untrans.res = 2^(data1$log2.imputed)) %>% # un-transform the results (for easier interpretation)
    select(c(subject.id, lod.flag, untrans.res))
  
  # fix column names
  name_untrans <- paste0(short_name, ".res")
  name_lod <- paste0(short_name, ".flag")
  
  colnames(data2)[2] <- name_lod
  colnames(data2)[3] <- name_untrans
  
  return(data2)
}


```


##Step 1: Initial set up
```{r}
# clean the data...
```

## Step 2: Flag below LOD results, remove the "<"
in the MIREC study, below LOD rows are written as "<" followed by whatever the LOD is. I want to remove the "<" and just keep the LOD values, and create a separate variable ("lod.flag") which indicates that a row is below the LOD. 
```{r}
# Create a new column in labres which acts as a flag: indicates whether or not a lab result value is is BELOW the LOD
labres <- mutate(labres, lod.flag = ifelse((substr(labres$test.result,1,1) == "<"), 1, 0))

# removing the "<" (this code actually removes all symbols except for ".") 
labres$test.result <- str_replace_all(labres$test.result, "[^[:alnum:]\\.\\s]","")

```

## Step 3: Dividing test results below the LOD by sqrt(2)
```{r}
labres$test.result <- as.numeric(labres$test.result)

labres$test.result <- ifelse(labres$lod.flag==1, labres$test.result/sqrt(2), labres$test.result) #IF a row is'flagged'(ie. lod.flag ==1), replace its test.result value with LOD/ sqrt(2). Otherwise, keep the test.result value as is. 

```
We are making our best guess as to what the chemical concentrations are for the below LOD participants by assuming that chemical exposure is normally distributed and assigning everyone with a below LOD value a test.result value of LOD/sqrt(2). Otherwise, we could exclude the below LOD people from the mean, sd calculation in step 5. Doing so would (probably) introduce a lot of bias. 


##Step 4: Make a new column in the labres dataframe representing the log (base 2) of the test results data
```{r step 4}
labres <-  mutate(labres, log2.test.result = log2(labres$test.result)) 
```



##Step 5: Making a new dataframe with columns for each chemicals' log2(mean) and log2(sd). 
Information used in imputation later
```{r}
LODmeansd.all <- labres %>% 
  group_by(test.name, test.time) %>% 
  summarise(meanlog2 = mean(log2.test.result, na.rm= TRUE), sdlog2 = sd(log2.test.result, na.rm=TRUE)) %>% 
  mutate(test.name = (paste0(test.name, ".t", test.time)))
```


##Step 6: Undo step 3
This step is necessary because I will need to have the correct test.result value in flagged rows. If not, there will be problems performing the imputation. As we will see, the test.result value is the upper limit for imputation, so if the test.result value is wrong, imputation will be wrong. 
```{r}
#fix test.result values for flagged rows:
labres$test.result <- ifelse(labres$lod.flag==1, labres$test.result*(2/sqrt(2)), labres$test.result) 

#IF a row is 'flagged' (ie lod.flag ==1), replace its test.result value with LOD*(2/sqrt(2)), which undoes step 3. If a row is not flagged, keep the test.result value as is. 

#fix log2.test.results for flagged rows:
labres$log2.test.result <- log2(labres$test.result) 

```

## Step 7: Complete the imputation for each chemical, and store the results in a dataset for each individual chemical/ timepoint.  
### cotinine
```{r}
set.seed(1010) # imputation is a random process. Setting a seed ensures that the exact same numbers will always be imputed. 

# cot
## t1
cot.t1 <- labres %>% # make a dataset that is a subset of labres, which only includes lab results on cotinine during visit #1. 
  filter(test.name == "COTISE") %>% 
  filter(test.time == 1)

# The impute function (which I created in step 0) carries out the single imputation process given 4 inputs
## first, a dataset of a chemical that has not been imputed yet. It that was created above
## second, a dataset with information on the distribution (log2-mean and SD) of all chemicals. This is "LODmeansd.all", which was created in step 5. 
## third, a character string which indicates which row of "LODmeansd.all" should be used. 
## forth, a character string representing the simplified name that *I* want to give this chemical and time point. It can be whatever you want
cot.t1 <- impute(cot.t1, LODmeansd.all, "COTISE.t1", "cot.t1") 

## t3
cot.t3 <- labres %>% 
  filter(test.name == "COTISE") %>% 
  filter(test.time == 1)
cot.t3 <- impute(cot.t3, LODmeansd.all, "COTISE.t3", "cot.t3")

## t4
cot.t4 <- labres %>% 
  filter(test.name == "COTISE") %>% 
  filter(test.time == 4)
cot.t4 <- impute(cot.t4, LODmeansd.all, "COTISE.t4", "cot.t4")

```


### repeat for all remaining chemicals...

## Then merge the cot.t1, cot.t3, cot.t4, ..., dataframes together. Change units, make lipid/ specific  gravity adjustements as needed. 