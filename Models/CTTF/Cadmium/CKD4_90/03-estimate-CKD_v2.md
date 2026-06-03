Global incidence of CKD4 Age 85-125 • Estimate incidence
================
LoVa3397
2025-10-16

- [Settings](#settings)
- [Model fit](#model-fit)
- [Predict all](#predict-all)
- [Summarize predictions:](#summarize-predictions)
  - [Global](#global)
  - [Regions](#regions)
  - [Subregions](#subregions)
  - [Countries](#countries)
- [Session info](#session-info)

# Settings

``` r
## required packages ----
library(bd)
library(brms)
library(FERG2)
library(ggplot2)
library(knitr)
library(rmarkdown)
library(sf)
library(tidyr)
library(dplyr)
library(DescTools)
```

    ## 
    ## Attaching package: 'DescTools'

    ## The following objects are masked from 'package:Hmisc':
    ## 
    ##     %nin%, Label, Mean, Quantile

``` r
library(readxl)
library(kableExtra)

## global options ----
knitr::opts_chunk$set(fig.width = 10)
do.call(file.remove, list(list.files(params$PlotDir, full.names = TRUE)))
```

    ##  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE

# Model fit

``` r
es <- readRDS(paste0(params$Dir, "/", params$es))
es <- subset(es, as.integer(FLAG) == 1)
fit_brms_reg_s <- readRDS(paste0(params$Dir, "/", params$fit))
summary(fit_brms_reg_s)
```

    ##  Family: gaussian 
    ##   Links: mu = identity; sigma = identity 
    ## Formula: yi | se(sei) ~ 1 + YEAR + (1 | REG2) + (1 | REG2:SUB2) + (1 | REG2:SUB2:COUNTRY) + (1 | REG2:SUB2:COUNTRY:ID) + (1 | REG2:SUB2:COUNTRY:ID:DTP_ID) 
    ##    Data: subset(es, as.integer(FLAG) == 1) (Number of observations: 200) 
    ##   Draws: 5 chains, each with iter = 5000; warmup = 3000; thin = 1;
    ##          total post-warmup draws = 10000
    ## 
    ## Multilevel Hyperparameters:
    ## ~REG2 (Number of levels: 4) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     1.53      0.81     0.12     3.13 1.00     5342     5709
    ## 
    ## ~REG2:SUB2 (Number of levels: 7) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.93      0.68     0.04     2.50 1.00     7990     6207
    ## 
    ## ~REG2:SUB2:COUNTRY (Number of levels: 17) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     1.99      0.76     0.34     3.41 1.00     2751     2406
    ## 
    ## ~REG2:SUB2:COUNTRY:ID (Number of levels: 110) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     3.81      0.37     3.12     4.58 1.00     5259     7413
    ## 
    ## ~REG2:SUB2:COUNTRY:ID:DTP_ID (Number of levels: 200) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     1.19      0.19     0.87     1.62 1.00     3454     5068
    ## 
    ## Regression Coefficients:
    ##           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept   219.30    111.25     0.39   438.06 1.00     5452     6701
    ## YEAR         -0.11      0.06    -0.22    -0.00 1.00     5455     6817
    ## 
    ## Further Distributional Parameters:
    ##       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sigma     0.00      0.00     0.00     0.00   NA       NA       NA
    ## 
    ## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

``` r
zero_cases<- read_xlsx("endemic_countries.xlsx")%>%
  select(REG2, SUB2, ISO3, Country, cttf_cadmium) %>% 
  rename(COUNTRY=ISO3, COUNTRY_LABEL = Country) %>%
  mutate(DISEASEFREE = cttf_cadmium)

kable(
  caption = "Disease-free countries",
  row.names = FALSE,
  subset(zero_cases, cttf_cadmium==0)[, 4])
```

| COUNTRY_LABEL |
|:--------------|

Disease-free countries

``` r
country_with_data <- es %>% select(ISO3) %>% distinct() %>% mutate(DATA=1, COUNTRY = ISO3)
Sub2_with_data <- es %>% select(SUB2) %>% distinct() %>% mutate(DATASUB2=1)
Reg2_with_data <- es %>% select(REG2) %>% distinct() %>% mutate(DATAREG2=1)
zero_cases <- left_join(zero_cases, country_with_data)
```

    ## Joining with `by = join_by(COUNTRY)`

``` r
zero_cases <- left_join(zero_cases, Sub2_with_data)
```

    ## Joining with `by = join_by(SUB2)`

``` r
zero_cases <- left_join(zero_cases, Reg2_with_data) %>%
  select(-c(ISO3)) %>%
  mutate(ESTIMATES = case_when(
    DATA == 1 ~ 1,
    DISEASEFREE == 0 ~ 2,
    is.na(DATA) & DISEASEFREE == 1 & DATASUB2 == 1 ~ 3,
    is.na(DATA) & DISEASEFREE == 1 & is.na(DATASUB2) & DATAREG2 == 1 ~ 4, 
    is.na(DATA) & DISEASEFREE == 1  & is.na(DATASUB2) & is.na(DATAREG2) ~5))
```

    ## Joining with `by = join_by(REG2)`

``` r
zero_cases$ESTIMATES <- factor(zero_cases$ESTIMATES, 
                               level = c(1,2,3,4,5),
                               labels = c("Data present", "Disease free", "Data in subregion", "Data in region", "Data in world"))
```

# Predict all

``` r
## set up dataframe
sim_all <-
  data.frame(
    sei = 0,
    REG2 = FERG2:::countries$REG2,
    SUB2 = FERG2:::countries$SUB2,
    COUNTRY = FERG2:::countries$ISO3,
    YEAR = rep(2000:2021, each = nrow(FERG2:::countries)))
sim_all <- sim_all %>% left_join(zero_cases) %>% select(sei, REG2, SUB2, COUNTRY, YEAR, ESTIMATES)
```

    ## Joining with `by = join_by(REG2, SUB2, COUNTRY)`

``` r
## draw from expected value of posterior predictive dist
set.seed(10)
# fit_all <- 
#   posterior_epred(
#     object = fit_brms_reg_s,
#     newdata = sim_all,
#     allow_new_levels = TRUE,
#     sample_new_levels = "uncertainty",
#     re_formula = ~ 1 + YEAR +
#       (1 | REG2) +
#       (1 | REG2:SUB2) +
#       (1 | REG2:SUB2:COUNTRY)
#   )

draws_fit <- as_draws_df(fit_brms_reg_s)
fit_all <- data.frame(1:10000)
for (x in 1:nrow(sim_all)){
  if (as.integer(sim_all[x, "ESTIMATES"]) == 1){
    # Data present for country
    fit_all[[paste0("V",x)]] <- draws_fit$b_Intercept +                                                                               # Global intercept
      sim_all[x, "YEAR"] * draws_fit$b_YEAR +                                                                                         # Year component
      draws_fit[[paste0("r_REG2[",sim_all[x,"REG2"],",Intercept]")]] +                                                                # Regional component
      draws_fit[[paste0("r_REG2:SUB2[",sim_all[x,"REG2"],"_",sim_all[x,"SUB2"],",Intercept]")]] +                                     # Sub regional component
      draws_fit[[paste0("r_REG2:SUB2:COUNTRY[",sim_all[x,"REG2"],"_",sim_all[x,"SUB2"],"_",sim_all[x,"COUNTRY"],",Intercept]")]]      # Country component
  } else if (as.integer(sim_all[x, "ESTIMATES"]) == 2) {
    # Disease-free country
    fit_all[[paste0("V",x)]] <- 0
  } else if (as.integer(sim_all[x, "ESTIMATES"]) == 3){
    # Data not present for country, but present in subregion
    fit_all[[paste0("V",x)]] <- draws_fit$b_Intercept +                                                                               # Global intercept
      sim_all[x, "YEAR"] * draws_fit$b_YEAR +                                                                                         # Year component
      draws_fit[[paste0("r_REG2[",sim_all[x,"REG2"],",Intercept]")]] +                                                                # Regional component
      draws_fit[[paste0("r_REG2:SUB2[",sim_all[x,"REG2"],"_",sim_all[x,"SUB2"],",Intercept]")]]                                       # Sub regional component
  } else if (as.integer(sim_all[x, "ESTIMATES"]) == 4){
    # Data not present for country, but present in region
    fit_all[[paste0("V",x)]] <- draws_fit$b_Intercept +                                                                               # Global intercept
      sim_all[x, "YEAR"] * draws_fit$b_YEAR +                                                                                         # Year component
      draws_fit[[paste0("r_REG2[",sim_all[x,"REG2"],",Intercept]")]]                                                                  # Regional component
  } else if (as.integer(sim_all[x, "ESTIMATES"]) == 5){
    # Data not present for country
    fit_all[[paste0("V",x)]] <- draws_fit$b_Intercept + 
      sim_all[x, "YEAR"] * draws_fit$b_YEAR
  } 
}

fit_all <- fit_all %>% select(-c(X1.10000))


## calculate cases
age_start <- as.numeric(gsub("-.*", "", params$Age))
age_end <- as.numeric(gsub(".*-", "", params$Age))
pop <- subset(FERG2:::pop, AGE >= age_start & AGE <= age_end)
sim_all$SIM <- t(fit_all)
pop_all <- aggregate(POP ~ ISO3 + YEAR, pop, sum)
sim_all <- merge(sim_all, pop_all,
                 by.x = c("COUNTRY", "YEAR"), by.y = c("ISO3", "YEAR"))
sim_all <- sim_all %>% left_join(zero_cases)
```

    ## Joining with `by = join_by(COUNTRY, REG2, SUB2, ESTIMATES)`

``` r
sim_all$CASES <- exp(sim_all$SIM) * sim_all$POP / 1e5
sim_all$CASES <- sim_all$CASES*sim_all$DISEASEFREE
sim_all$SIM<-sim_all$SIM*sim_all$DISEASEFREE
sim_all$sei<-sim_all$sei*sim_all$DISEASEFREE
saveRDS(sim_all, paste0(params$Dir, "/", params$sim_all))

## aggregate global
sim_all_glb <- with(sim_all, aggregate(CASES ~ YEAR, FUN = sum))
all_glb_id <- sim_all_glb[1]
all_glb_nr <-
  t(apply(sim_all_glb[, grepl("V", names(sim_all_glb))], 1, mean_median_ci))
all_glb_nr <- data.frame(all_glb_nr)
names(all_glb_nr) <- c("VAL_MEAN", "VAL_MEDIAN", "VAL_LWR", "VAL_UPR")
all_glb_nr <- cbind(all_glb_id, all_glb_nr)
all_glb_nr$LOCATION <- "Global"
all_glb_nr$LOCATION_NAME <- "Global"
all_glb_nr$METRIC <- "Number"
str(all_glb_nr)
```

    ## 'data.frame':    22 obs. of  8 variables:
    ##  $ YEAR         : int  2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 ...
    ##  $ VAL_MEAN     : num  33.9 31.1 28.5 26.3 24.2 ...
    ##  $ VAL_MEDIAN   : num  14.5 13.7 12.9 12 11.2 ...
    ##  $ VAL_LWR      : num  3.25 3.15 3.02 2.84 2.67 ...
    ##  $ VAL_UPR      : num  140 127 119 111 103 ...
    ##  $ LOCATION     : chr  "Global" "Global" "Global" "Global" ...
    ##  $ LOCATION_NAME: chr  "Global" "Global" "Global" "Global" ...
    ##  $ METRIC       : chr  "Number" "Number" "Number" "Number" ...

``` r
all_glb_rt <- all_glb_nr
all_glb_rt$POP <- with(sim_all, tapply(POP, YEAR, sum))
all_glb_rt$VAL_MEAN <- 1e5 * all_glb_rt$VAL_MEAN / all_glb_rt$POP
all_glb_rt$VAL_MEDIAN <- 1e5 * all_glb_rt$VAL_MEDIAN / all_glb_rt$POP
all_glb_rt$VAL_LWR <- 1e5 * all_glb_rt$VAL_LWR / all_glb_rt$POP
all_glb_rt$VAL_UPR <- 1e5 * all_glb_rt$VAL_UPR / all_glb_rt$POP
all_glb_rt$METRIC <- "Rate"
all_glb_rt$POP <- NULL
str(all_glb_rt)
```

    ## 'data.frame':    22 obs. of  8 variables:
    ##  $ YEAR         : int  2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 ...
    ##  $ VAL_MEAN     : num [1:22(1d)] 0.1159 0.1035 0.0933 0.0845 0.0765 ...
    ##  $ VAL_MEDIAN   : num [1:22(1d)] 0.0496 0.0456 0.0421 0.0385 0.0353 ...
    ##  $ VAL_LWR      : num [1:22(1d)] 0.01112 0.0105 0.00989 0.00914 0.00842 ...
    ##  $ VAL_UPR      : num [1:22(1d)] 0.479 0.423 0.388 0.355 0.326 ...
    ##  $ LOCATION     : chr  "Global" "Global" "Global" "Global" ...
    ##  $ LOCATION_NAME: chr  "Global" "Global" "Global" "Global" ...
    ##  $ METRIC       : chr  "Rate" "Rate" "Rate" "Rate" ...

``` r
## aggregate over regions
sim_all_reg <- with(sim_all, aggregate(CASES ~ REG2+YEAR, FUN = sum))
all_reg_id <- sim_all_reg[1:2]
all_reg_nr <-
  t(apply(sim_all_reg[, grepl("V", names(sim_all_reg))], 1, mean_median_ci))
all_reg_nr <- data.frame(all_reg_nr)
names(all_reg_nr) <- c("VAL_MEAN", "VAL_MEDIAN", "VAL_LWR", "VAL_UPR")
all_reg_nr <- cbind(all_reg_id, all_reg_nr)
all_reg_nr$LOCATION <- "Region"
all_reg_nr$LOCATION_NAME <- all_reg_nr$REG2
all_reg_nr$REG2 <- NULL
all_reg_nr$METRIC <- "Number"
str(all_reg_nr)
```

    ## 'data.frame':    132 obs. of  8 variables:
    ##  $ YEAR         : int  2000 2000 2000 2000 2000 2000 2001 2001 2001 2001 ...
    ##  $ VAL_MEAN     : num  0.0717 0.4225 0.0529 1.0774 18.751 ...
    ##  $ VAL_MEDIAN   : num  0.0275 0.0263 0.0203 0.2977 1.0622 ...
    ##  $ VAL_LWR      : num  0.001549 0.000234 0.001144 0.027297 0.020603 ...
    ##  $ VAL_UPR      : num  0.412 2.06 0.304 5.47 108.169 ...
    ##  $ LOCATION     : chr  "Region" "Region" "Region" "Region" ...
    ##  $ LOCATION_NAME: chr  "AFR" "AMR" "EMR" "EUR" ...
    ##  $ METRIC       : chr  "Number" "Number" "Number" "Number" ...

``` r
all_reg_rt <- all_reg_nr
all_reg_rt$POP <-
  with(sim_all, aggregate(POP ~ REG2 + YEAR, FUN = sum))$POP
all_reg_rt$VAL_MEAN <- 1e5 * all_reg_rt$VAL_MEAN / all_reg_rt$POP
all_reg_rt$VAL_MEDIAN <- 1e5 * all_reg_rt$VAL_MEDIAN / all_reg_rt$POP
all_reg_rt$VAL_LWR <- 1e5 * all_reg_rt$VAL_LWR / all_reg_rt$POP
all_reg_rt$VAL_UPR <- 1e5 * all_reg_rt$VAL_UPR / all_reg_rt$POP
all_reg_rt$METRIC <- "Rate"
all_reg_rt$POP <- NULL
str(all_reg_rt)
```

    ## 'data.frame':    132 obs. of  8 variables:
    ##  $ YEAR         : int  2000 2000 2000 2000 2000 2000 2001 2001 2001 2001 ...
    ##  $ VAL_MEAN     : num  0.00936 0.00662 0.00936 0.0098 0.61723 ...
    ##  $ VAL_MEDIAN   : num  0.003592 0.000411 0.003592 0.002708 0.034965 ...
    ##  $ VAL_LWR      : num  2.02e-04 3.67e-06 2.02e-04 2.48e-04 6.78e-04 ...
    ##  $ VAL_UPR      : num  0.0538 0.0323 0.0538 0.0497 3.5606 ...
    ##  $ LOCATION     : chr  "Region" "Region" "Region" "Region" ...
    ##  $ LOCATION_NAME: chr  "AFR" "AMR" "EMR" "EUR" ...
    ##  $ METRIC       : chr  "Rate" "Rate" "Rate" "Rate" ...

``` r
## aggregate over subregions
sim_all_sub <- with(sim_all, aggregate(CASES ~ SUB2+YEAR, FUN = sum))
all_sub_id <- sim_all_sub[1:2]
all_sub_nr <-
  t(apply(sim_all_sub[, grepl("V", names(sim_all_sub))], 1, mean_median_ci))
all_sub_nr <- data.frame(all_sub_nr)
names(all_sub_nr) <- c("VAL_MEAN", "VAL_MEDIAN", "VAL_LWR", "VAL_UPR")
all_sub_nr <- cbind(all_sub_id, all_sub_nr)
all_sub_nr$LOCATION <- "Subregion"
all_sub_nr$LOCATION_NAME <- all_sub_nr$SUB2
all_sub_nr$SUB2 <- NULL
all_sub_nr$METRIC <- "Number"
str(all_sub_nr)
```

    ## 'data.frame':    374 obs. of  8 variables:
    ##  $ YEAR         : int  2000 2000 2000 2000 2000 2000 2000 2000 2000 2000 ...
    ##  $ VAL_MEAN     : num  0.0132 0.0348 0.0237 0.1228 0.296 ...
    ##  $ VAL_MEDIAN   : num  0.00507 0.01333 0.00908 0.00992 0.00652 ...
    ##  $ VAL_LWR      : num  2.86e-04 7.52e-04 5.12e-04 7.76e-05 4.11e-05 ...
    ##  $ VAL_UPR      : num  0.076 0.2 0.136 0.734 0.995 ...
    ##  $ LOCATION     : chr  "Subregion" "Subregion" "Subregion" "Subregion" ...
    ##  $ LOCATION_NAME: chr  "AFRAB" "AFRC" "AFRD" "AMRA" ...
    ##  $ METRIC       : chr  "Number" "Number" "Number" "Number" ...

``` r
all_sub_rt <- all_sub_nr
all_sub_rt$POP <-
  with(sim_all, aggregate(POP ~ SUB2 + YEAR, FUN = sum))$POP
all_sub_rt$VAL_MEAN <- 1e5 * all_sub_rt$VAL_MEAN / all_sub_rt$POP
all_sub_rt$VAL_MEDIAN <- 1e5 * all_sub_rt$VAL_MEDIAN / all_sub_rt$POP
all_sub_rt$VAL_LWR <- 1e5 * all_sub_rt$VAL_LWR / all_sub_rt$POP
all_sub_rt$VAL_UPR <- 1e5 * all_sub_rt$VAL_UPR / all_sub_rt$POP
all_sub_rt$METRIC <- "Rate"
all_sub_rt$POP <- NULL
str(all_sub_rt)
```

    ## 'data.frame':    374 obs. of  8 variables:
    ##  $ YEAR         : int  2000 2000 2000 2000 2000 2000 2000 2000 2000 2000 ...
    ##  $ VAL_MEAN     : num  0.00936 0.00936 0.00936 0.00263 0.01832 ...
    ##  $ VAL_MEDIAN   : num  0.003592 0.003592 0.003592 0.000213 0.000403 ...
    ##  $ VAL_LWR      : num  2.02e-04 2.02e-04 2.02e-04 1.67e-06 2.54e-06 ...
    ##  $ VAL_UPR      : num  0.0538 0.0538 0.0538 0.0158 0.0616 ...
    ##  $ LOCATION     : chr  "Subregion" "Subregion" "Subregion" "Subregion" ...
    ##  $ LOCATION_NAME: chr  "AFRAB" "AFRC" "AFRD" "AMRA" ...
    ##  $ METRIC       : chr  "Rate" "Rate" "Rate" "Rate" ...

``` r
## aggregate over countries
all_cnt_nr <- t(apply(sim_all$CASES, 1, mean_median_ci))
all_cnt_nr <- data.frame(all_cnt_nr)
names(all_cnt_nr) <- c("VAL_MEAN", "VAL_MEDIAN", "VAL_LWR", "VAL_UPR")
all_cnt_nr <- cbind(sim_all[1:2], all_cnt_nr)
all_cnt_nr$LOCATION <- "Country"
all_cnt_nr$LOCATION_NAME <- all_cnt_nr$COUNTRY
all_cnt_nr$COUNTRY <- NULL
all_cnt_nr$METRIC <- "Number"
str(all_cnt_nr)
```

    ## 'data.frame':    4268 obs. of  8 variables:
    ##  $ YEAR         : int  2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 ...
    ##  $ VAL_MEAN     : num  0.000964 0.000851 0.000784 0.000761 0.000721 ...
    ##  $ VAL_MEDIAN   : num  0.00037 0.000332 0.000311 0.000305 0.000291 ...
    ##  $ VAL_LWR      : num  2.08e-05 1.90e-05 1.81e-05 1.81e-05 1.75e-05 ...
    ##  $ VAL_UPR      : num  0.00554 0.00487 0.0044 0.00427 0.00401 ...
    ##  $ LOCATION     : chr  "Country" "Country" "Country" "Country" ...
    ##  $ LOCATION_NAME: chr  "AFG" "AFG" "AFG" "AFG" ...
    ##  $ METRIC       : chr  "Number" "Number" "Number" "Number" ...

``` r
all_cnt_rt <- all_cnt_nr%>%left_join(pop_all, by=c("LOCATION_NAME"="ISO3","YEAR"="YEAR"))
all_cnt_rt$VAL_MEAN <- 1e5 * all_cnt_rt$VAL_MEAN / all_cnt_rt$POP
all_cnt_rt$VAL_MEDIAN <- 1e5 * all_cnt_rt$VAL_MEDIAN / all_cnt_rt$POP
all_cnt_rt$VAL_LWR <- 1e5 * all_cnt_rt$VAL_LWR / all_cnt_rt$POP
all_cnt_rt$VAL_UPR <- 1e5 * all_cnt_rt$VAL_UPR / all_cnt_rt$POP
all_cnt_rt$LOCATION <- "Country"
all_cnt_rt$METRIC <- "Rate"
all_cnt_rt$POP <- NULL
str(all_cnt_rt)
```

    ## 'data.frame':    4268 obs. of  8 variables:
    ##  $ YEAR         : int  2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 ...
    ##  $ VAL_MEAN     : num  0.00936 0.00822 0.00723 0.00638 0.00565 ...
    ##  $ VAL_MEDIAN   : num  0.00359 0.00321 0.00287 0.00256 0.00228 ...
    ##  $ VAL_LWR      : num  0.000202 0.000183 0.000167 0.000152 0.000137 ...
    ##  $ VAL_UPR      : num  0.0538 0.047 0.0406 0.0358 0.0314 ...
    ##  $ LOCATION     : chr  "Country" "Country" "Country" "Country" ...
    ##  $ LOCATION_NAME: chr  "AFG" "AFG" "AFG" "AFG" ...
    ##  $ METRIC       : chr  "Rate" "Rate" "Rate" "Rate" ...

``` r
## compile all
all_est <-
  rbind(all_glb_rt, all_glb_nr,
        all_reg_rt, all_reg_nr,
        all_sub_rt, all_sub_nr,
        all_cnt_rt, all_cnt_nr)
str(all_est)
```

    ## 'data.frame':    9592 obs. of  8 variables:
    ##  $ YEAR         : int  2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 ...
    ##  $ VAL_MEAN     : num  0.1159 0.1035 0.0933 0.0845 0.0765 ...
    ##  $ VAL_MEDIAN   : num  0.0496 0.0456 0.0421 0.0385 0.0353 ...
    ##  $ VAL_LWR      : num  0.01112 0.0105 0.00989 0.00914 0.00842 ...
    ##  $ VAL_UPR      : num  0.479 0.423 0.388 0.355 0.326 ...
    ##  $ LOCATION     : chr  "Global" "Global" "Global" "Global" ...
    ##  $ LOCATION_NAME: chr  "Global" "Global" "Global" "Global" ...
    ##  $ METRIC       : chr  "Rate" "Rate" "Rate" "Rate" ...

``` r
saveRDS(all_est, file = paste0(params$Dir, "/", params$all_est))
```

# Summarize predictions:

## Global

``` r
kable(
  caption = paste0("Global number of ",params$Pathogen," cases, 2010 vs 2020"),
  row.names = FALSE,
  subset(all_glb_nr, YEAR %in% c(2010, 2020))[, 1:5])
```

| YEAR | VAL_MEAN | VAL_MEDIAN |  VAL_LWR |  VAL_UPR |
|-----:|---------:|-----------:|---------:|---------:|
| 2010 | 17.16669 |   8.083843 | 1.657451 | 81.46895 |
| 2020 | 13.06033 |   4.473299 | 0.447186 | 68.42967 |

Global number of CKD4 Age 85-125 cases, 2010 vs 2020

## Regions

``` r
kbl(subset(all_reg_rt, YEAR == 2020)[,c(7,2:5)],
    align = c("l", "c", "c", "c"), row.names = FALSE,
    col.names = c("Region", "Mean", "Median", "Lower", "Upper"),
    caption=paste0("Incidence per 100.000 of ",params$Pathogen," by WHO region in 2020")) %>%
  kable_styling("striped", "hover")
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">

<caption>

Incidence per 100.000 of CKD4 Age 85-125 by WHO region in 2020
</caption>

<thead>

<tr>

<th style="text-align:left;">

Region
</th>

<th style="text-align:center;">

Mean
</th>

<th style="text-align:center;">

Median
</th>

<th style="text-align:center;">

Lower
</th>

<th style="text-align:left;">

Upper
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

AFR
</td>

<td style="text-align:center;">

0.0012106
</td>

<td style="text-align:center;">

0.0003633
</td>

<td style="text-align:center;">

0.0000158
</td>

<td style="text-align:left;">

0.0079237
</td>

</tr>

<tr>

<td style="text-align:left;">

AMR
</td>

<td style="text-align:center;">

0.0008635
</td>

<td style="text-align:center;">

0.0000436
</td>

<td style="text-align:center;">

0.0000003
</td>

<td style="text-align:left;">

0.0037594
</td>

</tr>

<tr>

<td style="text-align:left;">

EMR
</td>

<td style="text-align:center;">

0.0012106
</td>

<td style="text-align:center;">

0.0003633
</td>

<td style="text-align:center;">

0.0000158
</td>

<td style="text-align:left;">

0.0079237
</td>

</tr>

<tr>

<td style="text-align:left;">

EUR
</td>

<td style="text-align:center;">

0.0016046
</td>

<td style="text-align:center;">

0.0002777
</td>

<td style="text-align:center;">

0.0000128
</td>

<td style="text-align:left;">

0.0092320
</td>

</tr>

<tr>

<td style="text-align:left;">

SEAR
</td>

<td style="text-align:center;">

0.0683305
</td>

<td style="text-align:center;">

0.0042588
</td>

<td style="text-align:center;">

0.0000655
</td>

<td style="text-align:left;">

0.4344860
</td>

</tr>

<tr>

<td style="text-align:left;">

WPR
</td>

<td style="text-align:center;">

0.0322691
</td>

<td style="text-align:center;">

0.0138860
</td>

<td style="text-align:center;">

0.0013783
</td>

<td style="text-align:left;">

0.1658926
</td>

</tr>

</tbody>

</table>

``` r
kbl(subset(all_reg_nr, YEAR == 2020)[,c(7,2:5)],
    align = c("l", "c", "c", "c"), row.names = FALSE,
    col.names = c("Region", "Mean", "Median", "Lower", "Upper"),
    caption=paste0("Cases of ",params$Pathogen," by WHO region in 2020")) %>%
  kable_styling("striped", "hover")
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">

<caption>

Cases of CKD4 Age 85-125 by WHO region in 2020
</caption>

<thead>

<tr>

<th style="text-align:left;">

Region
</th>

<th style="text-align:center;">

Mean
</th>

<th style="text-align:center;">

Median
</th>

<th style="text-align:center;">

Lower
</th>

<th style="text-align:left;">

Upper
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

AFR
</td>

<td style="text-align:center;">

0.0215461
</td>

<td style="text-align:center;">

0.0064654
</td>

<td style="text-align:center;">

0.0002805
</td>

<td style="text-align:left;">

0.1410284
</td>

</tr>

<tr>

<td style="text-align:left;">

AMR
</td>

<td style="text-align:center;">

0.1000239
</td>

<td style="text-align:center;">

0.0050489
</td>

<td style="text-align:center;">

0.0000390
</td>

<td style="text-align:left;">

0.4354877
</td>

</tr>

<tr>

<td style="text-align:left;">

EMR
</td>

<td style="text-align:center;">

0.0177615
</td>

<td style="text-align:center;">

0.0053297
</td>

<td style="text-align:center;">

0.0002312
</td>

<td style="text-align:left;">

0.1162566
</td>

</tr>

<tr>

<td style="text-align:left;">

EUR
</td>

<td style="text-align:center;">

0.3102499
</td>

<td style="text-align:center;">

0.0536901
</td>

<td style="text-align:center;">

0.0024754
</td>

<td style="text-align:left;">

1.7850082
</td>

</tr>

<tr>

<td style="text-align:left;">

SEAR
</td>

<td style="text-align:center;">

5.5854938
</td>

<td style="text-align:center;">

0.3481215
</td>

<td style="text-align:center;">

0.0053522
</td>

<td style="text-align:left;">

35.5159071
</td>

</tr>

<tr>

<td style="text-align:left;">

WPR
</td>

<td style="text-align:center;">

7.0252572
</td>

<td style="text-align:center;">

3.0230937
</td>

<td style="text-align:center;">

0.3000620
</td>

<td style="text-align:left;">

36.1161796
</td>

</tr>

</tbody>

</table>

``` r
png(paste0(params$PlotDir, "/r_CI_2010.png"), width=480, height=480)
ggplot(subset(all_reg_rt, YEAR==2010),
       aes(y = VAL_MEAN, x = LOCATION_NAME)) +
  geom_pointrange(aes(ymin = VAL_LWR, ymax = VAL_UPR), size = 0.2) +
  coord_flip() +
  theme_bw() +
  scale_x_discrete(NULL, limits = rev(unique(all_reg_rt$LOCATION_NAME))) +
  scale_y_continuous(NULL) +
  ggtitle(paste0("Incidence per 100.000 of ", params$Pathogen, " by WHO region, 2010"))
dev.off()
```

    ## png 
    ##   2

``` r
setwd(params$Dir)
image <- paste0("03-estimate-CKD_v2_files/figure-gfm/r_CI_2010.png")
cat("![](",image,")")
```

![](03-estimate-CKD_v2_files/figure-gfm/r_CI_2010.png)

``` r
png(paste0(params$PlotDir, "/r_CI_2020.png"), width=480, height=480)
ggplot(subset(all_reg_rt, YEAR==2020),
       aes(y = VAL_MEAN, x = LOCATION_NAME)) +
  geom_pointrange(aes(ymin = VAL_LWR, ymax = VAL_UPR), size = 0.2) +
  coord_flip() +
  theme_bw() +
  scale_x_discrete(NULL, limits = rev(unique(all_reg_rt$LOCATION_NAME))) +
  scale_y_continuous(NULL) +
  ggtitle(paste0("Incidence per 100.000 of ", params$Pathogen, " by WHO region, 2020"))
dev.off()
```

    ## png 
    ##   2

``` r
setwd(params$Dir)
image <- paste0("03-estimate-CKD_v2_files/figure-gfm/r_CI_2020.png")
cat("![](",image,")")
```

![](03-estimate-CKD_v2_files/figure-gfm/r_CI_2020.png)

``` r
png(paste0(params$PlotDir, "/r_CASES_2010.png"), width=480, height=480)
ggplot(subset(all_reg_nr, YEAR==2010),
       aes(y = VAL_MEAN, x = LOCATION_NAME)) +
  geom_pointrange(aes(ymin = VAL_LWR, ymax = VAL_UPR), size = 0.2) +
  coord_flip() +
  theme_bw() +
  scale_x_discrete(NULL, limits = rev(unique(all_reg_nr$LOCATION_NAME))) +
  scale_y_continuous(NULL) +
  ggtitle(paste0("Cases of ", params$Pathogen, " by WHO region, 2010"))
dev.off()
```

    ## png 
    ##   2

``` r
setwd(params$Dir)
image <- paste0("03-estimate-CKD_v2_files/figure-gfm/r_CASES_2010.png")
cat("![](",image,")")
```

![](03-estimate-CKD_v2_files/figure-gfm/r_CASES_2010.png)

``` r
png(paste0(params$PlotDir, "/r_CASES_2020.png"), width=480, height=480)
ggplot(subset(all_reg_nr, YEAR==2020),
       aes(y = VAL_MEAN, x = LOCATION_NAME)) +
  geom_pointrange(aes(ymin = VAL_LWR, ymax = VAL_UPR), size = 0.2) +
  coord_flip() +
  theme_bw() +
  scale_x_discrete(NULL, limits = rev(unique(all_reg_nr$LOCATION_NAME))) +
  scale_y_continuous(NULL) +
  ggtitle(paste0("Cases of ", params$Pathogen, " by WHO region, 2020"))
dev.off()
```

    ## png 
    ##   2

``` r
setwd(params$Dir)
image <- paste0("03-estimate-CKD_v2_files/figure-gfm/r_CASES_2020.png")
cat("![](",image,")")
```

![](03-estimate-CKD_v2_files/figure-gfm/r_CASES_2020.png)

``` r
sim_all_reg <-
  merge(sim_all_reg,
        with(sim_all, aggregate(POP ~ REG2 + YEAR, FUN = sum)))
sim_all_reg_long <-
  pivot_longer(sim_all_reg, cols = starts_with("V"))
sim_all_reg_long$CASES <- sim_all_reg_long$value
```

``` r
png(paste0(params$PlotDir, "/r_hist_2010.png"), width=480, height=480)
ggplot(subset(sim_all_reg_long, YEAR==2010), aes(x = CASES)) +
  geom_density() +
  facet_wrap(~REG2) +
  theme_bw() +
  scale_x_log10() +
  ggtitle(paste0("Incidence per 100.000 of ", params$Pathogen, " by WHO region, 2010"))
dev.off()
```

    ## png 
    ##   2

``` r
setwd(params$Dir)
image <- paste0("03-estimate-CKD_v2_files/figure-gfm/r_hist_2010.png")
cat("![](",image,")")
```

![](03-estimate-CKD_v2_files/figure-gfm/r_hist_2010.png)

``` r
png(paste0(params$PlotDir, "/r_hist_2020_2010.png"), width=480, height=480)
ggplot(subset(sim_all_reg_long, YEAR==2010), aes(x = CASES)) +
  geom_density() +
  facet_wrap(~REG2) +
  theme_bw() +
  scale_x_log10() +
  ggtitle(paste0("Incidence per 100.000 of ", params$Pathogen, " by WHO region, 2010"))
dev.off()
```

    ## png 
    ##   2

``` r
setwd(params$Dir)
image <- paste0("03-estimate-CKD_v2_files/figure-gfm/r_hist_2020_2010.png")
cat("![](",image,")")
```

![](03-estimate-CKD_v2_files/figure-gfm/r_hist_2020_2010.png)

## Subregions

``` r
kbl(subset(all_sub_rt, YEAR == 2020)[,c(7,2:5)],
    align = c("l", "c", "c", "c"), row.names = FALSE,
    col.names = c("Region", "Mean", "Median", "Lower", "Upper"),
    caption=paste0("Incidence per 100.000 of ",params$Pathogen," by WHO subregion in 2020")) %>%
  kable_styling("striped", "hover")
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">

<caption>

Incidence per 100.000 of CKD4 Age 85-125 by WHO subregion in 2020
</caption>

<thead>

<tr>

<th style="text-align:left;">

Region
</th>

<th style="text-align:center;">

Mean
</th>

<th style="text-align:center;">

Median
</th>

<th style="text-align:center;">

Lower
</th>

<th style="text-align:left;">

Upper
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

AFRAB
</td>

<td style="text-align:center;">

0.0012106
</td>

<td style="text-align:center;">

0.0003633
</td>

<td style="text-align:center;">

0.0000158
</td>

<td style="text-align:left;">

0.0079237
</td>

</tr>

<tr>

<td style="text-align:left;">

AFRC
</td>

<td style="text-align:center;">

0.0012106
</td>

<td style="text-align:center;">

0.0003633
</td>

<td style="text-align:center;">

0.0000158
</td>

<td style="text-align:left;">

0.0079237
</td>

</tr>

<tr>

<td style="text-align:left;">

AFRD
</td>

<td style="text-align:center;">

0.0012106
</td>

<td style="text-align:center;">

0.0003633
</td>

<td style="text-align:center;">

0.0000158
</td>

<td style="text-align:left;">

0.0079237
</td>

</tr>

<tr>

<td style="text-align:left;">

AMRA
</td>

<td style="text-align:center;">

0.0003805
</td>

<td style="text-align:center;">

0.0000256
</td>

<td style="text-align:center;">

0.0000002
</td>

<td style="text-align:left;">

0.0022245
</td>

</tr>

<tr>

<td style="text-align:left;">

AMRB
</td>

<td style="text-align:center;">

0.0017811
</td>

<td style="text-align:center;">

0.0000376
</td>

<td style="text-align:center;">

0.0000002
</td>

<td style="text-align:left;">

0.0059634
</td>

</tr>

<tr>

<td style="text-align:left;">

AMRC
</td>

<td style="text-align:center;">

0.0004833
</td>

<td style="text-align:center;">

0.0000666
</td>

<td style="text-align:center;">

0.0000004
</td>

<td style="text-align:left;">

0.0035213
</td>

</tr>

<tr>

<td style="text-align:left;">

EMRA
</td>

<td style="text-align:center;">

0.0012106
</td>

<td style="text-align:center;">

0.0003633
</td>

<td style="text-align:center;">

0.0000158
</td>

<td style="text-align:left;">

0.0079237
</td>

</tr>

<tr>

<td style="text-align:left;">

EMRBC
</td>

<td style="text-align:center;">

0.0012106
</td>

<td style="text-align:center;">

0.0003633
</td>

<td style="text-align:center;">

0.0000158
</td>

<td style="text-align:left;">

0.0079237
</td>

</tr>

<tr>

<td style="text-align:left;">

EMRD
</td>

<td style="text-align:center;">

0.0012106
</td>

<td style="text-align:center;">

0.0003633
</td>

<td style="text-align:center;">

0.0000158
</td>

<td style="text-align:left;">

0.0079237
</td>

</tr>

<tr>

<td style="text-align:left;">

EURA
</td>

<td style="text-align:center;">

0.0019012
</td>

<td style="text-align:center;">

0.0002614
</td>

<td style="text-align:center;">

0.0000120
</td>

<td style="text-align:left;">

0.0109326
</td>

</tr>

<tr>

<td style="text-align:left;">

EURB
</td>

<td style="text-align:center;">

0.0006215
</td>

<td style="text-align:center;">

0.0001638
</td>

<td style="text-align:center;">

0.0000070
</td>

<td style="text-align:left;">

0.0039366
</td>

</tr>

<tr>

<td style="text-align:left;">

EURC
</td>

<td style="text-align:center;">

0.0006215
</td>

<td style="text-align:center;">

0.0001638
</td>

<td style="text-align:center;">

0.0000070
</td>

<td style="text-align:left;">

0.0039366
</td>

</tr>

<tr>

<td style="text-align:left;">

SEARB
</td>

<td style="text-align:center;">

0.1752916
</td>

<td style="text-align:center;">

0.0055419
</td>

<td style="text-align:center;">

0.0000666
</td>

<td style="text-align:left;">

0.8557718
</td>

</tr>

<tr>

<td style="text-align:left;">

SEARCD
</td>

<td style="text-align:center;">

0.0382153
</td>

<td style="text-align:center;">

0.0015824
</td>

<td style="text-align:center;">

0.0000204
</td>

<td style="text-align:left;">

0.2509158
</td>

</tr>

<tr>

<td style="text-align:left;">

WPRA
</td>

<td style="text-align:center;">

0.0673803
</td>

<td style="text-align:center;">

0.0195342
</td>

<td style="text-align:center;">

0.0010228
</td>

<td style="text-align:left;">

0.4041026
</td>

</tr>

<tr>

<td style="text-align:left;">

WPRB
</td>

<td style="text-align:center;">

0.0141880
</td>

<td style="text-align:center;">

0.0082626
</td>

<td style="text-align:center;">

0.0009491
</td>

<td style="text-align:left;">

0.0627363
</td>

</tr>

<tr>

<td style="text-align:left;">

WPRC
</td>

<td style="text-align:center;">

0.0071385
</td>

<td style="text-align:center;">

0.0018574
</td>

<td style="text-align:center;">

0.0000677
</td>

<td style="text-align:left;">

0.0459262
</td>

</tr>

</tbody>

</table>

``` r
kbl(subset(all_sub_nr, YEAR == 2020)[,c(7,2:5)],
    align = c("l", "c", "c", "c"), row.names = FALSE,
    col.names = c("Region", "Mean", "Median", "Lower", "Upper"),
    caption=paste0("Cases of ",params$Pathogen," by WHO sub region in 2020")) %>%
  kable_styling("striped", "hover")
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">

<caption>

Cases of CKD4 Age 85-125 by WHO sub region in 2020
</caption>

<thead>

<tr>

<th style="text-align:left;">

Region
</th>

<th style="text-align:center;">

Mean
</th>

<th style="text-align:center;">

Median
</th>

<th style="text-align:center;">

Lower
</th>

<th style="text-align:left;">

Upper
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

AFRAB
</td>

<td style="text-align:center;">

0.0046447
</td>

<td style="text-align:center;">

0.0013937
</td>

<td style="text-align:center;">

0.0000605
</td>

<td style="text-align:left;">

0.0304013
</td>

</tr>

<tr>

<td style="text-align:left;">

AFRC
</td>

<td style="text-align:center;">

0.0104871
</td>

<td style="text-align:center;">

0.0031469
</td>

<td style="text-align:center;">

0.0001365
</td>

<td style="text-align:left;">

0.0686427
</td>

</tr>

<tr>

<td style="text-align:left;">

AFRD
</td>

<td style="text-align:center;">

0.0064143
</td>

<td style="text-align:center;">

0.0019248
</td>

<td style="text-align:center;">

0.0000835
</td>

<td style="text-align:left;">

0.0419844
</td>

</tr>

<tr>

<td style="text-align:left;">

AMRA
</td>

<td style="text-align:center;">

0.0280414
</td>

<td style="text-align:center;">

0.0018893
</td>

<td style="text-align:center;">

0.0000126
</td>

<td style="text-align:left;">

0.1639213
</td>

</tr>

<tr>

<td style="text-align:left;">

AMRB
</td>

<td style="text-align:center;">

0.0708319
</td>

<td style="text-align:center;">

0.0014934
</td>

<td style="text-align:center;">

0.0000091
</td>

<td style="text-align:left;">

0.2371595
</td>

</tr>

<tr>

<td style="text-align:left;">

AMRC
</td>

<td style="text-align:center;">

0.0011505
</td>

<td style="text-align:center;">

0.0001586
</td>

<td style="text-align:center;">

0.0000009
</td>

<td style="text-align:left;">

0.0083836
</td>

</tr>

<tr>

<td style="text-align:left;">

EMRA
</td>

<td style="text-align:center;">

0.0008900
</td>

<td style="text-align:center;">

0.0002671
</td>

<td style="text-align:center;">

0.0000116
</td>

<td style="text-align:left;">

0.0058256
</td>

</tr>

<tr>

<td style="text-align:left;">

EMRBC
</td>

<td style="text-align:center;">

0.0145187
</td>

<td style="text-align:center;">

0.0043567
</td>

<td style="text-align:center;">

0.0001890
</td>

<td style="text-align:left;">

0.0950310
</td>

</tr>

<tr>

<td style="text-align:left;">

EMRD
</td>

<td style="text-align:center;">

0.0023528
</td>

<td style="text-align:center;">

0.0007060
</td>

<td style="text-align:center;">

0.0000306
</td>

<td style="text-align:left;">

0.0153999
</td>

</tr>

<tr>

<td style="text-align:left;">

EURA
</td>

<td style="text-align:center;">

0.2823981
</td>

<td style="text-align:center;">

0.0388307
</td>

<td style="text-align:center;">

0.0017842
</td>

<td style="text-align:left;">

1.6238723
</td>

</tr>

<tr>

<td style="text-align:left;">

EURB
</td>

<td style="text-align:center;">

0.0220885
</td>

<td style="text-align:center;">

0.0058209
</td>

<td style="text-align:center;">

0.0002491
</td>

<td style="text-align:left;">

0.1399115
</td>

</tr>

<tr>

<td style="text-align:left;">

EURC
</td>

<td style="text-align:center;">

0.0057633
</td>

<td style="text-align:center;">

0.0015188
</td>

<td style="text-align:center;">

0.0000650
</td>

<td style="text-align:left;">

0.0365056
</td>

</tr>

<tr>

<td style="text-align:left;">

SEARB
</td>

<td style="text-align:center;">

3.1479708
</td>

<td style="text-align:center;">

0.0995241
</td>

<td style="text-align:center;">

0.0011966
</td>

<td style="text-align:left;">

15.3683601
</td>

</tr>

<tr>

<td style="text-align:left;">

SEARCD
</td>

<td style="text-align:center;">

2.4375231
</td>

<td style="text-align:center;">

0.1009347
</td>

<td style="text-align:center;">

0.0013018
</td>

<td style="text-align:left;">

16.0043853
</td>

</tr>

<tr>

<td style="text-align:left;">

WPRA
</td>

<td style="text-align:center;">

5.0742989
</td>

<td style="text-align:center;">

1.4710872
</td>

<td style="text-align:center;">

0.0770262
</td>

<td style="text-align:left;">

30.4323052
</td>

</tr>

<tr>

<td style="text-align:left;">

WPRB
</td>

<td style="text-align:center;">

1.8806750
</td>

<td style="text-align:center;">

1.0952418
</td>

<td style="text-align:center;">

0.1258011
</td>

<td style="text-align:left;">

8.3159571
</td>

</tr>

<tr>

<td style="text-align:left;">

WPRC
</td>

<td style="text-align:center;">

0.0702833
</td>

<td style="text-align:center;">

0.0182873
</td>

<td style="text-align:center;">

0.0006668
</td>

<td style="text-align:left;">

0.4521727
</td>

</tr>

</tbody>

</table>

``` r
png(paste0(params$PlotDir, "/r_CI_SUB2_2010.png"), width=480, height=480)
ggplot(subset(all_sub_rt, YEAR==2010),
       aes(y = VAL_MEAN, x = LOCATION_NAME)) +
  geom_pointrange(aes(ymin = VAL_LWR, ymax = VAL_UPR), size = 0.2) +
  coord_flip() +
  theme_bw() +
  scale_x_discrete(NULL, limits = rev(unique(all_sub_rt$LOCATION_NAME))) +
  scale_y_continuous(NULL) +
  ggtitle(paste0("Incidence per 100.000 of ", params$Pathogen, " by WHO sub region, 2010"))
dev.off()
```

    ## png 
    ##   2

``` r
setwd(params$Dir)
image <- paste0("03-estimate-CKD_v2_files/figure-gfm/r_CI_SUB2_2010.png")
cat("![](",image,")")
```

![](03-estimate-CKD_v2_files/figure-gfm/r_CI_SUB2_2010.png)

``` r
png(paste0(params$PlotDir, "/r_CI_SUB2_2020.png"), width=480, height=480)
ggplot(subset(all_sub_rt, YEAR==2020),
       aes(y = VAL_MEAN, x = LOCATION_NAME)) +
  geom_pointrange(aes(ymin = VAL_LWR, ymax = VAL_UPR), size = 0.2) +
  coord_flip() +
  theme_bw() +
  scale_x_discrete(NULL, limits = rev(unique(all_sub_rt$LOCATION_NAME))) +
  scale_y_continuous(NULL) +
  ggtitle(paste0("Incidence per 100.000 of ", params$Pathogen, " by WHO sub region, 2020"))
dev.off()
```

    ## png 
    ##   2

``` r
setwd(params$Dir)
image <- paste0("03-estimate-CKD_v2_files/figure-gfm/r_CI_SUB2_2020.png")
cat("![](",image,")")
```

![](03-estimate-CKD_v2_files/figure-gfm/r_CI_SUB2_2020.png)

``` r
png(paste0(params$PlotDir, "/r_CASES_SUB2_2010.png"), width=480, height=480)
ggplot(subset(all_sub_nr, YEAR==2010),
       aes(y = VAL_MEAN, x = LOCATION_NAME)) +
  geom_pointrange(aes(ymin = VAL_LWR, ymax = VAL_UPR), size = 0.2) +
  coord_flip() +
  theme_bw() +
  scale_x_discrete(NULL, limits = rev(unique(all_sub_nr$LOCATION_NAME))) +
  scale_y_continuous(NULL) +
  ggtitle(paste0("Cases of ", params$Pathogen, " by WHO sub region, 2010"))
dev.off()
```

    ## png 
    ##   2

``` r
setwd(params$Dir)
image <- paste0("03-estimate-CKD_v2_files/figure-gfm/r_CASES_SUB2_2010.png")
cat("![](",image,")")
```

![](03-estimate-CKD_v2_files/figure-gfm/r_CASES_SUB2_2010.png)

``` r
png(paste0(params$PlotDir, "/r_CASES_SUB2_2020.png"), width=480, height=480)
ggplot(subset(all_sub_nr, YEAR==2020),
       aes(y = VAL_MEAN, x = LOCATION_NAME)) +
  geom_pointrange(aes(ymin = VAL_LWR, ymax = VAL_UPR), size = 0.2) +
  coord_flip() +
  theme_bw() +
  scale_x_discrete(NULL, limits = rev(unique(all_sub_nr$LOCATION_NAME))) +
  scale_y_continuous(NULL) +
  ggtitle(paste0("Cases of ", params$Pathogen, " by WHO sub region, 2020"))
dev.off()
```

    ## png 
    ##   2

``` r
setwd(params$Dir)
image <- paste0("03-estimate-CKD_v2_files/figure-gfm/r_CASES_SUB2_2020.png")
cat("![](",image,")")
```

![](03-estimate-CKD_v2_files/figure-gfm/r_CASES_SUB2_2020.png)

``` r
sim_all_sub <-
  merge(sim_all_sub,
        with(sim_all, aggregate(POP ~ SUB2 + YEAR, FUN = sum)))
sim_all_sub_long <-
  pivot_longer(sim_all_sub, cols = starts_with("V"))
sim_all_sub_long$CASES <- sim_all_sub_long$value
```

``` r
png(paste0(params$PlotDir, "/r_hist_SUB2_2010.png"), width=480, height=480)
ggplot(subset(sim_all_sub_long, YEAR==2010), aes(x = CASES)) +
  geom_density() +
  facet_wrap(~SUB2) +
  theme_bw() +
  scale_x_log10() +
  ggtitle(paste0("Incidence per 100.000 of ", params$Pathogen, "by WHO sub region, 2010"))
dev.off()
```

    ## png 
    ##   2

``` r
setwd(params$Dir)
image <- paste0("03-estimate-CKD_v2_files/figure-gfm/r_hist_SUB2_2010.png")
cat("![](",image,")")
```

![](03-estimate-CKD_v2_files/figure-gfm/r_hist_SUB2_2010.png)

``` r
png(paste0(params$PlotDir, "/r_hist_SUB2_2020_2010.png"), width=480, height=480)
ggplot(subset(sim_all_sub_long, YEAR==2010), aes(x = CASES)) +
  geom_density() +
  facet_wrap(~SUB2) +
  theme_bw() +
  scale_x_log10() +
  ggtitle(paste0("Incidence per 100.000 of ", params$Pathogen, "by WHO sub region, 2010"))
dev.off()
```

    ## png 
    ##   2

``` r
setwd(params$Dir)
image <- paste0("03-estimate-CKD_v2_files/figure-gfm/r_hist_SUB2_2020_2010.png")
cat("![](",image,")")
```

![](03-estimate-CKD_v2_files/figure-gfm/r_hist_SUB2_2020_2010.png)

## Countries

``` r
png(paste0(params$PlotDir, "/r_cnt_2010.png"), width=800, height=300)
plot_world(subset(all_cnt_rt, YEAR == 2010),
           "LOCATION_NAME", "VAL_MEAN", legend.title = "Incidence per 100k", diseasefree = zero_cases)
```

    ## [1] 0.0 0.2 0.4 0.6 0.8

``` r
dev.off()
```

    ## png 
    ##   2

``` r
setwd(params$Dir)
image <- paste0("03-estimate-CKD_v2_files/figure-gfm/r_cnt_2010.png")
cat("![](",image,")")
```

![](03-estimate-CKD_v2_files/figure-gfm/r_cnt_2010.png)

``` r
png(paste0(params$PlotDir, "/r_cnt_2020.png"), width=800, height=300)
plot_world(subset(all_cnt_rt, YEAR == 2020),
           "LOCATION_NAME", "VAL_MEAN", legend.title = "Incidence per 100k", diseasefree = zero_cases)
```

    ## [1] 0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35

``` r
dev.off()
```

    ## png 
    ##   2

``` r
setwd(params$Dir)
image <- paste0("03-estimate-CKD_v2_files/figure-gfm/r_cnt_2020.png")
cat("![](",image,")")
```

![](03-estimate-CKD_v2_files/figure-gfm/r_cnt_2020.png)

``` r
tab <-
  data.frame(subset(all_cnt_rt, YEAR == 2010)[,
                                              c("LOCATION_NAME", "VAL_MEAN", "VAL_MEDIAN", "VAL_LWR", "VAL_UPR")],
             subset(all_cnt_rt, YEAR == 2020)[,
                                              c("VAL_MEAN", "VAL_MEDIAN", "VAL_LWR", "VAL_UPR")])
tab$LOCATION_NAME <-
  FERG2:::countries$COUNTRY[match(tab$LOCATION_NAME, FERG2:::countries$ISO3)]
tab$LOCATION_NAME <- gsub(" \\(.*", "", tab$LOCATION_NAME)
names(tab) <-
  c("Country",
    "2010.mean", "2010.median", "2010.lwr", "2010.upr",
    "2020.mean", "2020.median", "2020.lwr", "2020.upr")

kable(tab, digits = 3, row.names = FALSE,
      caption = paste0("Estimated ", params$Pathogen, " incidence by country, 2010 vs 2020"))
```

| Country | 2010.mean | 2010.median | 2010.lwr | 2010.upr | 2020.mean | 2020.median | 2020.lwr | 2020.upr |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| Afghanistan | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Angola | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Albania | 0.001 | 0.001 | 0.000 | 0.008 | 0.001 | 0.000 | 0.000 | 0.004 |
| Andorra | 0.001 | 0.000 | 0.000 | 0.003 | 0.000 | 0.000 | 0.000 | 0.002 |
| United Arab Emirates | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Argentina | 0.001 | 0.000 | 0.000 | 0.006 | 0.000 | 0.000 | 0.000 | 0.002 |
| Armenia | 0.001 | 0.001 | 0.000 | 0.008 | 0.001 | 0.000 | 0.000 | 0.004 |
| Antigua and Barbuda | 0.001 | 0.000 | 0.000 | 0.008 | 0.000 | 0.000 | 0.000 | 0.004 |
| Australia | 0.024 | 0.005 | 0.000 | 0.156 | 0.009 | 0.002 | 0.000 | 0.064 |
| Austria | 0.001 | 0.000 | 0.000 | 0.003 | 0.000 | 0.000 | 0.000 | 0.002 |
| Azerbaijan | 0.001 | 0.001 | 0.000 | 0.008 | 0.001 | 0.000 | 0.000 | 0.004 |
| Burundi | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Belgium | 0.012 | 0.001 | 0.000 | 0.059 | 0.008 | 0.000 | 0.000 | 0.029 |
| Benin | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Burkina Faso | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Bangladesh | 0.052 | 0.003 | 0.000 | 0.286 | 0.018 | 0.001 | 0.000 | 0.100 |
| Bulgaria | 0.001 | 0.001 | 0.000 | 0.008 | 0.001 | 0.000 | 0.000 | 0.004 |
| Bahrain | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Bahamas | 0.001 | 0.000 | 0.000 | 0.008 | 0.000 | 0.000 | 0.000 | 0.004 |
| Bosnia and Herzegovina | 0.001 | 0.001 | 0.000 | 0.008 | 0.001 | 0.000 | 0.000 | 0.004 |
| Belarus | 0.001 | 0.001 | 0.000 | 0.008 | 0.001 | 0.000 | 0.000 | 0.004 |
| Belize | 0.001 | 0.000 | 0.000 | 0.006 | 0.000 | 0.000 | 0.000 | 0.002 |
| Bolivia | 0.001 | 0.000 | 0.000 | 0.008 | 0.000 | 0.000 | 0.000 | 0.004 |
| Brazil | 0.000 | 0.000 | 0.000 | 0.001 | 0.000 | 0.000 | 0.000 | 0.000 |
| Barbados | 0.001 | 0.000 | 0.000 | 0.008 | 0.000 | 0.000 | 0.000 | 0.004 |
| Brunei Darussalam | 0.026 | 0.010 | 0.000 | 0.157 | 0.011 | 0.003 | 0.000 | 0.071 |
| Bhutan | 0.052 | 0.003 | 0.000 | 0.286 | 0.018 | 0.001 | 0.000 | 0.100 |
| Botswana | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Central African Republic | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Canada | 0.001 | 0.000 | 0.000 | 0.008 | 0.000 | 0.000 | 0.000 | 0.004 |
| Switzerland | 0.001 | 0.000 | 0.000 | 0.003 | 0.000 | 0.000 | 0.000 | 0.002 |
| Chile | 0.001 | 0.000 | 0.000 | 0.008 | 0.000 | 0.000 | 0.000 | 0.004 |
| China | 0.034 | 0.025 | 0.005 | 0.117 | 0.014 | 0.008 | 0.001 | 0.062 |
| Côte d’Ivoire | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Cameroon | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Congo | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Congo | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Cook Islands | 0.026 | 0.010 | 0.000 | 0.157 | 0.011 | 0.003 | 0.000 | 0.071 |
| Colombia | 0.001 | 0.000 | 0.000 | 0.006 | 0.000 | 0.000 | 0.000 | 0.002 |
| Comoros | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Cabo Verde | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Costa Rica | 0.001 | 0.000 | 0.000 | 0.006 | 0.000 | 0.000 | 0.000 | 0.002 |
| Cuba | 0.001 | 0.000 | 0.000 | 0.006 | 0.000 | 0.000 | 0.000 | 0.002 |
| Cyprus | 0.001 | 0.000 | 0.000 | 0.003 | 0.000 | 0.000 | 0.000 | 0.002 |
| Czechia | 0.001 | 0.000 | 0.000 | 0.003 | 0.000 | 0.000 | 0.000 | 0.002 |
| Germany | 0.014 | 0.001 | 0.000 | 0.091 | 0.008 | 0.000 | 0.000 | 0.042 |
| Djibouti | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Dominica | 0.001 | 0.000 | 0.000 | 0.006 | 0.000 | 0.000 | 0.000 | 0.002 |
| Denmark | 0.001 | 0.000 | 0.000 | 0.003 | 0.000 | 0.000 | 0.000 | 0.002 |
| Dominican Republic | 0.001 | 0.000 | 0.000 | 0.006 | 0.000 | 0.000 | 0.000 | 0.002 |
| Algeria | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Ecuador | 0.001 | 0.000 | 0.000 | 0.006 | 0.000 | 0.000 | 0.000 | 0.002 |
| Egypt | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Eritrea | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Spain | 0.001 | 0.000 | 0.000 | 0.005 | 0.000 | 0.000 | 0.000 | 0.003 |
| Estonia | 0.001 | 0.000 | 0.000 | 0.003 | 0.000 | 0.000 | 0.000 | 0.002 |
| Ethiopia | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Finland | 0.001 | 0.000 | 0.000 | 0.003 | 0.000 | 0.000 | 0.000 | 0.002 |
| Fiji | 0.028 | 0.008 | 0.000 | 0.151 | 0.011 | 0.003 | 0.000 | 0.068 |
| France | 0.000 | 0.000 | 0.000 | 0.003 | 0.000 | 0.000 | 0.000 | 0.001 |
| Micronesia | 0.017 | 0.006 | 0.000 | 0.103 | 0.007 | 0.002 | 0.000 | 0.046 |
| Gabon | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| United Kingdom | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Georgia | 0.001 | 0.001 | 0.000 | 0.008 | 0.001 | 0.000 | 0.000 | 0.004 |
| Ghana | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Guinea | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Gambia | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Guinea-Bissau | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Equatorial Guinea | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Greece | 0.001 | 0.000 | 0.000 | 0.003 | 0.000 | 0.000 | 0.000 | 0.002 |
| Grenada | 0.001 | 0.000 | 0.000 | 0.006 | 0.000 | 0.000 | 0.000 | 0.002 |
| Guatemala | 0.001 | 0.000 | 0.000 | 0.006 | 0.000 | 0.000 | 0.000 | 0.002 |
| Guyana | 0.001 | 0.000 | 0.000 | 0.008 | 0.000 | 0.000 | 0.000 | 0.004 |
| Honduras | 0.001 | 0.000 | 0.000 | 0.008 | 0.000 | 0.000 | 0.000 | 0.004 |
| Croatia | 0.001 | 0.000 | 0.000 | 0.003 | 0.000 | 0.000 | 0.000 | 0.002 |
| Haiti | 0.001 | 0.000 | 0.000 | 0.008 | 0.000 | 0.000 | 0.000 | 0.004 |
| Hungary | 0.001 | 0.000 | 0.000 | 0.003 | 0.000 | 0.000 | 0.000 | 0.002 |
| Indonesia | 0.129 | 0.004 | 0.000 | 0.446 | 0.050 | 0.001 | 0.000 | 0.166 |
| India | 0.115 | 0.004 | 0.000 | 0.720 | 0.042 | 0.001 | 0.000 | 0.262 |
| Ireland | 0.001 | 0.000 | 0.000 | 0.003 | 0.000 | 0.000 | 0.000 | 0.002 |
| Iran | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Iraq | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Iceland | 0.001 | 0.000 | 0.000 | 0.003 | 0.000 | 0.000 | 0.000 | 0.002 |
| Israel | 0.001 | 0.000 | 0.000 | 0.003 | 0.000 | 0.000 | 0.000 | 0.002 |
| Italy | 0.001 | 0.000 | 0.000 | 0.009 | 0.001 | 0.000 | 0.000 | 0.004 |
| Jamaica | 0.001 | 0.000 | 0.000 | 0.006 | 0.000 | 0.000 | 0.000 | 0.002 |
| Jordan | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Japan | 0.135 | 0.062 | 0.006 | 0.695 | 0.072 | 0.020 | 0.001 | 0.458 |
| Kazakhstan | 0.001 | 0.001 | 0.000 | 0.008 | 0.001 | 0.000 | 0.000 | 0.004 |
| Kenya | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Kyrgyzstan | 0.001 | 0.001 | 0.000 | 0.008 | 0.001 | 0.000 | 0.000 | 0.004 |
| Cambodia | 0.017 | 0.006 | 0.000 | 0.103 | 0.007 | 0.002 | 0.000 | 0.046 |
| Kiribati | 0.017 | 0.006 | 0.000 | 0.103 | 0.007 | 0.002 | 0.000 | 0.046 |
| Saint Kitts and Nevis | 0.001 | 0.000 | 0.000 | 0.008 | 0.000 | 0.000 | 0.000 | 0.004 |
| Korea | 0.142 | 0.024 | 0.001 | 0.781 | 0.078 | 0.008 | 0.000 | 0.367 |
| Kuwait | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Lao People’s Dem. Republic | 0.017 | 0.006 | 0.000 | 0.103 | 0.007 | 0.002 | 0.000 | 0.046 |
| Lebanon | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Liberia | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Libya | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Saint Lucia | 0.001 | 0.000 | 0.000 | 0.006 | 0.000 | 0.000 | 0.000 | 0.002 |
| Sri Lanka | 0.052 | 0.003 | 0.000 | 0.286 | 0.018 | 0.001 | 0.000 | 0.100 |
| Lesotho | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Lithuania | 0.001 | 0.000 | 0.000 | 0.003 | 0.000 | 0.000 | 0.000 | 0.002 |
| Luxembourg | 0.001 | 0.000 | 0.000 | 0.003 | 0.000 | 0.000 | 0.000 | 0.002 |
| Latvia | 0.001 | 0.000 | 0.000 | 0.003 | 0.000 | 0.000 | 0.000 | 0.002 |
| Morocco | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Monaco | 0.001 | 0.000 | 0.000 | 0.003 | 0.000 | 0.000 | 0.000 | 0.002 |
| Republic of Moldova | 0.001 | 0.001 | 0.000 | 0.008 | 0.001 | 0.000 | 0.000 | 0.004 |
| Madagascar | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Maldives | 0.129 | 0.004 | 0.000 | 0.446 | 0.050 | 0.001 | 0.000 | 0.166 |
| Mexico | 0.019 | 0.000 | 0.000 | 0.057 | 0.008 | 0.000 | 0.000 | 0.022 |
| Marshall Islands | 0.028 | 0.008 | 0.000 | 0.151 | 0.011 | 0.003 | 0.000 | 0.068 |
| North Macedonia | 0.001 | 0.001 | 0.000 | 0.008 | 0.001 | 0.000 | 0.000 | 0.004 |
| Mali | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Malta | 0.001 | 0.000 | 0.000 | 0.003 | 0.000 | 0.000 | 0.000 | 0.002 |
| Myanmar | 0.052 | 0.003 | 0.000 | 0.286 | 0.018 | 0.001 | 0.000 | 0.100 |
| Montenegro | 0.001 | 0.001 | 0.000 | 0.008 | 0.001 | 0.000 | 0.000 | 0.004 |
| Mongolia | 0.017 | 0.006 | 0.000 | 0.103 | 0.007 | 0.002 | 0.000 | 0.046 |
| Mozambique | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Mauritania | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Mauritius | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Malawi | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Malaysia | 0.028 | 0.008 | 0.000 | 0.151 | 0.011 | 0.003 | 0.000 | 0.068 |
| Namibia | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Niger | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Nigeria | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Nicaragua | 0.001 | 0.000 | 0.000 | 0.008 | 0.000 | 0.000 | 0.000 | 0.004 |
| Niue | 0.026 | 0.010 | 0.000 | 0.157 | 0.011 | 0.003 | 0.000 | 0.071 |
| Netherlands | 0.020 | 0.001 | 0.000 | 0.107 | 0.008 | 0.000 | 0.000 | 0.045 |
| Norway | 0.001 | 0.000 | 0.000 | 0.003 | 0.000 | 0.000 | 0.000 | 0.002 |
| Nepal | 0.052 | 0.003 | 0.000 | 0.286 | 0.018 | 0.001 | 0.000 | 0.100 |
| Nauru | 0.026 | 0.010 | 0.000 | 0.157 | 0.011 | 0.003 | 0.000 | 0.071 |
| New Zealand | 0.026 | 0.010 | 0.000 | 0.157 | 0.011 | 0.003 | 0.000 | 0.071 |
| Oman | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Pakistan | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Panama | 0.001 | 0.000 | 0.000 | 0.008 | 0.000 | 0.000 | 0.000 | 0.004 |
| Peru | 0.001 | 0.000 | 0.000 | 0.006 | 0.000 | 0.000 | 0.000 | 0.002 |
| Philippines | 0.017 | 0.006 | 0.000 | 0.103 | 0.007 | 0.002 | 0.000 | 0.046 |
| Palau | 0.028 | 0.008 | 0.000 | 0.151 | 0.011 | 0.003 | 0.000 | 0.068 |
| Papua New Guinea | 0.017 | 0.006 | 0.000 | 0.103 | 0.007 | 0.002 | 0.000 | 0.046 |
| Poland | 0.001 | 0.000 | 0.000 | 0.003 | 0.000 | 0.000 | 0.000 | 0.002 |
| Korea | 0.052 | 0.003 | 0.000 | 0.286 | 0.018 | 0.001 | 0.000 | 0.100 |
| Portugal | 0.001 | 0.000 | 0.000 | 0.003 | 0.000 | 0.000 | 0.000 | 0.002 |
| Paraguay | 0.001 | 0.000 | 0.000 | 0.006 | 0.000 | 0.000 | 0.000 | 0.002 |
| Qatar | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Romania | 0.001 | 0.000 | 0.000 | 0.003 | 0.000 | 0.000 | 0.000 | 0.002 |
| Russian Federation | 0.001 | 0.001 | 0.000 | 0.008 | 0.001 | 0.000 | 0.000 | 0.004 |
| Rwanda | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Saudi Arabia | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Sudan | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Senegal | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Singapore | 0.026 | 0.010 | 0.000 | 0.157 | 0.011 | 0.003 | 0.000 | 0.071 |
| Solomon Islands | 0.017 | 0.006 | 0.000 | 0.103 | 0.007 | 0.002 | 0.000 | 0.046 |
| Sierra Leone | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| El Salvador | 0.001 | 0.000 | 0.000 | 0.006 | 0.000 | 0.000 | 0.000 | 0.002 |
| San Marino | 0.001 | 0.000 | 0.000 | 0.003 | 0.000 | 0.000 | 0.000 | 0.002 |
| Somalia | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Serbia | 0.001 | 0.001 | 0.000 | 0.008 | 0.001 | 0.000 | 0.000 | 0.004 |
| South Sudan | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Sao Tome and Principe | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Suriname | 0.001 | 0.000 | 0.000 | 0.006 | 0.000 | 0.000 | 0.000 | 0.002 |
| Slovakia | 0.001 | 0.000 | 0.000 | 0.003 | 0.000 | 0.000 | 0.000 | 0.002 |
| Slovenia | 0.001 | 0.000 | 0.000 | 0.003 | 0.000 | 0.000 | 0.000 | 0.002 |
| Sweden | 0.007 | 0.001 | 0.000 | 0.051 | 0.004 | 0.000 | 0.000 | 0.024 |
| Eswatini | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Seychelles | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Syrian Arab Republic | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Chad | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Togo | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Thailand | 0.769 | 0.022 | 0.000 | 3.801 | 0.304 | 0.007 | 0.000 | 1.439 |
| Tajikistan | 0.001 | 0.001 | 0.000 | 0.008 | 0.001 | 0.000 | 0.000 | 0.004 |
| Turkmenistan | 0.001 | 0.001 | 0.000 | 0.008 | 0.001 | 0.000 | 0.000 | 0.004 |
| Timor-Leste | 0.052 | 0.003 | 0.000 | 0.286 | 0.018 | 0.001 | 0.000 | 0.100 |
| Tonga | 0.028 | 0.008 | 0.000 | 0.151 | 0.011 | 0.003 | 0.000 | 0.068 |
| Trinidad and Tobago | 0.001 | 0.000 | 0.000 | 0.008 | 0.000 | 0.000 | 0.000 | 0.004 |
| Tunisia | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Turkiye | 0.001 | 0.001 | 0.000 | 0.008 | 0.001 | 0.000 | 0.000 | 0.004 |
| Tuvalu | 0.028 | 0.008 | 0.000 | 0.151 | 0.011 | 0.003 | 0.000 | 0.068 |
| United Republic of Tanzania | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Uganda | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Ukraine | 0.001 | 0.001 | 0.000 | 0.008 | 0.001 | 0.000 | 0.000 | 0.004 |
| Uruguay | 0.001 | 0.000 | 0.000 | 0.008 | 0.000 | 0.000 | 0.000 | 0.004 |
| United States of America | 0.001 | 0.000 | 0.000 | 0.005 | 0.000 | 0.000 | 0.000 | 0.002 |
| Uzbekistan | 0.001 | 0.001 | 0.000 | 0.008 | 0.001 | 0.000 | 0.000 | 0.004 |
| Saint Vincent and the Grenadines | 0.001 | 0.000 | 0.000 | 0.006 | 0.000 | 0.000 | 0.000 | 0.002 |
| Venezuela | 0.001 | 0.000 | 0.000 | 0.008 | 0.000 | 0.000 | 0.000 | 0.004 |
| Viet Nam | 0.017 | 0.006 | 0.000 | 0.103 | 0.007 | 0.002 | 0.000 | 0.046 |
| Vanuatu | 0.017 | 0.006 | 0.000 | 0.103 | 0.007 | 0.002 | 0.000 | 0.046 |
| Samoa | 0.017 | 0.006 | 0.000 | 0.103 | 0.007 | 0.002 | 0.000 | 0.046 |
| Yemen | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| South Africa | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Zambia | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |
| Zimbabwe | 0.003 | 0.001 | 0.000 | 0.016 | 0.001 | 0.000 | 0.000 | 0.008 |

Estimated CKD4 Age 85-125 incidence by country, 2010 vs 2020

``` r
tab2 <-
  data.frame(subset(all_cnt_nr, YEAR == 2010)[,
                                              c("LOCATION_NAME", "VAL_MEAN", "VAL_MEDIAN", "VAL_LWR", "VAL_UPR")],
             subset(all_cnt_nr, YEAR == 2020)[,
                                              c("VAL_MEAN", "VAL_MEDIAN", "VAL_LWR", "VAL_UPR")])
tab2$LOCATION_NAME <-
  FERG2:::countries$COUNTRY[match(tab2$LOCATION_NAME, FERG2:::countries$ISO3)]
tab2$LOCATION_NAME <- gsub(" \\(.*", "", tab2$LOCATION_NAME)
names(tab2) <-
  c("Country",
    "2010.mean", "2010.median", "2010.lwr", "2010.upr",
    "2020.mean", "2020.median", "2020.lwr", "2020.upr")

kable(tab2, digits = 1, row.names = FALSE,
      caption = paste0("Estimated ", params$Pathogen, " cases by country, 2010 vs 2020"))
```

| Country | 2010.mean | 2010.median | 2010.lwr | 2010.upr | 2020.mean | 2020.median | 2020.lwr | 2020.upr |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| Afghanistan | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Angola | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Albania | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Andorra | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| United Arab Emirates | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Argentina | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Armenia | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Antigua and Barbuda | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Australia | 0.1 | 0.0 | 0.0 | 0.6 | 0.0 | 0.0 | 0.0 | 0.3 |
| Austria | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Azerbaijan | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Burundi | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Belgium | 0.0 | 0.0 | 0.0 | 0.1 | 0.0 | 0.0 | 0.0 | 0.1 |
| Benin | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Burkina Faso | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Bangladesh | 0.1 | 0.0 | 0.0 | 0.5 | 0.1 | 0.0 | 0.0 | 0.5 |
| Bulgaria | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Bahrain | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Bahamas | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Bosnia and Herzegovina | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Belarus | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Belize | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Bolivia | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Brazil | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Barbados | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Brunei Darussalam | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Bhutan | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Botswana | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Central African Republic | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Canada | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Switzerland | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Chile | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| China | 2.7 | 2.0 | 0.4 | 9.4 | 1.9 | 1.1 | 0.1 | 8.2 |
| Côte d’Ivoire | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Cameroon | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Congo | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Congo | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Cook Islands | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Colombia | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Comoros | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Cabo Verde | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Costa Rica | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Cuba | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Cyprus | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Czechia | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Germany | 0.3 | 0.0 | 0.0 | 1.6 | 0.2 | 0.0 | 0.0 | 1.0 |
| Djibouti | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Dominica | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Denmark | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Dominican Republic | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Algeria | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Ecuador | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Egypt | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Eritrea | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Spain | 0.0 | 0.0 | 0.0 | 0.1 | 0.0 | 0.0 | 0.0 | 0.0 |
| Estonia | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Ethiopia | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Finland | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Fiji | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| France | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Micronesia | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Gabon | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| United Kingdom | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Georgia | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Ghana | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Guinea | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Gambia | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Guinea-Bissau | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Equatorial Guinea | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Greece | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Grenada | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Guatemala | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Guyana | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Honduras | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Croatia | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Haiti | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Hungary | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Indonesia | 0.7 | 0.0 | 0.0 | 2.3 | 0.5 | 0.0 | 0.0 | 1.5 |
| India | 3.7 | 0.1 | 0.0 | 22.8 | 2.3 | 0.1 | 0.0 | 14.1 |
| Ireland | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Iran | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Iraq | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Iceland | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Israel | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Italy | 0.0 | 0.0 | 0.0 | 0.1 | 0.0 | 0.0 | 0.0 | 0.1 |
| Jamaica | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Jordan | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Japan | 5.1 | 2.3 | 0.2 | 26.1 | 4.4 | 1.2 | 0.1 | 27.6 |
| Kazakhstan | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Kenya | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Kyrgyzstan | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Cambodia | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Kiribati | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Saint Kitts and Nevis | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Korea | 0.5 | 0.1 | 0.0 | 2.7 | 0.6 | 0.1 | 0.0 | 3.0 |
| Kuwait | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Lao People’s Dem. Republic | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Lebanon | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Liberia | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Libya | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Saint Lucia | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Sri Lanka | 0.0 | 0.0 | 0.0 | 0.3 | 0.0 | 0.0 | 0.0 | 0.1 |
| Lesotho | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Lithuania | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Luxembourg | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Latvia | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Morocco | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Monaco | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Republic of Moldova | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Madagascar | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Maldives | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Mexico | 0.1 | 0.0 | 0.0 | 0.4 | 0.1 | 0.0 | 0.0 | 0.2 |
| Marshall Islands | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| North Macedonia | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Mali | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Malta | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Myanmar | 0.1 | 0.0 | 0.0 | 0.3 | 0.0 | 0.0 | 0.0 | 0.1 |
| Montenegro | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Mongolia | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Mozambique | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Mauritania | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Mauritius | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Malawi | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Malaysia | 0.0 | 0.0 | 0.0 | 0.1 | 0.0 | 0.0 | 0.0 | 0.1 |
| Namibia | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Niger | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Nigeria | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Nicaragua | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Niue | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Netherlands | 0.1 | 0.0 | 0.0 | 0.3 | 0.0 | 0.0 | 0.0 | 0.2 |
| Norway | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Nepal | 0.0 | 0.0 | 0.0 | 0.1 | 0.0 | 0.0 | 0.0 | 0.1 |
| Nauru | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| New Zealand | 0.0 | 0.0 | 0.0 | 0.1 | 0.0 | 0.0 | 0.0 | 0.1 |
| Oman | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Pakistan | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Panama | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Peru | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Philippines | 0.0 | 0.0 | 0.0 | 0.1 | 0.0 | 0.0 | 0.0 | 0.1 |
| Palau | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Papua New Guinea | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Poland | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Korea | 0.0 | 0.0 | 0.0 | 0.2 | 0.0 | 0.0 | 0.0 | 0.2 |
| Portugal | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Paraguay | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Qatar | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Romania | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Russian Federation | 0.0 | 0.0 | 0.0 | 0.1 | 0.0 | 0.0 | 0.0 | 0.1 |
| Rwanda | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Saudi Arabia | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Sudan | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Senegal | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Singapore | 0.0 | 0.0 | 0.0 | 0.1 | 0.0 | 0.0 | 0.0 | 0.0 |
| Solomon Islands | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Sierra Leone | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| El Salvador | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| San Marino | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Somalia | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Serbia | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| South Sudan | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Sao Tome and Principe | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Suriname | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Slovakia | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Slovenia | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Sweden | 0.0 | 0.0 | 0.0 | 0.1 | 0.0 | 0.0 | 0.0 | 0.1 |
| Eswatini | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Seychelles | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Syrian Arab Republic | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Chad | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Togo | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Thailand | 3.3 | 0.1 | 0.0 | 16.2 | 2.7 | 0.1 | 0.0 | 12.7 |
| Tajikistan | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Turkmenistan | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Timor-Leste | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Tonga | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Trinidad and Tobago | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Tunisia | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Turkiye | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Tuvalu | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| United Republic of Tanzania | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Uganda | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Ukraine | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Uruguay | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| United States of America | 0.0 | 0.0 | 0.0 | 0.3 | 0.0 | 0.0 | 0.0 | 0.1 |
| Uzbekistan | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Saint Vincent and the Grenadines | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Venezuela | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Viet Nam | 0.1 | 0.0 | 0.0 | 0.5 | 0.0 | 0.0 | 0.0 | 0.3 |
| Vanuatu | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Samoa | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Yemen | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| South Africa | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Zambia | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Zimbabwe | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |

Estimated CKD4 Age 85-125 cases by country, 2010 vs 2020

# Session info

``` r
sessioninfo::session_info()
```

    ## Warning in system2("quarto", "-V", stdout = TRUE, env = paste0("TMPDIR=", : running command '"quarto"
    ## TMPDIR=C:/Users/LoVa3397/AppData/Local/Temp/RtmpOKNeDj/file24cc3ac01e2c -V' had status 1

    ## ─ Session info ────────────────────────────────────────────────────────────────────────────────────────────────────────
    ##  setting  value
    ##  version  R version 4.5.1 (2025-06-13 ucrt)
    ##  os       Windows 10 x64 (build 19045)
    ##  system   x86_64, mingw32
    ##  ui       RStudio
    ##  language (EN)
    ##  collate  English_United States.utf8
    ##  ctype    English_United States.utf8
    ##  tz       Europe/Brussels
    ##  date     2025-10-16
    ##  rstudio  2025.09.0+387 Cucumberleaf Sunflower (desktop)
    ##  pandoc   3.6.3 @ C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools/ (via rmarkdown)
    ##  quarto   ERROR: Unknown command "TMPDIR=C:/Users/LoVa3397/AppData/Local/Temp/RtmpOKNeDj/file24cc3ac01e2c". Did you mean command "install"? @ C:\\PROGRA~1\\RStudio\\RESOUR~1\\app\\bin\\quarto\\bin\\quarto.exe
    ## 
    ## ─ Packages ────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ##  ! package        * version    date (UTC) lib source
    ##    abind            1.4-8      2024-09-12 [1] CRAN (R 4.5.0)
    ##    backports        1.5.0      2024-05-23 [1] CRAN (R 4.5.0)
    ##    base64enc        0.1-3      2015-07-28 [1] CRAN (R 4.5.0)
    ##    bayesplot        1.13.0     2025-06-18 [1] CRAN (R 4.5.1)
    ##    bd             * 0.0.14     2025-07-14 [1] Github (brechtdv/bd@652191c)
    ##    boot             1.3-31     2024-08-28 [1] CRAN (R 4.5.1)
    ##    bridgesampling   1.1-2      2021-04-16 [1] CRAN (R 4.5.1)
    ##    brms           * 2.22.0     2024-09-23 [1] CRAN (R 4.5.1)
    ##    Brobdingnag      1.2-9      2022-10-19 [1] CRAN (R 4.5.1)
    ##    callr            3.7.6      2024-03-25 [1] CRAN (R 4.5.1)
    ##    cellranger       1.1.0      2016-07-27 [1] CRAN (R 4.5.1)
    ##    checkmate        2.3.2      2024-07-29 [1] CRAN (R 4.5.1)
    ##    class            7.3-23     2025-01-01 [1] CRAN (R 4.5.1)
    ##    classInt         0.4-11     2025-01-08 [1] CRAN (R 4.5.1)
    ##    cli              3.6.5      2025-04-23 [1] CRAN (R 4.5.1)
    ##    cluster          2.1.8.1    2025-03-12 [1] CRAN (R 4.5.1)
    ##    coda             0.19-4.1   2024-01-31 [1] CRAN (R 4.5.1)
    ##    codetools        0.2-20     2024-03-31 [1] CRAN (R 4.5.1)
    ##    colorspace       2.1-1      2024-07-26 [1] CRAN (R 4.5.1)
    ##    cowplot        * 1.2.0      2025-07-07 [1] CRAN (R 4.5.1)
    ##    curl             6.4.0      2025-06-22 [1] CRAN (R 4.5.1)
    ##    data.table       1.17.8     2025-07-10 [1] CRAN (R 4.5.1)
    ##    DBI              1.2.3      2024-06-02 [1] CRAN (R 4.5.1)
    ##    DescTools      * 0.99.60    2025-03-28 [1] CRAN (R 4.5.1)
    ##    digest           0.6.37     2024-08-19 [1] CRAN (R 4.5.1)
    ##    distributional   0.5.0      2024-09-17 [1] CRAN (R 4.5.1)
    ##    dplyr          * 1.1.4      2023-11-17 [1] CRAN (R 4.5.1)
    ##    e1071            1.7-16     2024-09-16 [1] CRAN (R 4.5.1)
    ##    evaluate         1.0.4      2025-06-18 [1] CRAN (R 4.5.1)
    ##    Exact            3.3        2024-07-21 [1] CRAN (R 4.5.0)
    ##    expm             1.0-0      2024-08-19 [1] CRAN (R 4.5.1)
    ##    farver           2.1.2      2024-05-13 [1] CRAN (R 4.5.1)
    ##    fastmap          1.2.0      2024-05-15 [1] CRAN (R 4.5.1)
    ##    FERG2          * 0.0.5      2025-07-15 [1] Github (brechtdv/FERG2@c2d4ac1)
    ##    forcats          1.0.0      2023-01-29 [1] CRAN (R 4.5.1)
    ##    foreign          0.8-90     2025-03-31 [1] CRAN (R 4.5.1)
    ##    Formula          1.2-5      2023-02-24 [1] CRAN (R 4.5.0)
    ##    fs               1.6.6      2025-04-12 [1] CRAN (R 4.5.1)
    ##    generics         0.1.4      2025-05-09 [1] CRAN (R 4.5.1)
    ##    ggplot2        * 3.5.2      2025-04-09 [1] CRAN (R 4.5.1)
    ##    gld              2.6.7      2025-01-17 [1] CRAN (R 4.5.1)
    ##    glue             1.8.0      2024-09-30 [1] CRAN (R 4.5.1)
    ##    gridExtra        2.3        2017-09-09 [1] CRAN (R 4.5.1)
    ##    gtable           0.3.6      2024-10-25 [1] CRAN (R 4.5.1)
    ##    haven            2.5.5      2025-05-30 [1] CRAN (R 4.5.1)
    ##    Hmisc          * 5.2-3      2025-03-16 [1] CRAN (R 4.5.1)
    ##    hms              1.1.3      2023-03-21 [1] CRAN (R 4.5.1)
    ##    htmlTable        2.4.3      2024-07-21 [1] CRAN (R 4.5.1)
    ##    htmltools        0.5.8.1    2024-04-04 [1] CRAN (R 4.5.1)
    ##    htmlwidgets      1.6.4      2023-12-06 [1] CRAN (R 4.5.1)
    ##    httr             1.4.7      2023-08-15 [1] CRAN (R 4.5.1)
    ##    inline           0.3.21     2025-01-09 [1] CRAN (R 4.5.1)
    ##    jsonlite         2.0.0      2025-03-27 [1] CRAN (R 4.5.1)
    ##    kableExtra     * 1.4.0      2024-01-24 [1] CRAN (R 4.5.1)
    ##    KernSmooth       2.23-26    2025-01-01 [1] CRAN (R 4.5.1)
    ##    knitr          * 1.50       2025-03-16 [1] CRAN (R 4.5.1)
    ##    labeling         0.4.3      2023-08-29 [1] CRAN (R 4.5.0)
    ##    lattice          0.22-7     2025-04-02 [1] CRAN (R 4.5.1)
    ##    lifecycle        1.0.4      2023-11-07 [1] CRAN (R 4.5.1)
    ##    lmom             3.2        2024-09-30 [1] CRAN (R 4.5.0)
    ##    loo              2.8.0      2024-07-03 [1] CRAN (R 4.5.1)
    ##    magrittr         2.0.3      2022-03-30 [1] CRAN (R 4.5.1)
    ##    MASS             7.3-65     2025-02-28 [1] CRAN (R 4.5.1)
    ##    mathjaxr         1.8-0      2025-04-30 [1] CRAN (R 4.5.1)
    ##    Matrix         * 1.7-3      2025-03-11 [1] CRAN (R 4.5.1)
    ##    MatrixModels     0.5-4      2025-03-26 [1] CRAN (R 4.5.1)
    ##    matrixStats      1.5.0      2025-01-07 [1] CRAN (R 4.5.1)
    ##    metadat        * 1.4-0      2025-02-04 [1] CRAN (R 4.5.1)
    ##    metafor        * 4.8-0      2025-01-28 [1] CRAN (R 4.5.1)
    ##    mgcv             1.9-3      2025-04-04 [1] CRAN (R 4.5.1)
    ##    multcomp         1.4-28     2025-01-29 [1] CRAN (R 4.5.1)
    ##    mvtnorm          1.3-3      2025-01-10 [1] CRAN (R 4.5.1)
    ##    nlme             3.1-168    2025-03-31 [1] CRAN (R 4.5.1)
    ##    nnet             7.3-20     2025-01-01 [1] CRAN (R 4.5.1)
    ##    numDeriv       * 2016.8-1.1 2019-06-06 [1] CRAN (R 4.5.0)
    ##    pillar           1.11.0     2025-07-04 [1] CRAN (R 4.5.1)
    ##    pkgbuild         1.4.8      2025-05-26 [1] CRAN (R 4.5.1)
    ##    pkgconfig        2.0.3      2019-09-22 [1] CRAN (R 4.5.1)
    ##    plyr             1.8.9      2023-10-02 [1] CRAN (R 4.5.1)
    ##    polspline        1.1.25     2024-05-10 [1] CRAN (R 4.5.0)
    ##    posterior        1.6.1      2025-02-27 [1] CRAN (R 4.5.1)
    ##    processx         3.8.6      2025-02-21 [1] CRAN (R 4.5.1)
    ##    proxy            0.4-27     2022-06-09 [1] CRAN (R 4.5.1)
    ##    ps               1.9.1      2025-04-12 [1] CRAN (R 4.5.1)
    ##    purrr            1.1.0      2025-07-10 [1] CRAN (R 4.5.1)
    ##    quantreg         6.1        2025-03-10 [1] CRAN (R 4.5.1)
    ##    QuickJSR         1.8.0      2025-06-09 [1] CRAN (R 4.5.1)
    ##    R6               2.6.1      2025-02-15 [1] CRAN (R 4.5.1)
    ##    RColorBrewer     1.1-3      2022-04-03 [1] CRAN (R 4.5.0)
    ##    Rcpp           * 1.1.0      2025-07-02 [1] CRAN (R 4.5.1)
    ##  D RcppParallel     5.1.10     2025-01-24 [1] CRAN (R 4.5.1)
    ##    readr            2.1.5      2024-01-10 [1] CRAN (R 4.5.1)
    ##    readxl         * 1.4.5      2025-03-07 [1] CRAN (R 4.5.1)
    ##    reshape2         1.4.4      2020-04-09 [1] CRAN (R 4.5.1)
    ##    rlang            1.1.6      2025-04-11 [1] CRAN (R 4.5.1)
    ##    rmarkdown      * 2.29       2024-11-04 [1] CRAN (R 4.5.1)
    ##    rms            * 8.0-0      2025-04-04 [1] CRAN (R 4.5.1)
    ##    rootSolve        1.8.2.4    2023-09-21 [1] CRAN (R 4.5.0)
    ##    rpart            4.1.24     2025-01-07 [1] CRAN (R 4.5.1)
    ##    rstan            2.32.7     2025-03-10 [1] CRAN (R 4.5.1)
    ##    rstantools       2.4.0      2024-01-31 [1] CRAN (R 4.5.1)
    ##    rstudioapi       0.17.1     2024-10-22 [1] CRAN (R 4.5.1)
    ##    sandwich         3.1-1      2024-09-15 [1] CRAN (R 4.5.1)
    ##    scales           1.4.0      2025-04-24 [1] CRAN (R 4.5.1)
    ##    sessioninfo      1.2.3      2025-02-05 [1] CRAN (R 4.5.1)
    ##    sf             * 1.0-21     2025-05-15 [1] CRAN (R 4.5.1)
    ##    SparseM          1.84-2     2024-07-17 [1] CRAN (R 4.5.1)
    ##    StanHeaders      2.32.10    2024-07-15 [1] CRAN (R 4.5.1)
    ##    stringi          1.8.7      2025-03-27 [1] CRAN (R 4.5.0)
    ##    stringr          1.5.1      2023-11-14 [1] CRAN (R 4.5.1)
    ##    survival         3.8-3      2024-12-17 [1] CRAN (R 4.5.1)
    ##    svglite          2.2.1      2025-05-12 [1] CRAN (R 4.5.1)
    ##    systemfonts      1.2.3      2025-04-30 [1] CRAN (R 4.5.1)
    ##    tensorA          0.36.2.1   2023-12-13 [1] CRAN (R 4.5.0)
    ##    textshaping      1.0.1      2025-05-01 [1] CRAN (R 4.5.1)
    ##    TH.data          1.1-3      2025-01-17 [1] CRAN (R 4.5.1)
    ##    tibble           3.3.0      2025-06-08 [1] CRAN (R 4.5.1)
    ##    tidyr          * 1.3.1      2024-01-24 [1] CRAN (R 4.5.1)
    ##    tidyselect       1.2.1      2024-03-11 [1] CRAN (R 4.5.1)
    ##    tzdb             0.5.0      2025-03-15 [1] CRAN (R 4.5.1)
    ##    units            0.8-7      2025-03-11 [1] CRAN (R 4.5.1)
    ##    V8               6.0.4      2025-06-04 [1] CRAN (R 4.5.1)
    ##    vctrs            0.6.5      2023-12-01 [1] CRAN (R 4.5.1)
    ##    viridisLite      0.4.2      2023-05-02 [1] CRAN (R 4.5.1)
    ##    withr            3.0.2      2024-10-28 [1] CRAN (R 4.5.1)
    ##    xfun             0.52       2025-04-02 [1] CRAN (R 4.5.1)
    ##    xml2             1.3.8      2025-03-14 [1] CRAN (R 4.5.1)
    ##    yaml             2.3.10     2024-07-26 [1] CRAN (R 4.5.0)
    ##    zoo              1.8-14     2025-04-10 [1] CRAN (R 4.5.1)
    ## 
    ##  [1] C:/Program Files/R/R-4.5.1/library
    ## 
    ##  * ── Packages attached to the search path.
    ##  D ── DLL MD5 mismatch, broken installation.
    ## 
    ## ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

``` r
##rmarkdown::render("03-estimate_v1.R")

# Save dataset for report created for expert to receive feedback
# save(all_cnt_rt, file="./00-Report_FB/all_cnt_rt.Rdata")
# save(all_glb_prop, file="./00-Report_FB/all_glb_prop.Rdata")
# save(all_reg_prop, file="./00-Report_FB/all_reg_prop.Rdata")
# save(all_reg_rt, file="./00-Report_FB/all_reg_rt.Rdata")
# save(all_sub_nr, file="./00-Report_FB/all_sub_nr.Rdata")
# save(all_sub_rt, file="./00-Report_FB/all_sub_rt.Rdata")
```
