Global incidence of CKD5 Age 85-125 • Estimate incidence
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
    ## sd(Intercept)     1.98      0.87     0.19     3.60 1.00     3901     3329
    ## 
    ## ~REG2:SUB2 (Number of levels: 7) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     1.05      0.78     0.04     2.84 1.00     5998     4840
    ## 
    ## ~REG2:SUB2:COUNTRY (Number of levels: 17) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     1.99      0.84     0.26     3.54 1.00     2832     3044
    ## 
    ## ~REG2:SUB2:COUNTRY:ID (Number of levels: 110) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     2.97      0.35     2.31     3.67 1.00     5051     7038
    ## 
    ## ~REG2:SUB2:COUNTRY:ID:DTP_ID (Number of levels: 200) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     1.23      0.20     0.90     1.69 1.00     4137     6107
    ## 
    ## Regression Coefficients:
    ##           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept   123.10     99.67   -70.25   321.64 1.00     6187     7499
    ## YEAR         -0.07      0.05    -0.17     0.03 1.00     6180     7453
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
    ##  $ VAL_MEAN     : num  11.8 11.4 11 10.7 10.3 ...
    ##  $ VAL_MEDIAN   : num  2.95 2.92 2.91 2.85 2.8 ...
    ##  $ VAL_LWR      : num  0.764 0.77 0.775 0.768 0.757 ...
    ##  $ VAL_UPR      : num  39.2 38.7 36.5 36.1 35.1 ...
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
    ##  $ VAL_MEAN     : num [1:22(1d)] 0.0402 0.0379 0.036 0.0343 0.0327 ...
    ##  $ VAL_MEDIAN   : num [1:22(1d)] 0.01006 0.00974 0.0095 0.00917 0.00886 ...
    ##  $ VAL_LWR      : num [1:22(1d)] 0.00261 0.00256 0.00253 0.00247 0.00239 ...
    ##  $ VAL_UPR      : num [1:22(1d)] 0.134 0.129 0.119 0.116 0.111 ...
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
    ##  $ VAL_MEAN     : num  0.00418 0.02028 0.00309 0.15935 8.89542 ...
    ##  $ VAL_MEDIAN   : num  0.00119 0.000413 0.000878 0.005561 0.201347 ...
    ##  $ VAL_LWR      : num  6.23e-05 1.62e-05 4.60e-05 2.27e-04 2.11e-03 ...
    ##  $ VAL_UPR      : num  0.0261 0.0385 0.0193 0.3662 34.5146 ...
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
    ##  $ VAL_MEAN     : num  0.000547 0.000318 0.000547 0.001449 0.292811 ...
    ##  $ VAL_MEDIAN   : num  1.55e-04 6.48e-06 1.55e-04 5.06e-05 6.63e-03 ...
    ##  $ VAL_LWR      : num  8.14e-06 2.54e-07 8.14e-06 2.07e-06 6.96e-05 ...
    ##  $ VAL_UPR      : num  0.003416 0.000603 0.003416 0.003331 1.136118 ...
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
    ##  $ VAL_MEAN     : num  0.000772 0.00203 0.001382 0.000632 0.019524 ...
    ##  $ VAL_MEDIAN   : num  0.00022 0.000577 0.000393 0.000134 0.000161 ...
    ##  $ VAL_LWR      : num  1.15e-05 3.02e-05 2.06e-05 6.11e-06 2.73e-06 ...
    ##  $ VAL_UPR      : num  0.00482 0.01268 0.00864 0.00336 0.03607 ...
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
    ##  $ VAL_MEAN     : num  5.47e-04 5.47e-04 5.47e-04 1.36e-05 1.21e-03 ...
    ##  $ VAL_MEDIAN   : num  1.55e-04 1.55e-04 1.55e-04 2.88e-06 9.95e-06 ...
    ##  $ VAL_LWR      : num  8.14e-06 8.14e-06 8.14e-06 1.31e-07 1.69e-07 ...
    ##  $ VAL_UPR      : num  3.42e-03 3.42e-03 3.42e-03 7.21e-05 2.23e-03 ...
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
    ##  $ VAL_MEAN     : num  5.63e-05 5.25e-05 5.10e-05 5.22e-05 5.21e-05 ...
    ##  $ VAL_MEDIAN   : num  1.60e-05 1.51e-05 1.48e-05 1.51e-05 1.52e-05 ...
    ##  $ VAL_LWR      : num  8.38e-07 7.90e-07 7.81e-07 8.01e-07 8.03e-07 ...
    ##  $ VAL_UPR      : num  0.000352 0.000328 0.000314 0.000322 0.000319 ...
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
    ##  $ VAL_MEAN     : num  0.000547 0.000507 0.00047 0.000438 0.000408 ...
    ##  $ VAL_MEDIAN   : num  0.000155 0.000146 0.000136 0.000127 0.00012 ...
    ##  $ VAL_LWR      : num  8.14e-06 7.63e-06 7.21e-06 6.72e-06 6.30e-06 ...
    ##  $ VAL_UPR      : num  0.00342 0.00317 0.0029 0.00271 0.0025 ...
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
    ##  $ VAL_MEAN     : num  0.0402 0.0379 0.036 0.0343 0.0327 ...
    ##  $ VAL_MEDIAN   : num  0.01006 0.00974 0.0095 0.00917 0.00886 ...
    ##  $ VAL_LWR      : num  0.00261 0.00256 0.00253 0.00247 0.00239 ...
    ##  $ VAL_UPR      : num  0.134 0.129 0.119 0.116 0.111 ...
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

| YEAR |  VAL_MEAN | VAL_MEDIAN |   VAL_LWR |  VAL_UPR |
|-----:|----------:|-----------:|----------:|---------:|
| 2010 |  9.249873 |   2.702478 | 0.6222929 | 33.50812 |
| 2020 | 10.064374 |   2.426620 | 0.2889789 | 41.72843 |

Global number of CKD5 Age 85-125 cases, 2010 vs 2020

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
Incidence per 100.000 of CKD5 Age 85-125 by WHO region in 2020
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
0.0001834
</td>
<td style="text-align:center;">
0.0000419
</td>
<td style="text-align:center;">
0.0000016
</td>
<td style="text-align:left;">
0.0011360
</td>
</tr>
<tr>
<td style="text-align:left;">
AMR
</td>
<td style="text-align:center;">
0.0001197
</td>
<td style="text-align:center;">
0.0000020
</td>
<td style="text-align:center;">
0.0000001
</td>
<td style="text-align:left;">
0.0001777
</td>
</tr>
<tr>
<td style="text-align:left;">
EMR
</td>
<td style="text-align:center;">
0.0001834
</td>
<td style="text-align:center;">
0.0000419
</td>
<td style="text-align:center;">
0.0000016
</td>
<td style="text-align:left;">
0.0011360
</td>
</tr>
<tr>
<td style="text-align:left;">
EUR
</td>
<td style="text-align:center;">
0.0005381
</td>
<td style="text-align:center;">
0.0000126
</td>
<td style="text-align:center;">
0.0000003
</td>
<td style="text-align:left;">
0.0013021
</td>
</tr>
<tr>
<td style="text-align:left;">
SEAR
</td>
<td style="text-align:center;">
0.0825341
</td>
<td style="text-align:center;">
0.0021190
</td>
<td style="text-align:center;">
0.0000196
</td>
<td style="text-align:left;">
0.3637722
</td>
</tr>
<tr>
<td style="text-align:left;">
WPR
</td>
<td style="text-align:center;">
0.0146710
</td>
<td style="text-align:center;">
0.0077908
</td>
<td style="text-align:center;">
0.0010117
</td>
<td style="text-align:left;">
0.0660877
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
Cases of CKD5 Age 85-125 by WHO region in 2020
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
0.0032637
</td>
<td style="text-align:center;">
0.0007465
</td>
<td style="text-align:center;">
0.0000282
</td>
<td style="text-align:left;">
0.0202183
</td>
</tr>
<tr>
<td style="text-align:left;">
AMR
</td>
<td style="text-align:center;">
0.0138629
</td>
<td style="text-align:center;">
0.0002294
</td>
<td style="text-align:center;">
0.0000066
</td>
<td style="text-align:left;">
0.0205877
</td>
</tr>
<tr>
<td style="text-align:left;">
EMR
</td>
<td style="text-align:center;">
0.0026904
</td>
<td style="text-align:center;">
0.0006154
</td>
<td style="text-align:center;">
0.0000232
</td>
<td style="text-align:left;">
0.0166669
</td>
</tr>
<tr>
<td style="text-align:left;">
EUR
</td>
<td style="text-align:center;">
0.1040323
</td>
<td style="text-align:center;">
0.0024300
</td>
<td style="text-align:center;">
0.0000610
</td>
<td style="text-align:left;">
0.2517646
</td>
</tr>
<tr>
<td style="text-align:left;">
SEAR
</td>
<td style="text-align:center;">
6.7465298
</td>
<td style="text-align:center;">
0.1732111
</td>
<td style="text-align:center;">
0.0016021
</td>
<td style="text-align:left;">
29.7355953
</td>
</tr>
<tr>
<td style="text-align:left;">
WPR
</td>
<td style="text-align:center;">
3.1939952
</td>
<td style="text-align:center;">
1.6961193
</td>
<td style="text-align:center;">
0.2202561
</td>
<td style="text-align:left;">
14.3878258
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
Incidence per 100.000 of CKD5 Age 85-125 by WHO subregion in 2020
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
0.0001834
</td>
<td style="text-align:center;">
0.0000419
</td>
<td style="text-align:center;">
0.0000016
</td>
<td style="text-align:left;">
0.0011360
</td>
</tr>
<tr>
<td style="text-align:left;">
AFRC
</td>
<td style="text-align:center;">
0.0001834
</td>
<td style="text-align:center;">
0.0000419
</td>
<td style="text-align:center;">
0.0000016
</td>
<td style="text-align:left;">
0.0011360
</td>
</tr>
<tr>
<td style="text-align:left;">
AFRD
</td>
<td style="text-align:center;">
0.0001834
</td>
<td style="text-align:center;">
0.0000419
</td>
<td style="text-align:center;">
0.0000016
</td>
<td style="text-align:left;">
0.0011360
</td>
</tr>
<tr>
<td style="text-align:left;">
AMRA
</td>
<td style="text-align:center;">
0.0000055
</td>
<td style="text-align:center;">
0.0000009
</td>
<td style="text-align:center;">
0.0000000
</td>
<td style="text-align:left;">
0.0000348
</td>
</tr>
<tr>
<td style="text-align:left;">
AMRB
</td>
<td style="text-align:center;">
0.0003359
</td>
<td style="text-align:center;">
0.0000025
</td>
<td style="text-align:center;">
0.0000000
</td>
<td style="text-align:left;">
0.0004590
</td>
</tr>
<tr>
<td style="text-align:left;">
AMRC
</td>
<td style="text-align:center;">
0.0000408
</td>
<td style="text-align:center;">
0.0000035
</td>
<td style="text-align:center;">
0.0000000
</td>
<td style="text-align:left;">
0.0002864
</td>
</tr>
<tr>
<td style="text-align:left;">
EMRA
</td>
<td style="text-align:center;">
0.0001834
</td>
<td style="text-align:center;">
0.0000419
</td>
<td style="text-align:center;">
0.0000016
</td>
<td style="text-align:left;">
0.0011360
</td>
</tr>
<tr>
<td style="text-align:left;">
EMRBC
</td>
<td style="text-align:center;">
0.0001834
</td>
<td style="text-align:center;">
0.0000419
</td>
<td style="text-align:center;">
0.0000016
</td>
<td style="text-align:left;">
0.0011360
</td>
</tr>
<tr>
<td style="text-align:left;">
EMRD
</td>
<td style="text-align:center;">
0.0001834
</td>
<td style="text-align:center;">
0.0000419
</td>
<td style="text-align:center;">
0.0000016
</td>
<td style="text-align:left;">
0.0011360
</td>
</tr>
<tr>
<td style="text-align:left;">
EURA
</td>
<td style="text-align:center;">
0.0006861
</td>
<td style="text-align:center;">
0.0000107
</td>
<td style="text-align:center;">
0.0000003
</td>
<td style="text-align:left;">
0.0016361
</td>
</tr>
<tr>
<td style="text-align:left;">
EURB
</td>
<td style="text-align:center;">
0.0000473
</td>
<td style="text-align:center;">
0.0000066
</td>
<td style="text-align:center;">
0.0000002
</td>
<td style="text-align:left;">
0.0003326
</td>
</tr>
<tr>
<td style="text-align:left;">
EURC
</td>
<td style="text-align:center;">
0.0000473
</td>
<td style="text-align:center;">
0.0000066
</td>
<td style="text-align:center;">
0.0000002
</td>
<td style="text-align:left;">
0.0003326
</td>
</tr>
<tr>
<td style="text-align:left;">
SEARB
</td>
<td style="text-align:center;">
0.1173329
</td>
<td style="text-align:center;">
0.0039492
</td>
<td style="text-align:center;">
0.0000316
</td>
<td style="text-align:left;">
0.6107408
</td>
</tr>
<tr>
<td style="text-align:left;">
SEARCD
</td>
<td style="text-align:center;">
0.0727364
</td>
<td style="text-align:center;">
0.0003904
</td>
<td style="text-align:center;">
0.0000010
</td>
<td style="text-align:left;">
0.2301078
</td>
</tr>
<tr>
<td style="text-align:left;">
WPRA
</td>
<td style="text-align:center;">
0.0253734
</td>
<td style="text-align:center;">
0.0085681
</td>
<td style="text-align:center;">
0.0005380
</td>
<td style="text-align:left;">
0.1402720
</td>
</tr>
<tr>
<td style="text-align:left;">
WPRB
</td>
<td style="text-align:center;">
0.0094354
</td>
<td style="text-align:center;">
0.0062371
</td>
<td style="text-align:center;">
0.0009533
</td>
<td style="text-align:left;">
0.0370606
</td>
</tr>
<tr>
<td style="text-align:left;">
WPRC
</td>
<td style="text-align:center;">
0.0032978
</td>
<td style="text-align:center;">
0.0007354
</td>
<td style="text-align:center;">
0.0000112
</td>
<td style="text-align:left;">
0.0231672
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
Cases of CKD5 Age 85-125 by WHO sub region in 2020
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
0.0007035
</td>
<td style="text-align:center;">
0.0001609
</td>
<td style="text-align:center;">
0.0000061
</td>
<td style="text-align:left;">
0.0043584
</td>
</tr>
<tr>
<td style="text-align:left;">
AFRC
</td>
<td style="text-align:center;">
0.0015885
</td>
<td style="text-align:center;">
0.0003634
</td>
<td style="text-align:center;">
0.0000137
</td>
<td style="text-align:left;">
0.0098408
</td>
</tr>
<tr>
<td style="text-align:left;">
AFRD
</td>
<td style="text-align:center;">
0.0009716
</td>
<td style="text-align:center;">
0.0002222
</td>
<td style="text-align:center;">
0.0000084
</td>
<td style="text-align:left;">
0.0060190
</td>
</tr>
<tr>
<td style="text-align:left;">
AMRA
</td>
<td style="text-align:center;">
0.0004053
</td>
<td style="text-align:center;">
0.0000641
</td>
<td style="text-align:center;">
0.0000021
</td>
<td style="text-align:left;">
0.0025616
</td>
</tr>
<tr>
<td style="text-align:left;">
AMRB
</td>
<td style="text-align:center;">
0.0133605
</td>
<td style="text-align:center;">
0.0000979
</td>
<td style="text-align:center;">
0.0000015
</td>
<td style="text-align:left;">
0.0182545
</td>
</tr>
<tr>
<td style="text-align:left;">
AMRC
</td>
<td style="text-align:center;">
0.0000972
</td>
<td style="text-align:center;">
0.0000083
</td>
<td style="text-align:center;">
0.0000001
</td>
<td style="text-align:left;">
0.0006818
</td>
</tr>
<tr>
<td style="text-align:left;">
EMRA
</td>
<td style="text-align:center;">
0.0001348
</td>
<td style="text-align:center;">
0.0000308
</td>
<td style="text-align:center;">
0.0000012
</td>
<td style="text-align:left;">
0.0008352
</td>
</tr>
<tr>
<td style="text-align:left;">
EMRBC
</td>
<td style="text-align:center;">
0.0021992
</td>
<td style="text-align:center;">
0.0005030
</td>
<td style="text-align:center;">
0.0000190
</td>
<td style="text-align:left;">
0.0136239
</td>
</tr>
<tr>
<td style="text-align:left;">
EMRD
</td>
<td style="text-align:center;">
0.0003564
</td>
<td style="text-align:center;">
0.0000815
</td>
<td style="text-align:center;">
0.0000031
</td>
<td style="text-align:left;">
0.0022078
</td>
</tr>
<tr>
<td style="text-align:left;">
EURA
</td>
<td style="text-align:center;">
0.1019131
</td>
<td style="text-align:center;">
0.0015872
</td>
<td style="text-align:center;">
0.0000383
</td>
<td style="text-align:left;">
0.2430253
</td>
</tr>
<tr>
<td style="text-align:left;">
EURB
</td>
<td style="text-align:center;">
0.0016807
</td>
<td style="text-align:center;">
0.0002348
</td>
<td style="text-align:center;">
0.0000066
</td>
<td style="text-align:left;">
0.0118194
</td>
</tr>
<tr>
<td style="text-align:left;">
EURC
</td>
<td style="text-align:center;">
0.0004385
</td>
<td style="text-align:center;">
0.0000613
</td>
<td style="text-align:center;">
0.0000017
</td>
<td style="text-align:left;">
0.0030839
</td>
</tr>
<tr>
<td style="text-align:left;">
SEARB
</td>
<td style="text-align:center;">
2.1071214
</td>
<td style="text-align:center;">
0.0709218
</td>
<td style="text-align:center;">
0.0005667
</td>
<td style="text-align:left;">
10.9679768
</td>
</tr>
<tr>
<td style="text-align:left;">
SEARCD
</td>
<td style="text-align:center;">
4.6394084
</td>
<td style="text-align:center;">
0.0249001
</td>
<td style="text-align:center;">
0.0000669
</td>
<td style="text-align:left;">
14.6771683
</td>
</tr>
<tr>
<td style="text-align:left;">
WPRA
</td>
<td style="text-align:center;">
1.9108291
</td>
<td style="text-align:center;">
0.6452527
</td>
<td style="text-align:center;">
0.0405168
</td>
<td style="text-align:left;">
10.5636533
</td>
</tr>
<tr>
<td style="text-align:left;">
WPRB
</td>
<td style="text-align:center;">
1.2506973
</td>
<td style="text-align:center;">
0.8267508
</td>
<td style="text-align:center;">
0.1263619
</td>
<td style="text-align:left;">
4.9125422
</td>
</tr>
<tr>
<td style="text-align:left;">
WPRC
</td>
<td style="text-align:center;">
0.0324688
</td>
<td style="text-align:center;">
0.0072404
</td>
<td style="text-align:center;">
0.0001098
</td>
<td style="text-align:left;">
0.2280963
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

    ## [1] 0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35

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

    ## [1] 0.00 0.05 0.10 0.15 0.20

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

| Country                          | 2010.mean | 2010.median | 2010.lwr | 2010.upr | 2020.mean | 2020.median | 2020.lwr | 2020.upr |
|:---------------------------------|----------:|------------:|---------:|---------:|----------:|------------:|---------:|---------:|
| Afghanistan                      |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Angola                           |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Albania                          |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Andorra                          |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| United Arab Emirates             |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Argentina                        |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Armenia                          |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Antigua and Barbuda              |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Australia                        |     0.004 |       0.001 |    0.000 |    0.026 |     0.002 |       0.000 |    0.000 |    0.016 |
| Austria                          |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Azerbaijan                       |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Burundi                          |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Belgium                          |     0.003 |       0.000 |    0.000 |    0.003 |     0.003 |       0.000 |    0.000 |    0.002 |
| Benin                            |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Burkina Faso                     |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Bangladesh                       |     0.044 |       0.001 |    0.000 |    0.179 |     0.025 |       0.000 |    0.000 |    0.098 |
| Bulgaria                         |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Bahrain                          |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Bahamas                          |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Bosnia and Herzegovina           |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Belarus                          |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Belize                           |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Bolivia                          |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Brazil                           |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Barbados                         |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Brunei Darussalam                |     0.007 |       0.003 |    0.000 |    0.040 |     0.005 |       0.001 |    0.000 |    0.030 |
| Bhutan                           |     0.044 |       0.001 |    0.000 |    0.179 |     0.025 |       0.000 |    0.000 |    0.098 |
| Botswana                         |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Central African Republic         |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Canada                           |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Switzerland                      |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Chile                            |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| China                            |     0.015 |       0.012 |    0.003 |    0.045 |     0.009 |       0.006 |    0.001 |    0.037 |
| Côte d’Ivoire                    |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Cameroon                         |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Congo                            |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Congo                            |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Cook Islands                     |     0.007 |       0.003 |    0.000 |    0.040 |     0.005 |       0.001 |    0.000 |    0.030 |
| Colombia                         |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Comoros                          |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Cabo Verde                       |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Costa Rica                       |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Cuba                             |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Cyprus                           |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Czechia                          |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Germany                          |     0.002 |       0.000 |    0.000 |    0.010 |     0.001 |       0.000 |    0.000 |    0.007 |
| Djibouti                         |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Dominica                         |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Denmark                          |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Dominican Republic               |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Algeria                          |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Ecuador                          |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Egypt                            |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Eritrea                          |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Spain                            |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Estonia                          |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Ethiopia                         |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Finland                          |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Fiji                             |     0.008 |       0.003 |    0.000 |    0.050 |     0.005 |       0.001 |    0.000 |    0.034 |
| France                           |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Micronesia                       |     0.005 |       0.001 |    0.000 |    0.032 |     0.003 |       0.001 |    0.000 |    0.023 |
| Gabon                            |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| United Kingdom                   |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Georgia                          |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Ghana                            |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Guinea                           |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Gambia                           |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Guinea-Bissau                    |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Equatorial Guinea                |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Greece                           |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Grenada                          |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Guatemala                        |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Guyana                           |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Honduras                         |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Croatia                          |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Haiti                            |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Hungary                          |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Indonesia                        |     0.075 |       0.001 |    0.000 |    0.285 |     0.039 |       0.001 |    0.000 |    0.161 |
| India                            |     0.138 |       0.001 |    0.000 |    0.431 |     0.082 |       0.000 |    0.000 |    0.250 |
| Ireland                          |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Iran                             |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Iraq                             |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Iceland                          |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Israel                           |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Italy                            |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Jamaica                          |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Jordan                           |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Japan                            |     0.031 |       0.016 |    0.002 |    0.146 |     0.024 |       0.008 |    0.001 |    0.143 |
| Kazakhstan                       |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Kenya                            |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Kyrgyzstan                       |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Cambodia                         |     0.005 |       0.001 |    0.000 |    0.032 |     0.003 |       0.001 |    0.000 |    0.023 |
| Kiribati                         |     0.005 |       0.001 |    0.000 |    0.032 |     0.003 |       0.001 |    0.000 |    0.023 |
| Saint Kitts and Nevis            |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Korea                            |     0.063 |       0.008 |    0.000 |    0.391 |     0.051 |       0.004 |    0.000 |    0.277 |
| Kuwait                           |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Lao People’s Dem. Republic       |     0.005 |       0.001 |    0.000 |    0.032 |     0.003 |       0.001 |    0.000 |    0.023 |
| Lebanon                          |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Liberia                          |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Libya                            |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Saint Lucia                      |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Sri Lanka                        |     0.044 |       0.001 |    0.000 |    0.179 |     0.025 |       0.000 |    0.000 |    0.098 |
| Lesotho                          |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Lithuania                        |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Luxembourg                       |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Latvia                           |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Morocco                          |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Monaco                           |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Republic of Moldova              |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Madagascar                       |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Maldives                         |     0.075 |       0.001 |    0.000 |    0.285 |     0.039 |       0.001 |    0.000 |    0.161 |
| Mexico                           |     0.002 |       0.000 |    0.000 |    0.003 |     0.002 |       0.000 |    0.000 |    0.002 |
| Marshall Islands                 |     0.008 |       0.003 |    0.000 |    0.050 |     0.005 |       0.001 |    0.000 |    0.034 |
| North Macedonia                  |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Mali                             |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Malta                            |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Myanmar                          |     0.044 |       0.001 |    0.000 |    0.179 |     0.025 |       0.000 |    0.000 |    0.098 |
| Montenegro                       |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Mongolia                         |     0.005 |       0.001 |    0.000 |    0.032 |     0.003 |       0.001 |    0.000 |    0.023 |
| Mozambique                       |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Mauritania                       |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Mauritius                        |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Malawi                           |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Malaysia                         |     0.008 |       0.003 |    0.000 |    0.050 |     0.005 |       0.001 |    0.000 |    0.034 |
| Namibia                          |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Niger                            |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Nigeria                          |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Nicaragua                        |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Niue                             |     0.007 |       0.003 |    0.000 |    0.040 |     0.005 |       0.001 |    0.000 |    0.030 |
| Netherlands                      |     0.005 |       0.000 |    0.000 |    0.003 |     0.003 |       0.000 |    0.000 |    0.002 |
| Norway                           |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Nepal                            |     0.044 |       0.001 |    0.000 |    0.179 |     0.025 |       0.000 |    0.000 |    0.098 |
| Nauru                            |     0.007 |       0.003 |    0.000 |    0.040 |     0.005 |       0.001 |    0.000 |    0.030 |
| New Zealand                      |     0.007 |       0.003 |    0.000 |    0.040 |     0.005 |       0.001 |    0.000 |    0.030 |
| Oman                             |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Pakistan                         |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Panama                           |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Peru                             |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Philippines                      |     0.005 |       0.001 |    0.000 |    0.032 |     0.003 |       0.001 |    0.000 |    0.023 |
| Palau                            |     0.008 |       0.003 |    0.000 |    0.050 |     0.005 |       0.001 |    0.000 |    0.034 |
| Papua New Guinea                 |     0.005 |       0.001 |    0.000 |    0.032 |     0.003 |       0.001 |    0.000 |    0.023 |
| Poland                           |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Korea                            |     0.044 |       0.001 |    0.000 |    0.179 |     0.025 |       0.000 |    0.000 |    0.098 |
| Portugal                         |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Paraguay                         |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Qatar                            |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Romania                          |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Russian Federation               |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Rwanda                           |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Saudi Arabia                     |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Sudan                            |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Senegal                          |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Singapore                        |     0.007 |       0.003 |    0.000 |    0.040 |     0.005 |       0.001 |    0.000 |    0.030 |
| Solomon Islands                  |     0.005 |       0.001 |    0.000 |    0.032 |     0.003 |       0.001 |    0.000 |    0.023 |
| Sierra Leone                     |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| El Salvador                      |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| San Marino                       |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Somalia                          |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Serbia                           |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| South Sudan                      |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Sao Tome and Principe            |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Suriname                         |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Slovakia                         |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Slovenia                         |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Sweden                           |     0.028 |       0.000 |    0.000 |    0.003 |     0.017 |       0.000 |    0.000 |    0.002 |
| Eswatini                         |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Seychelles                       |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Syrian Arab Republic             |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Chad                             |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Togo                             |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Thailand                         |     0.341 |       0.011 |    0.000 |    1.629 |     0.198 |       0.006 |    0.000 |    1.050 |
| Tajikistan                       |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Turkmenistan                     |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Timor-Leste                      |     0.044 |       0.001 |    0.000 |    0.179 |     0.025 |       0.000 |    0.000 |    0.098 |
| Tonga                            |     0.008 |       0.003 |    0.000 |    0.050 |     0.005 |       0.001 |    0.000 |    0.034 |
| Trinidad and Tobago              |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Tunisia                          |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Turkiye                          |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Tuvalu                           |     0.008 |       0.003 |    0.000 |    0.050 |     0.005 |       0.001 |    0.000 |    0.034 |
| United Republic of Tanzania      |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Uganda                           |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Ukraine                          |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Uruguay                          |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| United States of America         |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Uzbekistan                       |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Saint Vincent and the Grenadines |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Venezuela                        |     0.000 |       0.000 |    0.000 |    0.000 |     0.000 |       0.000 |    0.000 |    0.000 |
| Viet Nam                         |     0.005 |       0.001 |    0.000 |    0.032 |     0.003 |       0.001 |    0.000 |    0.023 |
| Vanuatu                          |     0.005 |       0.001 |    0.000 |    0.032 |     0.003 |       0.001 |    0.000 |    0.023 |
| Samoa                            |     0.005 |       0.001 |    0.000 |    0.032 |     0.003 |       0.001 |    0.000 |    0.023 |
| Yemen                            |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| South Africa                     |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Zambia                           |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |
| Zimbabwe                         |     0.000 |       0.000 |    0.000 |    0.002 |     0.000 |       0.000 |    0.000 |    0.001 |

Estimated CKD5 Age 85-125 incidence by country, 2010 vs 2020

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

| Country                          | 2010.mean | 2010.median | 2010.lwr | 2010.upr | 2020.mean | 2020.median | 2020.lwr | 2020.upr |
|:---------------------------------|----------:|------------:|---------:|---------:|----------:|------------:|---------:|---------:|
| Afghanistan                      |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Angola                           |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Albania                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Andorra                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| United Arab Emirates             |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Argentina                        |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Armenia                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Antigua and Barbuda              |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Australia                        |       0.0 |         0.0 |      0.0 |      0.1 |       0.0 |         0.0 |      0.0 |      0.1 |
| Austria                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Azerbaijan                       |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Burundi                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Belgium                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Benin                            |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Burkina Faso                     |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Bangladesh                       |       0.1 |         0.0 |      0.0 |      0.3 |       0.1 |         0.0 |      0.0 |      0.5 |
| Bulgaria                         |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Bahrain                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Bahamas                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Bosnia and Herzegovina           |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Belarus                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Belize                           |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Bolivia                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Brazil                           |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Barbados                         |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Brunei Darussalam                |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Bhutan                           |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Botswana                         |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Central African Republic         |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Canada                           |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Switzerland                      |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Chile                            |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| China                            |       1.2 |         1.0 |      0.2 |      3.6 |       1.2 |         0.8 |      0.1 |      4.9 |
| Côte d’Ivoire                    |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Cameroon                         |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Congo                            |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Congo                            |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Cook Islands                     |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Colombia                         |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Comoros                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Cabo Verde                       |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Costa Rica                       |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Cuba                             |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Cyprus                           |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Czechia                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Germany                          |       0.0 |         0.0 |      0.0 |      0.2 |       0.0 |         0.0 |      0.0 |      0.2 |
| Djibouti                         |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Dominica                         |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Denmark                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Dominican Republic               |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Algeria                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Ecuador                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Egypt                            |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Eritrea                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Spain                            |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Estonia                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Ethiopia                         |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Finland                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Fiji                             |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| France                           |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Micronesia                       |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Gabon                            |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| United Kingdom                   |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Georgia                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Ghana                            |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Guinea                           |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Gambia                           |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Guinea-Bissau                    |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Equatorial Guinea                |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Greece                           |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Grenada                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Guatemala                        |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Guyana                           |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Honduras                         |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Croatia                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Haiti                            |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Hungary                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Indonesia                        |       0.4 |         0.0 |      0.0 |      1.5 |       0.4 |         0.0 |      0.0 |      1.5 |
| India                            |       4.4 |         0.0 |      0.0 |     13.7 |       4.4 |         0.0 |      0.0 |     13.4 |
| Ireland                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Iran                             |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Iraq                             |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Iceland                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Israel                           |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Italy                            |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Jamaica                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Jordan                           |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Japan                            |       1.1 |         0.6 |      0.1 |      5.5 |       1.5 |         0.5 |      0.0 |      8.6 |
| Kazakhstan                       |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Kenya                            |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Kyrgyzstan                       |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Cambodia                         |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Kiribati                         |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Saint Kitts and Nevis            |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Korea                            |       0.2 |         0.0 |      0.0 |      1.4 |       0.4 |         0.0 |      0.0 |      2.3 |
| Kuwait                           |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Lao People’s Dem. Republic       |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Lebanon                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Liberia                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Libya                            |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Saint Lucia                      |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Sri Lanka                        |       0.0 |         0.0 |      0.0 |      0.2 |       0.0 |         0.0 |      0.0 |      0.1 |
| Lesotho                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Lithuania                        |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Luxembourg                       |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Latvia                           |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Morocco                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Monaco                           |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Republic of Moldova              |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Madagascar                       |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Maldives                         |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Mexico                           |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Marshall Islands                 |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| North Macedonia                  |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Mali                             |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Malta                            |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Myanmar                          |       0.0 |         0.0 |      0.0 |      0.2 |       0.0 |         0.0 |      0.0 |      0.1 |
| Montenegro                       |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Mongolia                         |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Mozambique                       |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Mauritania                       |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Mauritius                        |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Malawi                           |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Malaysia                         |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Namibia                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Niger                            |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Nigeria                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Nicaragua                        |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Niue                             |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Netherlands                      |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Norway                           |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Nepal                            |       0.0 |         0.0 |      0.0 |      0.1 |       0.0 |         0.0 |      0.0 |      0.1 |
| Nauru                            |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| New Zealand                      |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Oman                             |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Pakistan                         |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Panama                           |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Peru                             |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Philippines                      |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Palau                            |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Papua New Guinea                 |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Poland                           |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Korea                            |       0.0 |         0.0 |      0.0 |      0.1 |       0.0 |         0.0 |      0.0 |      0.2 |
| Portugal                         |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Paraguay                         |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Qatar                            |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Romania                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Russian Federation               |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Rwanda                           |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Saudi Arabia                     |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Sudan                            |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Senegal                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Singapore                        |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Solomon Islands                  |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Sierra Leone                     |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| El Salvador                      |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| San Marino                       |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Somalia                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Serbia                           |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| South Sudan                      |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Sao Tome and Principe            |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Suriname                         |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Slovakia                         |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Slovenia                         |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Sweden                           |       0.1 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Eswatini                         |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Seychelles                       |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Syrian Arab Republic             |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Chad                             |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Togo                             |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Thailand                         |       1.5 |         0.0 |      0.0 |      7.0 |       1.8 |         0.1 |      0.0 |      9.3 |
| Tajikistan                       |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Turkmenistan                     |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Timor-Leste                      |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Tonga                            |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Trinidad and Tobago              |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Tunisia                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Turkiye                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Tuvalu                           |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| United Republic of Tanzania      |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Uganda                           |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Ukraine                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Uruguay                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| United States of America         |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Uzbekistan                       |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Saint Vincent and the Grenadines |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Venezuela                        |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Viet Nam                         |       0.0 |         0.0 |      0.0 |      0.2 |       0.0 |         0.0 |      0.0 |      0.2 |
| Vanuatu                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Samoa                            |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Yemen                            |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| South Africa                     |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Zambia                           |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Zimbabwe                         |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |

Estimated CKD5 Age 85-125 cases by country, 2010 vs 2020

# Session info

``` r
sessioninfo::session_info()
```

    ## Warning in system2("quarto", "-V", stdout = TRUE, env = paste0("TMPDIR=", : running command '"quarto"
    ## TMPDIR=C:/Users/LoVa3397/AppData/Local/Temp/RtmpsVe1vM/file3e6c70f6c22 -V' had status 1

    ## ─ Session info ────────────────────────────────────────────────────────────────────────────────────────────────────────
    ##  setting  value
    ##  version  R version 4.5.0 (2025-04-11 ucrt)
    ##  os       Windows 10 x64 (build 19045)
    ##  system   x86_64, mingw32
    ##  ui       RStudio
    ##  language (EN)
    ##  collate  English_Belgium.utf8
    ##  ctype    English_Belgium.utf8
    ##  tz       Europe/Brussels
    ##  date     2025-10-16
    ##  rstudio  2024.04.2+764 Chocolate Cosmos (desktop)
    ##  pandoc   3.1.11 @ C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools/ (via rmarkdown)
    ##  quarto   ERROR: Unknown command "TMPDIR=C:/Users/LoVa3397/AppData/Local/Temp/RtmpsVe1vM/file3e6c70f6c22". Did you mean command "preview"? @ C:\\PROGRA~1\\RStudio\\RESOUR~1\\app\\bin\\quarto\\bin\\quarto.exe
    ## 
    ## ─ Packages ────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ##  ! package        * version    date (UTC) lib source
    ##    abind            1.4-8      2024-09-12 [1] CRAN (R 4.5.0)
    ##    backports        1.5.0      2024-05-23 [1] CRAN (R 4.5.0)
    ##    base64enc        0.1-3      2015-07-28 [1] CRAN (R 4.5.0)
    ##    bayesplot        1.12.0     2025-04-10 [1] CRAN (R 4.5.0)
    ##    bd             * 0.0.14     2025-04-26 [1] Github (brechtdv/bd@652191c)
    ##    boot             1.3-31     2024-08-28 [1] CRAN (R 4.5.0)
    ##    bridgesampling   1.1-2      2021-04-16 [1] CRAN (R 4.5.0)
    ##    brms           * 2.22.0     2024-09-23 [1] CRAN (R 4.5.0)
    ##    Brobdingnag      1.2-9      2022-10-19 [1] CRAN (R 4.5.0)
    ##    callr            3.7.6      2024-03-25 [1] CRAN (R 4.5.0)
    ##    cellranger       1.1.0      2016-07-27 [1] CRAN (R 4.5.0)
    ##    checkmate        2.3.2      2024-07-29 [1] CRAN (R 4.5.0)
    ##    class            7.3-23     2025-01-01 [1] CRAN (R 4.5.0)
    ##    classInt         0.4-11     2025-01-08 [1] CRAN (R 4.5.0)
    ##    cli              3.6.4      2025-02-13 [1] CRAN (R 4.5.0)
    ##    cluster          2.1.8.1    2025-03-12 [1] CRAN (R 4.5.0)
    ##    coda             0.19-4.1   2024-01-31 [1] CRAN (R 4.5.0)
    ##    codetools        0.2-20     2024-03-31 [1] CRAN (R 4.5.0)
    ##    colorspace       2.1-1      2024-07-26 [1] CRAN (R 4.5.0)
    ##    cowplot        * 1.1.3      2024-01-22 [1] CRAN (R 4.5.0)
    ##    curl             6.2.2      2025-03-24 [1] CRAN (R 4.5.0)
    ##    data.table       1.17.0     2025-02-22 [1] CRAN (R 4.5.0)
    ##    DBI              1.2.3      2024-06-02 [1] CRAN (R 4.5.0)
    ##    DescTools      * 0.99.60    2025-03-28 [1] CRAN (R 4.5.0)
    ##    digest           0.6.37     2024-08-19 [1] CRAN (R 4.5.0)
    ##    distributional   0.5.0      2024-09-17 [1] CRAN (R 4.5.0)
    ##    dplyr          * 1.1.4      2023-11-17 [1] CRAN (R 4.5.0)
    ##    e1071            1.7-16     2024-09-16 [1] CRAN (R 4.5.0)
    ##    evaluate         1.0.3      2025-01-10 [1] CRAN (R 4.5.0)
    ##    Exact            3.3        2024-07-21 [1] CRAN (R 4.5.0)
    ##    expm             1.0-0      2024-08-19 [1] CRAN (R 4.5.0)
    ##    farver           2.1.2      2024-05-13 [1] CRAN (R 4.5.0)
    ##    fastmap          1.2.0      2024-05-15 [1] CRAN (R 4.5.0)
    ##    FERG2          * 0.0.5      2025-07-28 [1] Github (brechtdv/FERG2@c2d4ac1)
    ##    forcats          1.0.0      2023-01-29 [1] CRAN (R 4.5.0)
    ##    foreign          0.8-90     2025-03-31 [1] CRAN (R 4.5.0)
    ##    Formula          1.2-5      2023-02-24 [1] CRAN (R 4.5.0)
    ##    fs               1.6.6      2025-04-12 [1] CRAN (R 4.5.0)
    ##    generics         0.1.3      2022-07-05 [1] CRAN (R 4.5.0)
    ##    ggplot2        * 3.5.2      2025-04-09 [1] CRAN (R 4.5.0)
    ##    gld              2.6.7      2025-01-17 [1] CRAN (R 4.5.0)
    ##    glue             1.8.0      2024-09-30 [1] CRAN (R 4.5.0)
    ##    gridExtra        2.3        2017-09-09 [1] CRAN (R 4.5.0)
    ##    gtable           0.3.6      2024-10-25 [1] CRAN (R 4.5.0)
    ##    haven            2.5.4      2023-11-30 [1] CRAN (R 4.5.0)
    ##    Hmisc          * 5.2-3      2025-03-16 [1] CRAN (R 4.5.0)
    ##    hms              1.1.3      2023-03-21 [1] CRAN (R 4.5.0)
    ##    htmlTable        2.4.3      2024-07-21 [1] CRAN (R 4.5.0)
    ##    htmltools        0.5.8.1    2024-04-04 [1] CRAN (R 4.5.0)
    ##    htmlwidgets      1.6.4      2023-12-06 [1] CRAN (R 4.5.0)
    ##    httr             1.4.7      2023-08-15 [1] CRAN (R 4.5.0)
    ##    inline           0.3.21     2025-01-09 [1] CRAN (R 4.5.0)
    ##    jsonlite         2.0.0      2025-03-27 [1] CRAN (R 4.5.0)
    ##    kableExtra     * 1.4.0      2024-01-24 [1] CRAN (R 4.5.0)
    ##    KernSmooth       2.23-26    2025-01-01 [1] CRAN (R 4.5.0)
    ##    knitr          * 1.50       2025-03-16 [1] CRAN (R 4.5.0)
    ##    labeling         0.4.3      2023-08-29 [1] CRAN (R 4.5.0)
    ##    lattice          0.22-6     2024-03-20 [1] CRAN (R 4.5.0)
    ##    lifecycle        1.0.4      2023-11-07 [1] CRAN (R 4.5.0)
    ##    lmom             3.2        2024-09-30 [1] CRAN (R 4.5.0)
    ##    loo              2.8.0      2024-07-03 [1] CRAN (R 4.5.0)
    ##    magrittr         2.0.3      2022-03-30 [1] CRAN (R 4.5.0)
    ##    MASS             7.3-65     2025-02-28 [1] CRAN (R 4.5.0)
    ##    mathjaxr         1.6-0      2022-02-28 [1] CRAN (R 4.5.0)
    ##    Matrix         * 1.7-3      2025-03-11 [1] CRAN (R 4.5.0)
    ##    MatrixModels     0.5-4      2025-03-26 [1] CRAN (R 4.5.0)
    ##    matrixStats      1.5.0      2025-01-07 [1] CRAN (R 4.5.0)
    ##    metadat        * 1.4-0      2025-02-04 [1] CRAN (R 4.5.0)
    ##    metafor        * 4.8-0      2025-01-28 [1] CRAN (R 4.5.0)
    ##    mgcv             1.9-1      2023-12-21 [1] CRAN (R 4.5.0)
    ##    multcomp         1.4-28     2025-01-29 [1] CRAN (R 4.5.0)
    ##    munsell          0.5.1      2024-04-01 [1] CRAN (R 4.5.0)
    ##    mvtnorm          1.3-3      2025-01-10 [1] CRAN (R 4.5.0)
    ##    nlme             3.1-168    2025-03-31 [1] CRAN (R 4.5.0)
    ##    nnet             7.3-20     2025-01-01 [1] CRAN (R 4.5.0)
    ##    numDeriv       * 2016.8-1.1 2019-06-06 [1] CRAN (R 4.5.0)
    ##    pillar           1.10.2     2025-04-05 [1] CRAN (R 4.5.0)
    ##    pkgbuild         1.4.7      2025-03-24 [1] CRAN (R 4.5.0)
    ##    pkgconfig        2.0.3      2019-09-22 [1] CRAN (R 4.5.0)
    ##    plyr             1.8.9      2023-10-02 [1] CRAN (R 4.5.0)
    ##    polspline        1.1.25     2024-05-10 [1] CRAN (R 4.5.0)
    ##    posterior        1.6.1      2025-02-27 [1] CRAN (R 4.5.0)
    ##    processx         3.8.6      2025-02-21 [1] CRAN (R 4.5.0)
    ##    proxy            0.4-27     2022-06-09 [1] CRAN (R 4.5.0)
    ##    ps               1.9.1      2025-04-12 [1] CRAN (R 4.5.0)
    ##    purrr            1.1.0      2025-07-10 [1] CRAN (R 4.5.1)
    ##    quantreg         6.1        2025-03-10 [1] CRAN (R 4.5.0)
    ##    QuickJSR         1.7.0      2025-03-31 [1] CRAN (R 4.5.0)
    ##    R6               2.6.1      2025-02-15 [1] CRAN (R 4.5.0)
    ##    RColorBrewer     1.1-3      2022-04-03 [1] CRAN (R 4.5.0)
    ##    Rcpp           * 1.0.14     2025-01-12 [1] CRAN (R 4.5.0)
    ##  D RcppParallel     5.1.10     2025-01-24 [1] CRAN (R 4.5.0)
    ##    readr            2.1.5      2024-01-10 [1] CRAN (R 4.5.0)
    ##    readxl         * 1.4.5      2025-03-07 [1] CRAN (R 4.5.0)
    ##    reshape2         1.4.4      2020-04-09 [1] CRAN (R 4.5.0)
    ##    rlang            1.1.6      2025-04-11 [1] CRAN (R 4.5.0)
    ##    rmarkdown      * 2.29       2024-11-04 [1] CRAN (R 4.5.0)
    ##    rms            * 8.0-0      2025-04-04 [1] CRAN (R 4.5.0)
    ##    rootSolve        1.8.2.4    2023-09-21 [1] CRAN (R 4.5.0)
    ##    rpart            4.1.24     2025-01-07 [1] CRAN (R 4.5.0)
    ##    rstan            2.32.7     2025-03-10 [1] CRAN (R 4.5.0)
    ##    rstantools       2.4.0      2024-01-31 [1] CRAN (R 4.5.0)
    ##    rstudioapi       0.17.1     2024-10-22 [1] CRAN (R 4.5.0)
    ##    sandwich         3.1-1      2024-09-15 [1] CRAN (R 4.5.0)
    ##    scales           1.3.0      2023-11-28 [1] CRAN (R 4.5.0)
    ##    sessioninfo      1.2.3      2025-02-05 [1] CRAN (R 4.5.0)
    ##    sf             * 1.0-20     2025-03-24 [1] CRAN (R 4.5.0)
    ##    SparseM          1.84-2     2024-07-17 [1] CRAN (R 4.5.0)
    ##    StanHeaders      2.32.10    2024-07-15 [1] CRAN (R 4.5.0)
    ##    stringi          1.8.7      2025-03-27 [1] CRAN (R 4.5.0)
    ##    stringr          1.5.1      2023-11-14 [1] CRAN (R 4.5.0)
    ##    survival         3.8-3      2024-12-17 [1] CRAN (R 4.5.0)
    ##    svglite          2.1.3      2023-12-08 [1] CRAN (R 4.5.0)
    ##    systemfonts      1.2.2      2025-04-04 [1] CRAN (R 4.5.0)
    ##    tensorA          0.36.2.1   2023-12-13 [1] CRAN (R 4.5.0)
    ##    TH.data          1.1-3      2025-01-17 [1] CRAN (R 4.5.0)
    ##    tibble           3.2.1      2023-03-20 [1] CRAN (R 4.5.0)
    ##    tidyr          * 1.3.1      2024-01-24 [1] CRAN (R 4.5.0)
    ##    tidyselect       1.2.1      2024-03-11 [1] CRAN (R 4.5.0)
    ##    tzdb             0.5.0      2025-03-15 [1] CRAN (R 4.5.0)
    ##    units            0.8-7      2025-03-11 [1] CRAN (R 4.5.0)
    ##    V8               6.0.3      2025-03-26 [1] CRAN (R 4.5.0)
    ##    vctrs            0.6.5      2023-12-01 [1] CRAN (R 4.5.0)
    ##    viridisLite      0.4.2      2023-05-02 [1] CRAN (R 4.5.0)
    ##    withr            3.0.2      2024-10-28 [1] CRAN (R 4.5.0)
    ##    xfun             0.52       2025-04-02 [1] CRAN (R 4.5.0)
    ##    xml2             1.3.8      2025-03-14 [1] CRAN (R 4.5.0)
    ##    yaml             2.3.10     2024-07-26 [1] CRAN (R 4.5.0)
    ##    zoo              1.8-14     2025-04-10 [1] CRAN (R 4.5.0)
    ## 
    ##  [1] C:/Users/LoVa3397/AppData/Local/Programs/R/R-4.5.0/library
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
