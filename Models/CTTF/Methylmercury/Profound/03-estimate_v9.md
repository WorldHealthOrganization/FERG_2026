Global incidence of Methyl mercury Profound • Estimate incidence
================
LoVa3397
2025-09-15

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

    ## Warning: There were 4 divergent transitions after warmup. Increasing adapt_delta above 0.8 may help. See
    ## http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup

    ##  Family: gaussian 
    ##   Links: mu = identity; sigma = identity 
    ## Formula: yi | se(sei) ~ 1 + (1 | REG2) + (1 | REG2:SUB2) + (1 | REG2:SUB2:COUNTRY) + (1 | REG2:SUB2:COUNTRY:ID) + (1 | REG2:SUB2:COUNTRY:ID:DTP_ID) 
    ##    Data: es (Number of observations: 1450) 
    ##   Draws: 5 chains, each with iter = 5000; warmup = 3000; thin = 1;
    ##          total post-warmup draws = 10000
    ## 
    ## Multilevel Hyperparameters:
    ## ~REG2 (Number of levels: 6) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.65      0.49     0.03     1.86 1.00     6272     5445
    ## 
    ## ~REG2:SUB2 (Number of levels: 17) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.78      0.51     0.03     1.90 1.00     4646     4385
    ## 
    ## ~REG2:SUB2:COUNTRY (Number of levels: 90) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     1.00      0.56     0.06     2.12 1.00     3940     4040
    ## 
    ## ~REG2:SUB2:COUNTRY:ID (Number of levels: 1034) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.32      0.24     0.01     0.89 1.00     6271     4290
    ## 
    ## ~REG2:SUB2:COUNTRY:ID:DTP_ID (Number of levels: 1450) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.29      0.22     0.01     0.81 1.00     7065     4745
    ## 
    ## Regression Coefficients:
    ##           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept    -6.81      0.76    -8.43    -5.39 1.00     7023     4766
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
  select(REG2, SUB2, ISO3, Country, cttf_mercury) %>% 
  rename(COUNTRY=ISO3, COUNTRY_LABEL = Country) %>%
  mutate(DISEASEFREE = cttf_mercury)

kable(
  caption = "Disease-free countries",
  row.names = FALSE,
  subset(zero_cases, cttf_mercury==0)[, 4])
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
    COUNTRY = FERG2:::countries$ISO3)
sim_all <- sim_all %>% left_join(zero_cases) %>% select(sei, REG2, SUB2, COUNTRY, ESTIMATES)
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
      # sim_all[x, "YEAR"] * draws_fit$b_YEAR +                                                                                         # Year component
      draws_fit[[paste0("r_REG2[",sim_all[x,"REG2"],",Intercept]")]] +                                                                # Regional component
      draws_fit[[paste0("r_REG2:SUB2[",sim_all[x,"REG2"],"_",sim_all[x,"SUB2"],",Intercept]")]] +                                     # Sub regional component
      draws_fit[[paste0("r_REG2:SUB2:COUNTRY[",sim_all[x,"REG2"],"_",sim_all[x,"SUB2"],"_",sim_all[x,"COUNTRY"],",Intercept]")]]      # Country component
  } else if (as.integer(sim_all[x, "ESTIMATES"]) == 2) {
    # Disease-free country
    fit_all[[paste0("V",x)]] <- 0
  } else if (as.integer(sim_all[x, "ESTIMATES"]) == 3){
    # Data not present for country, but present in subregion
    fit_all[[paste0("V",x)]] <- draws_fit$b_Intercept +                                                                               # Global intercept
      # sim_all[x, "YEAR"] * draws_fit$b_YEAR +                                                                                         # Year component
      draws_fit[[paste0("r_REG2[",sim_all[x,"REG2"],",Intercept]")]] +                                                                # Regional component
      draws_fit[[paste0("r_REG2:SUB2[",sim_all[x,"REG2"],"_",sim_all[x,"SUB2"],",Intercept]")]]                                       # Sub regional component
  } else if (as.integer(sim_all[x, "ESTIMATES"]) == 4){
    # Data not present for country, but present in region
    fit_all[[paste0("V",x)]] <- draws_fit$b_Intercept +                                                                               # Global intercept
      # sim_all[x, "YEAR"] * draws_fit$b_YEAR +                                                                                         # Year component
      draws_fit[[paste0("r_REG2[",sim_all[x,"REG2"],",Intercept]")]]                                                                  # Regional component
  } else if (as.integer(sim_all[x, "ESTIMATES"]) == 5){
    # Data not present for country
    fit_all[[paste0("V",x)]] <- draws_fit$b_Intercept  
      # sim_all[x, "YEAR"] * draws_fit$b_YEAR
  } 
}

fit_all <- fit_all %>% select(-c(X1.10000))


## calculate cases
sim_all$SIM <- t(fit_all)
sim_all <- do.call("rbind", replicate(2021-2000+1, sim_all, simplify = FALSE))
sim_all$YEAR <- sort(rep(2000:2021,194))
pop <- FERG2:::life_birth
pop$POP <- pop$LB
pop$YEAR <- as.integer(pop$YEAR)
pop_all <- aggregate(POP ~ ISO3 + YEAR, pop, sum)
sim_all <- merge(sim_all, pop_all,
                 by.x = c("COUNTRY", "YEAR"), by.y = c("ISO3", "YEAR"))
sim_all <- sim_all %>% left_join(zero_cases)
```

    ## Joining with `by = join_by(COUNTRY, REG2, SUB2, ESTIMATES)`

``` r
sim_all$CASES <- exp(sim_all$SIM) * sim_all$POP / 1e3
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
    ##  $ VAL_MEAN     : num  447 448 451 455 460 ...
    ##  $ VAL_MEDIAN   : num  263 262 262 263 265 ...
    ##  $ VAL_LWR      : num  91.7 91.2 91 91 91.3 ...
    ##  $ VAL_UPR      : num  1809 1809 1831 1861 1875 ...
    ##  $ LOCATION     : chr  "Global" "Global" "Global" "Global" ...
    ##  $ LOCATION_NAME: chr  "Global" "Global" "Global" "Global" ...
    ##  $ METRIC       : chr  "Number" "Number" "Number" "Number" ...

``` r
all_glb_rt <- all_glb_nr
all_glb_rt$POP <- with(sim_all, tapply(POP, YEAR, sum))
all_glb_rt$VAL_MEAN <- 1e3 * all_glb_rt$VAL_MEAN / all_glb_rt$POP
all_glb_rt$VAL_MEDIAN <- 1e3 * all_glb_rt$VAL_MEDIAN / all_glb_rt$POP
all_glb_rt$VAL_LWR <- 1e3 * all_glb_rt$VAL_LWR / all_glb_rt$POP
all_glb_rt$VAL_UPR <- 1e3 * all_glb_rt$VAL_UPR / all_glb_rt$POP
all_glb_rt$METRIC <- "Rate"
all_glb_rt$POP <- NULL
str(all_glb_rt)
```

    ## 'data.frame':    22 obs. of  8 variables:
    ##  $ YEAR         : int  2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 ...
    ##  $ VAL_MEAN     : num [1:22(1d)] 0.00333 0.00335 0.00337 0.00339 0.0034 ...
    ##  $ VAL_MEDIAN   : num [1:22(1d)] 0.00196 0.00196 0.00196 0.00196 0.00196 ...
    ##  $ VAL_LWR      : num [1:22(1d)] 0.000684 0.000682 0.00068 0.000678 0.000676 ...
    ##  $ VAL_UPR      : num [1:22(1d)] 0.0135 0.0135 0.0137 0.0139 0.0139 ...
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
    ##  $ VAL_MEAN     : num  180.8 55.2 47.3 42.1 52.2 ...
    ##  $ VAL_MEDIAN   : num  49.3 38.7 23.8 18.4 32 ...
    ##  $ VAL_LWR      : num  2.86 14.67 3.13 1.75 4.34 ...
    ##  $ VAL_UPR      : num  930 191 228 220 213 ...
    ##  $ LOCATION     : chr  "Region" "Region" "Region" "Region" ...
    ##  $ LOCATION_NAME: chr  "AFR" "AMR" "EMR" "EUR" ...
    ##  $ METRIC       : chr  "Number" "Number" "Number" "Number" ...

``` r
all_reg_rt <- all_reg_nr
all_reg_rt$POP <-
  with(sim_all, aggregate(POP ~ REG2 + YEAR, FUN = sum))$POP
all_reg_rt$VAL_MEAN <- 1e3 * all_reg_rt$VAL_MEAN / all_reg_rt$POP
all_reg_rt$VAL_MEDIAN <- 1e3 * all_reg_rt$VAL_MEDIAN / all_reg_rt$POP
all_reg_rt$VAL_LWR <- 1e3 * all_reg_rt$VAL_LWR / all_reg_rt$POP
all_reg_rt$VAL_UPR <- 1e3 * all_reg_rt$VAL_UPR / all_reg_rt$POP
all_reg_rt$METRIC <- "Rate"
all_reg_rt$POP <- NULL
str(all_reg_rt)
```

    ## 'data.frame':    132 obs. of  8 variables:
    ##  $ YEAR         : int  2000 2000 2000 2000 2000 2000 2001 2001 2001 2001 ...
    ##  $ VAL_MEAN     : num  0.00665 0.00352 0.0031 0.00412 0.00128 ...
    ##  $ VAL_MEDIAN   : num  0.001812 0.002465 0.001561 0.0018 0.000782 ...
    ##  $ VAL_LWR      : num  0.000105 0.000934 0.000205 0.000172 0.000106 ...
    ##  $ VAL_UPR      : num  0.03421 0.01213 0.01494 0.02156 0.00521 ...
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
    ##  $ VAL_MEAN     : num  8.21 103.39 69.18 8.59 39.82 ...
    ##  $ VAL_MEDIAN   : num  1.45 20.2 17.6 6.12 26.25 ...
    ##  $ VAL_LWR      : num  0.0426 0.7358 0.6092 1.1679 9.5779 ...
    ##  $ VAL_UPR      : num  47.8 464.9 371.2 31.2 143.1 ...
    ##  $ LOCATION     : chr  "Subregion" "Subregion" "Subregion" "Subregion" ...
    ##  $ LOCATION_NAME: chr  "AFRAB" "AFRC" "AFRD" "AMRA" ...
    ##  $ METRIC       : chr  "Number" "Number" "Number" "Number" ...

``` r
all_sub_rt <- all_sub_nr
all_sub_rt$POP <-
  with(sim_all, aggregate(POP ~ SUB2 + YEAR, FUN = sum))$POP
all_sub_rt$VAL_MEAN <- 1e3 * all_sub_rt$VAL_MEAN / all_sub_rt$POP
all_sub_rt$VAL_MEDIAN <- 1e3 * all_sub_rt$VAL_MEDIAN / all_sub_rt$POP
all_sub_rt$VAL_LWR <- 1e3 * all_sub_rt$VAL_LWR / all_sub_rt$POP
all_sub_rt$VAL_UPR <- 1e3 * all_sub_rt$VAL_UPR / all_sub_rt$POP
all_sub_rt$METRIC <- "Rate"
all_sub_rt$POP <- NULL
str(all_sub_rt)
```

    ## 'data.frame':    374 obs. of  8 variables:
    ##  $ YEAR         : int  2000 2000 2000 2000 2000 2000 2000 2000 2000 2000 ...
    ##  $ VAL_MEAN     : num  0.00679 0.00753 0.00564 0.00183 0.00417 ...
    ##  $ VAL_MEDIAN   : num  0.0012 0.00147 0.00144 0.0013 0.00275 ...
    ##  $ VAL_LWR      : num  3.53e-05 5.36e-05 4.97e-05 2.48e-04 1.00e-03 ...
    ##  $ VAL_UPR      : num  0.03958 0.03385 0.03029 0.00664 0.01499 ...
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
    ##  $ VAL_MEAN     : num  2.64 2.56 2.58 2.72 2.79 ...
    ##  $ VAL_MEDIAN   : num  1.14 1.1 1.11 1.18 1.2 ...
    ##  $ VAL_LWR      : num  0.0586 0.0569 0.0574 0.0605 0.062 ...
    ##  $ VAL_UPR      : num  13.3 12.9 13 13.7 14.1 ...
    ##  $ LOCATION     : chr  "Country" "Country" "Country" "Country" ...
    ##  $ LOCATION_NAME: chr  "AFG" "AFG" "AFG" "AFG" ...
    ##  $ METRIC       : chr  "Number" "Number" "Number" "Number" ...

``` r
all_cnt_rt <- all_cnt_nr%>%left_join(pop_all, by=c("LOCATION_NAME"="ISO3","YEAR"="YEAR"))
all_cnt_rt$VAL_MEAN <- 1e3 * all_cnt_rt$VAL_MEAN / all_cnt_rt$POP
all_cnt_rt$VAL_MEDIAN <- 1e3 * all_cnt_rt$VAL_MEDIAN / all_cnt_rt$POP
all_cnt_rt$VAL_LWR <- 1e3 * all_cnt_rt$VAL_LWR / all_cnt_rt$POP
all_cnt_rt$VAL_UPR <- 1e3 * all_cnt_rt$VAL_UPR / all_cnt_rt$POP
all_cnt_rt$LOCATION <- "Country"
all_cnt_rt$METRIC <- "Rate"
all_cnt_rt$POP <- NULL
str(all_cnt_rt)
```

    ## 'data.frame':    4268 obs. of  8 variables:
    ##  $ YEAR         : int  2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 ...
    ##  $ VAL_MEAN     : num  0.00258 0.00258 0.00258 0.00258 0.00258 ...
    ##  $ VAL_MEDIAN   : num  0.00111 0.00111 0.00111 0.00111 0.00111 ...
    ##  $ VAL_LWR      : num  5.73e-05 5.73e-05 5.73e-05 5.73e-05 5.73e-05 ...
    ##  $ VAL_UPR      : num  0.013 0.013 0.013 0.013 0.013 ...
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
    ##  $ VAL_MEAN     : num  0.00333 0.00335 0.00337 0.00339 0.0034 ...
    ##  $ VAL_MEDIAN   : num  0.00196 0.00196 0.00196 0.00196 0.00196 ...
    ##  $ VAL_LWR      : num  0.000684 0.000682 0.00068 0.000678 0.000676 ...
    ##  $ VAL_UPR      : num  0.0135 0.0135 0.0137 0.0139 0.0139 ...
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
| 2010 | 497.3967 |   281.7038 | 94.74635 | 2043.520 |
| 2020 | 500.7592 |   269.9118 | 86.55721 | 2125.185 |

Global number of Profound cases, 2010 vs 2020

## Regions

``` r
kbl(subset(all_reg_rt, YEAR == 2020)[,c(7,2:5)],
    align = c("l", "c", "c", "c"), row.names = FALSE,
    col.names = c("Region", "Mean", "Median", "Lower", "Upper"),
    caption=paste0("Incidence per 1000 life births of ",params$Pathogen," by WHO region in 2020")) %>%
  kable_styling("striped", "hover")
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<caption>
Incidence per 1000 life births of Profound by WHO region in 2020
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
0.0065498
</td>
<td style="text-align:center;">
0.0018117
</td>
<td style="text-align:center;">
0.0001031
</td>
<td style="text-align:left;">
0.0348828
</td>
</tr>
<tr>
<td style="text-align:left;">
AMR
</td>
<td style="text-align:center;">
0.0035049
</td>
<td style="text-align:center;">
0.0024566
</td>
<td style="text-align:center;">
0.0009204
</td>
<td style="text-align:left;">
0.0120867
</td>
</tr>
<tr>
<td style="text-align:left;">
EMR
</td>
<td style="text-align:center;">
0.0031076
</td>
<td style="text-align:center;">
0.0015647
</td>
<td style="text-align:center;">
0.0002019
</td>
<td style="text-align:left;">
0.0151112
</td>
</tr>
<tr>
<td style="text-align:left;">
EUR
</td>
<td style="text-align:center;">
0.0040162
</td>
<td style="text-align:center;">
0.0017870
</td>
<td style="text-align:center;">
0.0001710
</td>
<td style="text-align:left;">
0.0204064
</td>
</tr>
<tr>
<td style="text-align:left;">
SEAR
</td>
<td style="text-align:center;">
0.0012985
</td>
<td style="text-align:center;">
0.0007936
</td>
<td style="text-align:center;">
0.0001087
</td>
<td style="text-align:left;">
0.0052679
</td>
</tr>
<tr>
<td style="text-align:left;">
WPR
</td>
<td style="text-align:center;">
0.0031086
</td>
<td style="text-align:center;">
0.0018360
</td>
<td style="text-align:center;">
0.0005294
</td>
<td style="text-align:left;">
0.0123775
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
Cases of Profound by WHO region in 2020
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
251.11017
</td>
<td style="text-align:center;">
69.45743
</td>
<td style="text-align:center;">
3.951877
</td>
<td style="text-align:left;">
1337.3600
</td>
</tr>
<tr>
<td style="text-align:left;">
AMR
</td>
<td style="text-align:center;">
47.55424
</td>
<td style="text-align:center;">
33.33046
</td>
<td style="text-align:center;">
12.487681
</td>
<td style="text-align:left;">
163.9900
</td>
</tr>
<tr>
<td style="text-align:left;">
EMR
</td>
<td style="text-align:center;">
59.35215
</td>
<td style="text-align:center;">
29.88398
</td>
<td style="text-align:center;">
3.856823
</td>
<td style="text-align:left;">
288.6118
</td>
</tr>
<tr>
<td style="text-align:left;">
EUR
</td>
<td style="text-align:center;">
41.89595
</td>
<td style="text-align:center;">
18.64143
</td>
<td style="text-align:center;">
1.783340
</td>
<td style="text-align:left;">
212.8728
</td>
</tr>
<tr>
<td style="text-align:left;">
SEAR
</td>
<td style="text-align:center;">
44.48916
</td>
<td style="text-align:center;">
27.19202
</td>
<td style="text-align:center;">
3.725214
</td>
<td style="text-align:left;">
180.4923
</td>
</tr>
<tr>
<td style="text-align:left;">
WPR
</td>
<td style="text-align:center;">
56.35750
</td>
<td style="text-align:center;">
33.28633
</td>
<td style="text-align:center;">
9.598630
</td>
<td style="text-align:left;">
224.3990
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
  ggtitle(paste0("Incidence per 1000 life births of ", params$Pathogen, " by WHO region, 2010"))
dev.off()
```

    ## png 
    ##   2

``` r
setwd(params$Dir)
image <- paste0("03-estimate_v5_files/figure-gfm/r_CI_2010.png")
cat("![](",image,")")
```

![](03-estimate_v5_files/figure-gfm/r_CI_2010.png)

``` r
png(paste0(params$PlotDir, "/r_CI_2020.png"), width=480, height=480)
ggplot(subset(all_reg_rt, YEAR==2020),
       aes(y = VAL_MEAN, x = LOCATION_NAME)) +
  geom_pointrange(aes(ymin = VAL_LWR, ymax = VAL_UPR), size = 0.2) +
  coord_flip() +
  theme_bw() +
  scale_x_discrete(NULL, limits = rev(unique(all_reg_rt$LOCATION_NAME))) +
  scale_y_continuous(NULL) +
  ggtitle(paste0("Incidence per 1000 life births of ", params$Pathogen, " by WHO region, 2020"))
dev.off()
```

    ## png 
    ##   2

``` r
setwd(params$Dir)
image <- paste0("03-estimate_v5_files/figure-gfm/r_CI_2020.png")
cat("![](",image,")")
```

![](03-estimate_v5_files/figure-gfm/r_CI_2020.png)

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
image <- paste0("03-estimate_v5_files/figure-gfm/r_CASES_2010.png")
cat("![](",image,")")
```

![](03-estimate_v5_files/figure-gfm/r_CASES_2010.png)

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
image <- paste0("03-estimate_v5_files/figure-gfm/r_CASES_2020.png")
cat("![](",image,")")
```

![](03-estimate_v5_files/figure-gfm/r_CASES_2020.png)

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
  ggtitle(paste0("Incidence per 1000 life births of ", params$Pathogen, " by WHO region, 2010"))
dev.off()
```

    ## png 
    ##   2

``` r
setwd(params$Dir)
image <- paste0("03-estimate_v5_files/figure-gfm/r_hist_2010.png")
cat("![](",image,")")
```

![](03-estimate_v5_files/figure-gfm/r_hist_2010.png)

``` r
png(paste0(params$PlotDir, "/r_hist_2020_2010.png"), width=480, height=480)
ggplot(subset(sim_all_reg_long, YEAR==2010), aes(x = CASES)) +
  geom_density() +
  facet_wrap(~REG2) +
  theme_bw() +
  scale_x_log10() +
  ggtitle(paste0("Incidence per 1000 life births of ", params$Pathogen, " by WHO region, 2010"))
dev.off()
```

    ## png 
    ##   2

``` r
setwd(params$Dir)
image <- paste0("03-estimate_v5_files/figure-gfm/r_hist_2020_2010.png")
cat("![](",image,")")
```

![](03-estimate_v5_files/figure-gfm/r_hist_2020_2010.png)

## Subregions

``` r
kbl(subset(all_sub_rt, YEAR == 2020)[,c(7,2:5)],
    align = c("l", "c", "c", "c"), row.names = FALSE,
    col.names = c("Region", "Mean", "Median", "Lower", "Upper"),
    caption=paste0("Incidence per 1000 life births of ",params$Pathogen," by WHO subregion in 2020")) %>%
  kable_styling("striped", "hover")
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
<caption>
Incidence per 1000 life births of Profound by WHO subregion in 2020
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
0.0067052
</td>
<td style="text-align:center;">
0.0012046
</td>
<td style="text-align:center;">
0.0000363
</td>
<td style="text-align:left;">
0.0389227
</td>
</tr>
<tr>
<td style="text-align:left;">
AFRC
</td>
<td style="text-align:center;">
0.0073517
</td>
<td style="text-align:center;">
0.0014664
</td>
<td style="text-align:center;">
0.0000533
</td>
<td style="text-align:left;">
0.0337623
</td>
</tr>
<tr>
<td style="text-align:left;">
AFRD
</td>
<td style="text-align:center;">
0.0056952
</td>
<td style="text-align:center;">
0.0014331
</td>
<td style="text-align:center;">
0.0000501
</td>
<td style="text-align:left;">
0.0308913
</td>
</tr>
<tr>
<td style="text-align:left;">
AMRA
</td>
<td style="text-align:center;">
0.0018465
</td>
<td style="text-align:center;">
0.0013088
</td>
<td style="text-align:center;">
0.0002480
</td>
<td style="text-align:left;">
0.0067092
</td>
</tr>
<tr>
<td style="text-align:left;">
AMRB
</td>
<td style="text-align:center;">
0.0042322
</td>
<td style="text-align:center;">
0.0027703
</td>
<td style="text-align:center;">
0.0009966
</td>
<td style="text-align:left;">
0.0154943
</td>
</tr>
<tr>
<td style="text-align:left;">
AMRC
</td>
<td style="text-align:center;">
0.0046133
</td>
<td style="text-align:center;">
0.0016963
</td>
<td style="text-align:center;">
0.0001223
</td>
<td style="text-align:left;">
0.0222774
</td>
</tr>
<tr>
<td style="text-align:left;">
EMRA
</td>
<td style="text-align:center;">
0.0075182
</td>
<td style="text-align:center;">
0.0013542
</td>
<td style="text-align:center;">
0.0000678
</td>
<td style="text-align:left;">
0.0275128
</td>
</tr>
<tr>
<td style="text-align:left;">
EMRBC
</td>
<td style="text-align:center;">
0.0027567
</td>
<td style="text-align:center;">
0.0013373
</td>
<td style="text-align:center;">
0.0001462
</td>
<td style="text-align:left;">
0.0132337
</td>
</tr>
<tr>
<td style="text-align:left;">
EMRD
</td>
<td style="text-align:center;">
0.0032837
</td>
<td style="text-align:center;">
0.0012205
</td>
<td style="text-align:center;">
0.0000601
</td>
<td style="text-align:left;">
0.0184410
</td>
</tr>
<tr>
<td style="text-align:left;">
EURA
</td>
<td style="text-align:center;">
0.0044739
</td>
<td style="text-align:center;">
0.0016675
</td>
<td style="text-align:center;">
0.0001070
</td>
<td style="text-align:left;">
0.0243308
</td>
</tr>
<tr>
<td style="text-align:left;">
EURB
</td>
<td style="text-align:center;">
0.0038511
</td>
<td style="text-align:center;">
0.0012493
</td>
<td style="text-align:center;">
0.0000631
</td>
<td style="text-align:left;">
0.0223434
</td>
</tr>
<tr>
<td style="text-align:left;">
EURC
</td>
<td style="text-align:center;">
0.0030062
</td>
<td style="text-align:center;">
0.0011932
</td>
<td style="text-align:center;">
0.0000696
</td>
<td style="text-align:left;">
0.0155357
</td>
</tr>
<tr>
<td style="text-align:left;">
SEARB
</td>
<td style="text-align:center;">
0.0024739
</td>
<td style="text-align:center;">
0.0009489
</td>
<td style="text-align:center;">
0.0000599
</td>
<td style="text-align:left;">
0.0132990
</td>
</tr>
<tr>
<td style="text-align:left;">
SEARCD
</td>
<td style="text-align:center;">
0.0010876
</td>
<td style="text-align:center;">
0.0006420
</td>
<td style="text-align:center;">
0.0000720
</td>
<td style="text-align:left;">
0.0044033
</td>
</tr>
<tr>
<td style="text-align:left;">
WPRA
</td>
<td style="text-align:center;">
0.0110943
</td>
<td style="text-align:center;">
0.0053249
</td>
<td style="text-align:center;">
0.0008562
</td>
<td style="text-align:left;">
0.0557924
</td>
</tr>
<tr>
<td style="text-align:left;">
WPRB
</td>
<td style="text-align:center;">
0.0010287
</td>
<td style="text-align:center;">
0.0007814
</td>
<td style="text-align:center;">
0.0002004
</td>
<td style="text-align:left;">
0.0029784
</td>
</tr>
<tr>
<td style="text-align:left;">
WPRC
</td>
<td style="text-align:center;">
0.0062378
</td>
<td style="text-align:center;">
0.0020288
</td>
<td style="text-align:center;">
0.0001382
</td>
<td style="text-align:left;">
0.0338209
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
Cases of Profound by WHO sub region in 2020
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
9.726772
</td>
<td style="text-align:center;">
1.747388
</td>
<td style="text-align:center;">
0.0526918
</td>
<td style="text-align:left;">
56.46237
</td>
</tr>
<tr>
<td style="text-align:left;">
AFRC
</td>
<td style="text-align:center;">
138.903950
</td>
<td style="text-align:center;">
27.706423
</td>
<td style="text-align:center;">
1.0078404
</td>
<td style="text-align:left;">
637.90741
</td>
</tr>
<tr>
<td style="text-align:left;">
AFRD
</td>
<td style="text-align:center;">
102.479448
</td>
<td style="text-align:center;">
25.787627
</td>
<td style="text-align:center;">
0.9023812
</td>
<td style="text-align:left;">
555.85675
</td>
</tr>
<tr>
<td style="text-align:left;">
AMRA
</td>
<td style="text-align:center;">
8.026870
</td>
<td style="text-align:center;">
5.689435
</td>
<td style="text-align:center;">
1.0782901
</td>
<td style="text-align:left;">
29.16576
</td>
</tr>
<tr>
<td style="text-align:left;">
AMRB
</td>
<td style="text-align:center;">
33.434525
</td>
<td style="text-align:center;">
21.885219
</td>
<td style="text-align:center;">
7.8730552
</td>
<td style="text-align:left;">
122.40431
</td>
</tr>
<tr>
<td style="text-align:left;">
AMRC
</td>
<td style="text-align:center;">
6.092846
</td>
<td style="text-align:center;">
2.240387
</td>
<td style="text-align:center;">
0.1615803
</td>
<td style="text-align:left;">
29.42241
</td>
</tr>
<tr>
<td style="text-align:left;">
EMRA
</td>
<td style="text-align:center;">
6.005315
</td>
<td style="text-align:center;">
1.081709
</td>
<td style="text-align:center;">
0.0541515
</td>
<td style="text-align:left;">
21.97625
</td>
</tr>
<tr>
<td style="text-align:left;">
EMRBC
</td>
<td style="text-align:center;">
35.293373
</td>
<td style="text-align:center;">
17.120252
</td>
<td style="text-align:center;">
1.8714571
</td>
<td style="text-align:left;">
169.42535
</td>
</tr>
<tr>
<td style="text-align:left;">
EMRD
</td>
<td style="text-align:center;">
18.053465
</td>
<td style="text-align:center;">
6.709981
</td>
<td style="text-align:center;">
0.3303753
</td>
<td style="text-align:left;">
101.38632
</td>
</tr>
<tr>
<td style="text-align:left;">
EURA
</td>
<td style="text-align:center;">
22.297853
</td>
<td style="text-align:center;">
8.310999
</td>
<td style="text-align:center;">
0.5334114
</td>
<td style="text-align:left;">
121.26468
</td>
</tr>
<tr>
<td style="text-align:left;">
EURB
</td>
<td style="text-align:center;">
14.683828
</td>
<td style="text-align:center;">
4.763410
</td>
<td style="text-align:center;">
0.2404082
</td>
<td style="text-align:left;">
85.19354
</td>
</tr>
<tr>
<td style="text-align:left;">
EURC
</td>
<td style="text-align:center;">
4.914274
</td>
<td style="text-align:center;">
1.950644
</td>
<td style="text-align:center;">
0.1137824
</td>
<td style="text-align:left;">
25.39673
</td>
</tr>
<tr>
<td style="text-align:left;">
SEARB
</td>
<td style="text-align:center;">
12.893677
</td>
<td style="text-align:center;">
4.945608
</td>
<td style="text-align:center;">
0.3119512
</td>
<td style="text-align:left;">
69.31376
</td>
</tr>
<tr>
<td style="text-align:left;">
SEARCD
</td>
<td style="text-align:center;">
31.595480
</td>
<td style="text-align:center;">
18.651225
</td>
<td style="text-align:center;">
2.0901986
</td>
<td style="text-align:left;">
127.91907
</td>
</tr>
<tr>
<td style="text-align:left;">
WPRA
</td>
<td style="text-align:center;">
16.832463
</td>
<td style="text-align:center;">
8.079051
</td>
<td style="text-align:center;">
1.2991072
</td>
<td style="text-align:left;">
84.64937
</td>
</tr>
<tr>
<td style="text-align:left;">
WPRB
</td>
<td style="text-align:center;">
12.658970
</td>
<td style="text-align:center;">
9.614995
</td>
<td style="text-align:center;">
2.4657077
</td>
<td style="text-align:left;">
36.65065
</td>
</tr>
<tr>
<td style="text-align:left;">
WPRC
</td>
<td style="text-align:center;">
26.866065
</td>
<td style="text-align:center;">
8.738054
</td>
<td style="text-align:center;">
0.5953773
</td>
<td style="text-align:left;">
145.66595
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
  ggtitle(paste0("Incidence per 1000 life births of ", params$Pathogen, " by WHO sub region, 2010"))
dev.off()
```

    ## png 
    ##   2

``` r
setwd(params$Dir)
image <- paste0("03-estimate_v5_files/figure-gfm/r_CI_SUB2_2010.png")
cat("![](",image,")")
```

![](03-estimate_v5_files/figure-gfm/r_CI_SUB2_2010.png)

``` r
png(paste0(params$PlotDir, "/r_CI_SUB2_2020.png"), width=480, height=480)
ggplot(subset(all_sub_rt, YEAR==2020),
       aes(y = VAL_MEAN, x = LOCATION_NAME)) +
  geom_pointrange(aes(ymin = VAL_LWR, ymax = VAL_UPR), size = 0.2) +
  coord_flip() +
  theme_bw() +
  scale_x_discrete(NULL, limits = rev(unique(all_sub_rt$LOCATION_NAME))) +
  scale_y_continuous(NULL) +
  ggtitle(paste0("Incidence per 1000 life births of ", params$Pathogen, " by WHO sub region, 2020"))
dev.off()
```

    ## png 
    ##   2

``` r
setwd(params$Dir)
image <- paste0("03-estimate_v5_files/figure-gfm/r_CI_SUB2_2020.png")
cat("![](",image,")")
```

![](03-estimate_v5_files/figure-gfm/r_CI_SUB2_2020.png)

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
image <- paste0("03-estimate_v5_files/figure-gfm/r_CASES_SUB2_2010.png")
cat("![](",image,")")
```

![](03-estimate_v5_files/figure-gfm/r_CASES_SUB2_2010.png)

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
image <- paste0("03-estimate_v5_files/figure-gfm/r_CASES_SUB2_2020.png")
cat("![](",image,")")
```

![](03-estimate_v5_files/figure-gfm/r_CASES_SUB2_2020.png)

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
  ggtitle(paste0("Incidence per 1000 life births of ", params$Pathogen, "by WHO sub region, 2010"))
dev.off()
```

    ## png 
    ##   2

``` r
setwd(params$Dir)
image <- paste0("03-estimate_v5_files/figure-gfm/r_hist_SUB2_2010.png")
cat("![](",image,")")
```

![](03-estimate_v5_files/figure-gfm/r_hist_SUB2_2010.png)

``` r
png(paste0(params$PlotDir, "/r_hist_SUB2_2020_2010.png"), width=480, height=480)
ggplot(subset(sim_all_sub_long, YEAR==2010), aes(x = CASES)) +
  geom_density() +
  facet_wrap(~SUB2) +
  theme_bw() +
  scale_x_log10() +
  ggtitle(paste0("Incidence per 1000 life births of ", params$Pathogen, "by WHO sub region, 2010"))
dev.off()
```

    ## png 
    ##   2

``` r
setwd(params$Dir)
image <- paste0("03-estimate_v5_files/figure-gfm/r_hist_SUB2_2020_2010.png")
cat("![](",image,")")
```

![](03-estimate_v5_files/figure-gfm/r_hist_SUB2_2020_2010.png)

## Countries

``` r
png(paste0(params$PlotDir, "/r_cnt_2010.png"), width=800, height=300)
plot_world(subset(all_cnt_rt, YEAR == 2010),
           "LOCATION_NAME", "VAL_MEAN", legend.title = "Incidence per 1000", diseasefree = zero_cases)
```

    ## [1] 0.00 0.02 0.04 0.06 0.08

``` r
dev.off()
```

    ## png 
    ##   2

``` r
setwd(params$Dir)
image <- paste0("03-estimate_v5_files/figure-gfm/r_cnt_2010.png")
cat("![](",image,")")
```

![](03-estimate_v5_files/figure-gfm/r_cnt_2010.png)

``` r
png(paste0(params$PlotDir, "/r_cnt_2020.png"), width=800, height=300)
plot_world(subset(all_cnt_rt, YEAR == 2020),
           "LOCATION_NAME", "VAL_MEAN", legend.title = "Incidence per 1000", diseasefree = zero_cases)
```

    ## [1] 0.00 0.02 0.04 0.06 0.08

``` r
dev.off()
```

    ## png 
    ##   2

``` r
setwd(params$Dir)
image <- paste0("03-estimate_v5_files/figure-gfm/r_cnt_2020.png")
cat("![](",image,")")
```

![](03-estimate_v5_files/figure-gfm/r_cnt_2020.png)

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
| Afghanistan                      |     0.003 |       0.001 |    0.000 |    0.013 |     0.003 |       0.001 |    0.000 |    0.013 |
| Angola                           |     0.004 |       0.001 |    0.000 |    0.018 |     0.004 |       0.001 |    0.000 |    0.018 |
| Albania                          |     0.005 |       0.001 |    0.000 |    0.028 |     0.005 |       0.001 |    0.000 |    0.028 |
| Andorra                          |     0.002 |       0.001 |    0.000 |    0.010 |     0.002 |       0.001 |    0.000 |    0.010 |
| United Arab Emirates             |     0.006 |       0.001 |    0.000 |    0.033 |     0.006 |       0.001 |    0.000 |    0.033 |
| Argentina                        |     0.006 |       0.002 |    0.000 |    0.030 |     0.006 |       0.002 |    0.000 |    0.030 |
| Armenia                          |     0.002 |       0.001 |    0.000 |    0.010 |     0.002 |       0.001 |    0.000 |    0.010 |
| Antigua and Barbuda              |     0.002 |       0.001 |    0.000 |    0.007 |     0.002 |       0.001 |    0.000 |    0.007 |
| Australia                        |     0.014 |       0.004 |    0.000 |    0.078 |     0.014 |       0.004 |    0.000 |    0.078 |
| Austria                          |     0.005 |       0.001 |    0.000 |    0.026 |     0.005 |       0.001 |    0.000 |    0.026 |
| Azerbaijan                       |     0.002 |       0.001 |    0.000 |    0.010 |     0.002 |       0.001 |    0.000 |    0.010 |
| Burundi                          |     0.004 |       0.001 |    0.000 |    0.017 |     0.004 |       0.001 |    0.000 |    0.017 |
| Belgium                          |     0.005 |       0.001 |    0.000 |    0.029 |     0.005 |       0.001 |    0.000 |    0.029 |
| Benin                            |     0.008 |       0.001 |    0.000 |    0.040 |     0.008 |       0.001 |    0.000 |    0.040 |
| Burkina Faso                     |     0.004 |       0.001 |    0.000 |    0.017 |     0.004 |       0.001 |    0.000 |    0.017 |
| Bangladesh                       |     0.003 |       0.001 |    0.000 |    0.015 |     0.003 |       0.001 |    0.000 |    0.015 |
| Bulgaria                         |     0.002 |       0.001 |    0.000 |    0.010 |     0.002 |       0.001 |    0.000 |    0.010 |
| Bahrain                          |     0.003 |       0.001 |    0.000 |    0.013 |     0.003 |       0.001 |    0.000 |    0.013 |
| Bahamas                          |     0.002 |       0.001 |    0.000 |    0.007 |     0.002 |       0.001 |    0.000 |    0.007 |
| Bosnia and Herzegovina           |     0.002 |       0.001 |    0.000 |    0.010 |     0.002 |       0.001 |    0.000 |    0.010 |
| Belarus                          |     0.002 |       0.001 |    0.000 |    0.010 |     0.002 |       0.001 |    0.000 |    0.010 |
| Belize                           |     0.002 |       0.002 |    0.000 |    0.008 |     0.002 |       0.002 |    0.000 |    0.008 |
| Bolivia                          |     0.007 |       0.001 |    0.000 |    0.033 |     0.007 |       0.001 |    0.000 |    0.033 |
| Brazil                           |     0.002 |       0.002 |    0.001 |    0.004 |     0.002 |       0.002 |    0.001 |    0.004 |
| Barbados                         |     0.002 |       0.001 |    0.000 |    0.007 |     0.002 |       0.001 |    0.000 |    0.007 |
| Brunei Darussalam                |     0.007 |       0.003 |    0.001 |    0.040 |     0.007 |       0.003 |    0.001 |    0.040 |
| Bhutan                           |     0.001 |       0.001 |    0.000 |    0.005 |     0.001 |       0.001 |    0.000 |    0.005 |
| Botswana                         |     0.004 |       0.001 |    0.000 |    0.019 |     0.004 |       0.001 |    0.000 |    0.019 |
| Central African Republic         |     0.004 |       0.001 |    0.000 |    0.017 |     0.004 |       0.001 |    0.000 |    0.017 |
| Canada                           |     0.003 |       0.001 |    0.000 |    0.018 |     0.003 |       0.001 |    0.000 |    0.018 |
| Switzerland                      |     0.005 |       0.001 |    0.000 |    0.027 |     0.005 |       0.001 |    0.000 |    0.027 |
| Chile                            |     0.002 |       0.001 |    0.000 |    0.007 |     0.002 |       0.001 |    0.000 |    0.007 |
| China                            |     0.001 |       0.001 |    0.000 |    0.003 |     0.001 |       0.001 |    0.000 |    0.003 |
| Côte d’Ivoire                    |     0.004 |       0.001 |    0.000 |    0.018 |     0.004 |       0.001 |    0.000 |    0.018 |
| Cameroon                         |     0.004 |       0.001 |    0.000 |    0.018 |     0.004 |       0.001 |    0.000 |    0.018 |
| Congo                            |     0.008 |       0.001 |    0.000 |    0.042 |     0.008 |       0.001 |    0.000 |    0.042 |
| Congo                            |     0.006 |       0.001 |    0.000 |    0.040 |     0.006 |       0.001 |    0.000 |    0.040 |
| Cook Islands                     |     0.007 |       0.003 |    0.001 |    0.040 |     0.007 |       0.003 |    0.001 |    0.040 |
| Colombia                         |     0.005 |       0.002 |    0.000 |    0.024 |     0.005 |       0.002 |    0.000 |    0.024 |
| Comoros                          |     0.004 |       0.001 |    0.000 |    0.018 |     0.004 |       0.001 |    0.000 |    0.018 |
| Cabo Verde                       |     0.009 |       0.001 |    0.000 |    0.044 |     0.009 |       0.001 |    0.000 |    0.044 |
| Costa Rica                       |     0.002 |       0.002 |    0.000 |    0.008 |     0.002 |       0.002 |    0.000 |    0.008 |
| Cuba                             |     0.002 |       0.002 |    0.000 |    0.008 |     0.002 |       0.002 |    0.000 |    0.008 |
| Cyprus                           |     0.005 |       0.001 |    0.000 |    0.027 |     0.005 |       0.001 |    0.000 |    0.027 |
| Czechia                          |     0.005 |       0.001 |    0.000 |    0.025 |     0.005 |       0.001 |    0.000 |    0.025 |
| Germany                          |     0.005 |       0.001 |    0.000 |    0.026 |     0.005 |       0.001 |    0.000 |    0.026 |
| Djibouti                         |     0.002 |       0.001 |    0.000 |    0.008 |     0.002 |       0.001 |    0.000 |    0.008 |
| Dominica                         |     0.002 |       0.002 |    0.000 |    0.008 |     0.002 |       0.002 |    0.000 |    0.008 |
| Denmark                          |     0.004 |       0.001 |    0.000 |    0.025 |     0.004 |       0.001 |    0.000 |    0.025 |
| Dominican Republic               |     0.002 |       0.002 |    0.000 |    0.008 |     0.002 |       0.002 |    0.000 |    0.008 |
| Algeria                          |     0.004 |       0.001 |    0.000 |    0.018 |     0.004 |       0.001 |    0.000 |    0.018 |
| Ecuador                          |     0.008 |       0.002 |    0.000 |    0.031 |     0.008 |       0.002 |    0.000 |    0.031 |
| Egypt                            |     0.004 |       0.001 |    0.000 |    0.019 |     0.004 |       0.001 |    0.000 |    0.019 |
| Eritrea                          |     0.007 |       0.001 |    0.000 |    0.040 |     0.007 |       0.001 |    0.000 |    0.040 |
| Spain                            |     0.003 |       0.001 |    0.000 |    0.020 |     0.003 |       0.001 |    0.000 |    0.020 |
| Estonia                          |     0.002 |       0.001 |    0.000 |    0.010 |     0.002 |       0.001 |    0.000 |    0.010 |
| Ethiopia                         |     0.007 |       0.001 |    0.000 |    0.042 |     0.007 |       0.001 |    0.000 |    0.042 |
| Finland                          |     0.006 |       0.001 |    0.000 |    0.026 |     0.006 |       0.001 |    0.000 |    0.026 |
| Fiji                             |     0.002 |       0.001 |    0.000 |    0.007 |     0.002 |       0.001 |    0.000 |    0.007 |
| France                           |     0.004 |       0.001 |    0.000 |    0.024 |     0.004 |       0.001 |    0.000 |    0.024 |
| Micronesia                       |     0.003 |       0.002 |    0.000 |    0.016 |     0.003 |       0.002 |    0.000 |    0.016 |
| Gabon                            |     0.004 |       0.001 |    0.000 |    0.019 |     0.004 |       0.001 |    0.000 |    0.019 |
| United Kingdom                   |     0.004 |       0.001 |    0.000 |    0.025 |     0.004 |       0.001 |    0.000 |    0.025 |
| Georgia                          |     0.005 |       0.001 |    0.000 |    0.027 |     0.005 |       0.001 |    0.000 |    0.027 |
| Ghana                            |     0.027 |       0.001 |    0.000 |    0.039 |     0.027 |       0.001 |    0.000 |    0.039 |
| Guinea                           |     0.004 |       0.001 |    0.000 |    0.018 |     0.004 |       0.001 |    0.000 |    0.018 |
| Gambia                           |     0.004 |       0.001 |    0.000 |    0.017 |     0.004 |       0.001 |    0.000 |    0.017 |
| Guinea-Bissau                    |     0.004 |       0.001 |    0.000 |    0.017 |     0.004 |       0.001 |    0.000 |    0.017 |
| Equatorial Guinea                |     0.004 |       0.001 |    0.000 |    0.019 |     0.004 |       0.001 |    0.000 |    0.019 |
| Greece                           |     0.004 |       0.001 |    0.000 |    0.024 |     0.004 |       0.001 |    0.000 |    0.024 |
| Grenada                          |     0.002 |       0.002 |    0.000 |    0.008 |     0.002 |       0.002 |    0.000 |    0.008 |
| Guatemala                        |     0.002 |       0.002 |    0.000 |    0.008 |     0.002 |       0.002 |    0.000 |    0.008 |
| Guyana                           |     0.002 |       0.001 |    0.000 |    0.007 |     0.002 |       0.001 |    0.000 |    0.007 |
| Honduras                         |     0.003 |       0.001 |    0.000 |    0.014 |     0.003 |       0.001 |    0.000 |    0.014 |
| Croatia                          |     0.005 |       0.001 |    0.000 |    0.029 |     0.005 |       0.001 |    0.000 |    0.029 |
| Haiti                            |     0.003 |       0.001 |    0.000 |    0.014 |     0.003 |       0.001 |    0.000 |    0.014 |
| Hungary                          |     0.006 |       0.001 |    0.000 |    0.026 |     0.006 |       0.001 |    0.000 |    0.026 |
| Indonesia                        |     0.002 |       0.001 |    0.000 |    0.012 |     0.002 |       0.001 |    0.000 |    0.012 |
| India                            |     0.001 |       0.000 |    0.000 |    0.004 |     0.001 |       0.000 |    0.000 |    0.004 |
| Ireland                          |     0.005 |       0.001 |    0.000 |    0.028 |     0.005 |       0.001 |    0.000 |    0.028 |
| Iran                             |     0.003 |       0.001 |    0.000 |    0.017 |     0.003 |       0.001 |    0.000 |    0.017 |
| Iraq                             |     0.004 |       0.001 |    0.000 |    0.024 |     0.004 |       0.001 |    0.000 |    0.024 |
| Iceland                          |     0.004 |       0.001 |    0.000 |    0.027 |     0.004 |       0.001 |    0.000 |    0.027 |
| Israel                           |     0.002 |       0.001 |    0.000 |    0.010 |     0.002 |       0.001 |    0.000 |    0.010 |
| Italy                            |     0.004 |       0.001 |    0.000 |    0.021 |     0.004 |       0.001 |    0.000 |    0.021 |
| Jamaica                          |     0.006 |       0.002 |    0.000 |    0.032 |     0.006 |       0.002 |    0.000 |    0.032 |
| Jordan                           |     0.004 |       0.001 |    0.000 |    0.024 |     0.004 |       0.001 |    0.000 |    0.024 |
| Japan                            |     0.006 |       0.003 |    0.000 |    0.034 |     0.006 |       0.003 |    0.000 |    0.034 |
| Kazakhstan                       |     0.006 |       0.001 |    0.000 |    0.028 |     0.006 |       0.001 |    0.000 |    0.028 |
| Kenya                            |     0.008 |       0.001 |    0.000 |    0.040 |     0.008 |       0.001 |    0.000 |    0.040 |
| Kyrgyzstan                       |     0.002 |       0.001 |    0.000 |    0.012 |     0.002 |       0.001 |    0.000 |    0.012 |
| Cambodia                         |     0.008 |       0.002 |    0.000 |    0.040 |     0.008 |       0.002 |    0.000 |    0.040 |
| Kiribati                         |     0.003 |       0.002 |    0.000 |    0.016 |     0.003 |       0.002 |    0.000 |    0.016 |
| Saint Kitts and Nevis            |     0.002 |       0.001 |    0.000 |    0.007 |     0.002 |       0.001 |    0.000 |    0.007 |
| Korea                            |     0.012 |       0.004 |    0.000 |    0.069 |     0.012 |       0.004 |    0.000 |    0.069 |
| Kuwait                           |     0.006 |       0.001 |    0.000 |    0.035 |     0.006 |       0.001 |    0.000 |    0.035 |
| Lao People’s Dem. Republic       |     0.003 |       0.002 |    0.000 |    0.016 |     0.003 |       0.002 |    0.000 |    0.016 |
| Lebanon                          |     0.004 |       0.001 |    0.000 |    0.024 |     0.004 |       0.001 |    0.000 |    0.024 |
| Liberia                          |     0.004 |       0.001 |    0.000 |    0.017 |     0.004 |       0.001 |    0.000 |    0.017 |
| Libya                            |     0.004 |       0.001 |    0.000 |    0.023 |     0.004 |       0.001 |    0.000 |    0.023 |
| Saint Lucia                      |     0.002 |       0.002 |    0.000 |    0.008 |     0.002 |       0.002 |    0.000 |    0.008 |
| Sri Lanka                        |     0.003 |       0.001 |    0.000 |    0.017 |     0.003 |       0.001 |    0.000 |    0.017 |
| Lesotho                          |     0.004 |       0.001 |    0.000 |    0.018 |     0.004 |       0.001 |    0.000 |    0.018 |
| Lithuania                        |     0.005 |       0.001 |    0.000 |    0.025 |     0.005 |       0.001 |    0.000 |    0.025 |
| Luxembourg                       |     0.005 |       0.001 |    0.000 |    0.024 |     0.005 |       0.001 |    0.000 |    0.024 |
| Latvia                           |     0.002 |       0.001 |    0.000 |    0.010 |     0.002 |       0.001 |    0.000 |    0.010 |
| Morocco                          |     0.005 |       0.001 |    0.000 |    0.026 |     0.005 |       0.001 |    0.000 |    0.026 |
| Monaco                           |     0.002 |       0.001 |    0.000 |    0.010 |     0.002 |       0.001 |    0.000 |    0.010 |
| Republic of Moldova              |     0.002 |       0.001 |    0.000 |    0.010 |     0.002 |       0.001 |    0.000 |    0.010 |
| Madagascar                       |     0.004 |       0.001 |    0.000 |    0.017 |     0.004 |       0.001 |    0.000 |    0.017 |
| Maldives                         |     0.004 |       0.001 |    0.000 |    0.021 |     0.004 |       0.001 |    0.000 |    0.021 |
| Mexico                           |     0.007 |       0.003 |    0.001 |    0.035 |     0.007 |       0.003 |    0.001 |    0.035 |
| Marshall Islands                 |     0.002 |       0.001 |    0.000 |    0.007 |     0.002 |       0.001 |    0.000 |    0.007 |
| North Macedonia                  |     0.002 |       0.001 |    0.000 |    0.010 |     0.002 |       0.001 |    0.000 |    0.010 |
| Mali                             |     0.004 |       0.001 |    0.000 |    0.017 |     0.004 |       0.001 |    0.000 |    0.017 |
| Malta                            |     0.002 |       0.001 |    0.000 |    0.010 |     0.002 |       0.001 |    0.000 |    0.010 |
| Myanmar                          |     0.001 |       0.001 |    0.000 |    0.005 |     0.001 |       0.001 |    0.000 |    0.005 |
| Montenegro                       |     0.002 |       0.001 |    0.000 |    0.010 |     0.002 |       0.001 |    0.000 |    0.010 |
| Mongolia                         |     0.007 |       0.002 |    0.000 |    0.046 |     0.007 |       0.002 |    0.000 |    0.046 |
| Mozambique                       |     0.008 |       0.001 |    0.000 |    0.045 |     0.008 |       0.001 |    0.000 |    0.045 |
| Mauritania                       |     0.004 |       0.001 |    0.000 |    0.018 |     0.004 |       0.001 |    0.000 |    0.018 |
| Mauritius                        |     0.004 |       0.001 |    0.000 |    0.019 |     0.004 |       0.001 |    0.000 |    0.019 |
| Malawi                           |     0.004 |       0.001 |    0.000 |    0.017 |     0.004 |       0.001 |    0.000 |    0.017 |
| Malaysia                         |     0.005 |       0.001 |    0.000 |    0.021 |     0.005 |       0.001 |    0.000 |    0.021 |
| Namibia                          |     0.004 |       0.001 |    0.000 |    0.019 |     0.004 |       0.001 |    0.000 |    0.019 |
| Niger                            |     0.004 |       0.001 |    0.000 |    0.017 |     0.004 |       0.001 |    0.000 |    0.017 |
| Nigeria                          |     0.008 |       0.001 |    0.000 |    0.040 |     0.008 |       0.001 |    0.000 |    0.040 |
| Nicaragua                        |     0.003 |       0.001 |    0.000 |    0.014 |     0.003 |       0.001 |    0.000 |    0.014 |
| Niue                             |     0.007 |       0.003 |    0.001 |    0.040 |     0.007 |       0.003 |    0.001 |    0.040 |
| Netherlands                      |     0.006 |       0.001 |    0.000 |    0.025 |     0.006 |       0.001 |    0.000 |    0.025 |
| Norway                           |     0.005 |       0.001 |    0.000 |    0.026 |     0.005 |       0.001 |    0.000 |    0.026 |
| Nepal                            |     0.001 |       0.001 |    0.000 |    0.005 |     0.001 |       0.001 |    0.000 |    0.005 |
| Nauru                            |     0.007 |       0.003 |    0.001 |    0.040 |     0.007 |       0.003 |    0.001 |    0.040 |
| New Zealand                      |     0.014 |       0.004 |    0.000 |    0.087 |     0.014 |       0.004 |    0.000 |    0.087 |
| Oman                             |     0.003 |       0.001 |    0.000 |    0.013 |     0.003 |       0.001 |    0.000 |    0.013 |
| Pakistan                         |     0.002 |       0.001 |    0.000 |    0.009 |     0.002 |       0.001 |    0.000 |    0.009 |
| Panama                           |     0.002 |       0.001 |    0.000 |    0.007 |     0.002 |       0.001 |    0.000 |    0.007 |
| Peru                             |     0.005 |       0.002 |    0.000 |    0.026 |     0.005 |       0.002 |    0.000 |    0.026 |
| Philippines                      |     0.006 |       0.002 |    0.000 |    0.033 |     0.006 |       0.002 |    0.000 |    0.033 |
| Palau                            |     0.002 |       0.001 |    0.000 |    0.007 |     0.002 |       0.001 |    0.000 |    0.007 |
| Papua New Guinea                 |     0.003 |       0.002 |    0.000 |    0.016 |     0.003 |       0.002 |    0.000 |    0.016 |
| Poland                           |     0.005 |       0.001 |    0.000 |    0.026 |     0.005 |       0.001 |    0.000 |    0.026 |
| Korea                            |     0.001 |       0.001 |    0.000 |    0.005 |     0.001 |       0.001 |    0.000 |    0.005 |
| Portugal                         |     0.004 |       0.001 |    0.000 |    0.026 |     0.004 |       0.001 |    0.000 |    0.026 |
| Paraguay                         |     0.002 |       0.002 |    0.000 |    0.008 |     0.002 |       0.002 |    0.000 |    0.008 |
| Qatar                            |     0.003 |       0.001 |    0.000 |    0.013 |     0.003 |       0.001 |    0.000 |    0.013 |
| Romania                          |     0.005 |       0.001 |    0.000 |    0.027 |     0.005 |       0.001 |    0.000 |    0.027 |
| Russian Federation               |     0.003 |       0.001 |    0.000 |    0.020 |     0.003 |       0.001 |    0.000 |    0.020 |
| Rwanda                           |     0.004 |       0.001 |    0.000 |    0.017 |     0.004 |       0.001 |    0.000 |    0.017 |
| Saudi Arabia                     |     0.009 |       0.001 |    0.000 |    0.030 |     0.009 |       0.001 |    0.000 |    0.030 |
| Sudan                            |     0.005 |       0.001 |    0.000 |    0.033 |     0.005 |       0.001 |    0.000 |    0.033 |
| Senegal                          |     0.004 |       0.001 |    0.000 |    0.018 |     0.004 |       0.001 |    0.000 |    0.018 |
| Singapore                        |     0.071 |       0.020 |    0.001 |    0.448 |     0.071 |       0.020 |    0.001 |    0.448 |
| Solomon Islands                  |     0.003 |       0.002 |    0.000 |    0.016 |     0.003 |       0.002 |    0.000 |    0.016 |
| Sierra Leone                     |     0.004 |       0.001 |    0.000 |    0.017 |     0.004 |       0.001 |    0.000 |    0.017 |
| El Salvador                      |     0.002 |       0.002 |    0.000 |    0.008 |     0.002 |       0.002 |    0.000 |    0.008 |
| San Marino                       |     0.002 |       0.001 |    0.000 |    0.010 |     0.002 |       0.001 |    0.000 |    0.010 |
| Somalia                          |     0.003 |       0.001 |    0.000 |    0.013 |     0.003 |       0.001 |    0.000 |    0.013 |
| Serbia                           |     0.005 |       0.001 |    0.000 |    0.027 |     0.005 |       0.001 |    0.000 |    0.027 |
| South Sudan                      |     0.004 |       0.001 |    0.000 |    0.017 |     0.004 |       0.001 |    0.000 |    0.017 |
| Sao Tome and Principe            |     0.004 |       0.001 |    0.000 |    0.018 |     0.004 |       0.001 |    0.000 |    0.018 |
| Suriname                         |     0.006 |       0.002 |    0.000 |    0.033 |     0.006 |       0.002 |    0.000 |    0.033 |
| Slovakia                         |     0.005 |       0.001 |    0.000 |    0.026 |     0.005 |       0.001 |    0.000 |    0.026 |
| Slovenia                         |     0.006 |       0.001 |    0.000 |    0.029 |     0.006 |       0.001 |    0.000 |    0.029 |
| Sweden                           |     0.004 |       0.001 |    0.000 |    0.027 |     0.004 |       0.001 |    0.000 |    0.027 |
| Eswatini                         |     0.004 |       0.001 |    0.000 |    0.018 |     0.004 |       0.001 |    0.000 |    0.018 |
| Seychelles                       |     0.008 |       0.001 |    0.000 |    0.043 |     0.008 |       0.001 |    0.000 |    0.043 |
| Syrian Arab Republic             |     0.003 |       0.001 |    0.000 |    0.013 |     0.003 |       0.001 |    0.000 |    0.013 |
| Chad                             |     0.004 |       0.001 |    0.000 |    0.017 |     0.004 |       0.001 |    0.000 |    0.017 |
| Togo                             |     0.004 |       0.001 |    0.000 |    0.017 |     0.004 |       0.001 |    0.000 |    0.017 |
| Thailand                         |     0.004 |       0.001 |    0.000 |    0.022 |     0.004 |       0.001 |    0.000 |    0.022 |
| Tajikistan                       |     0.002 |       0.001 |    0.000 |    0.012 |     0.002 |       0.001 |    0.000 |    0.012 |
| Turkmenistan                     |     0.002 |       0.001 |    0.000 |    0.010 |     0.002 |       0.001 |    0.000 |    0.010 |
| Timor-Leste                      |     0.001 |       0.001 |    0.000 |    0.005 |     0.001 |       0.001 |    0.000 |    0.005 |
| Tonga                            |     0.002 |       0.001 |    0.000 |    0.007 |     0.002 |       0.001 |    0.000 |    0.007 |
| Trinidad and Tobago              |     0.002 |       0.001 |    0.000 |    0.007 |     0.002 |       0.001 |    0.000 |    0.007 |
| Tunisia                          |     0.004 |       0.001 |    0.000 |    0.021 |     0.004 |       0.001 |    0.000 |    0.021 |
| Turkiye                          |     0.005 |       0.001 |    0.000 |    0.026 |     0.005 |       0.001 |    0.000 |    0.026 |
| Tuvalu                           |     0.002 |       0.001 |    0.000 |    0.007 |     0.002 |       0.001 |    0.000 |    0.007 |
| United Republic of Tanzania      |     0.008 |       0.001 |    0.000 |    0.043 |     0.008 |       0.001 |    0.000 |    0.043 |
| Uganda                           |     0.004 |       0.001 |    0.000 |    0.017 |     0.004 |       0.001 |    0.000 |    0.017 |
| Ukraine                          |     0.005 |       0.001 |    0.000 |    0.029 |     0.005 |       0.001 |    0.000 |    0.029 |
| Uruguay                          |     0.002 |       0.001 |    0.000 |    0.007 |     0.002 |       0.001 |    0.000 |    0.007 |
| United States of America         |     0.002 |       0.001 |    0.000 |    0.006 |     0.002 |       0.001 |    0.000 |    0.006 |
| Uzbekistan                       |     0.002 |       0.001 |    0.000 |    0.012 |     0.002 |       0.001 |    0.000 |    0.012 |
| Saint Vincent and the Grenadines |     0.002 |       0.002 |    0.000 |    0.008 |     0.002 |       0.002 |    0.000 |    0.008 |
| Venezuela                        |     0.006 |       0.001 |    0.000 |    0.032 |     0.006 |       0.001 |    0.000 |    0.032 |
| Viet Nam                         |     0.007 |       0.002 |    0.000 |    0.039 |     0.007 |       0.002 |    0.000 |    0.039 |
| Vanuatu                          |     0.003 |       0.002 |    0.000 |    0.016 |     0.003 |       0.002 |    0.000 |    0.016 |
| Samoa                            |     0.003 |       0.002 |    0.000 |    0.016 |     0.003 |       0.002 |    0.000 |    0.016 |
| Yemen                            |     0.003 |       0.001 |    0.000 |    0.013 |     0.003 |       0.001 |    0.000 |    0.013 |
| South Africa                     |     0.007 |       0.001 |    0.000 |    0.044 |     0.007 |       0.001 |    0.000 |    0.044 |
| Zambia                           |     0.004 |       0.001 |    0.000 |    0.018 |     0.004 |       0.001 |    0.000 |    0.018 |
| Zimbabwe                         |     0.007 |       0.001 |    0.000 |    0.040 |     0.007 |       0.001 |    0.000 |    0.040 |

Estimated Profound incidence by country, 2010 vs 2020

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
| Afghanistan                      |       3.0 |         1.3 |      0.1 |     15.3 |       3.7 |         1.6 |      0.1 |     18.5 |
| Angola                           |       3.7 |         1.2 |      0.0 |     18.9 |       4.6 |         1.5 |      0.1 |     23.6 |
| Albania                          |       0.2 |         0.0 |      0.0 |      1.0 |       0.2 |         0.0 |      0.0 |      0.9 |
| Andorra                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| United Arab Emirates             |       0.5 |         0.1 |      0.0 |      2.6 |       0.6 |         0.1 |      0.0 |      3.2 |
| Argentina                        |       4.6 |         1.5 |      0.1 |     22.7 |       3.3 |         1.0 |      0.1 |     16.0 |
| Armenia                          |       0.1 |         0.0 |      0.0 |      0.4 |       0.1 |         0.0 |      0.0 |      0.3 |
| Antigua and Barbuda              |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Australia                        |       4.2 |         1.1 |      0.0 |     23.5 |       4.1 |         1.1 |      0.0 |     23.0 |
| Austria                          |       0.4 |         0.1 |      0.0 |      2.0 |       0.4 |         0.1 |      0.0 |      2.1 |
| Azerbaijan                       |       0.3 |         0.2 |      0.0 |      1.7 |       0.3 |         0.1 |      0.0 |      1.4 |
| Burundi                          |       1.6 |         0.5 |      0.0 |      7.7 |       1.7 |         0.5 |      0.0 |      7.9 |
| Belgium                          |       0.6 |         0.1 |      0.0 |      3.7 |       0.5 |         0.1 |      0.0 |      3.3 |
| Benin                            |       3.3 |         0.4 |      0.0 |     15.5 |       3.9 |         0.5 |      0.0 |     18.5 |
| Burkina Faso                     |       2.6 |         0.8 |      0.0 |     12.1 |       2.6 |         0.8 |      0.0 |     12.2 |
| Bangladesh                       |       8.4 |         2.3 |      0.1 |     49.4 |       8.6 |         2.4 |      0.1 |     50.9 |
| Bulgaria                         |       0.2 |         0.1 |      0.0 |      0.7 |       0.1 |         0.1 |      0.0 |      0.6 |
| Bahrain                          |       0.1 |         0.0 |      0.0 |      0.3 |       0.1 |         0.0 |      0.0 |      0.3 |
| Bahamas                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Bosnia and Herzegovina           |       0.1 |         0.0 |      0.0 |      0.4 |       0.1 |         0.0 |      0.0 |      0.3 |
| Belarus                          |       0.2 |         0.1 |      0.0 |      1.0 |       0.2 |         0.1 |      0.0 |      0.8 |
| Belize                           |       0.0 |         0.0 |      0.0 |      0.1 |       0.0 |         0.0 |      0.0 |      0.1 |
| Bolivia                          |       1.8 |         0.4 |      0.0 |      8.6 |       1.8 |         0.4 |      0.0 |      8.6 |
| Brazil                           |       5.9 |         5.4 |      2.4 |     11.9 |       5.4 |         5.0 |      2.2 |     10.9 |
| Barbados                         |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Brunei Darussalam                |       0.0 |         0.0 |      0.0 |      0.3 |       0.0 |         0.0 |      0.0 |      0.3 |
| Bhutan                           |       0.0 |         0.0 |      0.0 |      0.1 |       0.0 |         0.0 |      0.0 |      0.1 |
| Botswana                         |       0.2 |         0.1 |      0.0 |      1.1 |       0.2 |         0.1 |      0.0 |      1.2 |
| Central African Republic         |       0.7 |         0.2 |      0.0 |      3.5 |       0.8 |         0.3 |      0.0 |      3.8 |
| Canada                           |       1.3 |         0.5 |      0.0 |      6.9 |       1.2 |         0.5 |      0.0 |      6.6 |
| Switzerland                      |       0.4 |         0.1 |      0.0 |      2.2 |       0.4 |         0.1 |      0.0 |      2.3 |
| Chile                            |       0.5 |         0.3 |      0.1 |      1.7 |       0.4 |         0.3 |      0.0 |      1.3 |
| China                            |      15.9 |        12.7 |      3.2 |     45.7 |      10.5 |         8.4 |      2.1 |     30.2 |
| Côte d’Ivoire                    |       3.3 |         1.0 |      0.0 |     16.7 |       3.5 |         1.1 |      0.0 |     17.8 |
| Cameroon                         |       2.8 |         0.9 |      0.0 |     14.2 |       3.3 |         1.0 |      0.0 |     16.9 |
| Congo                            |      22.4 |         3.4 |      0.1 |    126.4 |      30.2 |         4.6 |      0.1 |    170.5 |
| Congo                            |       1.1 |         0.2 |      0.0 |      7.0 |       1.2 |         0.2 |      0.0 |      7.3 |
| Cook Islands                     |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Colombia                         |       3.4 |         1.4 |      0.1 |     18.1 |       3.2 |         1.3 |      0.1 |     17.1 |
| Comoros                          |       0.1 |         0.0 |      0.0 |      0.4 |       0.1 |         0.0 |      0.0 |      0.4 |
| Cabo Verde                       |       0.1 |         0.0 |      0.0 |      0.5 |       0.1 |         0.0 |      0.0 |      0.3 |
| Costa Rica                       |       0.2 |         0.1 |      0.0 |      0.6 |       0.1 |         0.1 |      0.0 |      0.5 |
| Cuba                             |       0.3 |         0.2 |      0.1 |      1.1 |       0.3 |         0.2 |      0.0 |      0.9 |
| Cyprus                           |       0.1 |         0.0 |      0.0 |      0.4 |       0.1 |         0.0 |      0.0 |      0.4 |
| Czechia                          |       0.5 |         0.1 |      0.0 |      2.9 |       0.5 |         0.1 |      0.0 |      2.7 |
| Germany                          |       3.4 |         0.7 |      0.0 |     17.9 |       3.9 |         0.8 |      0.0 |     20.5 |
| Djibouti                         |       0.0 |         0.0 |      0.0 |      0.2 |       0.0 |         0.0 |      0.0 |      0.2 |
| Dominica                         |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Denmark                          |       0.3 |         0.1 |      0.0 |      1.6 |       0.3 |         0.1 |      0.0 |      1.6 |
| Dominican Republic               |       0.5 |         0.4 |      0.1 |      1.8 |       0.5 |         0.4 |      0.1 |      1.7 |
| Algeria                          |       3.2 |         1.0 |      0.0 |     16.3 |       3.5 |         1.1 |      0.0 |     18.2 |
| Ecuador                          |       2.5 |         0.6 |      0.0 |     10.2 |       2.2 |         0.5 |      0.0 |      8.9 |
| Egypt                            |       8.6 |         2.3 |      0.1 |     47.5 |       8.5 |         2.3 |      0.1 |     46.5 |
| Eritrea                          |       0.7 |         0.1 |      0.0 |      3.9 |       0.6 |         0.1 |      0.0 |      3.8 |
| Spain                            |       1.6 |         0.5 |      0.0 |      9.8 |       1.1 |         0.3 |      0.0 |      7.0 |
| Estonia                          |       0.0 |         0.0 |      0.0 |      0.2 |       0.0 |         0.0 |      0.0 |      0.1 |
| Ethiopia                         |      24.5 |         3.7 |      0.1 |    138.7 |      29.5 |         4.4 |      0.1 |    167.2 |
| Finland                          |       0.4 |         0.1 |      0.0 |      1.6 |       0.3 |         0.0 |      0.0 |      1.2 |
| Fiji                             |       0.0 |         0.0 |      0.0 |      0.1 |       0.0 |         0.0 |      0.0 |      0.1 |
| France                           |       3.5 |         0.8 |      0.0 |     19.5 |       3.1 |         0.7 |      0.0 |     16.9 |
| Micronesia                       |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Gabon                            |       0.2 |         0.1 |      0.0 |      1.1 |       0.3 |         0.1 |      0.0 |      1.3 |
| United Kingdom                   |       3.6 |         0.8 |      0.0 |     20.2 |       3.0 |         0.7 |      0.0 |     17.1 |
| Georgia                          |       0.3 |         0.1 |      0.0 |      1.6 |       0.2 |         0.0 |      0.0 |      1.3 |
| Ghana                            |      22.2 |         0.9 |      0.0 |     31.9 |      23.4 |         0.9 |      0.0 |     33.5 |
| Guinea                           |       1.4 |         0.5 |      0.0 |      7.4 |       1.7 |         0.5 |      0.0 |      8.6 |
| Gambia                           |       0.3 |         0.1 |      0.0 |      1.3 |       0.3 |         0.1 |      0.0 |      1.4 |
| Guinea-Bissau                    |       0.2 |         0.1 |      0.0 |      1.0 |       0.2 |         0.1 |      0.0 |      1.1 |
| Equatorial Guinea                |       0.2 |         0.1 |      0.0 |      0.8 |       0.2 |         0.1 |      0.0 |      1.0 |
| Greece                           |       0.5 |         0.1 |      0.0 |      2.7 |       0.4 |         0.1 |      0.0 |      2.0 |
| Grenada                          |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Guatemala                        |       1.0 |         0.8 |      0.2 |      3.3 |       0.9 |         0.7 |      0.2 |      3.2 |
| Guyana                           |       0.0 |         0.0 |      0.0 |      0.1 |       0.0 |         0.0 |      0.0 |      0.1 |
| Honduras                         |       0.7 |         0.3 |      0.0 |      3.0 |       0.7 |         0.3 |      0.0 |      3.2 |
| Croatia                          |       0.2 |         0.0 |      0.0 |      1.2 |       0.2 |         0.0 |      0.0 |      1.0 |
| Haiti                            |       0.8 |         0.4 |      0.0 |      3.7 |       0.8 |         0.4 |      0.0 |      3.6 |
| Hungary                          |       0.5 |         0.1 |      0.0 |      2.3 |       0.5 |         0.1 |      0.0 |      2.4 |
| Indonesia                        |      11.2 |         4.2 |      0.2 |     62.5 |      10.2 |         3.9 |      0.2 |     57.1 |
| India                            |      22.4 |        13.1 |      1.4 |     98.8 |      19.6 |        11.5 |      1.2 |     86.5 |
| Ireland                          |       0.4 |         0.1 |      0.0 |      2.1 |       0.3 |         0.1 |      0.0 |      1.6 |
| Iran                             |       4.0 |         1.3 |      0.1 |     22.7 |       3.7 |         1.2 |      0.1 |     20.8 |
| Iraq                             |       4.7 |         1.1 |      0.0 |     25.5 |       5.0 |         1.1 |      0.0 |     26.9 |
| Iceland                          |       0.0 |         0.0 |      0.0 |      0.1 |       0.0 |         0.0 |      0.0 |      0.1 |
| Israel                           |       0.3 |         0.2 |      0.0 |      1.5 |       0.3 |         0.2 |      0.0 |      1.6 |
| Italy                            |       2.1 |         0.6 |      0.0 |     11.9 |       1.6 |         0.4 |      0.0 |      8.6 |
| Jamaica                          |       0.3 |         0.1 |      0.0 |      1.4 |       0.2 |         0.1 |      0.0 |      1.1 |
| Jordan                           |       0.9 |         0.2 |      0.0 |      4.9 |       1.0 |         0.2 |      0.0 |      5.6 |
| Japan                            |       6.6 |         3.0 |      0.3 |     36.6 |       5.2 |         2.3 |      0.2 |     28.8 |
| Kazakhstan                       |       2.1 |         0.4 |      0.0 |     10.7 |       2.4 |         0.4 |      0.0 |     12.4 |
| Kenya                            |      11.3 |         1.7 |      0.0 |     60.4 |      10.9 |         1.6 |      0.0 |     58.4 |
| Kyrgyzstan                       |       0.4 |         0.2 |      0.0 |      1.9 |       0.4 |         0.2 |      0.0 |      2.0 |
| Cambodia                         |       2.8 |         0.6 |      0.0 |     14.7 |       2.9 |         0.6 |      0.0 |     15.2 |
| Kiribati                         |       0.0 |         0.0 |      0.0 |      0.1 |       0.0 |         0.0 |      0.0 |      0.1 |
| Saint Kitts and Nevis            |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Korea                            |       5.5 |         1.6 |      0.1 |     30.9 |       3.3 |         1.0 |      0.0 |     18.5 |
| Kuwait                           |       0.3 |         0.1 |      0.0 |      1.9 |       0.3 |         0.1 |      0.0 |      1.8 |
| Lao People’s Dem. Republic       |       0.6 |         0.3 |      0.0 |      2.7 |       0.6 |         0.3 |      0.0 |      2.6 |
| Lebanon                          |       0.4 |         0.1 |      0.0 |      2.2 |       0.4 |         0.1 |      0.0 |      2.3 |
| Liberia                          |       0.6 |         0.2 |      0.0 |      2.7 |       0.6 |         0.2 |      0.0 |      2.8 |
| Libya                            |       0.7 |         0.2 |      0.0 |      3.5 |       0.6 |         0.1 |      0.0 |      3.0 |
| Saint Lucia                      |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Sri Lanka                        |       1.1 |         0.3 |      0.0 |      5.9 |       1.0 |         0.2 |      0.0 |      5.5 |
| Lesotho                          |       0.2 |         0.1 |      0.0 |      1.1 |       0.2 |         0.1 |      0.0 |      1.1 |
| Lithuania                        |       0.1 |         0.0 |      0.0 |      0.8 |       0.1 |         0.0 |      0.0 |      0.6 |
| Luxembourg                       |       0.0 |         0.0 |      0.0 |      0.1 |       0.0 |         0.0 |      0.0 |      0.2 |
| Latvia                           |       0.0 |         0.0 |      0.0 |      0.2 |       0.0 |         0.0 |      0.0 |      0.2 |
| Morocco                          |       3.3 |         0.7 |      0.0 |     18.6 |       3.0 |         0.7 |      0.0 |     17.1 |
| Monaco                           |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Republic of Moldova              |       0.1 |         0.1 |      0.0 |      0.5 |       0.1 |         0.0 |      0.0 |      0.4 |
| Madagascar                       |       2.9 |         0.9 |      0.0 |     13.8 |       3.5 |         1.1 |      0.0 |     16.5 |
| Maldives                         |       0.0 |         0.0 |      0.0 |      0.2 |       0.0 |         0.0 |      0.0 |      0.1 |
| Mexico                           |      15.3 |         7.7 |      1.5 |     80.1 |      13.9 |         7.0 |      1.4 |     72.9 |
| Marshall Islands                 |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| North Macedonia                  |       0.1 |         0.0 |      0.0 |      0.3 |       0.0 |         0.0 |      0.0 |      0.2 |
| Mali                             |       2.7 |         0.9 |      0.0 |     12.8 |       3.2 |         1.0 |      0.0 |     15.3 |
| Malta                            |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Myanmar                          |       1.2 |         0.7 |      0.1 |      4.8 |       1.1 |         0.7 |      0.1 |      4.7 |
| Montenegro                       |       0.0 |         0.0 |      0.0 |      0.1 |       0.0 |         0.0 |      0.0 |      0.1 |
| Mongolia                         |       0.4 |         0.1 |      0.0 |      3.0 |       0.5 |         0.1 |      0.0 |      3.4 |
| Mozambique                       |       8.2 |         1.1 |      0.0 |     43.8 |       9.9 |         1.3 |      0.0 |     52.9 |
| Mauritania                       |       0.5 |         0.1 |      0.0 |      2.4 |       0.6 |         0.2 |      0.0 |      3.0 |
| Mauritius                        |       0.1 |         0.0 |      0.0 |      0.3 |       0.1 |         0.0 |      0.0 |      0.3 |
| Malawi                           |       2.2 |         0.7 |      0.0 |     10.4 |       2.3 |         0.7 |      0.0 |     11.0 |
| Malaysia                         |       2.3 |         0.5 |      0.0 |     10.3 |       2.1 |         0.5 |      0.0 |      9.6 |
| Namibia                          |       0.3 |         0.1 |      0.0 |      1.3 |       0.3 |         0.1 |      0.0 |      1.4 |
| Niger                            |       3.0 |         0.9 |      0.0 |     14.0 |       3.7 |         1.2 |      0.0 |     17.5 |
| Nigeria                          |      55.7 |         7.5 |      0.1 |    281.2 |      57.6 |         7.7 |      0.2 |    290.3 |
| Nicaragua                        |       0.4 |         0.2 |      0.0 |      1.9 |       0.4 |         0.2 |      0.0 |      1.8 |
| Niue                             |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Netherlands                      |       1.2 |         0.2 |      0.0 |      4.7 |       1.1 |         0.2 |      0.0 |      4.4 |
| Norway                           |       0.3 |         0.1 |      0.0 |      1.6 |       0.3 |         0.1 |      0.0 |      1.4 |
| Nepal                            |       0.8 |         0.5 |      0.0 |      3.1 |       0.7 |         0.4 |      0.0 |      3.0 |
| Nauru                            |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| New Zealand                      |       0.9 |         0.2 |      0.0 |      5.5 |       0.8 |         0.2 |      0.0 |      5.0 |
| Oman                             |       0.2 |         0.1 |      0.0 |      0.9 |       0.3 |         0.1 |      0.0 |      1.2 |
| Pakistan                         |      12.2 |         6.3 |      0.6 |     61.1 |      12.3 |         6.3 |      0.6 |     61.6 |
| Panama                           |       0.1 |         0.1 |      0.0 |      0.5 |       0.1 |         0.1 |      0.0 |      0.5 |
| Peru                             |       2.8 |         1.1 |      0.1 |     15.1 |       2.7 |         1.0 |      0.1 |     14.2 |
| Philippines                      |      15.0 |         4.0 |      0.2 |     83.5 |      11.2 |         3.0 |      0.1 |     62.2 |
| Palau                            |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Papua New Guinea                 |       0.8 |         0.4 |      0.0 |      3.6 |       0.9 |         0.4 |      0.0 |      4.0 |
| Poland                           |       2.0 |         0.4 |      0.0 |     10.8 |       1.7 |         0.4 |      0.0 |      9.2 |
| Korea                            |       0.4 |         0.2 |      0.0 |      1.7 |       0.4 |         0.3 |      0.0 |      1.8 |
| Portugal                         |       0.4 |         0.1 |      0.0 |      2.7 |       0.4 |         0.1 |      0.0 |      2.2 |
| Paraguay                         |       0.3 |         0.2 |      0.1 |      1.1 |       0.3 |         0.3 |      0.1 |      1.2 |
| Qatar                            |       0.1 |         0.0 |      0.0 |      0.3 |       0.1 |         0.0 |      0.0 |      0.4 |
| Romania                          |       1.1 |         0.2 |      0.0 |      6.1 |       1.0 |         0.2 |      0.0 |      5.3 |
| Russian Federation               |       5.8 |         1.6 |      0.0 |     35.7 |       4.7 |         1.2 |      0.0 |     28.7 |
| Rwanda                           |       1.3 |         0.4 |      0.0 |      6.3 |       1.4 |         0.5 |      0.0 |      6.8 |
| Saudi Arabia                     |       4.6 |         0.6 |      0.0 |     15.2 |       4.7 |         0.6 |      0.0 |     15.4 |
| Sudan                            |       6.6 |         1.4 |      0.0 |     44.5 |       8.0 |         1.7 |      0.0 |     53.6 |
| Senegal                          |       1.7 |         0.5 |      0.0 |      8.6 |       1.8 |         0.6 |      0.0 |      9.2 |
| Singapore                        |       3.0 |         0.9 |      0.1 |     19.1 |       3.4 |         1.0 |      0.1 |     21.1 |
| Solomon Islands                  |       0.1 |         0.0 |      0.0 |      0.3 |       0.1 |         0.0 |      0.0 |      0.3 |
| Sierra Leone                     |       0.9 |         0.3 |      0.0 |      4.2 |       0.9 |         0.3 |      0.0 |      4.4 |
| El Salvador                      |       0.3 |         0.2 |      0.1 |      1.0 |       0.3 |         0.2 |      0.0 |      0.8 |
| San Marino                       |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Somalia                          |       1.5 |         0.7 |      0.0 |      7.6 |       1.9 |         0.8 |      0.0 |      9.6 |
| Serbia                           |       0.4 |         0.1 |      0.0 |      1.9 |       0.3 |         0.1 |      0.0 |      1.7 |
| South Sudan                      |       1.4 |         0.4 |      0.0 |      6.6 |       1.1 |         0.4 |      0.0 |      5.3 |
| Sao Tome and Principe            |       0.0 |         0.0 |      0.0 |      0.1 |       0.0 |         0.0 |      0.0 |      0.1 |
| Suriname                         |       0.1 |         0.0 |      0.0 |      0.4 |       0.1 |         0.0 |      0.0 |      0.4 |
| Slovakia                         |       0.3 |         0.1 |      0.0 |      1.6 |       0.3 |         0.1 |      0.0 |      1.5 |
| Slovenia                         |       0.1 |         0.0 |      0.0 |      0.7 |       0.1 |         0.0 |      0.0 |      0.5 |
| Sweden                           |       0.5 |         0.1 |      0.0 |      3.0 |       0.5 |         0.1 |      0.0 |      3.0 |
| Eswatini                         |       0.1 |         0.0 |      0.0 |      0.6 |       0.1 |         0.0 |      0.0 |      0.6 |
| Seychelles                       |       0.0 |         0.0 |      0.0 |      0.1 |       0.0 |         0.0 |      0.0 |      0.1 |
| Syrian Arab Republic             |       1.7 |         0.7 |      0.0 |      8.4 |       1.1 |         0.5 |      0.0 |      5.6 |
| Chad                             |       2.2 |         0.7 |      0.0 |     10.3 |       2.8 |         0.9 |      0.0 |     13.0 |
| Togo                             |       0.9 |         0.3 |      0.0 |      4.4 |       1.0 |         0.3 |      0.0 |      4.9 |
| Thailand                         |       3.4 |         0.7 |      0.0 |     17.9 |       2.6 |         0.6 |      0.0 |     13.7 |
| Tajikistan                       |       0.6 |         0.3 |      0.0 |      3.0 |       0.7 |         0.3 |      0.0 |      3.4 |
| Turkmenistan                     |       0.3 |         0.1 |      0.0 |      1.4 |       0.3 |         0.2 |      0.0 |      1.6 |
| Timor-Leste                      |       0.0 |         0.0 |      0.0 |      0.2 |       0.0 |         0.0 |      0.0 |      0.2 |
| Tonga                            |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Trinidad and Tobago              |       0.0 |         0.0 |      0.0 |      0.1 |       0.0 |         0.0 |      0.0 |      0.1 |
| Tunisia                          |       0.8 |         0.2 |      0.0 |      4.1 |       0.8 |         0.2 |      0.0 |      4.0 |
| Turkiye                          |       6.0 |         1.2 |      0.0 |     32.7 |       5.7 |         1.1 |      0.0 |     31.0 |
| Tuvalu                           |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| United Republic of Tanzania      |      13.4 |         1.9 |      0.0 |     73.5 |      17.1 |         2.4 |      0.1 |     94.3 |
| Uganda                           |       5.1 |         1.6 |      0.1 |     24.0 |       6.0 |         1.9 |      0.1 |     28.4 |
| Ukraine                          |       2.7 |         0.6 |      0.0 |     14.5 |       1.8 |         0.4 |      0.0 |      9.7 |
| Uruguay                          |       0.1 |         0.1 |      0.0 |      0.3 |       0.1 |         0.0 |      0.0 |      0.2 |
| United States of America         |       6.7 |         4.6 |      0.8 |     25.6 |       6.2 |         4.2 |      0.7 |     23.4 |
| Uzbekistan                       |       1.5 |         0.7 |      0.0 |      7.9 |       2.1 |         0.9 |      0.1 |     10.7 |
| Saint Vincent and the Grenadines |       0.0 |         0.0 |      0.0 |      0.0 |       0.0 |         0.0 |      0.0 |      0.0 |
| Venezuela                        |       3.2 |         0.8 |      0.0 |     18.7 |       2.5 |         0.6 |      0.0 |     14.2 |
| Viet Nam                         |      10.9 |         2.4 |      0.1 |     59.5 |      10.7 |         2.3 |      0.1 |     58.4 |
| Vanuatu                          |       0.0 |         0.0 |      0.0 |      0.1 |       0.0 |         0.0 |      0.0 |      0.1 |
| Samoa                            |       0.0 |         0.0 |      0.0 |      0.1 |       0.0 |         0.0 |      0.0 |      0.1 |
| Yemen                            |       2.4 |         1.1 |      0.1 |     12.3 |       3.3 |         1.4 |      0.1 |     16.8 |
| South Africa                     |       8.5 |         1.3 |      0.0 |     50.8 |       8.7 |         1.3 |      0.0 |     51.7 |
| Zambia                           |       2.1 |         0.7 |      0.0 |     10.6 |       2.3 |         0.7 |      0.0 |     12.0 |
| Zimbabwe                         |       3.3 |         0.5 |      0.0 |     20.2 |       3.2 |         0.5 |      0.0 |     19.2 |

Estimated Profound cases by country, 2010 vs 2020

# Session info

``` r
sessioninfo::session_info()
```

    ## Warning in system2("quarto", "-V", stdout = TRUE, env = paste0("TMPDIR=", : running command '"quarto"
    ## TMPDIR=C:/Users/LoVa3397/AppData/Local/Temp/RtmpUrHJ3m/file20a836d42de -V' had status 1

    ## ─ Session info ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ##  setting  value
    ##  version  R version 4.5.0 (2025-04-11 ucrt)
    ##  os       Windows 10 x64 (build 19045)
    ##  system   x86_64, mingw32
    ##  ui       RStudio
    ##  language (EN)
    ##  collate  English_Belgium.utf8
    ##  ctype    English_Belgium.utf8
    ##  tz       Europe/Brussels
    ##  date     2025-09-15
    ##  rstudio  2024.04.2+764 Chocolate Cosmos (desktop)
    ##  pandoc   3.1.11 @ C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools/ (via rmarkdown)
    ##  quarto   ERROR: Unknown command "TMPDIR=C:/Users/LoVa3397/AppData/Local/Temp/RtmpUrHJ3m/file20a836d42de". Did you mean command "create"? @ C:\\PROGRA~1\\RStudio\\RESOUR~1\\app\\bin\\quarto\\bin\\quarto.exe
    ## 
    ## ─ Packages ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
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
    ##    scales         * 1.3.0      2023-11-28 [1] CRAN (R 4.5.0)
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
    ## ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

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
