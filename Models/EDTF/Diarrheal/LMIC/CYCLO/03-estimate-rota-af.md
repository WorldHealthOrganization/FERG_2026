Global proportion of Cyclospora (Diarrheal) • ASYMPT+ROTA • Estimate
attributable fraction
================
fbbu6966
2025-11-27

- [Settings](#settings)
- [Model fit](#model-fit)
- [Predict all](#predict-all)
- [Summarize predictions: countries](#summarize-predictions-countries)
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

## global options ----
knitr::opts_chunk$set(fig.width = 10)
do.call(file.remove, list(list.files(params$PlotDir, full.names = TRUE)))
```

    ## logical(0)

``` r
Date <- format(Sys.Date(), "%Y%m%d")
```

# Model fit

``` r
es <- readRDS(paste0(params$Dir, "/es.rds"))
es <- subset(es, as.integer(FLAG) == 1)
```

``` r
png(paste0(params$PlotDir, "/imputation_map.png"), width=960, height=480)
plot_world_imputation(es)
```

NULL

``` r
dev.off()
```

png 2

``` r
setwd(params$Dir)
image <- paste0("03-estimate-rota-af_files/figure-gfm/imputation_map.png")
cat("![](", image, ")")
```

![](03-estimate-rota-af_files/figure-gfm/imputation_map.png)

``` r
fit_brms_reg_s <- readRDS(paste0(params$Dir, "/fit_brms_reg_s.rds"))
print(fit_brms_reg_s)
```

    ## Warning: There were 29 divergent transitions after warmup. Increasing adapt_delta above 0.95 may help. See
    ## http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup

    ##  Family: gaussian 
    ##   Links: mu = identity; sigma = identity 
    ## Formula: yi | se(sei) ~ 1 + YEAR + SYNDROMTYPE + AGE + REFERENCE + COVERAGE + (1 | REG2) + (1 | REG2:SUB2) + (1 | REG2:SUB2:COUNTRY) + (1 | REG2:SUB2:COUNTRY:ID) + (1 | REG2:SUB2:COUNTRY:ID:DTP_ID) 
    ##    Data: subset(es, as.integer(FLAG) == 1) (Number of observations: 168) 
    ##   Draws: 5 chains, each with iter = 5000; warmup = 3000; thin = 1;
    ##          total post-warmup draws = 10000
    ## 
    ## Multilevel Hyperparameters:
    ## ~REG2 (Number of levels: 6) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.62      0.47     0.03     1.74 1.00     4025     5273
    ## 
    ## ~REG2:SUB2 (Number of levels: 8) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.47      0.37     0.02     1.39 1.00     3640     4510
    ## 
    ## ~REG2:SUB2:COUNTRY (Number of levels: 19) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.83      0.31     0.23     1.49 1.00     2046     1860
    ## 
    ## ~REG2:SUB2:COUNTRY:ID (Number of levels: 30) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.33      0.23     0.01     0.86 1.00     1562     3006
    ## 
    ## ~REG2:SUB2:COUNTRY:ID:DTP_ID (Number of levels: 168) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.64      0.09     0.47     0.83 1.00     4661     6451
    ## 
    ## Regression Coefficients:
    ##                       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept               101.42    117.08  -136.00   326.81 1.00     2717     4345
    ## YEAR                     -0.05      0.06    -0.16     0.07 1.00     2720     4344
    ## SYNDROMTYPEInpatient      1.25      0.29     0.69     1.82 1.00     6440     7426
    ## SYNDROMTYPEOutpatient     0.39      0.23    -0.05     0.84 1.00     5989     6856
    ## AGEAgebelow5             -0.62      0.46    -1.53     0.28 1.00     7153     6568
    ## AGEMixedages              0.15      0.86    -1.48     1.84 1.00     6388     7394
    ## REFERENCEOther           -1.76      0.91    -3.42     0.16 1.00     2374     3666
    ## COVERAGE                 -0.01      0.01    -0.02     0.01 1.00     6000     6016
    ## 
    ## Further Distributional Parameters:
    ##       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sigma     0.00      0.00     0.00     0.00   NA       NA       NA
    ## 
    ## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

``` r
zero_cases <-
read_xlsx("endemic_countries.xlsx") %>%
  select(REG2, SUB2, ISO3, Country, edtf_diarrheal) %>% 
  rename(COUNTRY=ISO3, COUNTRY_LABEL = Country) %>%
  mutate(DISEASEFREE = edtf_diarrheal)

kable(
  caption = "Countries for which no estimates were generated",
  row.names = FALSE,
  subset(zero_cases, edtf_diarrheal==0)[, 4])
```

| COUNTRY_LABEL        |
|:---------------------|
| Andorra              |
| United Arab Emirates |
| Antigua and Barbuda  |
| Australia            |
| Austria              |
| Belgium              |
| Bahrain              |
| Bahamas, The         |
| Barbados             |
| Brunei Darussalam    |
| Canada               |
| Switzerland          |
| Chile                |
| Cook Islands         |
| Cyprus               |
| Czech Republic       |
| Germany              |
| Denmark              |
| Spain                |
| Estonia              |
| Finland              |
| France               |
| United Kingdom       |
| Greece               |
| Guyana               |
| Croatia              |
| Hungary              |
| Ireland              |
| Iceland              |
| Israel               |
| Italy                |
| Japan                |
| St. Kitts and Nevis  |
| Korea, Rep.          |
| Kuwait               |
| Lithuania            |
| Luxembourg           |
| Latvia               |
| Monaco               |
| Malta                |
| Niue                 |
| Netherlands          |
| Norway               |
| Nauru                |
| New Zealand          |
| Oman                 |
| Panama               |
| Poland               |
| Portugal             |
| Qatar                |
| Romania              |
| Saudi Arabia         |
| Singapore            |
| San Marino           |
| Slovak Republic      |
| Slovenia             |
| Sweden               |
| Trinidad and Tobago  |
| Uruguay              |
| United States        |

Countries for which no estimates were generated

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

VacCov <- readRDS(paste0(params$Dir, "/vaccine_coverage.rds"))
```

# Predict all

``` r
## set up dataframe
# sim_all <-
#   data.frame(
#     sei = 0,
#     REG2 = FERG2:::countries$REG2,
#     SUB2 = FERG2:::countries$SUB2,
#     COUNTRY = FERG2:::countries$ISO3,
#     YEAR = rep(2000:2021, each = nrow(FERG2:::countries)),
#     SYNDROMTYPE = rep(1:3, each=8536),
#     AGE = rep(c(1,2,1,2,1,2), each=4268)) %>%
#   distinct()
sim_all <-
  expand.grid(
    sei = 0,
    COUNTRY = FERG2:::countries$ISO3,
    YEAR = 2000:2021,
    SYNDROMTYPE = 1:3,
    AGE = 1:2)
sim_all$REG2 <-
  FERG2:::countries$REG2[match(sim_all$COUNTRY, FERG2:::countries$ISO3)]
sim_all$SUB2 <-
  FERG2:::countries$SUB2[match(sim_all$COUNTRY, FERG2:::countries$ISO3)]

sim_all <- sim_all %>%
  left_join(zero_cases) %>%
  select(sei, REG2, SUB2, COUNTRY, YEAR, SYNDROMTYPE, AGE, ESTIMATES)
```

    ## Joining with `by = join_by(COUNTRY, REG2, SUB2)`

``` r
sim_all$SYNDROMTYPE <-
  factor(sim_all$SYNDROMTYPE,
         levels = 1:3,
         labels = c("Asymptomatic", "Inpatient", "Outpatient"))
sim_all$AGE <-
  factor(sim_all$AGE, 
         levels = 1:2, 
         labels = c("Age above or equal 5", "Age below 5"))
sim_all <-
  merge(sim_all, VacCov,
        by.x = c("COUNTRY", "YEAR"), by.y = c("ISO3", "YEAR"), all.x = TRUE)

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

##

draws_fit <- as_draws_df(fit_brms_reg_s) %>% as.data.frame()
# fit_all <- data.frame(1:10000)

# str(sim_all)

## compile fixed effect
b0_ad_as <- draws_fit$b_Intercept
b0_ad_si <- b0_ad_as + draws_fit$b_SYNDROMTYPEInpatient
b0_ad_so <- b0_ad_as + draws_fit$b_SYNDROMTYPEOutpatient
b0_ch_as <- b0_ad_as + draws_fit$b_AGEAgebelow5
b0_ch_si <- b0_ch_as + draws_fit$b_SYNDROMTYPEInpatient
b0_ch_so <- b0_ch_as + draws_fit$b_SYNDROMTYPEOutpatient
b0 <-
  function(age, type) {
    b0 <-
      paste(
        "b0",
        c("ad", "ch")[age],
        c("as", "si", "so")[type],
        sep = "_")
    get(b0)
  }

bY <- draws_fit$b_YEAR
bCov <- draws_fit$b_COVERAGE

## compile random effects
rCOUNTRY <- data.frame(ranef(fit_brms_reg_s, summary = FALSE)$`REG2:SUB2:COUNTRY`)
names(rCOUNTRY) <- gsub(".*_(.*)\\..*", "\\1", names(rCOUNTRY))
rSUB <- data.frame(ranef(fit_brms_reg_s, summary = FALSE)$`REG2:SUB2`)
names(rSUB) <- gsub(".*_(.*)\\..*", "\\1", names(rSUB))
rREG <- data.frame(ranef(fit_brms_reg_s, summary = FALSE)$`REG2`)
names(rREG) <- gsub("\\..*", "", names(rREG))

## functions to predict estimate by year and location
est_cnt <-
  function(b0, .year, .coverage, .country, .sub, .reg) {
    est <-
      b0 + bY * .year + bCov * .coverage +
      rCOUNTRY[[.country]] +
      rSUB[[.sub]] +
      rREG[[.reg]]
    list(est)
  }

est_sub <-
  function(b0, .year, .coverage, .sub, .reg) {
    est <-
      b0 + bY * .year + bCov * .coverage +
      rSUB[[.sub]] +
      rREG[[.reg]]
    list(est)
  }

est_reg <-
  function(b0, .year, .coverage, .reg) {
    est <-
      b0 + bY * .year + bCov * .coverage +
      rREG[[.reg]]
    list(est)
  }

est_glb <-
  function(b0, .year, .coverage) {
    est <-
      b0 + bY * .year + bCov * .coverage
    list(est)
  }

## prepare 'SIM' as list
sim_all$SIM <- vector("list", nrow(sim_all))

## data present for country
id <- which(as.integer(sim_all$ESTIMATES) == 1); length(id)
```

    ## [1] 2508

``` r
for (x in id) {
  sim_all$SIM[x] <-
    est_cnt(
      b0(as.numeric(sim_all[x, "AGE"]), as.numeric(sim_all[x, "SYNDROMTYPE"])),
      sim_all[x, "YEAR"],
      sim_all[x, "COVERAGE"],
      sim_all[x, "COUNTRY"],
      sim_all[x, "SUB2"],
      sim_all[x, "REG2"])
}

## disease-free country
id <- which(as.integer(sim_all$ESTIMATES) == 2); length(id)
```

    ## [1] 7920

``` r
for (x in id) {
  sim_all$SIM[x] <- -Inf
}

## data in subregion
id <- which(as.integer(sim_all$ESTIMATES) == 3); length(id)
```

    ## [1] 12012

``` r
for (x in id) {
  sim_all$SIM[x] <-
    est_sub(
      b0(as.numeric(sim_all[x, "AGE"]), as.numeric(sim_all[x, "SYNDROMTYPE"])),
      sim_all[x, "YEAR"],
      sim_all[x, "COVERAGE"],
      sim_all[x, "SUB2"],
      sim_all[x, "REG2"])
}

## data in region
id <- which(as.integer(sim_all$ESTIMATES) == 4); length(id)
```

    ## [1] 3168

``` r
for (x in id) {
  sim_all$SIM[x] <-
    est_reg(
      b0(as.numeric(sim_all[x, "AGE"]), as.numeric(sim_all[x, "SYNDROMTYPE"])),
      sim_all[x, "YEAR"],
      sim_all[x, "COVERAGE"],
      sim_all[x, "REG2"])
}

## data in world
id <- which(as.integer(sim_all$ESTIMATES) == 5); length(id)
```

    ## [1] 0

``` r
for (x in id) {
  sim_all$SIM[x] <-
    est_glb(
      b0(as.numeric(sim_all[x, "AGE"]), as.numeric(sim_all[x, "SYNDROMTYPE"])),
      sim_all[x, "YEAR"],
      sim_all[x, "COVERAGE"])
}

# str(sim_all)

##

# draws_fit <- as_draws_df(fit_brms_reg_s) %>% as.data.frame()
# fit_all <- data.frame(1:10000)
# for (x in 1:nrow(sim_all)){
#   # Fixed effects
#   if (sim_all[x, "AGE"] == "Age above or equal 5" & sim_all[x, "SYNDROMTYPE"] == "Inpatient"){                                        # Inpatient and age above 5
#     draws_fit$beta <-  draws_fit$b_Intercept
#   } else if (sim_all[x, "AGE"] ==  "Age above or equal 5" & sim_all[x, "SYNDROMTYPE"] == "Outpatient"){                               # Outpatient and age above 5
#     draws_fit$beta <-  draws_fit$b_Intercept + draws_fit$b_SYNDROMTYPEOutpatient
#   } else if (sim_all[x, "AGE"] ==  "Age below 5" & sim_all[x, "SYNDROMTYPE"] == "Inpatient"){                                         # Inpatient and age below 5
#     draws_fit$beta <-  draws_fit$b_Intercept + draws_fit$b_AGEAgebelow5
#   } else if (sim_all[x, "AGE"] ==  "Age below 5" & sim_all[x, "SYNDROMTYPE"] == "Outpatient"){                                        # Outpatient and age below 5
#     draws_fit$beta <-  draws_fit$b_Intercept + draws_fit$b_SYNDROMTYPEOutpatient + draws_fit$b_AGEAgebelow5
#   }
#   # Data present for country
#   if (as.integer(sim_all[x, "ESTIMATES"]) == 1){
#     fit_all[[paste0("V",x)]] <- draws_fit$beta +                                                                                      # Global intercept
#       sim_all[x, "YEAR"] * draws_fit$b_YEAR +                                                                                         # Year component
#       sim_all[x, "COVERAGE"] * draws_fit$b_COVERAGE +                                                                                 # Coverage component
#       draws_fit[[paste0("r_REG2[",sim_all[x,"REG2"],",Intercept]")]] +                                                                # Regional component
#       draws_fit[[paste0("r_REG2:SUB2[",sim_all[x,"REG2"],"_",sim_all[x,"SUB2"],",Intercept]")]] +                                     # Sub regional component
#       draws_fit[[paste0("r_REG2:SUB2:COUNTRY[",sim_all[x,"REG2"],"_",sim_all[x,"SUB2"],"_",sim_all[x,"COUNTRY"],",Intercept]")]]      # Country component
#   } else if (as.integer(sim_all[x, "ESTIMATES"]) == 2) {
#     # Disease-free country
#     fit_all[[paste0("V",x)]] <- 0
#   } else if (as.integer(sim_all[x, "ESTIMATES"]) == 3){
#     # Data not present for country, but present in subregion
#     fit_all[[paste0("V",x)]] <- draws_fit$beta +                                                                                      # Global intercept
#       sim_all[x, "YEAR"] * draws_fit$b_YEAR +                                                                                         # Year component
#       sim_all[x, "COVERAGE"] * draws_fit$b_COVERAGE +                                                                                 # Coverage component
#       draws_fit[[paste0("r_REG2[",sim_all[x,"REG2"],",Intercept]")]] +                                                                # Regional component
#       draws_fit[[paste0("r_REG2:SUB2[",sim_all[x,"REG2"],"_",sim_all[x,"SUB2"],",Intercept]")]]                                       # Sub regional component
#   } else if (as.integer(sim_all[x, "ESTIMATES"]) == 4){
#     # Data not present for country, but present in region
#     fit_all[[paste0("V",x)]] <- draws_fit$beta +                                                                                      # Global intercept
#       sim_all[x, "YEAR"] * draws_fit$b_YEAR +                                                                                         # Year component
#       sim_all[x, "COVERAGE"] * draws_fit$b_COVERAGE +                                                                                 # Coverage component
#       draws_fit[[paste0("r_REG2[",sim_all[x,"REG2"],",Intercept]")]]                                                                  # Regional component
#   } else if (as.integer(sim_all[x, "ESTIMATES"]) == 5){
#     # Data not present for country
#     fit_all[[paste0("V",x)]] <- draws_fit$beta +                                                                                      # Global intercept
#       sim_all[x, "YEAR"] * draws_fit$b_YEAR +
#       sim_all[x, "COVERAGE"] * draws_fit$b_COVERAGE                                                                                   # Coverage component                                                                              # Year component
#   } 
# }
# 
# fit_all <- fit_all %>% select(-c(X1.10000))

## calculate proportions
sim_all <- sim_all %>% left_join(zero_cases)
```

    ## Joining with `by = join_by(COUNTRY, REG2, SUB2, ESTIMATES)`

``` r
sim_all$PROP <- lapply(sim_all$SIM, expit)
sim_all$PROP <- mapply(`*`, sim_all$PROP, sim_all$edtf_diarrheal)

sim_all <- sort_by(sim_all, ~COUNTRY+YEAR+AGE+SYNDROMTYPE)

sim_all_as <- subset(sim_all, SYNDROMTYPE == "Asymptomatic")
sim_all_si <- subset(sim_all, SYNDROMTYPE == "Inpatient")
sim_all_so <- subset(sim_all, SYNDROMTYPE == "Outpatient")

# sim_all_si$PROP <- mapply(`-`, sim_all_si$PROP, sim_all_as$PROP)
# sim_all_so$PROP <- mapply(`-`, sim_all_so$PROP, sim_all_as$PROP)

# odds
odds_as <- lapply(sim_all_as$PROP, function(p) p/(1-p))
odds_si <- lapply(sim_all_si$PROP, function(p) p/(1-p))
odds_so <- lapply(sim_all_so$PROP, function(p) p/(1-p))

# odds ratios
or_si <- mapply(function(x,y) x/y, odds_si, odds_as)
or_so <- mapply(function(x,y) x/y, odds_so, odds_as)

# attributable fraction among exposed
af_si <- lapply(or_si, function(or) 1-1/or)
af_so <- lapply(or_so, function(or) 1-1/or)

# attributable fraction
sim_all_si$PROP <- mapply(`*`, sim_all_si$PROP, af_si)
sim_all_so$PROP <- mapply(`*`, sim_all_so$PROP, af_so)

# replace NaN by 0
sim_all_si$PROP <- lapply(sim_all_si$PROP, function(x) {x[is.nan(x)] <- 0; x})
sim_all_so$PROP <- lapply(sim_all_so$PROP, function(x) {x[is.nan(x)] <- 0; x})

sim_all <- rbind(sim_all_si, sim_all_so)

# saveRDS(subset(sim_all %>% filter(SYNDROMTYPE == "Asymptomatic" & AGE == "Age above or equal 5")), 
#         file = paste0(params$Dir, "/sim_all_Asymptomatic_Older5_",Date,".rds"))
# saveRDS(subset(sim_all %>% filter(SYNDROMTYPE == "Asymptomatic" & AGE == "Age below 5")), 
#         file = paste0(params$Dir, "/sim_all_Asymptomatic_Younger5_",Date,".rds"))
saveRDS(subset(sim_all %>% filter(SYNDROMTYPE == "Inpatient" & AGE == "Age above or equal 5")),
        file = paste0(params$Dir, "/sim_all_Inpatient_Older5_",Date,".rds"))
saveRDS(subset(sim_all %>% filter(SYNDROMTYPE == "Inpatient" & AGE == "Age below 5")),
        file = paste0(params$Dir, "/sim_all_Inpatient_Younger5_",Date,".rds"))
saveRDS(subset(sim_all %>% filter(SYNDROMTYPE == "Outpatient" & AGE == "Age above or equal 5")),
        file = paste0(params$Dir, "/sim_all_Outpatient_Older5_",Date,".rds"))
saveRDS(subset(sim_all %>% filter(SYNDROMTYPE == "Outpatient" & AGE == "Age below 5")),
        file = paste0(params$Dir, "/sim_all_Outpatient_Younger5_",Date,".rds"))

## aggregate over countries
all_cnt_prop <- t(sapply(sim_all$PROP, mean_ci))
all_cnt_prop <- data.frame(all_cnt_prop)
names(all_cnt_prop) <- c("VAL_MEAN", "VAL_LWR", "VAL_UPR")
all_cnt_prop <- cbind(sim_all[c("YEAR","REG2", "SUB2", "COUNTRY", "SYNDROMTYPE","AGE","COUNTRY_LABEL")], all_cnt_prop)
all_cnt_prop$LOCATION <- "Country"
all_cnt_prop$LOCATION_NAME <- all_cnt_prop$COUNTRY_LABEL
all_cnt_prop$COUNTRY_LABEL <- NULL
all_cnt_prop$METRIC <- "Proportion"
str(all_cnt_prop)
```

    ## 'data.frame':    17072 obs. of  12 variables:
    ##  $ YEAR         : int  2000 2000 2001 2001 2002 2002 2003 2003 2004 2004 ...
    ##  $ REG2         : chr  "EMR" "EMR" "EMR" "EMR" ...
    ##  $ SUB2         : chr  "EMRD" "EMRD" "EMRD" "EMRD" ...
    ##  $ COUNTRY      : chr  "AFG" "AFG" "AFG" "AFG" ...
    ##  $ SYNDROMTYPE  : Factor w/ 3 levels "Asymptomatic",..: 2 2 2 2 2 2 2 2 2 2 ...
    ##  $ AGE          : Factor w/ 2 levels "Age above or equal 5",..: 1 2 1 2 1 2 1 2 1 2 ...
    ##  $ VAL_MEAN     : num  0.083 0.0477 0.0779 0.0442 0.0731 ...
    ##  $ VAL_LWR      : num  0.00588 0.00385 0.006 0.004 0.00621 ...
    ##  $ VAL_UPR      : num  0.326 0.199 0.305 0.179 0.282 ...
    ##  $ LOCATION     : chr  "Country" "Country" "Country" "Country" ...
    ##  $ LOCATION_NAME: chr  "Afghanistan" "Afghanistan" "Afghanistan" "Afghanistan" ...
    ##  $ METRIC       : chr  "Proportion" "Proportion" "Proportion" "Proportion" ...

``` r
## compile all
all_est <- all_cnt_prop
str(all_est)
```

    ## 'data.frame':    17072 obs. of  12 variables:
    ##  $ YEAR         : int  2000 2000 2001 2001 2002 2002 2003 2003 2004 2004 ...
    ##  $ REG2         : chr  "EMR" "EMR" "EMR" "EMR" ...
    ##  $ SUB2         : chr  "EMRD" "EMRD" "EMRD" "EMRD" ...
    ##  $ COUNTRY      : chr  "AFG" "AFG" "AFG" "AFG" ...
    ##  $ SYNDROMTYPE  : Factor w/ 3 levels "Asymptomatic",..: 2 2 2 2 2 2 2 2 2 2 ...
    ##  $ AGE          : Factor w/ 2 levels "Age above or equal 5",..: 1 2 1 2 1 2 1 2 1 2 ...
    ##  $ VAL_MEAN     : num  0.083 0.0477 0.0779 0.0442 0.0731 ...
    ##  $ VAL_LWR      : num  0.00588 0.00385 0.006 0.004 0.00621 ...
    ##  $ VAL_UPR      : num  0.326 0.199 0.305 0.179 0.282 ...
    ##  $ LOCATION     : chr  "Country" "Country" "Country" "Country" ...
    ##  $ LOCATION_NAME: chr  "Afghanistan" "Afghanistan" "Afghanistan" "Afghanistan" ...
    ##  $ METRIC       : chr  "Proportion" "Proportion" "Proportion" "Proportion" ...

``` r
saveRDS(all_est, file = paste0(params$Dir, "/all_estimates_",Date,".rds"))
```

# Summarize predictions: countries

``` r
# #+ fig.height=4
# ggplot(all_cnt_prop, aes(x = YEAR, y = VAL_MEAN, group=LOCATION_NAME)) +
#   # geom_line(data = all_cnt_prop, linewidth = 2) +
#   geom_line(aes(col = COUNTRY), linewidth = 0.5) +
# theme_bw() 
```

``` r
# breaks <- plot_world(subset(all_cnt_prop, YEAR == 2020 & AGE == "Age below 5" & SYNDROMTYPE == "Outpatient"), 
# "COUNTRY", "VAL_MEAN", legend.title = "Proportion", diseasefree = zero_cases, text.width = 23)
```

``` r
png(paste0(params$PlotDir, "/r_world1.png"), width=960, height=480)
plot_world(subset(all_cnt_prop, YEAR == 2010 & AGE == "Age above or equal 5" & SYNDROMTYPE == "Inpatient"),
           "COUNTRY", "VAL_MEAN", legend.title = "Proportion", diseasefree = zero_cases, text.width = 23)
```

\[1\] 0.00 0.02 0.04 0.06 0.08 0.10 0.12 0.14

``` r
dev.off()
```

png 2

``` r
setwd(params$Dir)
image <- paste0("03-estimate-rota-af_files/figure-gfm/r_world1.png")
cat("![](",image,")")
```

![](03-estimate-rota-af_files/figure-gfm/r_world1.png)

``` r
png(paste0(params$PlotDir, "/r_world2.png"), width=960, height=480)
plot_world(subset(all_cnt_prop, YEAR == 2010 & AGE == "Age above or equal 5" & SYNDROMTYPE == "Outpatient"),  
           "COUNTRY", "VAL_MEAN", legend.title = "Proportion", diseasefree = zero_cases, text.width = 23)
```

\[1\] 0.000 0.005 0.010 0.015 0.020 0.025 0.030

``` r
dev.off()
```

png 2

``` r
setwd(params$Dir)
image <- paste0("03-estimate-rota-af_files/figure-gfm/r_world2.png")
cat("![](",image,")")
```

![](03-estimate-rota-af_files/figure-gfm/r_world2.png)

``` r
png(paste0(params$PlotDir, "/r_world3.png"), width=960, height=480)
plot_world(subset(all_cnt_prop, YEAR == 2010 & AGE == "Age below 5" & SYNDROMTYPE == "Inpatient"),  
           "COUNTRY", "VAL_MEAN", legend.title = "Proportion", diseasefree = zero_cases, text.width = 23)
```

\[1\] 0.00 0.02 0.04 0.06 0.08

``` r
dev.off()
```

png 2

``` r
setwd(params$Dir)
image <- paste0("03-estimate-rota-af_files/figure-gfm/r_world3.png")
cat("![](",image,")")
```

![](03-estimate-rota-af_files/figure-gfm/r_world3.png)

``` r
png(paste0(params$PlotDir, "/r_world4.png"), width=960, height=480)
plot_world(subset(all_cnt_prop, YEAR == 2010 & AGE == "Age below 5" & SYNDROMTYPE == "Outpatient"), 
           "COUNTRY", "VAL_MEAN", legend.title = "Proportion", diseasefree = zero_cases, text.width = 23)
```

\[1\] 0.000 0.005 0.010 0.015

``` r
dev.off()
```

png 2

``` r
setwd(params$Dir)
image <- paste0("03-estimate-rota-af_files/figure-gfm/r_world4.png")
cat("![](",image,")")
```

![](03-estimate-rota-af_files/figure-gfm/r_world4.png)

``` r
png(paste0(params$PlotDir, "/r_world5.png"), width=960, height=480)
plot_world(subset(all_cnt_prop, YEAR == 2020 & AGE == "Age above or equal 5" & SYNDROMTYPE == "Inpatient"),  
           "COUNTRY", "VAL_MEAN", legend.title = "Proportion", diseasefree = zero_cases, text.width = 23)
```

\[1\] 0.00 0.02 0.04 0.06 0.08

``` r
dev.off()
```

png 2

``` r
setwd(params$Dir)
image <- paste0("03-estimate-rota-af_files/figure-gfm/r_world5.png")
cat("![](",image,")")
```

![](03-estimate-rota-af_files/figure-gfm/r_world5.png)

``` r
png(paste0(params$PlotDir, "/r_world6.png"), width=960, height=480)
plot_world(subset(all_cnt_prop, YEAR == 2020 & AGE == "Age above or equal 5" & SYNDROMTYPE == "Outpatient"),  
           "COUNTRY", "VAL_MEAN", legend.title = "Proportion", diseasefree = zero_cases, text.width = 23)
```

\[1\] 0.000 0.005 0.010 0.015 0.020

``` r
dev.off()
```

png 2

``` r
setwd(params$Dir)
image <- paste0("03-estimate-rota-af_files/figure-gfm/r_world6.png")
cat("![](",image,")")
```

![](03-estimate-rota-af_files/figure-gfm/r_world6.png)

``` r
png(paste0(params$PlotDir, "/r_world7.png"), width=960, height=480)
plot_world(subset(all_cnt_prop, YEAR == 2020 & AGE == "Age below 5" & SYNDROMTYPE == "Inpatient"),
           "COUNTRY", "VAL_MEAN", legend.title = "Proportion", diseasefree = zero_cases, text.width = 23)
```

\[1\] 0.00 0.01 0.02 0.03 0.04 0.05

``` r
dev.off()
```

png 2

``` r
setwd(params$Dir)
image <- paste0("03-estimate-rota-af_files/figure-gfm/r_world7.png")
cat("![](",image,")")
```

![](03-estimate-rota-af_files/figure-gfm/r_world7.png)

``` r
png(paste0(params$PlotDir, "/r_world8.png"), width=960, height=480)
plot_world(subset(all_cnt_prop, YEAR == 2020 & AGE == "Age below 5" & SYNDROMTYPE == "Outpatient"),  
           "COUNTRY", "VAL_MEAN", legend.title = "Proportion", diseasefree = zero_cases, text.width = 23)
```

\[1\] 0.000 0.002 0.004 0.006 0.008 0.010

``` r
dev.off()
```

png 2

``` r
setwd(params$Dir)
image <- paste0("03-estimate-rota-af_files/figure-gfm/r_world8.png")
cat("![](",image,")")
```

![](03-estimate-rota-af_files/figure-gfm/r_world8.png)

``` r
all_cntg_prop_sub <- subset(all_cnt_prop, LOCATION_NAME == "India" & YEAR %in% c(2010,2020))
all_cntg_prop_sub$TYPE <- paste0(all_cntg_prop_sub$YEAR, ": ", all_cntg_prop_sub$SYNDROMTYPE)
png(paste0(params$PlotDir, "/r_CI1.png"), width=960, height=480)
ggplot(all_cntg_prop_sub,
       aes(y = VAL_MEAN, x = AGE)) +
  geom_pointrange(aes(ymin = VAL_LWR, ymax = VAL_UPR, group=TYPE, colour = TYPE), size = 0.2,
                  position=position_dodge(width=0.40), show.legend=TRUE) +
  coord_flip() +
  theme_bw() +
  scale_x_discrete(NULL, limits = rev(unique(all_cnt_prop$AGE))) +
  scale_y_continuous(NULL) +
  ggtitle(paste0("Proportion of ", params$Pathogen," in India in 2010 vs 2020")) +
  scale_color_manual(breaks = c("2010: Inpatient", "2010: Outpatient", "2020: Inpatient", "2020: Outpatient"),
                     values=c("springgreen1", "orange","springgreen4", "orange3"))
dev.off()
```

png 2

``` r
setwd(params$Dir)
image <- paste0("03-estimate-rota-af_files/figure-gfm/r_CI1.png")
cat("![](",image,")")
```

![](03-estimate-rota-af_files/figure-gfm/r_CI1.png)

``` r
all_cntg_prop_sub <- subset(all_cnt_prop, LOCATION_NAME == "India") %>%
  mutate(AGE = case_when(
    AGE == "Age above or equal 5" ~ ">= 5",
    AGE == "Age below 5" ~ "< 5"))
all_cntg_prop_sub$TYPE <- paste0(all_cntg_prop_sub$SYNDROMTYPE, ": ", all_cntg_prop_sub$AGE)
png(paste0(params$PlotDir, "/r_CI2.png"), width=960, height=480)
ggplot(all_cntg_prop_sub,
       aes(y = VAL_MEAN, x = YEAR)) +
  geom_pointrange(aes(ymin = VAL_LWR, ymax = VAL_UPR, group=TYPE, colour = TYPE), size = 0.2,
                  position=position_dodge(width=0.40), show.legend=TRUE) +
  # coord_flip() +
  theme_bw() +
  # scale_x_discrete(NULL, limits = rev(unique(all_cnt_prop$YEAR))) +
  scale_y_continuous(NULL) +
  ggtitle(paste0("Proportion of ", params$Pathogen," in India from 2000 to 2020")) +
  scale_color_manual(breaks = c("Inpatient: < 5", "Inpatient: >= 5", "Outpatient: < 5", "Outpatient: >= 5"),
                     values=c("springgreen1", "orange","springgreen4", "orange3"))
dev.off()
```

png 2

``` r
setwd(params$Dir)
image <- paste0("03-estimate-rota-af_files/figure-gfm/r_CI2.png")
cat("![](",image,")")
```

![](03-estimate-rota-af_files/figure-gfm/r_CI2.png)

``` r
all_cntg_prop_sub <- subset(all_cnt_prop, LOCATION_NAME == "Tanzania" & YEAR %in% c(2010,2020))
all_cntg_prop_sub$TYPE <- paste0(all_cntg_prop_sub$YEAR, ": ", all_cntg_prop_sub$SYNDROMTYPE)
png(paste0(params$PlotDir, "/r_CI3.png"), width=960, height=480)
ggplot(all_cntg_prop_sub,
       aes(y = VAL_MEAN, x = AGE)) +
  geom_pointrange(aes(ymin = VAL_LWR, ymax = VAL_UPR, group=TYPE, colour = TYPE), size = 0.2,
                  position=position_dodge(width=0.40), show.legend=TRUE) +
  coord_flip() +
  theme_bw() +
  scale_x_discrete(NULL, limits = rev(unique(all_cnt_prop$AGE))) +
  scale_y_continuous(NULL) +
  ggtitle(paste0("Proportion of ", params$Pathogen," in Tanzania in 2010 vs 2020")) +
  scale_color_manual(breaks = c("2010: Inpatient", "2010: Outpatient", "2020: Inpatient", "2020: Outpatient"),
                     values=c("springgreen1", "orange","springgreen4", "orange3"))
dev.off()
```

png 2

``` r
setwd(params$Dir)
image <- paste0("03-estimate-rota-af_files/figure-gfm/r_CI3.png")
cat("![](",image,")")
```

![](03-estimate-rota-af_files/figure-gfm/r_CI3.png)

``` r
all_cntg_prop_sub <- subset(all_cnt_prop, LOCATION_NAME == "Tanzania") %>%
  mutate(AGE = case_when(
    AGE == "Age above or equal 5" ~ ">= 5",
    AGE == "Age below 5" ~ "< 5"))
all_cntg_prop_sub$TYPE <- paste0(all_cntg_prop_sub$SYNDROMTYPE, ": ", all_cntg_prop_sub$AGE)
png(paste0(params$PlotDir, "/r_CI4.png"), width=960, height=480)
ggplot(all_cntg_prop_sub,
       aes(y = VAL_MEAN, x = YEAR)) +
  geom_pointrange(aes(ymin = VAL_LWR, ymax = VAL_UPR, group=TYPE, colour = TYPE), size = 0.2,
                  position=position_dodge(width=0.40), show.legend=TRUE) +
  # coord_flip() +
  theme_bw() +
  # scale_x_discrete(NULL, limits = rev(unique(all_cnt_prop$YEAR))) +
  scale_y_continuous(NULL) +
  ggtitle(paste0("Proportion of ", params$Pathogen," in Tanzania from 2000 to 2020")) +
  scale_color_manual(breaks = c("Inpatient: < 5", "Inpatient: >= 5", "Outpatient: < 5", "Outpatient: >= 5"),
                     values=c("springgreen1", "orange","springgreen4", "orange3"))
dev.off()
```

png 2

``` r
setwd(params$Dir)
image <- paste0("03-estimate-rota-af_files/figure-gfm/r_CI4.png")
cat("![](",image,")")
```

![](03-estimate-rota-af_files/figure-gfm/r_CI4.png)

``` r
all_cntg_prop_sub <- subset(all_cnt_prop, LOCATION_NAME == "Peru" & YEAR %in% c(2010,2020))
all_cntg_prop_sub$TYPE <- paste0(all_cntg_prop_sub$YEAR, ": ", all_cntg_prop_sub$SYNDROMTYPE)
png(paste0(params$PlotDir, "/r_CI5.png"), width=960, height=480)
ggplot(all_cntg_prop_sub,
       aes(y = VAL_MEAN, x = AGE)) +
  geom_pointrange(aes(ymin = VAL_LWR, ymax = VAL_UPR, group=TYPE, colour = TYPE), size = 0.2,
                  position=position_dodge(width=0.40), show.legend=TRUE) +
  coord_flip() +
  theme_bw() +
  scale_x_discrete(NULL, limits = rev(unique(all_cnt_prop$AGE))) +
  scale_y_continuous(NULL) +
  ggtitle(paste0("Proportion of ", params$Pathogen," in Peru in 2010 vs 2020")) +
  scale_color_manual(breaks = c("2010: Inpatient", "2010: Outpatient", "2020: Inpatient", "2020: Outpatient"),
                     values=c("springgreen1", "orange","springgreen4", "orange3"))
dev.off()
```

png 2

``` r
setwd(params$Dir)
image <- paste0("03-estimate-rota-af_files/figure-gfm/r_CI5.png")
cat("![](",image,")")
```

![](03-estimate-rota-af_files/figure-gfm/r_CI5.png)

``` r
all_cntg_prop_sub <- subset(all_cnt_prop, LOCATION_NAME == "Peru") %>%
  mutate(AGE = case_when(
    AGE == "Age above or equal 5" ~ ">= 5",
    AGE == "Age below 5" ~ "< 5"))
all_cntg_prop_sub$TYPE <- paste0(all_cntg_prop_sub$SYNDROMTYPE, ": ", all_cntg_prop_sub$AGE)
png(paste0(params$PlotDir, "/r_CI6.png"), width=960, height=480)
ggplot(all_cntg_prop_sub,
       aes(y = VAL_MEAN, x = YEAR)) +
  geom_pointrange(aes(ymin = VAL_LWR, ymax = VAL_UPR, group=TYPE, colour = TYPE), size = 0.2,
                  position=position_dodge(width=0.40), show.legend=TRUE) +
  # coord_flip() +
  theme_bw() +
  # scale_x_discrete(NULL, limits = rev(unique(all_cnt_prop$YEAR))) +
  scale_y_continuous(NULL) +
  ggtitle(paste0("Proportion of ", params$Pathogen," in Peru from 2000 to 2020")) +
  scale_color_manual(breaks = c("Inpatient: < 5", "Inpatient: >= 5", "Outpatient: < 5", "Outpatient: >= 5"),
                     values=c("springgreen1", "orange","springgreen4", "orange3"))
dev.off()
```

png 2

``` r
setwd(params$Dir)
image <- paste0("03-estimate-rota-af_files/figure-gfm/r_CI6.png")
cat("![](",image,")")
```

![](03-estimate-rota-af_files/figure-gfm/r_CI6.png)

``` r
all_cnt_prop$AGE <- relevel(all_cnt_prop$AGE, "Age below 5")
all_cnt_prop <- sort_by(all_cnt_prop, ~LOCATION_NAME+SYNDROMTYPE+AGE)

tab <-
  data.frame(subset(all_cnt_prop, YEAR == 2010)[,c("LOCATION_NAME", "AGE", "SYNDROMTYPE", "VAL_MEAN", "VAL_LWR", "VAL_UPR")],
             subset(all_cnt_prop, YEAR == 2020)[,c("VAL_MEAN", "VAL_LWR", "VAL_UPR")])
tab$LOCATION_NAME <- gsub(" \\(.*", "", tab$LOCATION_NAME)
names(tab) <-
  c("Country", "Age", "Syndrome",
    "2010.mean", "2010.lwr", "2010.upr",
    "2020.mean", "2020.lwr", "2020.upr")
tab <- tab %>%
  mutate(Age = case_when(
    Age == "Age above or equal 5" ~ ">= 5",
    Age == "Age below 5" ~ "< 5"))

kable(tab, digits = 3, row.names = FALSE,
      caption = "Estimated `r params$Pathogen` proportion by country, 2010 vs 2020")
```

| Country | Age | Syndrome | 2010.mean | 2010.lwr | 2010.upr | 2020.mean | 2020.lwr | 2020.upr |
|:---|:---|:---|---:|---:|---:|---:|---:|---:|
| Afghanistan | \< 5 | Inpatient | 0.024 | 0.004 | 0.077 | 0.010 | 0.002 | 0.033 |
| Afghanistan | \>= 5 | Inpatient | 0.044 | 0.007 | 0.151 | 0.020 | 0.003 | 0.071 |
| Afghanistan | \< 5 | Outpatient | 0.004 | 0.000 | 0.016 | 0.002 | 0.000 | 0.008 |
| Afghanistan | \>= 5 | Outpatient | 0.008 | -0.001 | 0.033 | 0.004 | 0.000 | 0.015 |
| Albania | \< 5 | Inpatient | 0.019 | 0.002 | 0.067 | 0.007 | 0.001 | 0.028 |
| Albania | \>= 5 | Inpatient | 0.035 | 0.003 | 0.128 | 0.013 | 0.001 | 0.054 |
| Albania | \< 5 | Outpatient | 0.004 | 0.000 | 0.014 | 0.001 | 0.000 | 0.006 |
| Albania | \>= 5 | Outpatient | 0.007 | -0.001 | 0.027 | 0.003 | 0.000 | 0.012 |
| Algeria | \< 5 | Inpatient | 0.028 | 0.008 | 0.076 | 0.018 | 0.004 | 0.054 |
| Algeria | \>= 5 | Inpatient | 0.054 | 0.011 | 0.155 | 0.035 | 0.005 | 0.112 |
| Algeria | \< 5 | Outpatient | 0.005 | -0.001 | 0.016 | 0.003 | 0.000 | 0.012 |
| Algeria | \>= 5 | Outpatient | 0.010 | -0.001 | 0.035 | 0.007 | -0.001 | 0.025 |
| Andorra | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Andorra | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Andorra | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Andorra | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Angola | \< 5 | Inpatient | 0.028 | 0.008 | 0.076 | 0.014 | 0.003 | 0.036 |
| Angola | \>= 5 | Inpatient | 0.054 | 0.011 | 0.155 | 0.026 | 0.005 | 0.077 |
| Angola | \< 5 | Outpatient | 0.005 | -0.001 | 0.016 | 0.003 | 0.000 | 0.009 |
| Angola | \>= 5 | Outpatient | 0.010 | -0.001 | 0.035 | 0.005 | -0.001 | 0.018 |
| Antigua and Barbuda | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Antigua and Barbuda | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Antigua and Barbuda | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Antigua and Barbuda | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Argentina | \< 5 | Inpatient | 0.020 | 0.003 | 0.066 | 0.008 | 0.001 | 0.025 |
| Argentina | \>= 5 | Inpatient | 0.038 | 0.004 | 0.137 | 0.015 | 0.002 | 0.054 |
| Argentina | \< 5 | Outpatient | 0.004 | 0.000 | 0.013 | 0.001 | 0.000 | 0.005 |
| Argentina | \>= 5 | Outpatient | 0.007 | -0.001 | 0.028 | 0.003 | 0.000 | 0.011 |
| Armenia | \< 5 | Inpatient | 0.019 | 0.002 | 0.067 | 0.007 | 0.001 | 0.028 |
| Armenia | \>= 5 | Inpatient | 0.035 | 0.003 | 0.128 | 0.014 | 0.001 | 0.054 |
| Armenia | \< 5 | Outpatient | 0.004 | 0.000 | 0.014 | 0.001 | 0.000 | 0.006 |
| Armenia | \>= 5 | Outpatient | 0.007 | -0.001 | 0.027 | 0.003 | 0.000 | 0.012 |
| Australia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Australia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Australia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Australia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Austria | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Austria | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Austria | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Austria | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Azerbaijan | \< 5 | Inpatient | 0.019 | 0.002 | 0.067 | 0.012 | 0.001 | 0.045 |
| Azerbaijan | \>= 5 | Inpatient | 0.035 | 0.003 | 0.128 | 0.022 | 0.002 | 0.085 |
| Azerbaijan | \< 5 | Outpatient | 0.004 | 0.000 | 0.014 | 0.002 | 0.000 | 0.009 |
| Azerbaijan | \>= 5 | Outpatient | 0.007 | -0.001 | 0.027 | 0.004 | 0.000 | 0.018 |
| Bahamas, The | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahamas, The | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahamas, The | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahamas, The | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahrain | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahrain | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahrain | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahrain | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bangladesh | \< 5 | Inpatient | 0.014 | 0.004 | 0.030 | 0.009 | 0.002 | 0.028 |
| Bangladesh | \>= 5 | Inpatient | 0.027 | 0.005 | 0.079 | 0.018 | 0.002 | 0.063 |
| Bangladesh | \< 5 | Outpatient | 0.003 | 0.000 | 0.007 | 0.002 | 0.000 | 0.006 |
| Bangladesh | \>= 5 | Outpatient | 0.005 | -0.001 | 0.017 | 0.003 | 0.000 | 0.014 |
| Barbados | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Barbados | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Barbados | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Barbados | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belarus | \< 5 | Inpatient | 0.019 | 0.002 | 0.067 | 0.012 | 0.001 | 0.045 |
| Belarus | \>= 5 | Inpatient | 0.035 | 0.003 | 0.128 | 0.022 | 0.002 | 0.085 |
| Belarus | \< 5 | Outpatient | 0.004 | 0.000 | 0.014 | 0.002 | 0.000 | 0.009 |
| Belarus | \>= 5 | Outpatient | 0.007 | -0.001 | 0.027 | 0.004 | 0.000 | 0.018 |
| Belgium | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belgium | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belgium | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belgium | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belize | \< 5 | Inpatient | 0.020 | 0.003 | 0.066 | 0.013 | 0.001 | 0.050 |
| Belize | \>= 5 | Inpatient | 0.038 | 0.004 | 0.137 | 0.025 | 0.002 | 0.100 |
| Belize | \< 5 | Outpatient | 0.004 | 0.000 | 0.013 | 0.002 | 0.000 | 0.010 |
| Belize | \>= 5 | Outpatient | 0.007 | -0.001 | 0.028 | 0.005 | 0.000 | 0.020 |
| Benin | \< 5 | Inpatient | 0.028 | 0.008 | 0.076 | 0.012 | 0.003 | 0.033 |
| Benin | \>= 5 | Inpatient | 0.054 | 0.011 | 0.155 | 0.024 | 0.004 | 0.070 |
| Benin | \< 5 | Outpatient | 0.005 | -0.001 | 0.016 | 0.002 | 0.000 | 0.008 |
| Benin | \>= 5 | Outpatient | 0.010 | -0.001 | 0.035 | 0.005 | -0.001 | 0.016 |
| Bhutan | \< 5 | Inpatient | 0.018 | 0.004 | 0.049 | 0.011 | 0.002 | 0.034 |
| Bhutan | \>= 5 | Inpatient | 0.036 | 0.006 | 0.112 | 0.023 | 0.003 | 0.077 |
| Bhutan | \< 5 | Outpatient | 0.003 | 0.000 | 0.011 | 0.002 | 0.000 | 0.008 |
| Bhutan | \>= 5 | Outpatient | 0.007 | -0.001 | 0.024 | 0.004 | 0.000 | 0.016 |
| Bolivia | \< 5 | Inpatient | 0.012 | 0.002 | 0.038 | 0.008 | 0.001 | 0.024 |
| Bolivia | \>= 5 | Inpatient | 0.023 | 0.003 | 0.084 | 0.015 | 0.002 | 0.052 |
| Bolivia | \< 5 | Outpatient | 0.002 | 0.000 | 0.009 | 0.001 | 0.000 | 0.005 |
| Bolivia | \>= 5 | Outpatient | 0.004 | 0.000 | 0.018 | 0.003 | 0.000 | 0.011 |
| Bosnia and Herzegovina | \< 5 | Inpatient | 0.019 | 0.002 | 0.067 | 0.012 | 0.001 | 0.045 |
| Bosnia and Herzegovina | \>= 5 | Inpatient | 0.035 | 0.003 | 0.128 | 0.022 | 0.002 | 0.085 |
| Bosnia and Herzegovina | \< 5 | Outpatient | 0.004 | 0.000 | 0.014 | 0.002 | 0.000 | 0.009 |
| Bosnia and Herzegovina | \>= 5 | Outpatient | 0.007 | -0.001 | 0.027 | 0.004 | 0.000 | 0.018 |
| Botswana | \< 5 | Inpatient | 0.027 | 0.006 | 0.077 | 0.010 | 0.002 | 0.029 |
| Botswana | \>= 5 | Inpatient | 0.051 | 0.009 | 0.163 | 0.020 | 0.003 | 0.064 |
| Botswana | \< 5 | Outpatient | 0.005 | -0.001 | 0.016 | 0.002 | 0.000 | 0.007 |
| Botswana | \>= 5 | Outpatient | 0.010 | -0.001 | 0.035 | 0.004 | 0.000 | 0.014 |
| Brazil | \< 5 | Inpatient | 0.010 | 0.002 | 0.033 | 0.007 | 0.001 | 0.023 |
| Brazil | \>= 5 | Inpatient | 0.020 | 0.003 | 0.071 | 0.013 | 0.001 | 0.050 |
| Brazil | \< 5 | Outpatient | 0.002 | 0.000 | 0.007 | 0.001 | 0.000 | 0.005 |
| Brazil | \>= 5 | Outpatient | 0.004 | 0.000 | 0.015 | 0.002 | 0.000 | 0.010 |
| Brunei Darussalam | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Brunei Darussalam | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Brunei Darussalam | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Brunei Darussalam | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bulgaria | \< 5 | Inpatient | 0.019 | 0.002 | 0.067 | 0.009 | 0.001 | 0.033 |
| Bulgaria | \>= 5 | Inpatient | 0.035 | 0.003 | 0.128 | 0.017 | 0.001 | 0.064 |
| Bulgaria | \< 5 | Outpatient | 0.004 | 0.000 | 0.014 | 0.002 | 0.000 | 0.007 |
| Bulgaria | \>= 5 | Outpatient | 0.007 | -0.001 | 0.027 | 0.003 | 0.000 | 0.014 |
| Burkina Faso | \< 5 | Inpatient | 0.020 | 0.006 | 0.051 | 0.008 | 0.002 | 0.022 |
| Burkina Faso | \>= 5 | Inpatient | 0.040 | 0.008 | 0.118 | 0.016 | 0.002 | 0.049 |
| Burkina Faso | \< 5 | Outpatient | 0.004 | 0.000 | 0.011 | 0.002 | 0.000 | 0.006 |
| Burkina Faso | \>= 5 | Outpatient | 0.008 | -0.001 | 0.025 | 0.003 | 0.000 | 0.012 |
| Burundi | \< 5 | Inpatient | 0.020 | 0.006 | 0.051 | 0.008 | 0.001 | 0.022 |
| Burundi | \>= 5 | Inpatient | 0.040 | 0.008 | 0.118 | 0.015 | 0.002 | 0.048 |
| Burundi | \< 5 | Outpatient | 0.004 | 0.000 | 0.011 | 0.002 | 0.000 | 0.005 |
| Burundi | \>= 5 | Outpatient | 0.008 | -0.001 | 0.025 | 0.003 | 0.000 | 0.011 |
| Cabo Verde | \< 5 | Inpatient | 0.028 | 0.008 | 0.076 | 0.018 | 0.004 | 0.054 |
| Cabo Verde | \>= 5 | Inpatient | 0.054 | 0.011 | 0.155 | 0.035 | 0.005 | 0.112 |
| Cabo Verde | \< 5 | Outpatient | 0.005 | -0.001 | 0.016 | 0.003 | 0.000 | 0.012 |
| Cabo Verde | \>= 5 | Outpatient | 0.010 | -0.001 | 0.035 | 0.007 | -0.001 | 0.025 |
| Cambodia | \< 5 | Inpatient | 0.004 | 0.000 | 0.018 | 0.003 | 0.000 | 0.015 |
| Cambodia | \>= 5 | Inpatient | 0.007 | 0.000 | 0.032 | 0.005 | 0.000 | 0.026 |
| Cambodia | \< 5 | Outpatient | 0.001 | 0.000 | 0.004 | 0.001 | 0.000 | 0.003 |
| Cambodia | \>= 5 | Outpatient | 0.001 | 0.000 | 0.006 | 0.001 | 0.000 | 0.005 |
| Cameroon | \< 5 | Inpatient | 0.047 | 0.008 | 0.162 | 0.018 | 0.004 | 0.050 |
| Cameroon | \>= 5 | Inpatient | 0.083 | 0.012 | 0.275 | 0.034 | 0.006 | 0.100 |
| Cameroon | \< 5 | Outpatient | 0.009 | -0.001 | 0.032 | 0.003 | 0.000 | 0.011 |
| Cameroon | \>= 5 | Outpatient | 0.016 | -0.002 | 0.058 | 0.006 | -0.001 | 0.020 |
| Canada | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Canada | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Canada | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Canada | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Central African Republic | \< 5 | Inpatient | 0.020 | 0.006 | 0.051 | 0.013 | 0.002 | 0.040 |
| Central African Republic | \>= 5 | Inpatient | 0.040 | 0.008 | 0.118 | 0.026 | 0.004 | 0.085 |
| Central African Republic | \< 5 | Outpatient | 0.004 | 0.000 | 0.011 | 0.003 | 0.000 | 0.009 |
| Central African Republic | \>= 5 | Outpatient | 0.008 | -0.001 | 0.025 | 0.005 | -0.001 | 0.019 |
| Chad | \< 5 | Inpatient | 0.020 | 0.006 | 0.051 | 0.013 | 0.002 | 0.040 |
| Chad | \>= 5 | Inpatient | 0.040 | 0.008 | 0.118 | 0.026 | 0.004 | 0.085 |
| Chad | \< 5 | Outpatient | 0.004 | 0.000 | 0.011 | 0.003 | 0.000 | 0.009 |
| Chad | \>= 5 | Outpatient | 0.008 | -0.001 | 0.025 | 0.005 | -0.001 | 0.019 |
| Chile | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Chile | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Chile | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Chile | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| China | \< 5 | Inpatient | 0.013 | 0.000 | 0.044 | 0.008 | 0.000 | 0.030 |
| China | \>= 5 | Inpatient | 0.024 | 0.001 | 0.090 | 0.015 | 0.001 | 0.058 |
| China | \< 5 | Outpatient | 0.002 | 0.000 | 0.009 | 0.001 | 0.000 | 0.006 |
| China | \>= 5 | Outpatient | 0.004 | 0.000 | 0.018 | 0.003 | 0.000 | 0.012 |
| Colombia | \< 5 | Inpatient | 0.013 | 0.002 | 0.038 | 0.007 | 0.001 | 0.023 |
| Colombia | \>= 5 | Inpatient | 0.025 | 0.003 | 0.085 | 0.014 | 0.002 | 0.050 |
| Colombia | \< 5 | Outpatient | 0.002 | 0.000 | 0.008 | 0.001 | 0.000 | 0.005 |
| Colombia | \>= 5 | Outpatient | 0.005 | 0.000 | 0.018 | 0.003 | 0.000 | 0.010 |
| Comoros | \< 5 | Inpatient | 0.028 | 0.008 | 0.076 | 0.018 | 0.004 | 0.054 |
| Comoros | \>= 5 | Inpatient | 0.054 | 0.011 | 0.155 | 0.035 | 0.005 | 0.112 |
| Comoros | \< 5 | Outpatient | 0.005 | -0.001 | 0.016 | 0.003 | 0.000 | 0.012 |
| Comoros | \>= 5 | Outpatient | 0.010 | -0.001 | 0.035 | 0.007 | -0.001 | 0.025 |
| Congo, Dem. Rep. | \< 5 | Inpatient | 0.020 | 0.006 | 0.051 | 0.011 | 0.002 | 0.031 |
| Congo, Dem. Rep. | \>= 5 | Inpatient | 0.040 | 0.008 | 0.118 | 0.022 | 0.003 | 0.068 |
| Congo, Dem. Rep. | \< 5 | Outpatient | 0.004 | 0.000 | 0.011 | 0.002 | 0.000 | 0.007 |
| Congo, Dem. Rep. | \>= 5 | Outpatient | 0.008 | -0.001 | 0.025 | 0.004 | 0.000 | 0.015 |
| Congo, Rep. | \< 5 | Inpatient | 0.028 | 0.008 | 0.076 | 0.012 | 0.003 | 0.033 |
| Congo, Rep. | \>= 5 | Inpatient | 0.054 | 0.011 | 0.155 | 0.024 | 0.005 | 0.071 |
| Congo, Rep. | \< 5 | Outpatient | 0.005 | -0.001 | 0.016 | 0.002 | 0.000 | 0.008 |
| Congo, Rep. | \>= 5 | Outpatient | 0.010 | -0.001 | 0.035 | 0.005 | -0.001 | 0.016 |
| Cook Islands | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Cook Islands | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Cook Islands | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Cook Islands | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Costa Rica | \< 5 | Inpatient | 0.020 | 0.003 | 0.066 | 0.007 | 0.001 | 0.022 |
| Costa Rica | \>= 5 | Inpatient | 0.038 | 0.004 | 0.137 | 0.013 | 0.002 | 0.047 |
| Costa Rica | \< 5 | Outpatient | 0.004 | 0.000 | 0.013 | 0.001 | 0.000 | 0.005 |
| Costa Rica | \>= 5 | Outpatient | 0.007 | -0.001 | 0.028 | 0.002 | 0.000 | 0.010 |
| Côte d’Ivoire | \< 5 | Inpatient | 0.028 | 0.008 | 0.076 | 0.012 | 0.003 | 0.032 |
| Côte d’Ivoire | \>= 5 | Inpatient | 0.054 | 0.011 | 0.155 | 0.023 | 0.004 | 0.068 |
| Côte d’Ivoire | \< 5 | Outpatient | 0.005 | -0.001 | 0.016 | 0.002 | 0.000 | 0.008 |
| Côte d’Ivoire | \>= 5 | Outpatient | 0.010 | -0.001 | 0.035 | 0.004 | 0.000 | 0.016 |
| Croatia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Croatia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Croatia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Croatia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Cuba | \< 5 | Inpatient | 0.020 | 0.003 | 0.066 | 0.013 | 0.001 | 0.050 |
| Cuba | \>= 5 | Inpatient | 0.038 | 0.004 | 0.137 | 0.025 | 0.002 | 0.100 |
| Cuba | \< 5 | Outpatient | 0.004 | 0.000 | 0.013 | 0.002 | 0.000 | 0.010 |
| Cuba | \>= 5 | Outpatient | 0.007 | -0.001 | 0.028 | 0.005 | 0.000 | 0.020 |
| Cyprus | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Cyprus | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Cyprus | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Cyprus | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Czech Republic | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Czech Republic | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Czech Republic | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Czech Republic | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Denmark | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Denmark | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Denmark | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Denmark | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Djibouti | \< 5 | Inpatient | 0.027 | 0.005 | 0.085 | 0.011 | 0.002 | 0.036 |
| Djibouti | \>= 5 | Inpatient | 0.050 | 0.008 | 0.163 | 0.021 | 0.003 | 0.071 |
| Djibouti | \< 5 | Outpatient | 0.005 | 0.000 | 0.018 | 0.002 | 0.000 | 0.008 |
| Djibouti | \>= 5 | Outpatient | 0.010 | -0.001 | 0.036 | 0.004 | 0.000 | 0.016 |
| Dominica | \< 5 | Inpatient | 0.020 | 0.003 | 0.066 | 0.013 | 0.001 | 0.050 |
| Dominica | \>= 5 | Inpatient | 0.038 | 0.004 | 0.137 | 0.025 | 0.002 | 0.100 |
| Dominica | \< 5 | Outpatient | 0.004 | 0.000 | 0.013 | 0.002 | 0.000 | 0.010 |
| Dominica | \>= 5 | Outpatient | 0.007 | -0.001 | 0.028 | 0.005 | 0.000 | 0.020 |
| Dominican Republic | \< 5 | Inpatient | 0.020 | 0.003 | 0.066 | 0.008 | 0.001 | 0.024 |
| Dominican Republic | \>= 5 | Inpatient | 0.038 | 0.004 | 0.137 | 0.015 | 0.002 | 0.053 |
| Dominican Republic | \< 5 | Outpatient | 0.004 | 0.000 | 0.013 | 0.001 | 0.000 | 0.005 |
| Dominican Republic | \>= 5 | Outpatient | 0.007 | -0.001 | 0.028 | 0.003 | 0.000 | 0.011 |
| Ecuador | \< 5 | Inpatient | 0.011 | 0.002 | 0.036 | 0.008 | 0.001 | 0.025 |
| Ecuador | \>= 5 | Inpatient | 0.022 | 0.003 | 0.080 | 0.015 | 0.002 | 0.054 |
| Ecuador | \< 5 | Outpatient | 0.002 | 0.000 | 0.008 | 0.001 | 0.000 | 0.005 |
| Ecuador | \>= 5 | Outpatient | 0.004 | 0.000 | 0.017 | 0.003 | 0.000 | 0.011 |
| Egypt, Arab Rep. | \< 5 | Inpatient | 0.027 | 0.005 | 0.085 | 0.016 | 0.003 | 0.050 |
| Egypt, Arab Rep. | \>= 5 | Inpatient | 0.050 | 0.008 | 0.163 | 0.031 | 0.004 | 0.103 |
| Egypt, Arab Rep. | \< 5 | Outpatient | 0.005 | 0.000 | 0.018 | 0.003 | 0.000 | 0.011 |
| Egypt, Arab Rep. | \>= 5 | Outpatient | 0.010 | -0.001 | 0.036 | 0.006 | -0.001 | 0.022 |
| El Salvador | \< 5 | Inpatient | 0.012 | 0.002 | 0.036 | 0.008 | 0.001 | 0.024 |
| El Salvador | \>= 5 | Inpatient | 0.023 | 0.003 | 0.081 | 0.015 | 0.002 | 0.053 |
| El Salvador | \< 5 | Outpatient | 0.002 | 0.000 | 0.008 | 0.001 | 0.000 | 0.005 |
| El Salvador | \>= 5 | Outpatient | 0.004 | 0.000 | 0.017 | 0.003 | 0.000 | 0.011 |
| Equatorial Guinea | \< 5 | Inpatient | 0.027 | 0.006 | 0.077 | 0.017 | 0.003 | 0.053 |
| Equatorial Guinea | \>= 5 | Inpatient | 0.051 | 0.009 | 0.163 | 0.032 | 0.004 | 0.111 |
| Equatorial Guinea | \< 5 | Outpatient | 0.005 | -0.001 | 0.016 | 0.003 | 0.000 | 0.011 |
| Equatorial Guinea | \>= 5 | Outpatient | 0.010 | -0.001 | 0.035 | 0.006 | -0.001 | 0.023 |
| Eritrea | \< 5 | Inpatient | 0.020 | 0.006 | 0.051 | 0.008 | 0.001 | 0.022 |
| Eritrea | \>= 5 | Inpatient | 0.040 | 0.008 | 0.118 | 0.015 | 0.002 | 0.048 |
| Eritrea | \< 5 | Outpatient | 0.004 | 0.000 | 0.011 | 0.002 | 0.000 | 0.005 |
| Eritrea | \>= 5 | Outpatient | 0.008 | -0.001 | 0.025 | 0.003 | 0.000 | 0.011 |
| Estonia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Estonia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Estonia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Estonia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Eswatini | \< 5 | Inpatient | 0.028 | 0.008 | 0.076 | 0.012 | 0.003 | 0.033 |
| Eswatini | \>= 5 | Inpatient | 0.054 | 0.011 | 0.155 | 0.024 | 0.004 | 0.070 |
| Eswatini | \< 5 | Outpatient | 0.005 | -0.001 | 0.016 | 0.002 | 0.000 | 0.008 |
| Eswatini | \>= 5 | Outpatient | 0.010 | -0.001 | 0.035 | 0.005 | -0.001 | 0.016 |
| Ethiopia | \< 5 | Inpatient | 0.020 | 0.006 | 0.051 | 0.009 | 0.002 | 0.023 |
| Ethiopia | \>= 5 | Inpatient | 0.040 | 0.008 | 0.118 | 0.017 | 0.003 | 0.053 |
| Ethiopia | \< 5 | Outpatient | 0.004 | 0.000 | 0.011 | 0.002 | 0.000 | 0.006 |
| Ethiopia | \>= 5 | Outpatient | 0.008 | -0.001 | 0.025 | 0.003 | 0.000 | 0.012 |
| Fiji | \< 5 | Inpatient | 0.013 | 0.000 | 0.044 | 0.005 | 0.000 | 0.018 |
| Fiji | \>= 5 | Inpatient | 0.024 | 0.001 | 0.090 | 0.009 | 0.000 | 0.035 |
| Fiji | \< 5 | Outpatient | 0.002 | 0.000 | 0.009 | 0.001 | 0.000 | 0.004 |
| Fiji | \>= 5 | Outpatient | 0.004 | 0.000 | 0.018 | 0.002 | 0.000 | 0.008 |
| Finland | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Finland | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Finland | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Finland | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| France | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| France | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| France | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| France | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Gabon | \< 5 | Inpatient | 0.060 | 0.009 | 0.189 | 0.035 | 0.007 | 0.100 |
| Gabon | \>= 5 | Inpatient | 0.106 | 0.013 | 0.341 | 0.065 | 0.009 | 0.208 |
| Gabon | \< 5 | Outpatient | 0.011 | -0.001 | 0.041 | 0.006 | -0.001 | 0.021 |
| Gabon | \>= 5 | Outpatient | 0.021 | -0.003 | 0.084 | 0.012 | -0.002 | 0.045 |
| Gambia, The | \< 5 | Inpatient | 0.013 | 0.004 | 0.035 | 0.006 | 0.001 | 0.020 |
| Gambia, The | \>= 5 | Inpatient | 0.027 | 0.005 | 0.083 | 0.011 | 0.001 | 0.043 |
| Gambia, The | \< 5 | Outpatient | 0.003 | 0.000 | 0.008 | 0.001 | 0.000 | 0.005 |
| Gambia, The | \>= 5 | Outpatient | 0.005 | -0.001 | 0.018 | 0.002 | 0.000 | 0.010 |
| Georgia | \< 5 | Inpatient | 0.019 | 0.002 | 0.067 | 0.008 | 0.001 | 0.029 |
| Georgia | \>= 5 | Inpatient | 0.035 | 0.003 | 0.128 | 0.015 | 0.001 | 0.056 |
| Georgia | \< 5 | Outpatient | 0.004 | 0.000 | 0.014 | 0.001 | 0.000 | 0.006 |
| Georgia | \>= 5 | Outpatient | 0.007 | -0.001 | 0.027 | 0.003 | 0.000 | 0.012 |
| Germany | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Germany | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Germany | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Germany | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Ghana | \< 5 | Inpatient | 0.028 | 0.008 | 0.076 | 0.011 | 0.002 | 0.030 |
| Ghana | \>= 5 | Inpatient | 0.054 | 0.011 | 0.155 | 0.021 | 0.004 | 0.064 |
| Ghana | \< 5 | Outpatient | 0.005 | -0.001 | 0.016 | 0.002 | 0.000 | 0.007 |
| Ghana | \>= 5 | Outpatient | 0.010 | -0.001 | 0.035 | 0.004 | 0.000 | 0.015 |
| Greece | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Greece | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Greece | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Greece | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Grenada | \< 5 | Inpatient | 0.020 | 0.003 | 0.066 | 0.013 | 0.001 | 0.050 |
| Grenada | \>= 5 | Inpatient | 0.038 | 0.004 | 0.137 | 0.025 | 0.002 | 0.100 |
| Grenada | \< 5 | Outpatient | 0.004 | 0.000 | 0.013 | 0.002 | 0.000 | 0.010 |
| Grenada | \>= 5 | Outpatient | 0.007 | -0.001 | 0.028 | 0.005 | 0.000 | 0.020 |
| Guatemala | \< 5 | Inpatient | 0.015 | 0.003 | 0.047 | 0.007 | 0.001 | 0.023 |
| Guatemala | \>= 5 | Inpatient | 0.030 | 0.004 | 0.105 | 0.014 | 0.002 | 0.051 |
| Guatemala | \< 5 | Outpatient | 0.003 | 0.000 | 0.010 | 0.001 | 0.000 | 0.005 |
| Guatemala | \>= 5 | Outpatient | 0.006 | -0.001 | 0.021 | 0.003 | 0.000 | 0.010 |
| Guinea | \< 5 | Inpatient | 0.028 | 0.008 | 0.076 | 0.018 | 0.004 | 0.054 |
| Guinea | \>= 5 | Inpatient | 0.054 | 0.011 | 0.155 | 0.035 | 0.005 | 0.112 |
| Guinea | \< 5 | Outpatient | 0.005 | -0.001 | 0.016 | 0.003 | 0.000 | 0.012 |
| Guinea | \>= 5 | Outpatient | 0.010 | -0.001 | 0.035 | 0.007 | -0.001 | 0.025 |
| Guinea-Bissau | \< 5 | Inpatient | 0.020 | 0.006 | 0.051 | 0.008 | 0.002 | 0.022 |
| Guinea-Bissau | \>= 5 | Inpatient | 0.040 | 0.008 | 0.118 | 0.016 | 0.003 | 0.050 |
| Guinea-Bissau | \< 5 | Outpatient | 0.004 | 0.000 | 0.011 | 0.002 | 0.000 | 0.006 |
| Guinea-Bissau | \>= 5 | Outpatient | 0.008 | -0.001 | 0.025 | 0.003 | 0.000 | 0.012 |
| Guyana | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Guyana | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Guyana | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Guyana | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Haiti | \< 5 | Inpatient | 0.020 | 0.003 | 0.063 | 0.009 | 0.001 | 0.029 |
| Haiti | \>= 5 | Inpatient | 0.037 | 0.005 | 0.131 | 0.017 | 0.002 | 0.060 |
| Haiti | \< 5 | Outpatient | 0.004 | 0.000 | 0.013 | 0.002 | 0.000 | 0.006 |
| Haiti | \>= 5 | Outpatient | 0.007 | -0.001 | 0.027 | 0.003 | 0.000 | 0.013 |
| Honduras | \< 5 | Inpatient | 0.012 | 0.002 | 0.038 | 0.007 | 0.001 | 0.023 |
| Honduras | \>= 5 | Inpatient | 0.023 | 0.003 | 0.083 | 0.014 | 0.002 | 0.050 |
| Honduras | \< 5 | Outpatient | 0.002 | 0.000 | 0.009 | 0.001 | 0.000 | 0.005 |
| Honduras | \>= 5 | Outpatient | 0.004 | 0.000 | 0.018 | 0.003 | 0.000 | 0.011 |
| Hungary | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Hungary | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Hungary | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Hungary | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Iceland | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Iceland | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Iceland | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Iceland | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| India | \< 5 | Inpatient | 0.007 | 0.002 | 0.017 | 0.003 | 0.001 | 0.011 |
| India | \>= 5 | Inpatient | 0.015 | 0.003 | 0.041 | 0.006 | 0.001 | 0.023 |
| India | \< 5 | Outpatient | 0.001 | 0.000 | 0.004 | 0.001 | 0.000 | 0.003 |
| India | \>= 5 | Outpatient | 0.003 | 0.000 | 0.008 | 0.001 | 0.000 | 0.005 |
| Indonesia | \< 5 | Inpatient | 0.019 | 0.004 | 0.053 | 0.012 | 0.002 | 0.037 |
| Indonesia | \>= 5 | Inpatient | 0.036 | 0.005 | 0.115 | 0.023 | 0.003 | 0.080 |
| Indonesia | \< 5 | Outpatient | 0.003 | 0.000 | 0.011 | 0.002 | 0.000 | 0.008 |
| Indonesia | \>= 5 | Outpatient | 0.007 | -0.001 | 0.025 | 0.004 | 0.000 | 0.017 |
| Iran, Islamic Rep. | \< 5 | Inpatient | 0.027 | 0.005 | 0.085 | 0.016 | 0.003 | 0.050 |
| Iran, Islamic Rep. | \>= 5 | Inpatient | 0.050 | 0.008 | 0.163 | 0.031 | 0.004 | 0.103 |
| Iran, Islamic Rep. | \< 5 | Outpatient | 0.005 | 0.000 | 0.018 | 0.003 | 0.000 | 0.011 |
| Iran, Islamic Rep. | \>= 5 | Outpatient | 0.010 | -0.001 | 0.036 | 0.006 | -0.001 | 0.022 |
| Iraq | \< 5 | Inpatient | 0.027 | 0.005 | 0.085 | 0.013 | 0.002 | 0.039 |
| Iraq | \>= 5 | Inpatient | 0.050 | 0.008 | 0.163 | 0.024 | 0.004 | 0.080 |
| Iraq | \< 5 | Outpatient | 0.005 | 0.000 | 0.018 | 0.002 | 0.000 | 0.009 |
| Iraq | \>= 5 | Outpatient | 0.010 | -0.001 | 0.036 | 0.005 | 0.000 | 0.018 |
| Ireland | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Ireland | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Ireland | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Ireland | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Israel | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Israel | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Israel | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Israel | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Italy | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Italy | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Italy | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Italy | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Jamaica | \< 5 | Inpatient | 0.020 | 0.003 | 0.066 | 0.013 | 0.001 | 0.050 |
| Jamaica | \>= 5 | Inpatient | 0.038 | 0.004 | 0.137 | 0.025 | 0.002 | 0.100 |
| Jamaica | \< 5 | Outpatient | 0.004 | 0.000 | 0.013 | 0.002 | 0.000 | 0.010 |
| Jamaica | \>= 5 | Outpatient | 0.007 | -0.001 | 0.028 | 0.005 | 0.000 | 0.020 |
| Japan | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Japan | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Japan | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Japan | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Jordan | \< 5 | Inpatient | 0.027 | 0.005 | 0.085 | 0.011 | 0.002 | 0.035 |
| Jordan | \>= 5 | Inpatient | 0.050 | 0.008 | 0.163 | 0.021 | 0.003 | 0.071 |
| Jordan | \< 5 | Outpatient | 0.005 | 0.000 | 0.018 | 0.002 | 0.000 | 0.008 |
| Jordan | \>= 5 | Outpatient | 0.010 | -0.001 | 0.036 | 0.004 | 0.000 | 0.016 |
| Kazakhstan | \< 5 | Inpatient | 0.019 | 0.002 | 0.067 | 0.012 | 0.001 | 0.045 |
| Kazakhstan | \>= 5 | Inpatient | 0.035 | 0.003 | 0.128 | 0.022 | 0.002 | 0.085 |
| Kazakhstan | \< 5 | Outpatient | 0.004 | 0.000 | 0.014 | 0.002 | 0.000 | 0.009 |
| Kazakhstan | \>= 5 | Outpatient | 0.007 | -0.001 | 0.027 | 0.004 | 0.000 | 0.018 |
| Kenya | \< 5 | Inpatient | 0.016 | 0.005 | 0.035 | 0.006 | 0.001 | 0.020 |
| Kenya | \>= 5 | Inpatient | 0.031 | 0.007 | 0.081 | 0.012 | 0.002 | 0.042 |
| Kenya | \< 5 | Outpatient | 0.003 | 0.000 | 0.009 | 0.001 | 0.000 | 0.005 |
| Kenya | \>= 5 | Outpatient | 0.006 | -0.001 | 0.020 | 0.002 | 0.000 | 0.010 |
| Kiribati | \< 5 | Inpatient | 0.010 | 0.000 | 0.038 | 0.004 | 0.000 | 0.015 |
| Kiribati | \>= 5 | Inpatient | 0.018 | 0.001 | 0.077 | 0.007 | 0.000 | 0.029 |
| Kiribati | \< 5 | Outpatient | 0.002 | 0.000 | 0.008 | 0.001 | 0.000 | 0.003 |
| Kiribati | \>= 5 | Outpatient | 0.003 | 0.000 | 0.016 | 0.001 | 0.000 | 0.006 |
| Korea, Dem. People’s Rep. | \< 5 | Inpatient | 0.018 | 0.004 | 0.049 | 0.011 | 0.002 | 0.034 |
| Korea, Dem. People’s Rep. | \>= 5 | Inpatient | 0.036 | 0.006 | 0.112 | 0.023 | 0.003 | 0.077 |
| Korea, Dem. People’s Rep. | \< 5 | Outpatient | 0.003 | 0.000 | 0.011 | 0.002 | 0.000 | 0.008 |
| Korea, Dem. People’s Rep. | \>= 5 | Outpatient | 0.007 | -0.001 | 0.024 | 0.004 | 0.000 | 0.016 |
| Korea, Rep. | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Korea, Rep. | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Korea, Rep. | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Korea, Rep. | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kuwait | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kuwait | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kuwait | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kuwait | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kyrgyz Republic | \< 5 | Inpatient | 0.019 | 0.002 | 0.062 | 0.009 | 0.001 | 0.029 |
| Kyrgyz Republic | \>= 5 | Inpatient | 0.035 | 0.004 | 0.123 | 0.016 | 0.002 | 0.059 |
| Kyrgyz Republic | \< 5 | Outpatient | 0.004 | 0.000 | 0.013 | 0.002 | 0.000 | 0.006 |
| Kyrgyz Republic | \>= 5 | Outpatient | 0.007 | -0.001 | 0.025 | 0.003 | 0.000 | 0.012 |
| Lao PDR | \< 5 | Inpatient | 0.010 | 0.000 | 0.038 | 0.006 | 0.000 | 0.025 |
| Lao PDR | \>= 5 | Inpatient | 0.018 | 0.001 | 0.077 | 0.011 | 0.000 | 0.048 |
| Lao PDR | \< 5 | Outpatient | 0.002 | 0.000 | 0.008 | 0.001 | 0.000 | 0.005 |
| Lao PDR | \>= 5 | Outpatient | 0.003 | 0.000 | 0.016 | 0.002 | 0.000 | 0.010 |
| Latvia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Latvia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Latvia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Latvia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Lebanon | \< 5 | Inpatient | 0.056 | 0.007 | 0.198 | 0.030 | 0.006 | 0.087 |
| Lebanon | \>= 5 | Inpatient | 0.096 | 0.012 | 0.309 | 0.055 | 0.010 | 0.158 |
| Lebanon | \< 5 | Outpatient | 0.011 | -0.001 | 0.043 | 0.006 | -0.001 | 0.018 |
| Lebanon | \>= 5 | Outpatient | 0.019 | -0.002 | 0.072 | 0.010 | -0.001 | 0.032 |
| Lesotho | \< 5 | Inpatient | 0.028 | 0.008 | 0.076 | 0.011 | 0.002 | 0.030 |
| Lesotho | \>= 5 | Inpatient | 0.054 | 0.011 | 0.155 | 0.021 | 0.004 | 0.064 |
| Lesotho | \< 5 | Outpatient | 0.005 | -0.001 | 0.016 | 0.002 | 0.000 | 0.007 |
| Lesotho | \>= 5 | Outpatient | 0.010 | -0.001 | 0.035 | 0.004 | 0.000 | 0.015 |
| Liberia | \< 5 | Inpatient | 0.020 | 0.006 | 0.051 | 0.009 | 0.002 | 0.024 |
| Liberia | \>= 5 | Inpatient | 0.040 | 0.008 | 0.118 | 0.017 | 0.003 | 0.054 |
| Liberia | \< 5 | Outpatient | 0.004 | 0.000 | 0.011 | 0.002 | 0.000 | 0.006 |
| Liberia | \>= 5 | Outpatient | 0.008 | -0.001 | 0.025 | 0.003 | 0.000 | 0.013 |
| Libya | \< 5 | Inpatient | 0.027 | 0.005 | 0.085 | 0.011 | 0.002 | 0.036 |
| Libya | \>= 5 | Inpatient | 0.050 | 0.008 | 0.163 | 0.021 | 0.003 | 0.071 |
| Libya | \< 5 | Outpatient | 0.005 | 0.000 | 0.018 | 0.002 | 0.000 | 0.008 |
| Libya | \>= 5 | Outpatient | 0.010 | -0.001 | 0.036 | 0.004 | 0.000 | 0.016 |
| Lithuania | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Lithuania | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Lithuania | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Lithuania | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Luxembourg | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Luxembourg | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Luxembourg | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Luxembourg | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Madagascar | \< 5 | Inpatient | 0.020 | 0.006 | 0.051 | 0.009 | 0.002 | 0.024 |
| Madagascar | \>= 5 | Inpatient | 0.040 | 0.008 | 0.118 | 0.017 | 0.003 | 0.054 |
| Madagascar | \< 5 | Outpatient | 0.004 | 0.000 | 0.011 | 0.002 | 0.000 | 0.006 |
| Madagascar | \>= 5 | Outpatient | 0.008 | -0.001 | 0.025 | 0.003 | 0.000 | 0.013 |
| Malawi | \< 5 | Inpatient | 0.032 | 0.006 | 0.095 | 0.011 | 0.003 | 0.028 |
| Malawi | \>= 5 | Inpatient | 0.061 | 0.009 | 0.199 | 0.021 | 0.004 | 0.067 |
| Malawi | \< 5 | Outpatient | 0.006 | -0.001 | 0.023 | 0.002 | 0.000 | 0.008 |
| Malawi | \>= 5 | Outpatient | 0.012 | -0.001 | 0.049 | 0.004 | 0.000 | 0.017 |
| Malaysia | \< 5 | Inpatient | 0.013 | 0.000 | 0.044 | 0.008 | 0.000 | 0.030 |
| Malaysia | \>= 5 | Inpatient | 0.024 | 0.001 | 0.090 | 0.015 | 0.001 | 0.058 |
| Malaysia | \< 5 | Outpatient | 0.002 | 0.000 | 0.009 | 0.001 | 0.000 | 0.006 |
| Malaysia | \>= 5 | Outpatient | 0.004 | 0.000 | 0.018 | 0.003 | 0.000 | 0.012 |
| Maldives | \< 5 | Inpatient | 0.019 | 0.004 | 0.053 | 0.012 | 0.002 | 0.037 |
| Maldives | \>= 5 | Inpatient | 0.036 | 0.005 | 0.115 | 0.023 | 0.003 | 0.080 |
| Maldives | \< 5 | Outpatient | 0.003 | 0.000 | 0.011 | 0.002 | 0.000 | 0.008 |
| Maldives | \>= 5 | Outpatient | 0.007 | -0.001 | 0.025 | 0.004 | 0.000 | 0.017 |
| Mali | \< 5 | Inpatient | 0.015 | 0.004 | 0.039 | 0.007 | 0.001 | 0.023 |
| Mali | \>= 5 | Inpatient | 0.030 | 0.005 | 0.093 | 0.014 | 0.002 | 0.050 |
| Mali | \< 5 | Outpatient | 0.003 | 0.000 | 0.009 | 0.001 | 0.000 | 0.005 |
| Mali | \>= 5 | Outpatient | 0.006 | -0.001 | 0.020 | 0.003 | 0.000 | 0.012 |
| Malta | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Malta | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Malta | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Malta | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Marshall Islands | \< 5 | Inpatient | 0.011 | 0.000 | 0.038 | 0.006 | 0.000 | 0.020 |
| Marshall Islands | \>= 5 | Inpatient | 0.021 | 0.001 | 0.077 | 0.011 | 0.000 | 0.039 |
| Marshall Islands | \< 5 | Outpatient | 0.002 | 0.000 | 0.008 | 0.001 | 0.000 | 0.004 |
| Marshall Islands | \>= 5 | Outpatient | 0.004 | 0.000 | 0.016 | 0.002 | 0.000 | 0.008 |
| Mauritania | \< 5 | Inpatient | 0.028 | 0.008 | 0.076 | 0.011 | 0.002 | 0.030 |
| Mauritania | \>= 5 | Inpatient | 0.054 | 0.011 | 0.155 | 0.021 | 0.004 | 0.065 |
| Mauritania | \< 5 | Outpatient | 0.005 | -0.001 | 0.016 | 0.002 | 0.000 | 0.007 |
| Mauritania | \>= 5 | Outpatient | 0.010 | -0.001 | 0.035 | 0.004 | 0.000 | 0.015 |
| Mauritius | \< 5 | Inpatient | 0.027 | 0.006 | 0.077 | 0.010 | 0.002 | 0.029 |
| Mauritius | \>= 5 | Inpatient | 0.051 | 0.009 | 0.163 | 0.019 | 0.003 | 0.063 |
| Mauritius | \< 5 | Outpatient | 0.005 | -0.001 | 0.016 | 0.002 | 0.000 | 0.007 |
| Mauritius | \>= 5 | Outpatient | 0.010 | -0.001 | 0.035 | 0.004 | 0.000 | 0.014 |
| Mexico | \< 5 | Inpatient | 0.012 | 0.002 | 0.036 | 0.008 | 0.001 | 0.025 |
| Mexico | \>= 5 | Inpatient | 0.023 | 0.003 | 0.081 | 0.015 | 0.002 | 0.055 |
| Mexico | \< 5 | Outpatient | 0.002 | 0.000 | 0.008 | 0.001 | 0.000 | 0.006 |
| Mexico | \>= 5 | Outpatient | 0.004 | 0.000 | 0.017 | 0.003 | 0.000 | 0.011 |
| Micronesia, Fed. Sts. | \< 5 | Inpatient | 0.007 | 0.000 | 0.028 | 0.005 | 0.000 | 0.018 |
| Micronesia, Fed. Sts. | \>= 5 | Inpatient | 0.013 | 0.000 | 0.055 | 0.009 | 0.000 | 0.036 |
| Micronesia, Fed. Sts. | \< 5 | Outpatient | 0.001 | 0.000 | 0.006 | 0.001 | 0.000 | 0.004 |
| Micronesia, Fed. Sts. | \>= 5 | Outpatient | 0.002 | 0.000 | 0.012 | 0.002 | 0.000 | 0.007 |
| Moldova | \< 5 | Inpatient | 0.019 | 0.002 | 0.067 | 0.009 | 0.001 | 0.031 |
| Moldova | \>= 5 | Inpatient | 0.035 | 0.003 | 0.128 | 0.016 | 0.001 | 0.061 |
| Moldova | \< 5 | Outpatient | 0.004 | 0.000 | 0.014 | 0.002 | 0.000 | 0.007 |
| Moldova | \>= 5 | Outpatient | 0.007 | -0.001 | 0.027 | 0.003 | 0.000 | 0.013 |
| Monaco | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Monaco | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Monaco | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Monaco | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Mongolia | \< 5 | Inpatient | 0.010 | 0.000 | 0.038 | 0.006 | 0.000 | 0.025 |
| Mongolia | \>= 5 | Inpatient | 0.018 | 0.001 | 0.077 | 0.011 | 0.000 | 0.048 |
| Mongolia | \< 5 | Outpatient | 0.002 | 0.000 | 0.008 | 0.001 | 0.000 | 0.005 |
| Mongolia | \>= 5 | Outpatient | 0.003 | 0.000 | 0.016 | 0.002 | 0.000 | 0.010 |
| Montenegro | \< 5 | Inpatient | 0.019 | 0.002 | 0.067 | 0.012 | 0.001 | 0.045 |
| Montenegro | \>= 5 | Inpatient | 0.035 | 0.003 | 0.128 | 0.022 | 0.002 | 0.085 |
| Montenegro | \< 5 | Outpatient | 0.004 | 0.000 | 0.014 | 0.002 | 0.000 | 0.009 |
| Montenegro | \>= 5 | Outpatient | 0.007 | -0.001 | 0.027 | 0.004 | 0.000 | 0.018 |
| Morocco | \< 5 | Inpatient | 0.026 | 0.005 | 0.084 | 0.010 | 0.001 | 0.035 |
| Morocco | \>= 5 | Inpatient | 0.049 | 0.008 | 0.162 | 0.019 | 0.002 | 0.069 |
| Morocco | \< 5 | Outpatient | 0.005 | 0.000 | 0.018 | 0.002 | 0.000 | 0.008 |
| Morocco | \>= 5 | Outpatient | 0.009 | -0.001 | 0.036 | 0.004 | 0.000 | 0.016 |
| Mozambique | \< 5 | Inpatient | 0.018 | 0.005 | 0.044 | 0.008 | 0.002 | 0.026 |
| Mozambique | \>= 5 | Inpatient | 0.035 | 0.007 | 0.100 | 0.017 | 0.002 | 0.057 |
| Mozambique | \< 5 | Outpatient | 0.003 | 0.000 | 0.011 | 0.002 | 0.000 | 0.006 |
| Mozambique | \>= 5 | Outpatient | 0.007 | -0.001 | 0.023 | 0.003 | 0.000 | 0.014 |
| Myanmar | \< 5 | Inpatient | 0.018 | 0.004 | 0.049 | 0.009 | 0.002 | 0.025 |
| Myanmar | \>= 5 | Inpatient | 0.036 | 0.006 | 0.112 | 0.017 | 0.002 | 0.056 |
| Myanmar | \< 5 | Outpatient | 0.003 | 0.000 | 0.011 | 0.002 | 0.000 | 0.006 |
| Myanmar | \>= 5 | Outpatient | 0.007 | -0.001 | 0.024 | 0.003 | 0.000 | 0.012 |
| Namibia | \< 5 | Inpatient | 0.027 | 0.006 | 0.077 | 0.010 | 0.002 | 0.030 |
| Namibia | \>= 5 | Inpatient | 0.051 | 0.009 | 0.163 | 0.020 | 0.003 | 0.066 |
| Namibia | \< 5 | Outpatient | 0.005 | -0.001 | 0.016 | 0.002 | 0.000 | 0.007 |
| Namibia | \>= 5 | Outpatient | 0.010 | -0.001 | 0.035 | 0.004 | 0.000 | 0.014 |
| Nauru | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nauru | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nauru | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nauru | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nepal | \< 5 | Inpatient | 0.051 | 0.008 | 0.138 | 0.029 | 0.007 | 0.071 |
| Nepal | \>= 5 | Inpatient | 0.095 | 0.011 | 0.283 | 0.057 | 0.009 | 0.168 |
| Nepal | \< 5 | Outpatient | 0.010 | -0.001 | 0.035 | 0.006 | -0.001 | 0.018 |
| Nepal | \>= 5 | Outpatient | 0.020 | -0.002 | 0.076 | 0.012 | -0.001 | 0.043 |
| Netherlands | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Netherlands | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Netherlands | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Netherlands | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| New Zealand | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| New Zealand | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| New Zealand | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| New Zealand | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nicaragua | \< 5 | Inpatient | 0.012 | 0.002 | 0.038 | 0.007 | 0.001 | 0.022 |
| Nicaragua | \>= 5 | Inpatient | 0.022 | 0.003 | 0.083 | 0.014 | 0.002 | 0.048 |
| Nicaragua | \< 5 | Outpatient | 0.002 | 0.000 | 0.009 | 0.001 | 0.000 | 0.005 |
| Nicaragua | \>= 5 | Outpatient | 0.004 | 0.000 | 0.018 | 0.003 | 0.000 | 0.010 |
| Niger | \< 5 | Inpatient | 0.020 | 0.006 | 0.051 | 0.008 | 0.002 | 0.022 |
| Niger | \>= 5 | Inpatient | 0.040 | 0.008 | 0.118 | 0.016 | 0.002 | 0.049 |
| Niger | \< 5 | Outpatient | 0.004 | 0.000 | 0.011 | 0.002 | 0.000 | 0.006 |
| Niger | \>= 5 | Outpatient | 0.008 | -0.001 | 0.025 | 0.003 | 0.000 | 0.012 |
| Nigeria | \< 5 | Inpatient | 0.074 | 0.010 | 0.260 | 0.046 | 0.006 | 0.164 |
| Nigeria | \>= 5 | Inpatient | 0.121 | 0.016 | 0.380 | 0.079 | 0.010 | 0.264 |
| Nigeria | \< 5 | Outpatient | 0.015 | -0.001 | 0.061 | 0.009 | -0.001 | 0.036 |
| Nigeria | \>= 5 | Outpatient | 0.025 | -0.003 | 0.097 | 0.016 | -0.002 | 0.061 |
| Niue | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Niue | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Niue | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Niue | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| North Macedonia | \< 5 | Inpatient | 0.019 | 0.002 | 0.067 | 0.008 | 0.001 | 0.030 |
| North Macedonia | \>= 5 | Inpatient | 0.035 | 0.003 | 0.128 | 0.016 | 0.001 | 0.059 |
| North Macedonia | \< 5 | Outpatient | 0.004 | 0.000 | 0.014 | 0.002 | 0.000 | 0.007 |
| North Macedonia | \>= 5 | Outpatient | 0.007 | -0.001 | 0.027 | 0.003 | 0.000 | 0.013 |
| Norway | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Norway | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Norway | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Norway | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Oman | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Oman | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Oman | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Oman | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Pakistan | \< 5 | Inpatient | 0.019 | 0.006 | 0.044 | 0.008 | 0.002 | 0.025 |
| Pakistan | \>= 5 | Inpatient | 0.037 | 0.007 | 0.106 | 0.016 | 0.002 | 0.053 |
| Pakistan | \< 5 | Outpatient | 0.004 | 0.000 | 0.010 | 0.002 | 0.000 | 0.006 |
| Pakistan | \>= 5 | Outpatient | 0.007 | -0.001 | 0.023 | 0.003 | 0.000 | 0.012 |
| Palau | \< 5 | Inpatient | 0.009 | 0.000 | 0.032 | 0.005 | 0.000 | 0.018 |
| Palau | \>= 5 | Inpatient | 0.017 | 0.001 | 0.066 | 0.009 | 0.000 | 0.036 |
| Palau | \< 5 | Outpatient | 0.002 | 0.000 | 0.007 | 0.001 | 0.000 | 0.004 |
| Palau | \>= 5 | Outpatient | 0.003 | 0.000 | 0.014 | 0.002 | 0.000 | 0.008 |
| Panama | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Panama | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Panama | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Panama | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Papua New Guinea | \< 5 | Inpatient | 0.010 | 0.000 | 0.038 | 0.006 | 0.000 | 0.025 |
| Papua New Guinea | \>= 5 | Inpatient | 0.018 | 0.001 | 0.077 | 0.011 | 0.000 | 0.048 |
| Papua New Guinea | \< 5 | Outpatient | 0.002 | 0.000 | 0.008 | 0.001 | 0.000 | 0.005 |
| Papua New Guinea | \>= 5 | Outpatient | 0.003 | 0.000 | 0.016 | 0.002 | 0.000 | 0.010 |
| Paraguay | \< 5 | Inpatient | 0.014 | 0.003 | 0.041 | 0.007 | 0.001 | 0.024 |
| Paraguay | \>= 5 | Inpatient | 0.027 | 0.004 | 0.094 | 0.015 | 0.002 | 0.052 |
| Paraguay | \< 5 | Outpatient | 0.003 | 0.000 | 0.009 | 0.001 | 0.000 | 0.005 |
| Paraguay | \>= 5 | Outpatient | 0.005 | -0.001 | 0.019 | 0.003 | 0.000 | 0.011 |
| Peru | \< 5 | Inpatient | 0.014 | 0.003 | 0.037 | 0.008 | 0.002 | 0.025 |
| Peru | \>= 5 | Inpatient | 0.027 | 0.004 | 0.087 | 0.017 | 0.002 | 0.057 |
| Peru | \< 5 | Outpatient | 0.003 | 0.000 | 0.008 | 0.002 | 0.000 | 0.005 |
| Peru | \>= 5 | Outpatient | 0.005 | -0.001 | 0.019 | 0.003 | 0.000 | 0.012 |
| Philippines | \< 5 | Inpatient | 0.010 | 0.000 | 0.038 | 0.006 | 0.000 | 0.025 |
| Philippines | \>= 5 | Inpatient | 0.018 | 0.001 | 0.077 | 0.011 | 0.000 | 0.048 |
| Philippines | \< 5 | Outpatient | 0.002 | 0.000 | 0.008 | 0.001 | 0.000 | 0.005 |
| Philippines | \>= 5 | Outpatient | 0.003 | 0.000 | 0.016 | 0.002 | 0.000 | 0.010 |
| Poland | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Poland | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Poland | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Poland | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Portugal | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Portugal | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Portugal | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Portugal | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Qatar | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Qatar | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Qatar | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Qatar | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Romania | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Romania | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Romania | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Romania | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Russian Federation | \< 5 | Inpatient | 0.019 | 0.002 | 0.067 | 0.012 | 0.001 | 0.045 |
| Russian Federation | \>= 5 | Inpatient | 0.035 | 0.003 | 0.128 | 0.022 | 0.002 | 0.085 |
| Russian Federation | \< 5 | Outpatient | 0.004 | 0.000 | 0.014 | 0.002 | 0.000 | 0.009 |
| Russian Federation | \>= 5 | Outpatient | 0.007 | -0.001 | 0.027 | 0.004 | 0.000 | 0.018 |
| Rwanda | \< 5 | Inpatient | 0.020 | 0.006 | 0.051 | 0.008 | 0.002 | 0.022 |
| Rwanda | \>= 5 | Inpatient | 0.040 | 0.008 | 0.118 | 0.015 | 0.002 | 0.049 |
| Rwanda | \< 5 | Outpatient | 0.004 | 0.000 | 0.011 | 0.002 | 0.000 | 0.006 |
| Rwanda | \>= 5 | Outpatient | 0.008 | -0.001 | 0.025 | 0.003 | 0.000 | 0.012 |
| Samoa | \< 5 | Inpatient | 0.010 | 0.000 | 0.038 | 0.006 | 0.000 | 0.025 |
| Samoa | \>= 5 | Inpatient | 0.018 | 0.001 | 0.077 | 0.011 | 0.000 | 0.048 |
| Samoa | \< 5 | Outpatient | 0.002 | 0.000 | 0.008 | 0.001 | 0.000 | 0.005 |
| Samoa | \>= 5 | Outpatient | 0.003 | 0.000 | 0.016 | 0.002 | 0.000 | 0.010 |
| San Marino | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| San Marino | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| San Marino | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| San Marino | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| São Tomé and Principe | \< 5 | Inpatient | 0.028 | 0.008 | 0.076 | 0.011 | 0.002 | 0.030 |
| São Tomé and Principe | \>= 5 | Inpatient | 0.054 | 0.011 | 0.155 | 0.021 | 0.004 | 0.064 |
| São Tomé and Principe | \< 5 | Outpatient | 0.005 | -0.001 | 0.016 | 0.002 | 0.000 | 0.007 |
| São Tomé and Principe | \>= 5 | Outpatient | 0.010 | -0.001 | 0.035 | 0.004 | 0.000 | 0.015 |
| Saudi Arabia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Saudi Arabia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Saudi Arabia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Saudi Arabia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Senegal | \< 5 | Inpatient | 0.028 | 0.008 | 0.076 | 0.011 | 0.002 | 0.030 |
| Senegal | \>= 5 | Inpatient | 0.054 | 0.011 | 0.155 | 0.020 | 0.004 | 0.064 |
| Senegal | \< 5 | Outpatient | 0.005 | -0.001 | 0.016 | 0.002 | 0.000 | 0.007 |
| Senegal | \>= 5 | Outpatient | 0.010 | -0.001 | 0.035 | 0.004 | 0.000 | 0.015 |
| Serbia | \< 5 | Inpatient | 0.019 | 0.002 | 0.067 | 0.012 | 0.001 | 0.045 |
| Serbia | \>= 5 | Inpatient | 0.035 | 0.003 | 0.128 | 0.022 | 0.002 | 0.085 |
| Serbia | \< 5 | Outpatient | 0.004 | 0.000 | 0.014 | 0.002 | 0.000 | 0.009 |
| Serbia | \>= 5 | Outpatient | 0.007 | -0.001 | 0.027 | 0.004 | 0.000 | 0.018 |
| Seychelles | \< 5 | Inpatient | 0.027 | 0.006 | 0.077 | 0.009 | 0.002 | 0.029 |
| Seychelles | \>= 5 | Inpatient | 0.051 | 0.009 | 0.163 | 0.019 | 0.003 | 0.062 |
| Seychelles | \< 5 | Outpatient | 0.005 | -0.001 | 0.016 | 0.002 | 0.000 | 0.007 |
| Seychelles | \>= 5 | Outpatient | 0.010 | -0.001 | 0.035 | 0.004 | 0.000 | 0.014 |
| Sierra Leone | \< 5 | Inpatient | 0.020 | 0.006 | 0.051 | 0.008 | 0.002 | 0.022 |
| Sierra Leone | \>= 5 | Inpatient | 0.040 | 0.008 | 0.118 | 0.015 | 0.002 | 0.049 |
| Sierra Leone | \< 5 | Outpatient | 0.004 | 0.000 | 0.011 | 0.002 | 0.000 | 0.005 |
| Sierra Leone | \>= 5 | Outpatient | 0.008 | -0.001 | 0.025 | 0.003 | 0.000 | 0.012 |
| Singapore | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Singapore | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Singapore | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Singapore | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Slovak Republic | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Slovak Republic | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Slovak Republic | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Slovak Republic | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Slovenia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Slovenia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Slovenia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Slovenia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Solomon Islands | \< 5 | Inpatient | 0.010 | 0.000 | 0.038 | 0.005 | 0.000 | 0.018 |
| Solomon Islands | \>= 5 | Inpatient | 0.018 | 0.001 | 0.077 | 0.009 | 0.000 | 0.036 |
| Solomon Islands | \< 5 | Outpatient | 0.002 | 0.000 | 0.008 | 0.001 | 0.000 | 0.004 |
| Solomon Islands | \>= 5 | Outpatient | 0.003 | 0.000 | 0.016 | 0.002 | 0.000 | 0.007 |
| Somalia | \< 5 | Inpatient | 0.024 | 0.004 | 0.077 | 0.015 | 0.002 | 0.048 |
| Somalia | \>= 5 | Inpatient | 0.044 | 0.007 | 0.151 | 0.028 | 0.004 | 0.098 |
| Somalia | \< 5 | Outpatient | 0.004 | 0.000 | 0.016 | 0.003 | 0.000 | 0.011 |
| Somalia | \>= 5 | Outpatient | 0.008 | -0.001 | 0.033 | 0.005 | -0.001 | 0.021 |
| South Africa | \< 5 | Inpatient | 0.013 | 0.003 | 0.035 | 0.007 | 0.001 | 0.022 |
| South Africa | \>= 5 | Inpatient | 0.025 | 0.004 | 0.082 | 0.014 | 0.002 | 0.048 |
| South Africa | \< 5 | Outpatient | 0.002 | 0.000 | 0.007 | 0.001 | 0.000 | 0.005 |
| South Africa | \>= 5 | Outpatient | 0.005 | -0.001 | 0.016 | 0.003 | 0.000 | 0.010 |
| South Sudan | \< 5 | Inpatient | 0.020 | 0.006 | 0.051 | 0.013 | 0.002 | 0.040 |
| South Sudan | \>= 5 | Inpatient | 0.040 | 0.008 | 0.118 | 0.026 | 0.004 | 0.085 |
| South Sudan | \< 5 | Outpatient | 0.004 | 0.000 | 0.011 | 0.003 | 0.000 | 0.009 |
| South Sudan | \>= 5 | Outpatient | 0.008 | -0.001 | 0.025 | 0.005 | -0.001 | 0.019 |
| Spain | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Spain | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Spain | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Spain | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Sri Lanka | \< 5 | Inpatient | 0.018 | 0.004 | 0.049 | 0.011 | 0.002 | 0.034 |
| Sri Lanka | \>= 5 | Inpatient | 0.036 | 0.006 | 0.112 | 0.023 | 0.003 | 0.077 |
| Sri Lanka | \< 5 | Outpatient | 0.003 | 0.000 | 0.011 | 0.002 | 0.000 | 0.008 |
| Sri Lanka | \>= 5 | Outpatient | 0.007 | -0.001 | 0.024 | 0.004 | 0.000 | 0.016 |
| St. Kitts and Nevis | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| St. Kitts and Nevis | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| St. Kitts and Nevis | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| St. Kitts and Nevis | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| St. Lucia | \< 5 | Inpatient | 0.020 | 0.003 | 0.066 | 0.013 | 0.001 | 0.050 |
| St. Lucia | \>= 5 | Inpatient | 0.038 | 0.004 | 0.137 | 0.025 | 0.002 | 0.100 |
| St. Lucia | \< 5 | Outpatient | 0.004 | 0.000 | 0.013 | 0.002 | 0.000 | 0.010 |
| St. Lucia | \>= 5 | Outpatient | 0.007 | -0.001 | 0.028 | 0.005 | 0.000 | 0.020 |
| St. Vincent and the Grenadines | \< 5 | Inpatient | 0.020 | 0.003 | 0.066 | 0.013 | 0.001 | 0.050 |
| St. Vincent and the Grenadines | \>= 5 | Inpatient | 0.038 | 0.004 | 0.137 | 0.025 | 0.002 | 0.100 |
| St. Vincent and the Grenadines | \< 5 | Outpatient | 0.004 | 0.000 | 0.013 | 0.002 | 0.000 | 0.010 |
| St. Vincent and the Grenadines | \>= 5 | Outpatient | 0.007 | -0.001 | 0.028 | 0.005 | 0.000 | 0.020 |
| Sudan | \< 5 | Inpatient | 0.024 | 0.004 | 0.077 | 0.009 | 0.001 | 0.031 |
| Sudan | \>= 5 | Inpatient | 0.044 | 0.007 | 0.151 | 0.017 | 0.002 | 0.064 |
| Sudan | \< 5 | Outpatient | 0.004 | 0.000 | 0.016 | 0.002 | 0.000 | 0.007 |
| Sudan | \>= 5 | Outpatient | 0.008 | -0.001 | 0.033 | 0.003 | 0.000 | 0.014 |
| Suriname | \< 5 | Inpatient | 0.020 | 0.003 | 0.066 | 0.013 | 0.001 | 0.050 |
| Suriname | \>= 5 | Inpatient | 0.038 | 0.004 | 0.137 | 0.025 | 0.002 | 0.100 |
| Suriname | \< 5 | Outpatient | 0.004 | 0.000 | 0.013 | 0.002 | 0.000 | 0.010 |
| Suriname | \>= 5 | Outpatient | 0.007 | -0.001 | 0.028 | 0.005 | 0.000 | 0.020 |
| Sweden | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Sweden | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Sweden | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Sweden | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Switzerland | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Switzerland | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Switzerland | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Switzerland | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Syrian Arab Republic | \< 5 | Inpatient | 0.024 | 0.004 | 0.077 | 0.015 | 0.002 | 0.048 |
| Syrian Arab Republic | \>= 5 | Inpatient | 0.044 | 0.007 | 0.151 | 0.028 | 0.004 | 0.098 |
| Syrian Arab Republic | \< 5 | Outpatient | 0.004 | 0.000 | 0.016 | 0.003 | 0.000 | 0.011 |
| Syrian Arab Republic | \>= 5 | Outpatient | 0.008 | -0.001 | 0.033 | 0.005 | -0.001 | 0.021 |
| Tajikistan | \< 5 | Inpatient | 0.019 | 0.002 | 0.062 | 0.007 | 0.001 | 0.026 |
| Tajikistan | \>= 5 | Inpatient | 0.035 | 0.004 | 0.123 | 0.013 | 0.001 | 0.050 |
| Tajikistan | \< 5 | Outpatient | 0.004 | 0.000 | 0.013 | 0.001 | 0.000 | 0.006 |
| Tajikistan | \>= 5 | Outpatient | 0.007 | -0.001 | 0.025 | 0.003 | 0.000 | 0.011 |
| Tanzania | \< 5 | Inpatient | 0.033 | 0.009 | 0.081 | 0.013 | 0.003 | 0.034 |
| Tanzania | \>= 5 | Inpatient | 0.063 | 0.012 | 0.177 | 0.025 | 0.004 | 0.075 |
| Tanzania | \< 5 | Outpatient | 0.007 | -0.001 | 0.020 | 0.003 | 0.000 | 0.009 |
| Tanzania | \>= 5 | Outpatient | 0.013 | -0.001 | 0.044 | 0.005 | 0.000 | 0.020 |
| Thailand | \< 5 | Inpatient | 0.019 | 0.004 | 0.053 | 0.012 | 0.002 | 0.037 |
| Thailand | \>= 5 | Inpatient | 0.036 | 0.005 | 0.115 | 0.023 | 0.003 | 0.080 |
| Thailand | \< 5 | Outpatient | 0.003 | 0.000 | 0.011 | 0.002 | 0.000 | 0.008 |
| Thailand | \>= 5 | Outpatient | 0.007 | -0.001 | 0.025 | 0.004 | 0.000 | 0.017 |
| Timor-Leste | \< 5 | Inpatient | 0.018 | 0.004 | 0.049 | 0.011 | 0.002 | 0.034 |
| Timor-Leste | \>= 5 | Inpatient | 0.036 | 0.006 | 0.112 | 0.023 | 0.003 | 0.077 |
| Timor-Leste | \< 5 | Outpatient | 0.003 | 0.000 | 0.011 | 0.002 | 0.000 | 0.008 |
| Timor-Leste | \>= 5 | Outpatient | 0.007 | -0.001 | 0.024 | 0.004 | 0.000 | 0.016 |
| Togo | \< 5 | Inpatient | 0.020 | 0.006 | 0.051 | 0.008 | 0.002 | 0.023 |
| Togo | \>= 5 | Inpatient | 0.040 | 0.008 | 0.118 | 0.016 | 0.003 | 0.050 |
| Togo | \< 5 | Outpatient | 0.004 | 0.000 | 0.011 | 0.002 | 0.000 | 0.006 |
| Togo | \>= 5 | Outpatient | 0.008 | -0.001 | 0.025 | 0.003 | 0.000 | 0.012 |
| Tonga | \< 5 | Inpatient | 0.013 | 0.000 | 0.044 | 0.008 | 0.000 | 0.030 |
| Tonga | \>= 5 | Inpatient | 0.024 | 0.001 | 0.090 | 0.015 | 0.001 | 0.058 |
| Tonga | \< 5 | Outpatient | 0.002 | 0.000 | 0.009 | 0.001 | 0.000 | 0.006 |
| Tonga | \>= 5 | Outpatient | 0.004 | 0.000 | 0.018 | 0.003 | 0.000 | 0.012 |
| Trinidad and Tobago | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Trinidad and Tobago | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Trinidad and Tobago | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Trinidad and Tobago | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Tunisia | \< 5 | Inpatient | 0.027 | 0.005 | 0.085 | 0.016 | 0.003 | 0.050 |
| Tunisia | \>= 5 | Inpatient | 0.050 | 0.008 | 0.163 | 0.031 | 0.004 | 0.103 |
| Tunisia | \< 5 | Outpatient | 0.005 | 0.000 | 0.018 | 0.003 | 0.000 | 0.011 |
| Tunisia | \>= 5 | Outpatient | 0.010 | -0.001 | 0.036 | 0.006 | -0.001 | 0.022 |
| Türkiye | \< 5 | Inpatient | 0.017 | 0.002 | 0.065 | 0.011 | 0.001 | 0.043 |
| Türkiye | \>= 5 | Inpatient | 0.030 | 0.003 | 0.113 | 0.020 | 0.002 | 0.080 |
| Türkiye | \< 5 | Outpatient | 0.003 | 0.000 | 0.013 | 0.002 | 0.000 | 0.009 |
| Türkiye | \>= 5 | Outpatient | 0.006 | -0.001 | 0.023 | 0.004 | 0.000 | 0.016 |
| Turkmenistan | \< 5 | Inpatient | 0.019 | 0.002 | 0.067 | 0.007 | 0.001 | 0.028 |
| Turkmenistan | \>= 5 | Inpatient | 0.035 | 0.003 | 0.128 | 0.013 | 0.001 | 0.053 |
| Turkmenistan | \< 5 | Outpatient | 0.004 | 0.000 | 0.014 | 0.001 | 0.000 | 0.006 |
| Turkmenistan | \>= 5 | Outpatient | 0.007 | -0.001 | 0.027 | 0.003 | 0.000 | 0.012 |
| Tuvalu | \< 5 | Inpatient | 0.013 | 0.000 | 0.044 | 0.008 | 0.000 | 0.030 |
| Tuvalu | \>= 5 | Inpatient | 0.024 | 0.001 | 0.090 | 0.015 | 0.001 | 0.058 |
| Tuvalu | \< 5 | Outpatient | 0.002 | 0.000 | 0.009 | 0.001 | 0.000 | 0.006 |
| Tuvalu | \>= 5 | Outpatient | 0.004 | 0.000 | 0.018 | 0.003 | 0.000 | 0.012 |
| Uganda | \< 5 | Inpatient | 0.020 | 0.006 | 0.051 | 0.008 | 0.002 | 0.022 |
| Uganda | \>= 5 | Inpatient | 0.040 | 0.008 | 0.118 | 0.016 | 0.002 | 0.049 |
| Uganda | \< 5 | Outpatient | 0.004 | 0.000 | 0.011 | 0.002 | 0.000 | 0.006 |
| Uganda | \>= 5 | Outpatient | 0.008 | -0.001 | 0.025 | 0.003 | 0.000 | 0.012 |
| Ukraine | \< 5 | Inpatient | 0.019 | 0.002 | 0.062 | 0.012 | 0.001 | 0.043 |
| Ukraine | \>= 5 | Inpatient | 0.035 | 0.004 | 0.123 | 0.022 | 0.002 | 0.084 |
| Ukraine | \< 5 | Outpatient | 0.004 | 0.000 | 0.013 | 0.002 | 0.000 | 0.009 |
| Ukraine | \>= 5 | Outpatient | 0.007 | -0.001 | 0.025 | 0.004 | 0.000 | 0.017 |
| United Arab Emirates | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| United Arab Emirates | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| United Arab Emirates | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| United Arab Emirates | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| United Kingdom | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| United Kingdom | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| United Kingdom | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| United Kingdom | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| United States | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| United States | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| United States | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| United States | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Uruguay | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Uruguay | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Uruguay | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Uruguay | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Uzbekistan | \< 5 | Inpatient | 0.019 | 0.002 | 0.062 | 0.007 | 0.001 | 0.026 |
| Uzbekistan | \>= 5 | Inpatient | 0.035 | 0.004 | 0.123 | 0.014 | 0.001 | 0.050 |
| Uzbekistan | \< 5 | Outpatient | 0.004 | 0.000 | 0.013 | 0.001 | 0.000 | 0.006 |
| Uzbekistan | \>= 5 | Outpatient | 0.007 | -0.001 | 0.025 | 0.003 | 0.000 | 0.011 |
| Vanuatu | \< 5 | Inpatient | 0.010 | 0.000 | 0.038 | 0.006 | 0.000 | 0.025 |
| Vanuatu | \>= 5 | Inpatient | 0.018 | 0.001 | 0.077 | 0.011 | 0.000 | 0.048 |
| Vanuatu | \< 5 | Outpatient | 0.002 | 0.000 | 0.008 | 0.001 | 0.000 | 0.005 |
| Vanuatu | \>= 5 | Outpatient | 0.003 | 0.000 | 0.016 | 0.002 | 0.000 | 0.010 |
| Venezuela | \< 5 | Inpatient | 0.014 | 0.002 | 0.043 | 0.013 | 0.001 | 0.046 |
| Venezuela | \>= 5 | Inpatient | 0.028 | 0.004 | 0.094 | 0.024 | 0.002 | 0.092 |
| Venezuela | \< 5 | Outpatient | 0.003 | 0.000 | 0.009 | 0.002 | 0.000 | 0.010 |
| Venezuela | \>= 5 | Outpatient | 0.005 | -0.001 | 0.020 | 0.004 | 0.000 | 0.019 |
| Vietnam | \< 5 | Inpatient | 0.010 | 0.000 | 0.038 | 0.006 | 0.000 | 0.025 |
| Vietnam | \>= 5 | Inpatient | 0.018 | 0.001 | 0.077 | 0.011 | 0.000 | 0.048 |
| Vietnam | \< 5 | Outpatient | 0.002 | 0.000 | 0.008 | 0.001 | 0.000 | 0.005 |
| Vietnam | \>= 5 | Outpatient | 0.003 | 0.000 | 0.016 | 0.002 | 0.000 | 0.010 |
| Yemen, Rep. | \< 5 | Inpatient | 0.024 | 0.004 | 0.077 | 0.011 | 0.002 | 0.033 |
| Yemen, Rep. | \>= 5 | Inpatient | 0.044 | 0.007 | 0.151 | 0.020 | 0.003 | 0.071 |
| Yemen, Rep. | \< 5 | Outpatient | 0.004 | 0.000 | 0.016 | 0.002 | 0.000 | 0.008 |
| Yemen, Rep. | \>= 5 | Outpatient | 0.008 | -0.001 | 0.033 | 0.004 | 0.000 | 0.016 |
| Zambia | \< 5 | Inpatient | 0.028 | 0.008 | 0.076 | 0.011 | 0.002 | 0.030 |
| Zambia | \>= 5 | Inpatient | 0.054 | 0.011 | 0.155 | 0.021 | 0.004 | 0.064 |
| Zambia | \< 5 | Outpatient | 0.005 | -0.001 | 0.016 | 0.002 | 0.000 | 0.007 |
| Zambia | \>= 5 | Outpatient | 0.010 | -0.001 | 0.035 | 0.004 | 0.000 | 0.015 |
| Zimbabwe | \< 5 | Inpatient | 0.028 | 0.008 | 0.076 | 0.011 | 0.002 | 0.030 |
| Zimbabwe | \>= 5 | Inpatient | 0.054 | 0.011 | 0.155 | 0.021 | 0.004 | 0.064 |
| Zimbabwe | \< 5 | Outpatient | 0.005 | -0.001 | 0.016 | 0.002 | 0.000 | 0.007 |
| Zimbabwe | \>= 5 | Outpatient | 0.010 | -0.001 | 0.035 | 0.004 | 0.000 | 0.015 |

Estimated `r params$Pathogen` proportion by country, 2010 vs 2020

# Session info

``` r
sessioninfo::session_info()
```

    ## Warning in system2("quarto", "-V", stdout = TRUE, env = paste0("TMPDIR=", : running command '"quarto"
    ## TMPDIR=C:/Users/fbbu6966/AppData/Local/Temp/RtmpKSA1q2/file3d4c6e7f7ef3 -V' had status 1

    ## ─ Session info ─────────────────────────────────────────────────────────────────────────────────────────────────────────
    ##  setting  value
    ##  version  R version 4.5.0 (2025-04-11 ucrt)
    ##  os       Windows 10 x64 (build 19045)
    ##  system   x86_64, mingw32
    ##  ui       RStudio
    ##  language (EN)
    ##  collate  English_United States.utf8
    ##  ctype    English_United States.utf8
    ##  tz       Europe/Brussels
    ##  date     2025-11-27
    ##  rstudio  2025.05.0+496 Mariposa Orchid (desktop)
    ##  pandoc   3.4 @ C:/Users/fbbu6966/AppData/Local/Programs/RStudio/resources/app/bin/quarto/bin/tools/ (via rmarkdown)
    ##  quarto   ERROR: Unknown command "TMPDIR=C:/Users/fbbu6966/AppData/Local/Temp/RtmpKSA1q2/file3d4c6e7f7ef3". Did you mean command "update"? @ C:\\Users\\fbbu6966\\AppData\\Local\\Programs\\RStudio\\RESOUR~1\\app\\bin\\quarto\\bin\\quarto.exe
    ## 
    ## ─ Packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ##  ! package        * version    date (UTC) lib source
    ##    abind            1.4-8      2024-09-12 [1] CRAN (R 4.5.0)
    ##    backports        1.5.0      2024-05-23 [1] CRAN (R 4.5.0)
    ##    base64enc        0.1-3      2015-07-28 [1] CRAN (R 4.5.0)
    ##    bayesplot        1.12.0     2025-04-10 [1] CRAN (R 4.5.0)
    ##    bd             * 0.0.14     2025-04-14 [1] Github (brechtdv/bd@652191c)
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
    ##    FERG2          * 0.0.5      2025-08-07 [1] Github (brechtdv/FERG2@c2d4ac1)
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
    ##    lubridate      * 1.9.4      2024-12-08 [1] CRAN (R 4.5.0)
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
    ##    pillar           1.11.0     2025-07-04 [1] CRAN (R 4.5.1)
    ##    pkgbuild         1.4.7      2025-03-24 [1] CRAN (R 4.5.0)
    ##    pkgconfig        2.0.3      2019-09-22 [1] CRAN (R 4.5.0)
    ##    plyr             1.8.9      2023-10-02 [1] CRAN (R 4.5.0)
    ##    polspline        1.1.25     2024-05-10 [1] CRAN (R 4.5.0)
    ##    posterior        1.6.1      2025-02-27 [1] CRAN (R 4.5.0)
    ##    processx         3.8.6      2025-02-21 [1] CRAN (R 4.5.0)
    ##    proxy            0.4-27     2022-06-09 [1] CRAN (R 4.5.0)
    ##    ps               1.9.1      2025-04-12 [1] CRAN (R 4.5.0)
    ##    purrr            1.0.4      2025-02-05 [1] CRAN (R 4.5.0)
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
    ##    stringr        * 1.5.1      2023-11-14 [1] CRAN (R 4.5.0)
    ##    survival         3.8-3      2024-12-17 [1] CRAN (R 4.5.0)
    ##    svglite          2.1.3      2023-12-08 [1] CRAN (R 4.5.0)
    ##    systemfonts      1.2.2      2025-04-04 [1] CRAN (R 4.5.0)
    ##    tensorA          0.36.2.1   2023-12-13 [1] CRAN (R 4.5.0)
    ##    TH.data          1.1-3      2025-01-17 [1] CRAN (R 4.5.0)
    ##    tibble           3.3.0      2025-06-08 [1] CRAN (R 4.5.1)
    ##    tidyr          * 1.3.1      2024-01-24 [1] CRAN (R 4.5.0)
    ##    tidyselect       1.2.1      2024-03-11 [1] CRAN (R 4.5.0)
    ##    timechange       0.3.0      2024-01-18 [1] CRAN (R 4.5.0)
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
    ##  [1] C:/Users/fbbu6966/AppData/Local/Programs/R/R-4.5.0/library
    ## 
    ##  * ── Packages attached to the search path.
    ##  D ── DLL MD5 mismatch, broken installation.
    ## 
    ## ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

``` r
# Save dataset for report created for expert to receive feedback
# save(all_cnt_rt, file="./00-Report_FB/all_cnt_rt.Rdata")
# save(all_glb_prop, file="./00-Report_FB/all_glb_prop.Rdata")
# save(all_reg_prop, file="./00-Report_FB/all_reg_prop.Rdata")
# save(all_reg_rt, file="./00-Report_FB/all_reg_rt.Rdata")
# save(all_sub_nr, file="./00-Report_FB/all_sub_nr.Rdata")
# save(all_sub_rt, file="./00-Report_FB/all_sub_rt.Rdata")
```
