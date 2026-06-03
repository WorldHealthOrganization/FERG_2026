Global proportion of Rotavirus (Diarrheal) • ASYMPT+ROTA • Estimate
attributable fraction
================
fbbu6966
2025-09-28

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

    ## Warning: There were 27 divergent transitions after warmup. Increasing adapt_delta above 0.95 may help. See
    ## http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup

    ##  Family: gaussian 
    ##   Links: mu = identity; sigma = identity 
    ## Formula: yi | se(sei) ~ 1 + YEAR + SYNDROMTYPE + AGE + REFERENCE + COVERAGE + STRAIN + (1 | REG2) + (1 | REG2:SUB2) + (1 | REG2:SUB2:COUNTRY) + (1 | REG2:SUB2:COUNTRY:ID) + (1 | REG2:SUB2:COUNTRY:ID:DTP_ID) 
    ##    Data: subset(es, as.integer(FLAG) == 1) (Number of observations: 589) 
    ##   Draws: 5 chains, each with iter = 5000; warmup = 3000; thin = 1;
    ##          total post-warmup draws = 10000
    ## 
    ## Multilevel Hyperparameters:
    ## ~REG2 (Number of levels: 6) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.10      0.10     0.00     0.36 1.00     6256     6627
    ## 
    ## ~REG2:SUB2 (Number of levels: 12) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.08      0.07     0.00     0.25 1.00     5794     6048
    ## 
    ## ~REG2:SUB2:COUNTRY (Number of levels: 53) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.13      0.09     0.01     0.34 1.00     1700     2891
    ## 
    ## ~REG2:SUB2:COUNTRY:ID (Number of levels: 177) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.63      0.06     0.52     0.75 1.00     2932     5389
    ## 
    ## ~REG2:SUB2:COUNTRY:ID:DTP_ID (Number of levels: 589) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.60      0.03     0.55     0.66 1.00     3159     5386
    ## 
    ## Regression Coefficients:
    ##                       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept               -25.31     23.81   -71.98    21.30 1.00     4999     6164
    ## YEAR                      0.01      0.01    -0.01     0.03 1.00     5001     6196
    ## SYNDROMTYPEInpatient      2.35      0.10     2.15     2.56 1.00     4447     6426
    ## SYNDROMTYPEOutpatient     1.53      0.09     1.37     1.71 1.00     4244     5628
    ## AGEAgebelow5              1.26      0.13     1.00     1.51 1.00     4434     6486
    ## AGEMixedages              0.67      0.18     0.32     1.02 1.00     5131     6375
    ## REFERENCEOther           -0.24      0.14    -0.52     0.03 1.00     4540     5020
    ## COVERAGE                 -0.00      0.00    -0.01     0.00 1.00     5827     6159
    ## STRAINGroup_A             0.27      0.15    -0.02     0.57 1.00     5420     6193
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

    ## [1] 6996

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

    ## [1] 10164

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

    ## [1] 528

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
    ##  $ VAL_MEAN     : num  0.104 0.28 0.105 0.282 0.106 ...
    ##  $ VAL_LWR      : num  0.0638 0.1976 0.065 0.2012 0.0662 ...
    ##  $ VAL_UPR      : num  0.156 0.374 0.156 0.374 0.155 ...
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
    ##  $ VAL_MEAN     : num  0.104 0.28 0.105 0.282 0.106 ...
    ##  $ VAL_LWR      : num  0.0638 0.1976 0.065 0.2012 0.0662 ...
    ##  $ VAL_UPR      : num  0.156 0.374 0.156 0.374 0.155 ...
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

\[1\] 0.00 0.01 0.02 0.03 0.04 0.05 0.06

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

\[1\] 0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35

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

\[1\] 0.00 0.05 0.10 0.15 0.20

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

\[1\] 0.00 0.05 0.10 0.15

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

\[1\] 0.00 0.01 0.02 0.03 0.04 0.05 0.06 0.07

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

\[1\] 0.0 0.1 0.2 0.3 0.4

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

\[1\] 0.00 0.05 0.10 0.15 0.20

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
| Afghanistan | \< 5 | Inpatient | 0.300 | 0.226 | 0.382 | 0.283 | 0.204 | 0.369 |
| Afghanistan | \>= 5 | Inpatient | 0.113 | 0.075 | 0.161 | 0.105 | 0.066 | 0.156 |
| Afghanistan | \< 5 | Outpatient | 0.141 | 0.100 | 0.188 | 0.131 | 0.090 | 0.181 |
| Afghanistan | \>= 5 | Outpatient | 0.046 | 0.030 | 0.067 | 0.043 | 0.027 | 0.065 |
| Albania | \< 5 | Inpatient | 0.293 | 0.208 | 0.383 | 0.252 | 0.160 | 0.355 |
| Albania | \>= 5 | Inpatient | 0.110 | 0.070 | 0.162 | 0.091 | 0.050 | 0.148 |
| Albania | \< 5 | Outpatient | 0.137 | 0.092 | 0.190 | 0.115 | 0.068 | 0.174 |
| Albania | \>= 5 | Outpatient | 0.045 | 0.028 | 0.068 | 0.037 | 0.020 | 0.062 |
| Algeria | \< 5 | Inpatient | 0.309 | 0.250 | 0.376 | 0.331 | 0.258 | 0.410 |
| Algeria | \>= 5 | Inpatient | 0.117 | 0.083 | 0.160 | 0.129 | 0.088 | 0.180 |
| Algeria | \< 5 | Outpatient | 0.146 | 0.114 | 0.184 | 0.159 | 0.119 | 0.207 |
| Algeria | \>= 5 | Outpatient | 0.048 | 0.033 | 0.067 | 0.053 | 0.036 | 0.076 |
| Andorra | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Andorra | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Andorra | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Andorra | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Angola | \< 5 | Inpatient | 0.313 | 0.232 | 0.415 | 0.306 | 0.223 | 0.407 |
| Angola | \>= 5 | Inpatient | 0.120 | 0.077 | 0.183 | 0.116 | 0.073 | 0.179 |
| Angola | \< 5 | Outpatient | 0.149 | 0.104 | 0.211 | 0.144 | 0.099 | 0.205 |
| Angola | \>= 5 | Outpatient | 0.050 | 0.031 | 0.078 | 0.048 | 0.030 | 0.075 |
| Antigua and Barbuda | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Antigua and Barbuda | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Antigua and Barbuda | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Antigua and Barbuda | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Argentina | \< 5 | Inpatient | 0.309 | 0.224 | 0.413 | 0.283 | 0.197 | 0.388 |
| Argentina | \>= 5 | Inpatient | 0.118 | 0.073 | 0.180 | 0.105 | 0.063 | 0.165 |
| Argentina | \< 5 | Outpatient | 0.147 | 0.100 | 0.211 | 0.131 | 0.086 | 0.194 |
| Argentina | \>= 5 | Outpatient | 0.049 | 0.030 | 0.077 | 0.043 | 0.025 | 0.069 |
| Armenia | \< 5 | Inpatient | 0.293 | 0.223 | 0.365 | 0.255 | 0.175 | 0.344 |
| Armenia | \>= 5 | Inpatient | 0.109 | 0.074 | 0.152 | 0.092 | 0.055 | 0.141 |
| Armenia | \< 5 | Outpatient | 0.137 | 0.098 | 0.178 | 0.116 | 0.075 | 0.166 |
| Armenia | \>= 5 | Outpatient | 0.045 | 0.030 | 0.063 | 0.037 | 0.022 | 0.058 |
| Australia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Australia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Australia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Australia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Austria | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Austria | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Austria | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Austria | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Azerbaijan | \< 5 | Inpatient | 0.293 | 0.223 | 0.365 | 0.314 | 0.230 | 0.403 |
| Azerbaijan | \>= 5 | Inpatient | 0.109 | 0.074 | 0.152 | 0.120 | 0.078 | 0.173 |
| Azerbaijan | \< 5 | Outpatient | 0.137 | 0.098 | 0.178 | 0.149 | 0.103 | 0.201 |
| Azerbaijan | \>= 5 | Outpatient | 0.045 | 0.030 | 0.063 | 0.050 | 0.031 | 0.073 |
| Bahamas, The | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahamas, The | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahamas, The | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahamas, The | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahrain | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahrain | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahrain | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahrain | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bangladesh | \< 5 | Inpatient | 0.310 | 0.243 | 0.387 | 0.332 | 0.249 | 0.425 |
| Bangladesh | \>= 5 | Inpatient | 0.118 | 0.081 | 0.167 | 0.129 | 0.085 | 0.189 |
| Bangladesh | \< 5 | Outpatient | 0.146 | 0.110 | 0.192 | 0.160 | 0.113 | 0.216 |
| Bangladesh | \>= 5 | Outpatient | 0.048 | 0.033 | 0.070 | 0.054 | 0.034 | 0.080 |
| Barbados | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Barbados | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Barbados | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Barbados | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belarus | \< 5 | Inpatient | 0.293 | 0.223 | 0.365 | 0.314 | 0.230 | 0.403 |
| Belarus | \>= 5 | Inpatient | 0.109 | 0.074 | 0.152 | 0.120 | 0.078 | 0.173 |
| Belarus | \< 5 | Outpatient | 0.137 | 0.098 | 0.178 | 0.149 | 0.103 | 0.201 |
| Belarus | \>= 5 | Outpatient | 0.045 | 0.030 | 0.063 | 0.050 | 0.031 | 0.073 |
| Belgium | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belgium | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belgium | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belgium | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belize | \< 5 | Inpatient | 0.304 | 0.236 | 0.378 | 0.326 | 0.242 | 0.420 |
| Belize | \>= 5 | Inpatient | 0.115 | 0.078 | 0.160 | 0.126 | 0.082 | 0.185 |
| Belize | \< 5 | Outpatient | 0.143 | 0.106 | 0.186 | 0.156 | 0.110 | 0.214 |
| Belize | \>= 5 | Outpatient | 0.047 | 0.032 | 0.067 | 0.052 | 0.033 | 0.078 |
| Benin | \< 5 | Inpatient | 0.309 | 0.250 | 0.376 | 0.288 | 0.225 | 0.360 |
| Benin | \>= 5 | Inpatient | 0.117 | 0.083 | 0.160 | 0.107 | 0.073 | 0.151 |
| Benin | \< 5 | Outpatient | 0.146 | 0.114 | 0.184 | 0.134 | 0.100 | 0.174 |
| Benin | \>= 5 | Outpatient | 0.048 | 0.033 | 0.067 | 0.044 | 0.029 | 0.063 |
| Bhutan | \< 5 | Inpatient | 0.308 | 0.249 | 0.375 | 0.330 | 0.253 | 0.412 |
| Bhutan | \>= 5 | Inpatient | 0.117 | 0.082 | 0.160 | 0.128 | 0.086 | 0.182 |
| Bhutan | \< 5 | Outpatient | 0.145 | 0.113 | 0.185 | 0.159 | 0.116 | 0.210 |
| Bhutan | \>= 5 | Outpatient | 0.048 | 0.034 | 0.067 | 0.053 | 0.035 | 0.077 |
| Bolivia | \< 5 | Inpatient | 0.252 | 0.173 | 0.349 | 0.280 | 0.193 | 0.381 |
| Bolivia | \>= 5 | Inpatient | 0.091 | 0.054 | 0.145 | 0.104 | 0.062 | 0.162 |
| Bolivia | \< 5 | Outpatient | 0.115 | 0.074 | 0.168 | 0.130 | 0.084 | 0.189 |
| Bolivia | \>= 5 | Outpatient | 0.037 | 0.021 | 0.060 | 0.042 | 0.025 | 0.068 |
| Bosnia and Herzegovina | \< 5 | Inpatient | 0.293 | 0.223 | 0.365 | 0.314 | 0.230 | 0.403 |
| Bosnia and Herzegovina | \>= 5 | Inpatient | 0.109 | 0.074 | 0.152 | 0.120 | 0.078 | 0.173 |
| Bosnia and Herzegovina | \< 5 | Outpatient | 0.137 | 0.098 | 0.178 | 0.149 | 0.103 | 0.201 |
| Bosnia and Herzegovina | \>= 5 | Outpatient | 0.045 | 0.030 | 0.063 | 0.050 | 0.031 | 0.073 |
| Botswana | \< 5 | Inpatient | 0.300 | 0.213 | 0.391 | 0.266 | 0.178 | 0.364 |
| Botswana | \>= 5 | Inpatient | 0.113 | 0.070 | 0.167 | 0.097 | 0.056 | 0.152 |
| Botswana | \< 5 | Outpatient | 0.141 | 0.094 | 0.195 | 0.122 | 0.077 | 0.177 |
| Botswana | \>= 5 | Outpatient | 0.046 | 0.028 | 0.071 | 0.040 | 0.022 | 0.063 |
| Brazil | \< 5 | Inpatient | 0.256 | 0.184 | 0.346 | 0.280 | 0.201 | 0.375 |
| Brazil | \>= 5 | Inpatient | 0.093 | 0.057 | 0.143 | 0.104 | 0.064 | 0.158 |
| Brazil | \< 5 | Outpatient | 0.117 | 0.080 | 0.166 | 0.130 | 0.088 | 0.185 |
| Brazil | \>= 5 | Outpatient | 0.038 | 0.023 | 0.058 | 0.042 | 0.025 | 0.066 |
| Brunei Darussalam | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Brunei Darussalam | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Brunei Darussalam | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Brunei Darussalam | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bulgaria | \< 5 | Inpatient | 0.293 | 0.223 | 0.365 | 0.284 | 0.208 | 0.365 |
| Bulgaria | \>= 5 | Inpatient | 0.109 | 0.074 | 0.152 | 0.105 | 0.068 | 0.152 |
| Bulgaria | \< 5 | Outpatient | 0.137 | 0.098 | 0.178 | 0.132 | 0.092 | 0.177 |
| Bulgaria | \>= 5 | Outpatient | 0.045 | 0.030 | 0.063 | 0.043 | 0.027 | 0.063 |
| Burkina Faso | \< 5 | Inpatient | 0.297 | 0.215 | 0.382 | 0.261 | 0.176 | 0.357 |
| Burkina Faso | \>= 5 | Inpatient | 0.112 | 0.070 | 0.162 | 0.095 | 0.055 | 0.148 |
| Burkina Faso | \< 5 | Outpatient | 0.139 | 0.094 | 0.189 | 0.120 | 0.075 | 0.173 |
| Burkina Faso | \>= 5 | Outpatient | 0.046 | 0.028 | 0.069 | 0.039 | 0.022 | 0.062 |
| Burundi | \< 5 | Inpatient | 0.303 | 0.242 | 0.367 | 0.262 | 0.188 | 0.347 |
| Burundi | \>= 5 | Inpatient | 0.114 | 0.079 | 0.155 | 0.095 | 0.059 | 0.143 |
| Burundi | \< 5 | Outpatient | 0.142 | 0.109 | 0.179 | 0.120 | 0.081 | 0.166 |
| Burundi | \>= 5 | Outpatient | 0.047 | 0.032 | 0.065 | 0.039 | 0.023 | 0.060 |
| Cabo Verde | \< 5 | Inpatient | 0.314 | 0.235 | 0.416 | 0.336 | 0.245 | 0.444 |
| Cabo Verde | \>= 5 | Inpatient | 0.121 | 0.077 | 0.182 | 0.132 | 0.084 | 0.201 |
| Cabo Verde | \< 5 | Outpatient | 0.149 | 0.105 | 0.212 | 0.163 | 0.111 | 0.232 |
| Cabo Verde | \>= 5 | Outpatient | 0.050 | 0.031 | 0.077 | 0.055 | 0.034 | 0.085 |
| Cambodia | \< 5 | Inpatient | 0.292 | 0.202 | 0.382 | 0.313 | 0.212 | 0.417 |
| Cambodia | \>= 5 | Inpatient | 0.109 | 0.066 | 0.161 | 0.120 | 0.070 | 0.182 |
| Cambodia | \< 5 | Outpatient | 0.136 | 0.088 | 0.189 | 0.149 | 0.092 | 0.212 |
| Cambodia | \>= 5 | Outpatient | 0.045 | 0.026 | 0.067 | 0.050 | 0.028 | 0.077 |
| Cameroon | \< 5 | Inpatient | 0.315 | 0.238 | 0.417 | 0.291 | 0.210 | 0.396 |
| Cameroon | \>= 5 | Inpatient | 0.121 | 0.079 | 0.183 | 0.109 | 0.068 | 0.170 |
| Cameroon | \< 5 | Outpatient | 0.150 | 0.107 | 0.211 | 0.136 | 0.093 | 0.199 |
| Cameroon | \>= 5 | Outpatient | 0.050 | 0.032 | 0.077 | 0.045 | 0.028 | 0.071 |
| Canada | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Canada | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Canada | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Canada | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Central African Republic | \< 5 | Inpatient | 0.301 | 0.217 | 0.390 | 0.323 | 0.227 | 0.421 |
| Central African Republic | \>= 5 | Inpatient | 0.114 | 0.071 | 0.167 | 0.125 | 0.076 | 0.185 |
| Central African Republic | \< 5 | Outpatient | 0.142 | 0.096 | 0.195 | 0.155 | 0.101 | 0.216 |
| Central African Republic | \>= 5 | Outpatient | 0.047 | 0.029 | 0.071 | 0.052 | 0.030 | 0.079 |
| Chad | \< 5 | Inpatient | 0.303 | 0.242 | 0.367 | 0.324 | 0.250 | 0.405 |
| Chad | \>= 5 | Inpatient | 0.114 | 0.079 | 0.155 | 0.125 | 0.084 | 0.176 |
| Chad | \< 5 | Outpatient | 0.142 | 0.109 | 0.179 | 0.155 | 0.114 | 0.204 |
| Chad | \>= 5 | Outpatient | 0.047 | 0.032 | 0.065 | 0.052 | 0.034 | 0.074 |
| Chile | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Chile | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Chile | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Chile | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| China | \< 5 | Inpatient | 0.305 | 0.240 | 0.379 | 0.327 | 0.249 | 0.410 |
| China | \>= 5 | Inpatient | 0.115 | 0.081 | 0.159 | 0.126 | 0.085 | 0.178 |
| China | \< 5 | Outpatient | 0.144 | 0.109 | 0.185 | 0.156 | 0.114 | 0.206 |
| China | \>= 5 | Outpatient | 0.047 | 0.033 | 0.066 | 0.052 | 0.035 | 0.075 |
| Colombia | \< 5 | Inpatient | 0.257 | 0.182 | 0.347 | 0.269 | 0.186 | 0.364 |
| Colombia | \>= 5 | Inpatient | 0.093 | 0.057 | 0.141 | 0.098 | 0.058 | 0.154 |
| Colombia | \< 5 | Outpatient | 0.117 | 0.079 | 0.165 | 0.123 | 0.081 | 0.178 |
| Colombia | \>= 5 | Outpatient | 0.038 | 0.022 | 0.058 | 0.040 | 0.023 | 0.063 |
| Comoros | \< 5 | Inpatient | 0.309 | 0.250 | 0.376 | 0.331 | 0.258 | 0.410 |
| Comoros | \>= 5 | Inpatient | 0.117 | 0.083 | 0.160 | 0.129 | 0.088 | 0.180 |
| Comoros | \< 5 | Outpatient | 0.146 | 0.114 | 0.184 | 0.159 | 0.119 | 0.207 |
| Comoros | \>= 5 | Outpatient | 0.048 | 0.033 | 0.067 | 0.053 | 0.036 | 0.076 |
| Congo, Dem. Rep. | \< 5 | Inpatient | 0.303 | 0.242 | 0.367 | 0.306 | 0.239 | 0.377 |
| Congo, Dem. Rep. | \>= 5 | Inpatient | 0.114 | 0.079 | 0.155 | 0.116 | 0.079 | 0.160 |
| Congo, Dem. Rep. | \< 5 | Outpatient | 0.142 | 0.109 | 0.179 | 0.144 | 0.107 | 0.186 |
| Congo, Dem. Rep. | \>= 5 | Outpatient | 0.047 | 0.032 | 0.065 | 0.048 | 0.032 | 0.067 |
| Congo, Rep. | \< 5 | Inpatient | 0.314 | 0.235 | 0.414 | 0.295 | 0.215 | 0.395 |
| Congo, Rep. | \>= 5 | Inpatient | 0.120 | 0.078 | 0.181 | 0.111 | 0.069 | 0.172 |
| Congo, Rep. | \< 5 | Outpatient | 0.149 | 0.106 | 0.210 | 0.138 | 0.094 | 0.198 |
| Congo, Rep. | \>= 5 | Outpatient | 0.050 | 0.031 | 0.076 | 0.045 | 0.028 | 0.073 |
| Cook Islands | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Cook Islands | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Cook Islands | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Cook Islands | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Costa Rica | \< 5 | Inpatient | 0.304 | 0.236 | 0.378 | 0.262 | 0.190 | 0.346 |
| Costa Rica | \>= 5 | Inpatient | 0.115 | 0.078 | 0.160 | 0.095 | 0.059 | 0.143 |
| Costa Rica | \< 5 | Outpatient | 0.143 | 0.106 | 0.186 | 0.120 | 0.083 | 0.167 |
| Costa Rica | \>= 5 | Outpatient | 0.047 | 0.032 | 0.067 | 0.039 | 0.024 | 0.059 |
| Côte d’Ivoire | \< 5 | Inpatient | 0.309 | 0.250 | 0.376 | 0.284 | 0.220 | 0.359 |
| Côte d’Ivoire | \>= 5 | Inpatient | 0.117 | 0.083 | 0.160 | 0.105 | 0.071 | 0.150 |
| Côte d’Ivoire | \< 5 | Outpatient | 0.146 | 0.114 | 0.184 | 0.132 | 0.097 | 0.173 |
| Côte d’Ivoire | \>= 5 | Outpatient | 0.048 | 0.033 | 0.067 | 0.043 | 0.028 | 0.063 |
| Croatia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Croatia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Croatia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Croatia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Cuba | \< 5 | Inpatient | 0.302 | 0.211 | 0.399 | 0.323 | 0.216 | 0.438 |
| Cuba | \>= 5 | Inpatient | 0.114 | 0.069 | 0.173 | 0.126 | 0.072 | 0.197 |
| Cuba | \< 5 | Outpatient | 0.142 | 0.094 | 0.201 | 0.155 | 0.096 | 0.227 |
| Cuba | \>= 5 | Outpatient | 0.047 | 0.028 | 0.072 | 0.052 | 0.029 | 0.084 |
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
| Djibouti | \< 5 | Inpatient | 0.295 | 0.228 | 0.364 | 0.270 | 0.197 | 0.351 |
| Djibouti | \>= 5 | Inpatient | 0.111 | 0.075 | 0.153 | 0.099 | 0.063 | 0.144 |
| Djibouti | \< 5 | Outpatient | 0.138 | 0.102 | 0.177 | 0.124 | 0.086 | 0.169 |
| Djibouti | \>= 5 | Outpatient | 0.045 | 0.031 | 0.063 | 0.040 | 0.025 | 0.060 |
| Dominica | \< 5 | Inpatient | 0.304 | 0.236 | 0.378 | 0.326 | 0.242 | 0.420 |
| Dominica | \>= 5 | Inpatient | 0.115 | 0.078 | 0.160 | 0.126 | 0.082 | 0.185 |
| Dominica | \< 5 | Outpatient | 0.143 | 0.106 | 0.186 | 0.156 | 0.110 | 0.214 |
| Dominica | \>= 5 | Outpatient | 0.047 | 0.032 | 0.067 | 0.052 | 0.033 | 0.078 |
| Dominican Republic | \< 5 | Inpatient | 0.304 | 0.236 | 0.378 | 0.274 | 0.206 | 0.354 |
| Dominican Republic | \>= 5 | Inpatient | 0.115 | 0.078 | 0.160 | 0.101 | 0.065 | 0.147 |
| Dominican Republic | \< 5 | Outpatient | 0.143 | 0.106 | 0.186 | 0.126 | 0.090 | 0.170 |
| Dominican Republic | \>= 5 | Outpatient | 0.047 | 0.032 | 0.067 | 0.041 | 0.026 | 0.061 |
| Ecuador | \< 5 | Inpatient | 0.245 | 0.179 | 0.322 | 0.276 | 0.208 | 0.355 |
| Ecuador | \>= 5 | Inpatient | 0.087 | 0.055 | 0.130 | 0.101 | 0.066 | 0.148 |
| Ecuador | \< 5 | Outpatient | 0.110 | 0.077 | 0.152 | 0.127 | 0.091 | 0.171 |
| Ecuador | \>= 5 | Outpatient | 0.035 | 0.022 | 0.053 | 0.041 | 0.027 | 0.061 |
| Egypt, Arab Rep. | \< 5 | Inpatient | 0.276 | 0.190 | 0.358 | 0.296 | 0.195 | 0.395 |
| Egypt, Arab Rep. | \>= 5 | Inpatient | 0.101 | 0.061 | 0.147 | 0.112 | 0.064 | 0.169 |
| Egypt, Arab Rep. | \< 5 | Outpatient | 0.127 | 0.083 | 0.172 | 0.139 | 0.086 | 0.196 |
| Egypt, Arab Rep. | \>= 5 | Outpatient | 0.041 | 0.025 | 0.060 | 0.046 | 0.026 | 0.071 |
| El Salvador | \< 5 | Inpatient | 0.247 | 0.182 | 0.323 | 0.274 | 0.206 | 0.354 |
| El Salvador | \>= 5 | Inpatient | 0.088 | 0.056 | 0.130 | 0.101 | 0.065 | 0.147 |
| El Salvador | \< 5 | Outpatient | 0.112 | 0.079 | 0.153 | 0.126 | 0.090 | 0.170 |
| El Salvador | \>= 5 | Outpatient | 0.036 | 0.022 | 0.053 | 0.041 | 0.026 | 0.061 |
| Equatorial Guinea | \< 5 | Inpatient | 0.302 | 0.235 | 0.370 | 0.323 | 0.246 | 0.407 |
| Equatorial Guinea | \>= 5 | Inpatient | 0.114 | 0.078 | 0.156 | 0.125 | 0.082 | 0.176 |
| Equatorial Guinea | \< 5 | Outpatient | 0.142 | 0.106 | 0.181 | 0.154 | 0.112 | 0.204 |
| Equatorial Guinea | \>= 5 | Outpatient | 0.047 | 0.032 | 0.065 | 0.052 | 0.034 | 0.075 |
| Eritrea | \< 5 | Inpatient | 0.303 | 0.242 | 0.367 | 0.262 | 0.188 | 0.347 |
| Eritrea | \>= 5 | Inpatient | 0.114 | 0.079 | 0.155 | 0.095 | 0.059 | 0.143 |
| Eritrea | \< 5 | Outpatient | 0.142 | 0.109 | 0.179 | 0.120 | 0.081 | 0.166 |
| Eritrea | \>= 5 | Outpatient | 0.047 | 0.032 | 0.065 | 0.039 | 0.023 | 0.060 |
| Estonia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Estonia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Estonia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Estonia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Eswatini | \< 5 | Inpatient | 0.309 | 0.250 | 0.376 | 0.288 | 0.225 | 0.360 |
| Eswatini | \>= 5 | Inpatient | 0.117 | 0.083 | 0.160 | 0.107 | 0.073 | 0.151 |
| Eswatini | \< 5 | Outpatient | 0.146 | 0.114 | 0.184 | 0.134 | 0.100 | 0.174 |
| Eswatini | \>= 5 | Outpatient | 0.048 | 0.033 | 0.067 | 0.044 | 0.029 | 0.063 |
| Ethiopia | \< 5 | Inpatient | 0.303 | 0.242 | 0.367 | 0.278 | 0.210 | 0.353 |
| Ethiopia | \>= 5 | Inpatient | 0.114 | 0.079 | 0.155 | 0.102 | 0.067 | 0.146 |
| Ethiopia | \< 5 | Outpatient | 0.142 | 0.109 | 0.179 | 0.128 | 0.092 | 0.170 |
| Ethiopia | \>= 5 | Outpatient | 0.047 | 0.032 | 0.065 | 0.042 | 0.027 | 0.061 |
| Fiji | \< 5 | Inpatient | 0.306 | 0.242 | 0.378 | 0.264 | 0.186 | 0.357 |
| Fiji | \>= 5 | Inpatient | 0.116 | 0.081 | 0.159 | 0.096 | 0.059 | 0.147 |
| Fiji | \< 5 | Outpatient | 0.144 | 0.109 | 0.186 | 0.121 | 0.080 | 0.174 |
| Fiji | \>= 5 | Outpatient | 0.048 | 0.033 | 0.066 | 0.039 | 0.023 | 0.062 |
| Finland | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Finland | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Finland | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Finland | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| France | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| France | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| France | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| France | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Gabon | \< 5 | Inpatient | 0.311 | 0.233 | 0.410 | 0.333 | 0.244 | 0.439 |
| Gabon | \>= 5 | Inpatient | 0.119 | 0.077 | 0.179 | 0.130 | 0.082 | 0.198 |
| Gabon | \< 5 | Outpatient | 0.148 | 0.105 | 0.205 | 0.161 | 0.110 | 0.226 |
| Gabon | \>= 5 | Outpatient | 0.049 | 0.031 | 0.075 | 0.054 | 0.033 | 0.084 |
| Gambia, The | \< 5 | Inpatient | 0.301 | 0.222 | 0.386 | 0.267 | 0.183 | 0.363 |
| Gambia, The | \>= 5 | Inpatient | 0.113 | 0.072 | 0.163 | 0.098 | 0.057 | 0.151 |
| Gambia, The | \< 5 | Outpatient | 0.141 | 0.099 | 0.194 | 0.122 | 0.079 | 0.178 |
| Gambia, The | \>= 5 | Outpatient | 0.047 | 0.029 | 0.069 | 0.040 | 0.023 | 0.064 |
| Georgia | \< 5 | Inpatient | 0.293 | 0.223 | 0.365 | 0.266 | 0.188 | 0.350 |
| Georgia | \>= 5 | Inpatient | 0.109 | 0.074 | 0.152 | 0.097 | 0.060 | 0.143 |
| Georgia | \< 5 | Outpatient | 0.137 | 0.098 | 0.178 | 0.122 | 0.081 | 0.169 |
| Georgia | \>= 5 | Outpatient | 0.045 | 0.030 | 0.063 | 0.039 | 0.024 | 0.060 |
| Germany | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Germany | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Germany | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Germany | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Ghana | \< 5 | Inpatient | 0.296 | 0.211 | 0.378 | 0.260 | 0.175 | 0.351 |
| Ghana | \>= 5 | Inpatient | 0.111 | 0.069 | 0.160 | 0.095 | 0.055 | 0.145 |
| Ghana | \< 5 | Outpatient | 0.139 | 0.092 | 0.186 | 0.119 | 0.074 | 0.169 |
| Ghana | \>= 5 | Outpatient | 0.046 | 0.028 | 0.067 | 0.038 | 0.022 | 0.061 |
| Greece | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Greece | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Greece | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Greece | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Grenada | \< 5 | Inpatient | 0.304 | 0.236 | 0.378 | 0.326 | 0.242 | 0.420 |
| Grenada | \>= 5 | Inpatient | 0.115 | 0.078 | 0.160 | 0.126 | 0.082 | 0.185 |
| Grenada | \< 5 | Outpatient | 0.143 | 0.106 | 0.186 | 0.156 | 0.110 | 0.214 |
| Grenada | \>= 5 | Outpatient | 0.047 | 0.032 | 0.067 | 0.052 | 0.033 | 0.078 |
| Guatemala | \< 5 | Inpatient | 0.277 | 0.204 | 0.357 | 0.268 | 0.187 | 0.361 |
| Guatemala | \>= 5 | Inpatient | 0.102 | 0.066 | 0.148 | 0.098 | 0.059 | 0.150 |
| Guatemala | \< 5 | Outpatient | 0.128 | 0.090 | 0.174 | 0.123 | 0.081 | 0.176 |
| Guatemala | \>= 5 | Outpatient | 0.042 | 0.027 | 0.062 | 0.040 | 0.024 | 0.062 |
| Guinea | \< 5 | Inpatient | 0.309 | 0.250 | 0.376 | 0.331 | 0.258 | 0.410 |
| Guinea | \>= 5 | Inpatient | 0.117 | 0.083 | 0.160 | 0.129 | 0.088 | 0.180 |
| Guinea | \< 5 | Outpatient | 0.146 | 0.114 | 0.184 | 0.159 | 0.119 | 0.207 |
| Guinea | \>= 5 | Outpatient | 0.048 | 0.033 | 0.067 | 0.053 | 0.036 | 0.076 |
| Guinea-Bissau | \< 5 | Inpatient | 0.309 | 0.231 | 0.405 | 0.279 | 0.195 | 0.382 |
| Guinea-Bissau | \>= 5 | Inpatient | 0.118 | 0.076 | 0.176 | 0.103 | 0.062 | 0.162 |
| Guinea-Bissau | \< 5 | Outpatient | 0.146 | 0.103 | 0.206 | 0.129 | 0.085 | 0.190 |
| Guinea-Bissau | \>= 5 | Outpatient | 0.049 | 0.031 | 0.075 | 0.042 | 0.025 | 0.069 |
| Guyana | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Guyana | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Guyana | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Guyana | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Haiti | \< 5 | Inpatient | 0.305 | 0.235 | 0.383 | 0.294 | 0.221 | 0.375 |
| Haiti | \>= 5 | Inpatient | 0.115 | 0.078 | 0.163 | 0.110 | 0.072 | 0.161 |
| Haiti | \< 5 | Outpatient | 0.143 | 0.105 | 0.189 | 0.137 | 0.098 | 0.185 |
| Haiti | \>= 5 | Outpatient | 0.047 | 0.031 | 0.068 | 0.045 | 0.029 | 0.067 |
| Honduras | \< 5 | Inpatient | 0.244 | 0.174 | 0.325 | 0.271 | 0.197 | 0.356 |
| Honduras | \>= 5 | Inpatient | 0.087 | 0.054 | 0.132 | 0.100 | 0.062 | 0.147 |
| Honduras | \< 5 | Outpatient | 0.110 | 0.074 | 0.155 | 0.125 | 0.086 | 0.172 |
| Honduras | \>= 5 | Outpatient | 0.035 | 0.021 | 0.054 | 0.040 | 0.025 | 0.062 |
| Hungary | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Hungary | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Hungary | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Hungary | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Iceland | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Iceland | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Iceland | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Iceland | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| India | \< 5 | Inpatient | 0.310 | 0.244 | 0.384 | 0.277 | 0.199 | 0.365 |
| India | \>= 5 | Inpatient | 0.118 | 0.081 | 0.165 | 0.102 | 0.064 | 0.154 |
| India | \< 5 | Outpatient | 0.147 | 0.110 | 0.191 | 0.128 | 0.086 | 0.180 |
| India | \>= 5 | Outpatient | 0.048 | 0.033 | 0.069 | 0.042 | 0.025 | 0.065 |
| Indonesia | \< 5 | Inpatient | 0.338 | 0.251 | 0.478 | 0.360 | 0.260 | 0.505 |
| Indonesia | \>= 5 | Inpatient | 0.133 | 0.085 | 0.222 | 0.146 | 0.090 | 0.245 |
| Indonesia | \< 5 | Outpatient | 0.164 | 0.114 | 0.257 | 0.178 | 0.119 | 0.280 |
| Indonesia | \>= 5 | Outpatient | 0.056 | 0.035 | 0.098 | 0.062 | 0.036 | 0.110 |
| Iran, Islamic Rep. | \< 5 | Inpatient | 0.287 | 0.205 | 0.371 | 0.309 | 0.213 | 0.407 |
| Iran, Islamic Rep. | \>= 5 | Inpatient | 0.107 | 0.067 | 0.155 | 0.118 | 0.070 | 0.175 |
| Iran, Islamic Rep. | \< 5 | Outpatient | 0.134 | 0.090 | 0.181 | 0.146 | 0.094 | 0.206 |
| Iran, Islamic Rep. | \>= 5 | Outpatient | 0.044 | 0.027 | 0.064 | 0.048 | 0.028 | 0.074 |
| Iraq | \< 5 | Inpatient | 0.295 | 0.228 | 0.364 | 0.288 | 0.217 | 0.366 |
| Iraq | \>= 5 | Inpatient | 0.111 | 0.075 | 0.153 | 0.107 | 0.071 | 0.153 |
| Iraq | \< 5 | Outpatient | 0.138 | 0.102 | 0.177 | 0.134 | 0.097 | 0.178 |
| Iraq | \>= 5 | Outpatient | 0.045 | 0.031 | 0.063 | 0.044 | 0.029 | 0.064 |
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
| Jamaica | \< 5 | Inpatient | 0.304 | 0.236 | 0.378 | 0.326 | 0.242 | 0.420 |
| Jamaica | \>= 5 | Inpatient | 0.115 | 0.078 | 0.160 | 0.126 | 0.082 | 0.185 |
| Jamaica | \< 5 | Outpatient | 0.143 | 0.106 | 0.186 | 0.156 | 0.110 | 0.214 |
| Jamaica | \>= 5 | Outpatient | 0.047 | 0.032 | 0.067 | 0.052 | 0.033 | 0.078 |
| Japan | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Japan | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Japan | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Japan | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Jordan | \< 5 | Inpatient | 0.302 | 0.218 | 0.403 | 0.274 | 0.185 | 0.383 |
| Jordan | \>= 5 | Inpatient | 0.114 | 0.071 | 0.173 | 0.101 | 0.059 | 0.162 |
| Jordan | \< 5 | Outpatient | 0.142 | 0.097 | 0.204 | 0.127 | 0.081 | 0.192 |
| Jordan | \>= 5 | Outpatient | 0.047 | 0.029 | 0.073 | 0.041 | 0.024 | 0.068 |
| Kazakhstan | \< 5 | Inpatient | 0.293 | 0.223 | 0.365 | 0.314 | 0.230 | 0.403 |
| Kazakhstan | \>= 5 | Inpatient | 0.109 | 0.074 | 0.152 | 0.120 | 0.078 | 0.173 |
| Kazakhstan | \< 5 | Outpatient | 0.137 | 0.098 | 0.178 | 0.149 | 0.103 | 0.201 |
| Kazakhstan | \>= 5 | Outpatient | 0.045 | 0.030 | 0.063 | 0.050 | 0.031 | 0.073 |
| Kenya | \< 5 | Inpatient | 0.301 | 0.231 | 0.374 | 0.263 | 0.187 | 0.351 |
| Kenya | \>= 5 | Inpatient | 0.114 | 0.077 | 0.159 | 0.096 | 0.059 | 0.144 |
| Kenya | \< 5 | Outpatient | 0.142 | 0.102 | 0.185 | 0.120 | 0.080 | 0.168 |
| Kenya | \>= 5 | Outpatient | 0.047 | 0.031 | 0.066 | 0.039 | 0.023 | 0.060 |
| Kiribati | \< 5 | Inpatient | 0.301 | 0.233 | 0.375 | 0.265 | 0.186 | 0.356 |
| Kiribati | \>= 5 | Inpatient | 0.113 | 0.077 | 0.158 | 0.097 | 0.059 | 0.147 |
| Kiribati | \< 5 | Outpatient | 0.141 | 0.104 | 0.184 | 0.121 | 0.080 | 0.173 |
| Kiribati | \>= 5 | Outpatient | 0.047 | 0.031 | 0.066 | 0.039 | 0.023 | 0.061 |
| Korea, Dem. People’s Rep. | \< 5 | Inpatient | 0.308 | 0.249 | 0.375 | 0.330 | 0.253 | 0.412 |
| Korea, Dem. People’s Rep. | \>= 5 | Inpatient | 0.117 | 0.082 | 0.160 | 0.128 | 0.086 | 0.182 |
| Korea, Dem. People’s Rep. | \< 5 | Outpatient | 0.145 | 0.113 | 0.185 | 0.159 | 0.116 | 0.210 |
| Korea, Dem. People’s Rep. | \>= 5 | Outpatient | 0.048 | 0.034 | 0.067 | 0.053 | 0.035 | 0.077 |
| Korea, Rep. | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Korea, Rep. | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Korea, Rep. | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Korea, Rep. | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kuwait | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kuwait | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kuwait | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kuwait | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kyrgyz Republic | \< 5 | Inpatient | 0.297 | 0.228 | 0.366 | 0.284 | 0.210 | 0.362 |
| Kyrgyz Republic | \>= 5 | Inpatient | 0.111 | 0.076 | 0.153 | 0.105 | 0.069 | 0.150 |
| Kyrgyz Republic | \< 5 | Outpatient | 0.139 | 0.102 | 0.177 | 0.131 | 0.092 | 0.175 |
| Kyrgyz Republic | \>= 5 | Outpatient | 0.046 | 0.031 | 0.063 | 0.043 | 0.028 | 0.063 |
| Lao PDR | \< 5 | Inpatient | 0.301 | 0.233 | 0.375 | 0.323 | 0.240 | 0.412 |
| Lao PDR | \>= 5 | Inpatient | 0.113 | 0.077 | 0.158 | 0.125 | 0.081 | 0.179 |
| Lao PDR | \< 5 | Outpatient | 0.141 | 0.104 | 0.184 | 0.154 | 0.107 | 0.208 |
| Lao PDR | \>= 5 | Outpatient | 0.047 | 0.031 | 0.066 | 0.052 | 0.033 | 0.075 |
| Latvia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Latvia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Latvia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Latvia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Lebanon | \< 5 | Inpatient | 0.296 | 0.212 | 0.389 | 0.317 | 0.220 | 0.422 |
| Lebanon | \>= 5 | Inpatient | 0.111 | 0.070 | 0.166 | 0.122 | 0.073 | 0.185 |
| Lebanon | \< 5 | Outpatient | 0.139 | 0.095 | 0.193 | 0.151 | 0.098 | 0.215 |
| Lebanon | \>= 5 | Outpatient | 0.046 | 0.028 | 0.069 | 0.050 | 0.030 | 0.078 |
| Lesotho | \< 5 | Inpatient | 0.309 | 0.250 | 0.376 | 0.273 | 0.204 | 0.353 |
| Lesotho | \>= 5 | Inpatient | 0.117 | 0.083 | 0.160 | 0.100 | 0.065 | 0.146 |
| Lesotho | \< 5 | Outpatient | 0.146 | 0.114 | 0.184 | 0.126 | 0.090 | 0.170 |
| Lesotho | \>= 5 | Outpatient | 0.048 | 0.033 | 0.067 | 0.041 | 0.026 | 0.061 |
| Liberia | \< 5 | Inpatient | 0.303 | 0.242 | 0.367 | 0.281 | 0.214 | 0.355 |
| Liberia | \>= 5 | Inpatient | 0.114 | 0.079 | 0.155 | 0.104 | 0.068 | 0.147 |
| Liberia | \< 5 | Outpatient | 0.142 | 0.109 | 0.179 | 0.130 | 0.094 | 0.171 |
| Liberia | \>= 5 | Outpatient | 0.047 | 0.032 | 0.065 | 0.042 | 0.027 | 0.061 |
| Libya | \< 5 | Inpatient | 0.295 | 0.211 | 0.388 | 0.269 | 0.181 | 0.371 |
| Libya | \>= 5 | Inpatient | 0.111 | 0.069 | 0.164 | 0.099 | 0.058 | 0.154 |
| Libya | \< 5 | Outpatient | 0.138 | 0.092 | 0.193 | 0.124 | 0.078 | 0.183 |
| Libya | \>= 5 | Outpatient | 0.045 | 0.027 | 0.069 | 0.040 | 0.023 | 0.065 |
| Lithuania | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Lithuania | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Lithuania | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Lithuania | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Luxembourg | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Luxembourg | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Luxembourg | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Luxembourg | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Madagascar | \< 5 | Inpatient | 0.310 | 0.229 | 0.407 | 0.288 | 0.207 | 0.390 |
| Madagascar | \>= 5 | Inpatient | 0.118 | 0.075 | 0.180 | 0.108 | 0.066 | 0.169 |
| Madagascar | \< 5 | Outpatient | 0.147 | 0.103 | 0.207 | 0.134 | 0.091 | 0.196 |
| Madagascar | \>= 5 | Outpatient | 0.049 | 0.031 | 0.076 | 0.044 | 0.026 | 0.071 |
| Malawi | \< 5 | Inpatient | 0.292 | 0.207 | 0.374 | 0.255 | 0.167 | 0.350 |
| Malawi | \>= 5 | Inpatient | 0.109 | 0.068 | 0.158 | 0.092 | 0.052 | 0.144 |
| Malawi | \< 5 | Outpatient | 0.136 | 0.090 | 0.184 | 0.116 | 0.071 | 0.169 |
| Malawi | \>= 5 | Outpatient | 0.045 | 0.027 | 0.066 | 0.037 | 0.021 | 0.060 |
| Malaysia | \< 5 | Inpatient | 0.314 | 0.232 | 0.418 | 0.336 | 0.241 | 0.447 |
| Malaysia | \>= 5 | Inpatient | 0.120 | 0.078 | 0.182 | 0.132 | 0.082 | 0.202 |
| Malaysia | \< 5 | Outpatient | 0.149 | 0.104 | 0.216 | 0.163 | 0.109 | 0.236 |
| Malaysia | \>= 5 | Outpatient | 0.050 | 0.031 | 0.078 | 0.055 | 0.033 | 0.088 |
| Maldives | \< 5 | Inpatient | 0.309 | 0.245 | 0.383 | 0.331 | 0.251 | 0.420 |
| Maldives | \>= 5 | Inpatient | 0.118 | 0.083 | 0.164 | 0.129 | 0.086 | 0.186 |
| Maldives | \< 5 | Outpatient | 0.146 | 0.111 | 0.191 | 0.159 | 0.115 | 0.216 |
| Maldives | \>= 5 | Outpatient | 0.048 | 0.033 | 0.069 | 0.054 | 0.035 | 0.079 |
| Mali | \< 5 | Inpatient | 0.299 | 0.220 | 0.381 | 0.276 | 0.194 | 0.366 |
| Mali | \>= 5 | Inpatient | 0.113 | 0.072 | 0.164 | 0.102 | 0.061 | 0.153 |
| Mali | \< 5 | Outpatient | 0.140 | 0.097 | 0.189 | 0.127 | 0.084 | 0.177 |
| Mali | \>= 5 | Outpatient | 0.046 | 0.029 | 0.069 | 0.041 | 0.025 | 0.064 |
| Malta | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Malta | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Malta | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Malta | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Marshall Islands | \< 5 | Inpatient | 0.289 | 0.230 | 0.359 | 0.289 | 0.219 | 0.370 |
| Marshall Islands | \>= 5 | Inpatient | 0.108 | 0.075 | 0.148 | 0.108 | 0.072 | 0.155 |
| Marshall Islands | \< 5 | Outpatient | 0.135 | 0.102 | 0.174 | 0.135 | 0.097 | 0.181 |
| Marshall Islands | \>= 5 | Outpatient | 0.044 | 0.030 | 0.061 | 0.044 | 0.029 | 0.065 |
| Mauritania | \< 5 | Inpatient | 0.309 | 0.250 | 0.376 | 0.275 | 0.207 | 0.354 |
| Mauritania | \>= 5 | Inpatient | 0.117 | 0.083 | 0.160 | 0.101 | 0.066 | 0.147 |
| Mauritania | \< 5 | Outpatient | 0.146 | 0.114 | 0.184 | 0.127 | 0.091 | 0.170 |
| Mauritania | \>= 5 | Outpatient | 0.048 | 0.033 | 0.067 | 0.041 | 0.026 | 0.062 |
| Mauritius | \< 5 | Inpatient | 0.302 | 0.235 | 0.370 | 0.263 | 0.191 | 0.345 |
| Mauritius | \>= 5 | Inpatient | 0.114 | 0.078 | 0.156 | 0.096 | 0.060 | 0.142 |
| Mauritius | \< 5 | Outpatient | 0.142 | 0.106 | 0.181 | 0.120 | 0.083 | 0.165 |
| Mauritius | \>= 5 | Outpatient | 0.047 | 0.032 | 0.065 | 0.039 | 0.024 | 0.059 |
| Mexico | \< 5 | Inpatient | 0.243 | 0.162 | 0.336 | 0.274 | 0.186 | 0.371 |
| Mexico | \>= 5 | Inpatient | 0.087 | 0.050 | 0.135 | 0.101 | 0.060 | 0.155 |
| Mexico | \< 5 | Outpatient | 0.109 | 0.069 | 0.161 | 0.127 | 0.081 | 0.181 |
| Mexico | \>= 5 | Outpatient | 0.035 | 0.020 | 0.056 | 0.041 | 0.024 | 0.065 |
| Micronesia, Fed. Sts. | \< 5 | Inpatient | 0.257 | 0.187 | 0.337 | 0.296 | 0.222 | 0.377 |
| Micronesia, Fed. Sts. | \>= 5 | Inpatient | 0.093 | 0.059 | 0.136 | 0.111 | 0.073 | 0.158 |
| Micronesia, Fed. Sts. | \< 5 | Outpatient | 0.117 | 0.081 | 0.161 | 0.138 | 0.098 | 0.185 |
| Micronesia, Fed. Sts. | \>= 5 | Outpatient | 0.038 | 0.024 | 0.057 | 0.045 | 0.029 | 0.066 |
| Moldova | \< 5 | Inpatient | 0.293 | 0.223 | 0.365 | 0.277 | 0.202 | 0.359 |
| Moldova | \>= 5 | Inpatient | 0.109 | 0.074 | 0.152 | 0.102 | 0.065 | 0.148 |
| Moldova | \< 5 | Outpatient | 0.137 | 0.098 | 0.178 | 0.128 | 0.088 | 0.173 |
| Moldova | \>= 5 | Outpatient | 0.045 | 0.030 | 0.063 | 0.042 | 0.026 | 0.061 |
| Monaco | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Monaco | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Monaco | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Monaco | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Mongolia | \< 5 | Inpatient | 0.301 | 0.233 | 0.375 | 0.323 | 0.240 | 0.412 |
| Mongolia | \>= 5 | Inpatient | 0.113 | 0.077 | 0.158 | 0.125 | 0.081 | 0.179 |
| Mongolia | \< 5 | Outpatient | 0.141 | 0.104 | 0.184 | 0.154 | 0.107 | 0.208 |
| Mongolia | \>= 5 | Outpatient | 0.047 | 0.031 | 0.066 | 0.052 | 0.033 | 0.075 |
| Montenegro | \< 5 | Inpatient | 0.293 | 0.223 | 0.365 | 0.314 | 0.230 | 0.403 |
| Montenegro | \>= 5 | Inpatient | 0.109 | 0.074 | 0.152 | 0.120 | 0.078 | 0.173 |
| Montenegro | \< 5 | Outpatient | 0.137 | 0.098 | 0.178 | 0.149 | 0.103 | 0.201 |
| Montenegro | \>= 5 | Outpatient | 0.045 | 0.030 | 0.063 | 0.050 | 0.031 | 0.073 |
| Morocco | \< 5 | Inpatient | 0.303 | 0.220 | 0.409 | 0.263 | 0.175 | 0.378 |
| Morocco | \>= 5 | Inpatient | 0.115 | 0.073 | 0.178 | 0.096 | 0.055 | 0.159 |
| Morocco | \< 5 | Outpatient | 0.143 | 0.099 | 0.209 | 0.121 | 0.075 | 0.188 |
| Morocco | \>= 5 | Outpatient | 0.047 | 0.029 | 0.075 | 0.039 | 0.022 | 0.067 |
| Mozambique | \< 5 | Inpatient | 0.312 | 0.237 | 0.409 | 0.294 | 0.215 | 0.395 |
| Mozambique | \>= 5 | Inpatient | 0.119 | 0.078 | 0.178 | 0.111 | 0.069 | 0.170 |
| Mozambique | \< 5 | Outpatient | 0.148 | 0.105 | 0.208 | 0.138 | 0.094 | 0.200 |
| Mozambique | \>= 5 | Outpatient | 0.049 | 0.031 | 0.076 | 0.045 | 0.027 | 0.072 |
| Myanmar | \< 5 | Inpatient | 0.308 | 0.249 | 0.375 | 0.297 | 0.229 | 0.372 |
| Myanmar | \>= 5 | Inpatient | 0.117 | 0.082 | 0.160 | 0.111 | 0.075 | 0.158 |
| Myanmar | \< 5 | Outpatient | 0.145 | 0.113 | 0.185 | 0.139 | 0.102 | 0.183 |
| Myanmar | \>= 5 | Outpatient | 0.048 | 0.034 | 0.067 | 0.046 | 0.030 | 0.066 |
| Namibia | \< 5 | Inpatient | 0.302 | 0.235 | 0.370 | 0.272 | 0.204 | 0.349 |
| Namibia | \>= 5 | Inpatient | 0.114 | 0.078 | 0.156 | 0.100 | 0.064 | 0.144 |
| Namibia | \< 5 | Outpatient | 0.142 | 0.106 | 0.181 | 0.125 | 0.089 | 0.168 |
| Namibia | \>= 5 | Outpatient | 0.047 | 0.032 | 0.065 | 0.041 | 0.026 | 0.060 |
| Nauru | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nauru | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nauru | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nauru | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nepal | \< 5 | Inpatient | 0.308 | 0.235 | 0.387 | 0.330 | 0.242 | 0.424 |
| Nepal | \>= 5 | Inpatient | 0.117 | 0.078 | 0.166 | 0.129 | 0.082 | 0.188 |
| Nepal | \< 5 | Outpatient | 0.146 | 0.106 | 0.193 | 0.159 | 0.109 | 0.218 |
| Nepal | \>= 5 | Outpatient | 0.048 | 0.032 | 0.070 | 0.053 | 0.033 | 0.080 |
| Netherlands | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Netherlands | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Netherlands | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Netherlands | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| New Zealand | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| New Zealand | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| New Zealand | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| New Zealand | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nicaragua | \< 5 | Inpatient | 0.243 | 0.161 | 0.338 | 0.265 | 0.179 | 0.364 |
| Nicaragua | \>= 5 | Inpatient | 0.087 | 0.050 | 0.138 | 0.097 | 0.056 | 0.154 |
| Nicaragua | \< 5 | Outpatient | 0.109 | 0.070 | 0.162 | 0.122 | 0.077 | 0.179 |
| Nicaragua | \>= 5 | Outpatient | 0.035 | 0.020 | 0.057 | 0.039 | 0.022 | 0.065 |
| Niger | \< 5 | Inpatient | 0.303 | 0.242 | 0.367 | 0.270 | 0.200 | 0.349 |
| Niger | \>= 5 | Inpatient | 0.114 | 0.079 | 0.155 | 0.099 | 0.063 | 0.144 |
| Niger | \< 5 | Outpatient | 0.142 | 0.109 | 0.179 | 0.124 | 0.087 | 0.168 |
| Niger | \>= 5 | Outpatient | 0.047 | 0.032 | 0.065 | 0.040 | 0.025 | 0.060 |
| Nigeria | \< 5 | Inpatient | 0.312 | 0.233 | 0.404 | 0.334 | 0.244 | 0.437 |
| Nigeria | \>= 5 | Inpatient | 0.119 | 0.077 | 0.177 | 0.131 | 0.082 | 0.197 |
| Nigeria | \< 5 | Outpatient | 0.148 | 0.105 | 0.205 | 0.161 | 0.110 | 0.226 |
| Nigeria | \>= 5 | Outpatient | 0.049 | 0.031 | 0.074 | 0.054 | 0.033 | 0.084 |
| Niue | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Niue | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Niue | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Niue | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| North Macedonia | \< 5 | Inpatient | 0.293 | 0.223 | 0.365 | 0.274 | 0.198 | 0.356 |
| North Macedonia | \>= 5 | Inpatient | 0.109 | 0.074 | 0.152 | 0.101 | 0.064 | 0.147 |
| North Macedonia | \< 5 | Outpatient | 0.137 | 0.098 | 0.178 | 0.126 | 0.086 | 0.172 |
| North Macedonia | \>= 5 | Outpatient | 0.045 | 0.030 | 0.063 | 0.041 | 0.026 | 0.061 |
| Norway | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Norway | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Norway | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Norway | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Oman | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Oman | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Oman | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Oman | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Pakistan | \< 5 | Inpatient | 0.301 | 0.224 | 0.391 | 0.268 | 0.184 | 0.370 |
| Pakistan | \>= 5 | Inpatient | 0.114 | 0.073 | 0.168 | 0.098 | 0.058 | 0.155 |
| Pakistan | \< 5 | Outpatient | 0.142 | 0.100 | 0.196 | 0.123 | 0.080 | 0.183 |
| Pakistan | \>= 5 | Outpatient | 0.047 | 0.030 | 0.071 | 0.040 | 0.023 | 0.065 |
| Palau | \< 5 | Inpatient | 0.258 | 0.189 | 0.338 | 0.274 | 0.199 | 0.361 |
| Palau | \>= 5 | Inpatient | 0.093 | 0.060 | 0.137 | 0.101 | 0.064 | 0.149 |
| Palau | \< 5 | Outpatient | 0.117 | 0.082 | 0.163 | 0.126 | 0.086 | 0.175 |
| Palau | \>= 5 | Outpatient | 0.038 | 0.024 | 0.057 | 0.041 | 0.025 | 0.062 |
| Panama | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Panama | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Panama | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Panama | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Papua New Guinea | \< 5 | Inpatient | 0.300 | 0.211 | 0.395 | 0.321 | 0.218 | 0.429 |
| Papua New Guinea | \>= 5 | Inpatient | 0.113 | 0.069 | 0.170 | 0.124 | 0.073 | 0.191 |
| Papua New Guinea | \< 5 | Outpatient | 0.141 | 0.092 | 0.200 | 0.154 | 0.096 | 0.221 |
| Papua New Guinea | \>= 5 | Outpatient | 0.046 | 0.028 | 0.072 | 0.051 | 0.029 | 0.081 |
| Paraguay | \< 5 | Inpatient | 0.268 | 0.208 | 0.333 | 0.273 | 0.204 | 0.352 |
| Paraguay | \>= 5 | Inpatient | 0.097 | 0.066 | 0.136 | 0.100 | 0.064 | 0.146 |
| Paraguay | \< 5 | Outpatient | 0.122 | 0.093 | 0.158 | 0.125 | 0.089 | 0.170 |
| Paraguay | \>= 5 | Outpatient | 0.039 | 0.027 | 0.056 | 0.041 | 0.026 | 0.061 |
| Peru | \< 5 | Inpatient | 0.258 | 0.185 | 0.347 | 0.277 | 0.198 | 0.373 |
| Peru | \>= 5 | Inpatient | 0.094 | 0.058 | 0.142 | 0.103 | 0.063 | 0.157 |
| Peru | \< 5 | Outpatient | 0.118 | 0.081 | 0.167 | 0.128 | 0.086 | 0.183 |
| Peru | \>= 5 | Outpatient | 0.038 | 0.023 | 0.059 | 0.042 | 0.025 | 0.066 |
| Philippines | \< 5 | Inpatient | 0.301 | 0.233 | 0.375 | 0.323 | 0.240 | 0.412 |
| Philippines | \>= 5 | Inpatient | 0.113 | 0.077 | 0.158 | 0.125 | 0.081 | 0.179 |
| Philippines | \< 5 | Outpatient | 0.141 | 0.104 | 0.184 | 0.154 | 0.107 | 0.208 |
| Philippines | \>= 5 | Outpatient | 0.047 | 0.031 | 0.066 | 0.052 | 0.033 | 0.075 |
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
| Russian Federation | \< 5 | Inpatient | 0.289 | 0.201 | 0.380 | 0.310 | 0.210 | 0.415 |
| Russian Federation | \>= 5 | Inpatient | 0.108 | 0.066 | 0.159 | 0.119 | 0.070 | 0.180 |
| Russian Federation | \< 5 | Outpatient | 0.135 | 0.088 | 0.189 | 0.147 | 0.092 | 0.211 |
| Russian Federation | \>= 5 | Outpatient | 0.044 | 0.027 | 0.066 | 0.049 | 0.028 | 0.076 |
| Rwanda | \< 5 | Inpatient | 0.303 | 0.242 | 0.367 | 0.265 | 0.193 | 0.348 |
| Rwanda | \>= 5 | Inpatient | 0.114 | 0.079 | 0.155 | 0.097 | 0.061 | 0.143 |
| Rwanda | \< 5 | Outpatient | 0.142 | 0.109 | 0.179 | 0.121 | 0.083 | 0.167 |
| Rwanda | \>= 5 | Outpatient | 0.047 | 0.032 | 0.065 | 0.039 | 0.024 | 0.060 |
| Samoa | \< 5 | Inpatient | 0.301 | 0.233 | 0.375 | 0.323 | 0.240 | 0.412 |
| Samoa | \>= 5 | Inpatient | 0.113 | 0.077 | 0.158 | 0.125 | 0.081 | 0.179 |
| Samoa | \< 5 | Outpatient | 0.141 | 0.104 | 0.184 | 0.154 | 0.107 | 0.208 |
| Samoa | \>= 5 | Outpatient | 0.047 | 0.031 | 0.066 | 0.052 | 0.033 | 0.075 |
| San Marino | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| San Marino | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| San Marino | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| San Marino | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| São Tomé and Principe | \< 5 | Inpatient | 0.309 | 0.250 | 0.376 | 0.273 | 0.204 | 0.353 |
| São Tomé and Principe | \>= 5 | Inpatient | 0.117 | 0.083 | 0.160 | 0.100 | 0.065 | 0.146 |
| São Tomé and Principe | \< 5 | Outpatient | 0.146 | 0.114 | 0.184 | 0.126 | 0.090 | 0.170 |
| São Tomé and Principe | \>= 5 | Outpatient | 0.048 | 0.033 | 0.067 | 0.041 | 0.026 | 0.061 |
| Saudi Arabia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Saudi Arabia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Saudi Arabia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Saudi Arabia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Senegal | \< 5 | Inpatient | 0.310 | 0.230 | 0.404 | 0.271 | 0.188 | 0.372 |
| Senegal | \>= 5 | Inpatient | 0.118 | 0.076 | 0.177 | 0.100 | 0.060 | 0.158 |
| Senegal | \< 5 | Outpatient | 0.147 | 0.102 | 0.205 | 0.125 | 0.081 | 0.184 |
| Senegal | \>= 5 | Outpatient | 0.049 | 0.031 | 0.075 | 0.041 | 0.024 | 0.066 |
| Serbia | \< 5 | Inpatient | 0.293 | 0.223 | 0.365 | 0.314 | 0.230 | 0.403 |
| Serbia | \>= 5 | Inpatient | 0.109 | 0.074 | 0.152 | 0.120 | 0.078 | 0.173 |
| Serbia | \< 5 | Outpatient | 0.137 | 0.098 | 0.178 | 0.149 | 0.103 | 0.201 |
| Serbia | \>= 5 | Outpatient | 0.045 | 0.030 | 0.063 | 0.050 | 0.031 | 0.073 |
| Seychelles | \< 5 | Inpatient | 0.302 | 0.235 | 0.370 | 0.260 | 0.185 | 0.344 |
| Seychelles | \>= 5 | Inpatient | 0.114 | 0.078 | 0.156 | 0.094 | 0.058 | 0.141 |
| Seychelles | \< 5 | Outpatient | 0.142 | 0.106 | 0.181 | 0.118 | 0.080 | 0.165 |
| Seychelles | \>= 5 | Outpatient | 0.047 | 0.032 | 0.065 | 0.038 | 0.023 | 0.059 |
| Sierra Leone | \< 5 | Inpatient | 0.303 | 0.242 | 0.367 | 0.265 | 0.194 | 0.348 |
| Sierra Leone | \>= 5 | Inpatient | 0.114 | 0.079 | 0.155 | 0.097 | 0.061 | 0.143 |
| Sierra Leone | \< 5 | Outpatient | 0.142 | 0.109 | 0.179 | 0.121 | 0.084 | 0.167 |
| Sierra Leone | \>= 5 | Outpatient | 0.047 | 0.032 | 0.065 | 0.039 | 0.024 | 0.060 |
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
| Solomon Islands | \< 5 | Inpatient | 0.301 | 0.233 | 0.375 | 0.295 | 0.222 | 0.377 |
| Solomon Islands | \>= 5 | Inpatient | 0.113 | 0.077 | 0.158 | 0.111 | 0.072 | 0.158 |
| Solomon Islands | \< 5 | Outpatient | 0.141 | 0.104 | 0.184 | 0.138 | 0.097 | 0.185 |
| Solomon Islands | \>= 5 | Outpatient | 0.047 | 0.031 | 0.066 | 0.045 | 0.029 | 0.066 |
| Somalia | \< 5 | Inpatient | 0.300 | 0.226 | 0.382 | 0.322 | 0.232 | 0.419 |
| Somalia | \>= 5 | Inpatient | 0.113 | 0.075 | 0.161 | 0.124 | 0.078 | 0.182 |
| Somalia | \< 5 | Outpatient | 0.141 | 0.100 | 0.188 | 0.154 | 0.104 | 0.214 |
| Somalia | \>= 5 | Outpatient | 0.046 | 0.030 | 0.067 | 0.051 | 0.032 | 0.077 |
| South Africa | \< 5 | Inpatient | 0.252 | 0.180 | 0.328 | 0.258 | 0.181 | 0.341 |
| South Africa | \>= 5 | Inpatient | 0.091 | 0.057 | 0.133 | 0.093 | 0.057 | 0.140 |
| South Africa | \< 5 | Outpatient | 0.114 | 0.078 | 0.155 | 0.117 | 0.078 | 0.164 |
| South Africa | \>= 5 | Outpatient | 0.037 | 0.023 | 0.055 | 0.038 | 0.023 | 0.058 |
| South Sudan | \< 5 | Inpatient | 0.303 | 0.242 | 0.367 | 0.324 | 0.250 | 0.405 |
| South Sudan | \>= 5 | Inpatient | 0.114 | 0.079 | 0.155 | 0.125 | 0.084 | 0.176 |
| South Sudan | \< 5 | Outpatient | 0.142 | 0.109 | 0.179 | 0.155 | 0.114 | 0.204 |
| South Sudan | \>= 5 | Outpatient | 0.047 | 0.032 | 0.065 | 0.052 | 0.034 | 0.074 |
| Spain | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Spain | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Spain | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Spain | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Sri Lanka | \< 5 | Inpatient | 0.308 | 0.249 | 0.375 | 0.330 | 0.253 | 0.412 |
| Sri Lanka | \>= 5 | Inpatient | 0.117 | 0.082 | 0.160 | 0.128 | 0.086 | 0.182 |
| Sri Lanka | \< 5 | Outpatient | 0.145 | 0.113 | 0.185 | 0.159 | 0.116 | 0.210 |
| Sri Lanka | \>= 5 | Outpatient | 0.048 | 0.034 | 0.067 | 0.053 | 0.035 | 0.077 |
| St. Kitts and Nevis | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| St. Kitts and Nevis | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| St. Kitts and Nevis | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| St. Kitts and Nevis | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| St. Lucia | \< 5 | Inpatient | 0.304 | 0.236 | 0.378 | 0.326 | 0.242 | 0.420 |
| St. Lucia | \>= 5 | Inpatient | 0.115 | 0.078 | 0.160 | 0.126 | 0.082 | 0.185 |
| St. Lucia | \< 5 | Outpatient | 0.143 | 0.106 | 0.186 | 0.156 | 0.110 | 0.214 |
| St. Lucia | \>= 5 | Outpatient | 0.047 | 0.032 | 0.067 | 0.052 | 0.033 | 0.078 |
| St. Vincent and the Grenadines | \< 5 | Inpatient | 0.304 | 0.236 | 0.378 | 0.326 | 0.242 | 0.420 |
| St. Vincent and the Grenadines | \>= 5 | Inpatient | 0.115 | 0.078 | 0.160 | 0.126 | 0.082 | 0.185 |
| St. Vincent and the Grenadines | \< 5 | Outpatient | 0.143 | 0.106 | 0.186 | 0.156 | 0.110 | 0.214 |
| St. Vincent and the Grenadines | \>= 5 | Outpatient | 0.047 | 0.032 | 0.067 | 0.052 | 0.033 | 0.078 |
| Sudan | \< 5 | Inpatient | 0.303 | 0.211 | 0.416 | 0.265 | 0.173 | 0.383 |
| Sudan | \>= 5 | Inpatient | 0.115 | 0.069 | 0.180 | 0.097 | 0.054 | 0.161 |
| Sudan | \< 5 | Outpatient | 0.143 | 0.093 | 0.211 | 0.122 | 0.074 | 0.189 |
| Sudan | \>= 5 | Outpatient | 0.047 | 0.028 | 0.076 | 0.040 | 0.021 | 0.068 |
| Suriname | \< 5 | Inpatient | 0.304 | 0.236 | 0.378 | 0.326 | 0.242 | 0.420 |
| Suriname | \>= 5 | Inpatient | 0.115 | 0.078 | 0.160 | 0.126 | 0.082 | 0.185 |
| Suriname | \< 5 | Outpatient | 0.143 | 0.106 | 0.186 | 0.156 | 0.110 | 0.214 |
| Suriname | \>= 5 | Outpatient | 0.047 | 0.032 | 0.067 | 0.052 | 0.033 | 0.078 |
| Sweden | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Sweden | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Sweden | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Sweden | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Switzerland | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Switzerland | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Switzerland | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Switzerland | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Syrian Arab Republic | \< 5 | Inpatient | 0.300 | 0.226 | 0.382 | 0.322 | 0.232 | 0.419 |
| Syrian Arab Republic | \>= 5 | Inpatient | 0.113 | 0.075 | 0.161 | 0.124 | 0.078 | 0.182 |
| Syrian Arab Republic | \< 5 | Outpatient | 0.141 | 0.100 | 0.188 | 0.154 | 0.104 | 0.214 |
| Syrian Arab Republic | \>= 5 | Outpatient | 0.046 | 0.030 | 0.067 | 0.051 | 0.032 | 0.077 |
| Tajikistan | \< 5 | Inpatient | 0.297 | 0.228 | 0.366 | 0.257 | 0.177 | 0.344 |
| Tajikistan | \>= 5 | Inpatient | 0.111 | 0.076 | 0.153 | 0.093 | 0.056 | 0.141 |
| Tajikistan | \< 5 | Outpatient | 0.139 | 0.102 | 0.177 | 0.117 | 0.077 | 0.166 |
| Tajikistan | \>= 5 | Outpatient | 0.046 | 0.031 | 0.063 | 0.038 | 0.022 | 0.059 |
| Tanzania | \< 5 | Inpatient | 0.318 | 0.244 | 0.413 | 0.281 | 0.201 | 0.379 |
| Tanzania | \>= 5 | Inpatient | 0.123 | 0.081 | 0.182 | 0.104 | 0.064 | 0.163 |
| Tanzania | \< 5 | Outpatient | 0.152 | 0.110 | 0.210 | 0.130 | 0.088 | 0.188 |
| Tanzania | \>= 5 | Outpatient | 0.051 | 0.033 | 0.077 | 0.043 | 0.026 | 0.068 |
| Thailand | \< 5 | Inpatient | 0.291 | 0.211 | 0.373 | 0.312 | 0.217 | 0.409 |
| Thailand | \>= 5 | Inpatient | 0.109 | 0.069 | 0.156 | 0.119 | 0.072 | 0.177 |
| Thailand | \< 5 | Outpatient | 0.136 | 0.091 | 0.184 | 0.148 | 0.094 | 0.207 |
| Thailand | \>= 5 | Outpatient | 0.044 | 0.028 | 0.065 | 0.049 | 0.028 | 0.075 |
| Timor-Leste | \< 5 | Inpatient | 0.308 | 0.249 | 0.375 | 0.330 | 0.253 | 0.412 |
| Timor-Leste | \>= 5 | Inpatient | 0.117 | 0.082 | 0.160 | 0.128 | 0.086 | 0.182 |
| Timor-Leste | \< 5 | Outpatient | 0.145 | 0.113 | 0.185 | 0.159 | 0.116 | 0.210 |
| Timor-Leste | \>= 5 | Outpatient | 0.048 | 0.034 | 0.067 | 0.053 | 0.035 | 0.077 |
| Togo | \< 5 | Inpatient | 0.303 | 0.242 | 0.367 | 0.273 | 0.204 | 0.351 |
| Togo | \>= 5 | Inpatient | 0.114 | 0.079 | 0.155 | 0.100 | 0.065 | 0.145 |
| Togo | \< 5 | Outpatient | 0.142 | 0.109 | 0.179 | 0.125 | 0.089 | 0.168 |
| Togo | \>= 5 | Outpatient | 0.047 | 0.032 | 0.065 | 0.041 | 0.026 | 0.061 |
| Tonga | \< 5 | Inpatient | 0.306 | 0.242 | 0.378 | 0.328 | 0.249 | 0.411 |
| Tonga | \>= 5 | Inpatient | 0.116 | 0.081 | 0.159 | 0.127 | 0.085 | 0.180 |
| Tonga | \< 5 | Outpatient | 0.144 | 0.109 | 0.186 | 0.157 | 0.113 | 0.209 |
| Tonga | \>= 5 | Outpatient | 0.048 | 0.033 | 0.066 | 0.053 | 0.035 | 0.076 |
| Trinidad and Tobago | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Trinidad and Tobago | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Trinidad and Tobago | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Trinidad and Tobago | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Tunisia | \< 5 | Inpatient | 0.295 | 0.228 | 0.364 | 0.317 | 0.234 | 0.403 |
| Tunisia | \>= 5 | Inpatient | 0.111 | 0.075 | 0.153 | 0.121 | 0.078 | 0.174 |
| Tunisia | \< 5 | Outpatient | 0.138 | 0.102 | 0.177 | 0.151 | 0.106 | 0.203 |
| Tunisia | \>= 5 | Outpatient | 0.045 | 0.031 | 0.063 | 0.050 | 0.032 | 0.073 |
| Türkiye | \< 5 | Inpatient | 0.288 | 0.211 | 0.367 | 0.309 | 0.221 | 0.401 |
| Türkiye | \>= 5 | Inpatient | 0.107 | 0.070 | 0.154 | 0.118 | 0.074 | 0.172 |
| Türkiye | \< 5 | Outpatient | 0.134 | 0.094 | 0.178 | 0.146 | 0.099 | 0.201 |
| Türkiye | \>= 5 | Outpatient | 0.044 | 0.028 | 0.063 | 0.048 | 0.030 | 0.072 |
| Turkmenistan | \< 5 | Inpatient | 0.293 | 0.223 | 0.365 | 0.254 | 0.172 | 0.344 |
| Turkmenistan | \>= 5 | Inpatient | 0.109 | 0.074 | 0.152 | 0.092 | 0.054 | 0.141 |
| Turkmenistan | \< 5 | Outpatient | 0.137 | 0.098 | 0.178 | 0.115 | 0.074 | 0.166 |
| Turkmenistan | \>= 5 | Outpatient | 0.045 | 0.030 | 0.063 | 0.037 | 0.021 | 0.058 |
| Tuvalu | \< 5 | Inpatient | 0.306 | 0.242 | 0.378 | 0.328 | 0.249 | 0.411 |
| Tuvalu | \>= 5 | Inpatient | 0.116 | 0.081 | 0.159 | 0.127 | 0.085 | 0.180 |
| Tuvalu | \< 5 | Outpatient | 0.144 | 0.109 | 0.186 | 0.157 | 0.113 | 0.209 |
| Tuvalu | \>= 5 | Outpatient | 0.048 | 0.033 | 0.066 | 0.053 | 0.035 | 0.076 |
| Uganda | \< 5 | Inpatient | 0.303 | 0.242 | 0.367 | 0.267 | 0.196 | 0.349 |
| Uganda | \>= 5 | Inpatient | 0.114 | 0.079 | 0.155 | 0.097 | 0.062 | 0.143 |
| Uganda | \< 5 | Outpatient | 0.142 | 0.109 | 0.179 | 0.122 | 0.084 | 0.167 |
| Uganda | \>= 5 | Outpatient | 0.047 | 0.032 | 0.065 | 0.040 | 0.025 | 0.060 |
| Ukraine | \< 5 | Inpatient | 0.297 | 0.228 | 0.366 | 0.318 | 0.235 | 0.405 |
| Ukraine | \>= 5 | Inpatient | 0.111 | 0.076 | 0.153 | 0.122 | 0.080 | 0.175 |
| Ukraine | \< 5 | Outpatient | 0.139 | 0.102 | 0.177 | 0.152 | 0.106 | 0.202 |
| Ukraine | \>= 5 | Outpatient | 0.046 | 0.031 | 0.063 | 0.050 | 0.032 | 0.074 |
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
| Uzbekistan | \< 5 | Inpatient | 0.297 | 0.228 | 0.366 | 0.260 | 0.181 | 0.346 |
| Uzbekistan | \>= 5 | Inpatient | 0.111 | 0.076 | 0.153 | 0.094 | 0.057 | 0.142 |
| Uzbekistan | \< 5 | Outpatient | 0.139 | 0.102 | 0.177 | 0.119 | 0.079 | 0.167 |
| Uzbekistan | \>= 5 | Outpatient | 0.046 | 0.031 | 0.063 | 0.038 | 0.023 | 0.059 |
| Vanuatu | \< 5 | Inpatient | 0.301 | 0.233 | 0.375 | 0.323 | 0.240 | 0.412 |
| Vanuatu | \>= 5 | Inpatient | 0.113 | 0.077 | 0.158 | 0.125 | 0.081 | 0.179 |
| Vanuatu | \< 5 | Outpatient | 0.141 | 0.104 | 0.184 | 0.154 | 0.107 | 0.208 |
| Vanuatu | \>= 5 | Outpatient | 0.047 | 0.031 | 0.066 | 0.052 | 0.033 | 0.075 |
| Venezuela | \< 5 | Inpatient | 0.273 | 0.197 | 0.359 | 0.327 | 0.226 | 0.437 |
| Venezuela | \>= 5 | Inpatient | 0.100 | 0.063 | 0.150 | 0.127 | 0.076 | 0.196 |
| Venezuela | \< 5 | Outpatient | 0.126 | 0.087 | 0.176 | 0.157 | 0.101 | 0.228 |
| Venezuela | \>= 5 | Outpatient | 0.041 | 0.025 | 0.063 | 0.053 | 0.031 | 0.084 |
| Vietnam | \< 5 | Inpatient | 0.308 | 0.228 | 0.405 | 0.329 | 0.236 | 0.438 |
| Vietnam | \>= 5 | Inpatient | 0.117 | 0.075 | 0.174 | 0.129 | 0.080 | 0.196 |
| Vietnam | \< 5 | Outpatient | 0.146 | 0.101 | 0.205 | 0.159 | 0.105 | 0.228 |
| Vietnam | \>= 5 | Outpatient | 0.048 | 0.030 | 0.073 | 0.053 | 0.032 | 0.084 |
| Yemen, Rep. | \< 5 | Inpatient | 0.300 | 0.226 | 0.382 | 0.285 | 0.206 | 0.372 |
| Yemen, Rep. | \>= 5 | Inpatient | 0.113 | 0.075 | 0.161 | 0.106 | 0.067 | 0.157 |
| Yemen, Rep. | \< 5 | Outpatient | 0.141 | 0.100 | 0.188 | 0.132 | 0.091 | 0.183 |
| Yemen, Rep. | \>= 5 | Outpatient | 0.046 | 0.030 | 0.067 | 0.043 | 0.027 | 0.066 |
| Zambia | \< 5 | Inpatient | 0.321 | 0.242 | 0.435 | 0.285 | 0.201 | 0.402 |
| Zambia | \>= 5 | Inpatient | 0.124 | 0.081 | 0.193 | 0.106 | 0.064 | 0.175 |
| Zambia | \< 5 | Outpatient | 0.154 | 0.109 | 0.225 | 0.133 | 0.087 | 0.204 |
| Zambia | \>= 5 | Outpatient | 0.051 | 0.032 | 0.083 | 0.044 | 0.025 | 0.074 |
| Zimbabwe | \< 5 | Inpatient | 0.307 | 0.227 | 0.401 | 0.271 | 0.187 | 0.370 |
| Zimbabwe | \>= 5 | Inpatient | 0.117 | 0.074 | 0.175 | 0.100 | 0.058 | 0.157 |
| Zimbabwe | \< 5 | Outpatient | 0.145 | 0.101 | 0.203 | 0.125 | 0.081 | 0.180 |
| Zimbabwe | \>= 5 | Outpatient | 0.048 | 0.030 | 0.073 | 0.041 | 0.024 | 0.066 |

Estimated `r params$Pathogen` proportion by country, 2010 vs 2020

# Session info

``` r
sessioninfo::session_info()
```

    ## Warning in system2("quarto", "-V", stdout = TRUE, env = paste0("TMPDIR=", : running command '"quarto"
    ## TMPDIR=C:/Users/fbbu6966/AppData/Local/Temp/RtmpOgtm2C/file37e86b335d37 -V' had status 1

    ## ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
    ##  setting  value
    ##  version  R version 4.5.0 (2025-04-11 ucrt)
    ##  os       Windows 10 x64 (build 19045)
    ##  system   x86_64, mingw32
    ##  ui       RStudio
    ##  language (EN)
    ##  collate  English_United States.utf8
    ##  ctype    English_United States.utf8
    ##  tz       Europe/Brussels
    ##  date     2025-09-28
    ##  rstudio  2025.05.0+496 Mariposa Orchid (desktop)
    ##  pandoc   3.4 @ C:/Users/fbbu6966/AppData/Local/Programs/RStudio/resources/app/bin/quarto/bin/tools/ (via rmarkdown)
    ##  quarto   ERROR: Unknown command "TMPDIR=C:/Users/fbbu6966/AppData/Local/Temp/RtmpOgtm2C/file37e86b335d37". Did you mean command "update"? @ C:\\Users\\fbbu6966\\AppData\\Local\\Programs\\RStudio\\RESOUR~1\\app\\bin\\quarto\\bin\\quarto.exe
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
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
    ## ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

``` r
# Save dataset for report created for expert to receive feedback
# save(all_cnt_rt, file="./00-Report_FB/all_cnt_rt.Rdata")
# save(all_glb_prop, file="./00-Report_FB/all_glb_prop.Rdata")
# save(all_reg_prop, file="./00-Report_FB/all_reg_prop.Rdata")
# save(all_reg_rt, file="./00-Report_FB/all_reg_rt.Rdata")
# save(all_sub_nr, file="./00-Report_FB/all_sub_nr.Rdata")
# save(all_sub_rt, file="./00-Report_FB/all_sub_rt.Rdata")
```
