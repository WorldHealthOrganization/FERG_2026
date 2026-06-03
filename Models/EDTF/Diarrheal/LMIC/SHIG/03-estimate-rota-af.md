Global proportion of Shigella (Diarrheal) • ASYMPT+ROTA • Estimate
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

    ## Warning: There were 3 divergent transitions after warmup. Increasing adapt_delta above 0.95 may help. See
    ## http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup

    ##  Family: gaussian 
    ##   Links: mu = identity; sigma = identity 
    ## Formula: yi | se(sei) ~ 1 + YEAR + SYNDROMTYPE + AGE + REFERENCE + COVERAGE + (1 | REG2) + (1 | REG2:SUB2) + (1 | REG2:SUB2:COUNTRY) + (1 | REG2:SUB2:COUNTRY:ID) + (1 | REG2:SUB2:COUNTRY:ID:DTP_ID) 
    ##    Data: subset(es, as.integer(FLAG) == 1) (Number of observations: 533) 
    ##   Draws: 5 chains, each with iter = 5000; warmup = 3000; thin = 1;
    ##          total post-warmup draws = 10000
    ## 
    ## Multilevel Hyperparameters:
    ## ~REG2 (Number of levels: 6) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.55      0.31     0.06     1.31 1.00     3855     3983
    ## 
    ## ~REG2:SUB2 (Number of levels: 12) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.25      0.20     0.01     0.75 1.00     3890     5407
    ## 
    ## ~REG2:SUB2:COUNTRY (Number of levels: 50) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.20      0.15     0.01     0.55 1.00     2002     3493
    ## 
    ## ~REG2:SUB2:COUNTRY:ID (Number of levels: 165) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.97      0.09     0.81     1.16 1.00     3569     5221
    ## 
    ## ~REG2:SUB2:COUNTRY:ID:DTP_ID (Number of levels: 533) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.82      0.04     0.75     0.90 1.00     3258     5865
    ## 
    ## Regression Coefficients:
    ##                       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept                88.76     35.98    16.60   158.46 1.00     5641     6679
    ## YEAR                     -0.05      0.02    -0.08    -0.01 1.00     5637     6699
    ## SYNDROMTYPEInpatient      1.25      0.14     0.96     1.53 1.00     4519     5993
    ## SYNDROMTYPEOutpatient     1.06      0.12     0.83     1.29 1.00     4339     5409
    ## AGEAgebelow5             -0.13      0.18    -0.48     0.22 1.00     4325     5881
    ## AGEMixedages              0.16      0.25    -0.31     0.64 1.00     6222     7228
    ## REFERENCEOther           -1.93      0.24    -2.40    -1.45 1.00     4616     5874
    ## COVERAGE                  0.00      0.00    -0.00     0.01 1.00     5778     6750
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

    ## [1] 6600

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

    ## [1] 10560

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
    ##  $ VAL_MEAN     : num  0.253 0.232 0.246 0.225 0.239 ...
    ##  $ VAL_LWR      : num  0.106 0.0961 0.103 0.0934 0.1004 ...
    ##  $ VAL_UPR      : num  0.448 0.415 0.438 0.405 0.429 ...
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
    ##  $ VAL_MEAN     : num  0.253 0.232 0.246 0.225 0.239 ...
    ##  $ VAL_LWR      : num  0.106 0.0961 0.103 0.0934 0.1004 ...
    ##  $ VAL_UPR      : num  0.448 0.415 0.438 0.405 0.429 ...
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

\[1\] 0.00 0.05 0.10 0.15 0.20 0.25

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

\[1\] 0.00 0.05 0.10 0.15 0.20

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

\[1\] 0.00 0.05 0.10 0.15 0.20

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

\[1\] 0.00 0.05 0.10 0.15

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

\[1\] 0.00 0.05 0.10 0.15 0.20

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

\[1\] 0.00 0.02 0.04 0.06 0.08 0.10 0.12 0.14

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

\[1\] 0.00 0.05 0.10 0.15 0.20

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

\[1\] 0.00 0.02 0.04 0.06 0.08 0.10 0.12 0.14

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
| Afghanistan | \< 5 | Inpatient | 0.169 | 0.070 | 0.314 | 0.136 | 0.054 | 0.262 |
| Afghanistan | \>= 5 | Inpatient | 0.186 | 0.078 | 0.351 | 0.151 | 0.059 | 0.295 |
| Afghanistan | \< 5 | Outpatient | 0.133 | 0.055 | 0.251 | 0.106 | 0.042 | 0.208 |
| Afghanistan | \>= 5 | Outpatient | 0.148 | 0.061 | 0.285 | 0.119 | 0.047 | 0.238 |
| Albania | \< 5 | Inpatient | 0.095 | 0.032 | 0.206 | 0.083 | 0.024 | 0.195 |
| Albania | \>= 5 | Inpatient | 0.106 | 0.036 | 0.229 | 0.094 | 0.027 | 0.222 |
| Albania | \< 5 | Outpatient | 0.073 | 0.025 | 0.159 | 0.064 | 0.019 | 0.152 |
| Albania | \>= 5 | Outpatient | 0.082 | 0.028 | 0.180 | 0.072 | 0.021 | 0.171 |
| Algeria | \< 5 | Inpatient | 0.158 | 0.093 | 0.240 | 0.110 | 0.060 | 0.177 |
| Algeria | \>= 5 | Inpatient | 0.175 | 0.100 | 0.272 | 0.123 | 0.066 | 0.204 |
| Algeria | \< 5 | Outpatient | 0.124 | 0.075 | 0.185 | 0.085 | 0.048 | 0.136 |
| Algeria | \>= 5 | Outpatient | 0.138 | 0.080 | 0.214 | 0.096 | 0.052 | 0.158 |
| Andorra | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Andorra | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Andorra | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Andorra | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Angola | \< 5 | Inpatient | 0.154 | 0.071 | 0.256 | 0.118 | 0.055 | 0.199 |
| Angola | \>= 5 | Inpatient | 0.171 | 0.077 | 0.287 | 0.132 | 0.059 | 0.228 |
| Angola | \< 5 | Outpatient | 0.120 | 0.056 | 0.199 | 0.091 | 0.043 | 0.155 |
| Angola | \>= 5 | Outpatient | 0.135 | 0.061 | 0.225 | 0.103 | 0.047 | 0.180 |
| Antigua and Barbuda | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Antigua and Barbuda | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Antigua and Barbuda | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Antigua and Barbuda | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Argentina | \< 5 | Inpatient | 0.150 | 0.067 | 0.297 | 0.123 | 0.055 | 0.260 |
| Argentina | \>= 5 | Inpatient | 0.166 | 0.073 | 0.328 | 0.138 | 0.059 | 0.291 |
| Argentina | \< 5 | Outpatient | 0.118 | 0.052 | 0.236 | 0.096 | 0.043 | 0.208 |
| Argentina | \>= 5 | Outpatient | 0.131 | 0.057 | 0.262 | 0.108 | 0.046 | 0.231 |
| Armenia | \< 5 | Inpatient | 0.092 | 0.033 | 0.191 | 0.080 | 0.025 | 0.175 |
| Armenia | \>= 5 | Inpatient | 0.104 | 0.037 | 0.211 | 0.090 | 0.027 | 0.195 |
| Armenia | \< 5 | Outpatient | 0.071 | 0.026 | 0.144 | 0.061 | 0.019 | 0.136 |
| Armenia | \>= 5 | Outpatient | 0.080 | 0.029 | 0.162 | 0.069 | 0.021 | 0.153 |
| Australia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Australia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Australia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Australia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Austria | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Austria | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Austria | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Austria | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Azerbaijan | \< 5 | Inpatient | 0.092 | 0.033 | 0.191 | 0.062 | 0.020 | 0.134 |
| Azerbaijan | \>= 5 | Inpatient | 0.104 | 0.037 | 0.211 | 0.070 | 0.023 | 0.151 |
| Azerbaijan | \< 5 | Outpatient | 0.071 | 0.026 | 0.144 | 0.048 | 0.016 | 0.102 |
| Azerbaijan | \>= 5 | Outpatient | 0.080 | 0.029 | 0.162 | 0.054 | 0.018 | 0.116 |
| Bahamas, The | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahamas, The | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahamas, The | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahamas, The | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahrain | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahrain | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahrain | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahrain | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bangladesh | \< 5 | Inpatient | 0.140 | 0.080 | 0.224 | 0.097 | 0.049 | 0.166 |
| Bangladesh | \>= 5 | Inpatient | 0.156 | 0.087 | 0.252 | 0.109 | 0.055 | 0.191 |
| Bangladesh | \< 5 | Outpatient | 0.109 | 0.063 | 0.170 | 0.075 | 0.039 | 0.127 |
| Bangladesh | \>= 5 | Outpatient | 0.122 | 0.069 | 0.197 | 0.084 | 0.043 | 0.146 |
| Barbados | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Barbados | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Barbados | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Barbados | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belarus | \< 5 | Inpatient | 0.092 | 0.033 | 0.191 | 0.062 | 0.020 | 0.134 |
| Belarus | \>= 5 | Inpatient | 0.104 | 0.037 | 0.211 | 0.070 | 0.023 | 0.151 |
| Belarus | \< 5 | Outpatient | 0.071 | 0.026 | 0.144 | 0.048 | 0.016 | 0.102 |
| Belarus | \>= 5 | Outpatient | 0.080 | 0.029 | 0.162 | 0.054 | 0.018 | 0.116 |
| Belgium | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belgium | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belgium | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belgium | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belize | \< 5 | Inpatient | 0.130 | 0.066 | 0.218 | 0.090 | 0.040 | 0.167 |
| Belize | \>= 5 | Inpatient | 0.145 | 0.072 | 0.247 | 0.101 | 0.044 | 0.194 |
| Belize | \< 5 | Outpatient | 0.101 | 0.052 | 0.170 | 0.069 | 0.031 | 0.128 |
| Belize | \>= 5 | Outpatient | 0.114 | 0.056 | 0.195 | 0.078 | 0.034 | 0.148 |
| Benin | \< 5 | Inpatient | 0.158 | 0.093 | 0.240 | 0.127 | 0.075 | 0.197 |
| Benin | \>= 5 | Inpatient | 0.175 | 0.100 | 0.272 | 0.142 | 0.079 | 0.227 |
| Benin | \< 5 | Outpatient | 0.124 | 0.075 | 0.185 | 0.098 | 0.059 | 0.153 |
| Benin | \>= 5 | Outpatient | 0.138 | 0.080 | 0.214 | 0.111 | 0.063 | 0.178 |
| Bhutan | \< 5 | Inpatient | 0.133 | 0.075 | 0.208 | 0.091 | 0.046 | 0.156 |
| Bhutan | \>= 5 | Inpatient | 0.148 | 0.082 | 0.237 | 0.103 | 0.051 | 0.179 |
| Bhutan | \< 5 | Outpatient | 0.103 | 0.060 | 0.159 | 0.070 | 0.037 | 0.120 |
| Bhutan | \>= 5 | Outpatient | 0.116 | 0.064 | 0.184 | 0.079 | 0.040 | 0.138 |
| Bolivia | \< 5 | Inpatient | 0.161 | 0.074 | 0.288 | 0.109 | 0.047 | 0.204 |
| Bolivia | \>= 5 | Inpatient | 0.179 | 0.078 | 0.324 | 0.123 | 0.051 | 0.239 |
| Bolivia | \< 5 | Outpatient | 0.127 | 0.057 | 0.232 | 0.085 | 0.036 | 0.162 |
| Bolivia | \>= 5 | Outpatient | 0.142 | 0.062 | 0.264 | 0.096 | 0.039 | 0.188 |
| Bosnia and Herzegovina | \< 5 | Inpatient | 0.092 | 0.033 | 0.191 | 0.062 | 0.020 | 0.134 |
| Bosnia and Herzegovina | \>= 5 | Inpatient | 0.104 | 0.037 | 0.211 | 0.070 | 0.023 | 0.151 |
| Bosnia and Herzegovina | \< 5 | Outpatient | 0.071 | 0.026 | 0.144 | 0.048 | 0.016 | 0.102 |
| Bosnia and Herzegovina | \>= 5 | Outpatient | 0.080 | 0.029 | 0.162 | 0.054 | 0.018 | 0.116 |
| Botswana | \< 5 | Inpatient | 0.160 | 0.083 | 0.259 | 0.135 | 0.070 | 0.226 |
| Botswana | \>= 5 | Inpatient | 0.178 | 0.089 | 0.291 | 0.151 | 0.073 | 0.261 |
| Botswana | \< 5 | Outpatient | 0.125 | 0.066 | 0.203 | 0.106 | 0.054 | 0.179 |
| Botswana | \>= 5 | Outpatient | 0.140 | 0.072 | 0.232 | 0.119 | 0.058 | 0.205 |
| Brazil | \< 5 | Inpatient | 0.153 | 0.080 | 0.250 | 0.105 | 0.052 | 0.180 |
| Brazil | \>= 5 | Inpatient | 0.171 | 0.086 | 0.284 | 0.118 | 0.056 | 0.207 |
| Brazil | \< 5 | Outpatient | 0.120 | 0.064 | 0.199 | 0.081 | 0.041 | 0.141 |
| Brazil | \>= 5 | Outpatient | 0.134 | 0.068 | 0.229 | 0.091 | 0.043 | 0.164 |
| Brunei Darussalam | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Brunei Darussalam | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Brunei Darussalam | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Brunei Darussalam | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bulgaria | \< 5 | Inpatient | 0.095 | 0.031 | 0.208 | 0.072 | 0.022 | 0.165 |
| Bulgaria | \>= 5 | Inpatient | 0.106 | 0.035 | 0.233 | 0.081 | 0.025 | 0.185 |
| Bulgaria | \< 5 | Outpatient | 0.073 | 0.024 | 0.160 | 0.055 | 0.018 | 0.127 |
| Bulgaria | \>= 5 | Outpatient | 0.082 | 0.027 | 0.181 | 0.062 | 0.020 | 0.144 |
| Burkina Faso | \< 5 | Inpatient | 0.157 | 0.074 | 0.257 | 0.134 | 0.063 | 0.231 |
| Burkina Faso | \>= 5 | Inpatient | 0.175 | 0.080 | 0.290 | 0.150 | 0.067 | 0.266 |
| Burkina Faso | \< 5 | Outpatient | 0.124 | 0.058 | 0.202 | 0.105 | 0.048 | 0.186 |
| Burkina Faso | \>= 5 | Outpatient | 0.138 | 0.063 | 0.231 | 0.118 | 0.052 | 0.216 |
| Burundi | \< 5 | Inpatient | 0.168 | 0.099 | 0.257 | 0.147 | 0.078 | 0.242 |
| Burundi | \>= 5 | Inpatient | 0.186 | 0.105 | 0.292 | 0.164 | 0.082 | 0.278 |
| Burundi | \< 5 | Outpatient | 0.132 | 0.079 | 0.201 | 0.115 | 0.060 | 0.195 |
| Burundi | \>= 5 | Outpatient | 0.147 | 0.083 | 0.232 | 0.129 | 0.064 | 0.223 |
| Cabo Verde | \< 5 | Inpatient | 0.163 | 0.084 | 0.271 | 0.114 | 0.054 | 0.205 |
| Cabo Verde | \>= 5 | Inpatient | 0.180 | 0.091 | 0.304 | 0.128 | 0.060 | 0.230 |
| Cabo Verde | \< 5 | Outpatient | 0.128 | 0.067 | 0.215 | 0.088 | 0.043 | 0.159 |
| Cabo Verde | \>= 5 | Outpatient | 0.142 | 0.072 | 0.240 | 0.099 | 0.047 | 0.182 |
| Cambodia | \< 5 | Inpatient | 0.081 | 0.032 | 0.170 | 0.054 | 0.020 | 0.121 |
| Cambodia | \>= 5 | Inpatient | 0.091 | 0.036 | 0.191 | 0.061 | 0.022 | 0.137 |
| Cambodia | \< 5 | Outpatient | 0.062 | 0.026 | 0.130 | 0.041 | 0.016 | 0.093 |
| Cambodia | \>= 5 | Outpatient | 0.070 | 0.029 | 0.147 | 0.047 | 0.018 | 0.107 |
| Cameroon | \< 5 | Inpatient | 0.163 | 0.084 | 0.272 | 0.133 | 0.068 | 0.230 |
| Cameroon | \>= 5 | Inpatient | 0.181 | 0.091 | 0.304 | 0.148 | 0.072 | 0.262 |
| Cameroon | \< 5 | Outpatient | 0.128 | 0.066 | 0.212 | 0.103 | 0.053 | 0.182 |
| Cameroon | \>= 5 | Outpatient | 0.142 | 0.073 | 0.240 | 0.116 | 0.057 | 0.207 |
| Canada | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Canada | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Canada | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Canada | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Central African Republic | \< 5 | Inpatient | 0.169 | 0.087 | 0.276 | 0.119 | 0.056 | 0.211 |
| Central African Republic | \>= 5 | Inpatient | 0.187 | 0.093 | 0.312 | 0.133 | 0.061 | 0.239 |
| Central African Republic | \< 5 | Outpatient | 0.133 | 0.068 | 0.219 | 0.092 | 0.043 | 0.164 |
| Central African Republic | \>= 5 | Outpatient | 0.148 | 0.073 | 0.249 | 0.104 | 0.048 | 0.190 |
| Chad | \< 5 | Inpatient | 0.168 | 0.099 | 0.257 | 0.118 | 0.063 | 0.194 |
| Chad | \>= 5 | Inpatient | 0.186 | 0.105 | 0.292 | 0.132 | 0.068 | 0.223 |
| Chad | \< 5 | Outpatient | 0.132 | 0.079 | 0.201 | 0.091 | 0.050 | 0.150 |
| Chad | \>= 5 | Outpatient | 0.147 | 0.083 | 0.232 | 0.103 | 0.053 | 0.176 |
| Chile | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Chile | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Chile | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Chile | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| China | \< 5 | Inpatient | 0.053 | 0.023 | 0.099 | 0.035 | 0.015 | 0.066 |
| China | \>= 5 | Inpatient | 0.060 | 0.027 | 0.113 | 0.039 | 0.017 | 0.076 |
| China | \< 5 | Outpatient | 0.040 | 0.019 | 0.072 | 0.026 | 0.012 | 0.048 |
| China | \>= 5 | Outpatient | 0.046 | 0.021 | 0.083 | 0.030 | 0.014 | 0.055 |
| Colombia | \< 5 | Inpatient | 0.150 | 0.074 | 0.253 | 0.107 | 0.049 | 0.192 |
| Colombia | \>= 5 | Inpatient | 0.166 | 0.078 | 0.286 | 0.120 | 0.053 | 0.221 |
| Colombia | \< 5 | Outpatient | 0.117 | 0.057 | 0.202 | 0.083 | 0.038 | 0.153 |
| Colombia | \>= 5 | Outpatient | 0.131 | 0.061 | 0.230 | 0.094 | 0.041 | 0.178 |
| Comoros | \< 5 | Inpatient | 0.158 | 0.093 | 0.240 | 0.110 | 0.060 | 0.177 |
| Comoros | \>= 5 | Inpatient | 0.175 | 0.100 | 0.272 | 0.123 | 0.066 | 0.204 |
| Comoros | \< 5 | Outpatient | 0.124 | 0.075 | 0.185 | 0.085 | 0.048 | 0.136 |
| Comoros | \>= 5 | Outpatient | 0.138 | 0.080 | 0.214 | 0.096 | 0.052 | 0.158 |
| Congo, Dem. Rep. | \< 5 | Inpatient | 0.168 | 0.099 | 0.257 | 0.125 | 0.071 | 0.197 |
| Congo, Dem. Rep. | \>= 5 | Inpatient | 0.186 | 0.105 | 0.292 | 0.140 | 0.076 | 0.228 |
| Congo, Dem. Rep. | \< 5 | Outpatient | 0.132 | 0.079 | 0.201 | 0.097 | 0.056 | 0.154 |
| Congo, Dem. Rep. | \>= 5 | Outpatient | 0.147 | 0.083 | 0.232 | 0.109 | 0.060 | 0.180 |
| Congo, Rep. | \< 5 | Inpatient | 0.158 | 0.093 | 0.240 | 0.126 | 0.074 | 0.195 |
| Congo, Rep. | \>= 5 | Inpatient | 0.175 | 0.100 | 0.272 | 0.141 | 0.079 | 0.226 |
| Congo, Rep. | \< 5 | Outpatient | 0.124 | 0.075 | 0.185 | 0.098 | 0.059 | 0.151 |
| Congo, Rep. | \>= 5 | Outpatient | 0.138 | 0.080 | 0.214 | 0.110 | 0.063 | 0.176 |
| Cook Islands | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Cook Islands | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Cook Islands | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Cook Islands | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Costa Rica | \< 5 | Inpatient | 0.130 | 0.066 | 0.218 | 0.112 | 0.056 | 0.191 |
| Costa Rica | \>= 5 | Inpatient | 0.145 | 0.072 | 0.247 | 0.126 | 0.060 | 0.222 |
| Costa Rica | \< 5 | Outpatient | 0.101 | 0.052 | 0.170 | 0.087 | 0.044 | 0.150 |
| Costa Rica | \>= 5 | Outpatient | 0.114 | 0.056 | 0.195 | 0.098 | 0.046 | 0.176 |
| Côte d’Ivoire | \< 5 | Inpatient | 0.158 | 0.093 | 0.240 | 0.129 | 0.075 | 0.202 |
| Côte d’Ivoire | \>= 5 | Inpatient | 0.175 | 0.100 | 0.272 | 0.144 | 0.079 | 0.232 |
| Côte d’Ivoire | \< 5 | Outpatient | 0.124 | 0.075 | 0.185 | 0.100 | 0.059 | 0.157 |
| Côte d’Ivoire | \>= 5 | Outpatient | 0.138 | 0.080 | 0.214 | 0.112 | 0.063 | 0.181 |
| Croatia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Croatia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Croatia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Croatia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Cuba | \< 5 | Inpatient | 0.127 | 0.055 | 0.232 | 0.088 | 0.033 | 0.176 |
| Cuba | \>= 5 | Inpatient | 0.142 | 0.060 | 0.261 | 0.099 | 0.037 | 0.201 |
| Cuba | \< 5 | Outpatient | 0.099 | 0.042 | 0.180 | 0.068 | 0.026 | 0.136 |
| Cuba | \>= 5 | Outpatient | 0.111 | 0.047 | 0.208 | 0.077 | 0.029 | 0.158 |
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
| Djibouti | \< 5 | Inpatient | 0.169 | 0.091 | 0.275 | 0.140 | 0.070 | 0.238 |
| Djibouti | \>= 5 | Inpatient | 0.187 | 0.098 | 0.308 | 0.156 | 0.076 | 0.272 |
| Djibouti | \< 5 | Outpatient | 0.132 | 0.072 | 0.213 | 0.109 | 0.055 | 0.186 |
| Djibouti | \>= 5 | Outpatient | 0.147 | 0.079 | 0.243 | 0.122 | 0.059 | 0.215 |
| Dominica | \< 5 | Inpatient | 0.130 | 0.066 | 0.218 | 0.090 | 0.040 | 0.167 |
| Dominica | \>= 5 | Inpatient | 0.145 | 0.072 | 0.247 | 0.101 | 0.044 | 0.194 |
| Dominica | \< 5 | Outpatient | 0.101 | 0.052 | 0.170 | 0.069 | 0.031 | 0.128 |
| Dominica | \>= 5 | Outpatient | 0.114 | 0.056 | 0.195 | 0.078 | 0.034 | 0.148 |
| Dominican Republic | \< 5 | Inpatient | 0.130 | 0.066 | 0.218 | 0.107 | 0.055 | 0.178 |
| Dominican Republic | \>= 5 | Inpatient | 0.145 | 0.072 | 0.247 | 0.120 | 0.059 | 0.207 |
| Dominican Republic | \< 5 | Outpatient | 0.101 | 0.052 | 0.170 | 0.083 | 0.043 | 0.140 |
| Dominican Republic | \>= 5 | Outpatient | 0.114 | 0.056 | 0.195 | 0.093 | 0.046 | 0.165 |
| Ecuador | \< 5 | Inpatient | 0.160 | 0.084 | 0.259 | 0.106 | 0.055 | 0.177 |
| Ecuador | \>= 5 | Inpatient | 0.178 | 0.090 | 0.296 | 0.120 | 0.059 | 0.206 |
| Ecuador | \< 5 | Outpatient | 0.126 | 0.066 | 0.207 | 0.082 | 0.043 | 0.139 |
| Ecuador | \>= 5 | Outpatient | 0.141 | 0.071 | 0.238 | 0.093 | 0.046 | 0.164 |
| Egypt, Arab Rep. | \< 5 | Inpatient | 0.156 | 0.073 | 0.265 | 0.109 | 0.047 | 0.200 |
| Egypt, Arab Rep. | \>= 5 | Inpatient | 0.173 | 0.078 | 0.299 | 0.122 | 0.051 | 0.226 |
| Egypt, Arab Rep. | \< 5 | Outpatient | 0.122 | 0.058 | 0.208 | 0.084 | 0.037 | 0.154 |
| Egypt, Arab Rep. | \>= 5 | Outpatient | 0.136 | 0.063 | 0.236 | 0.095 | 0.041 | 0.177 |
| El Salvador | \< 5 | Inpatient | 0.159 | 0.084 | 0.255 | 0.107 | 0.055 | 0.178 |
| El Salvador | \>= 5 | Inpatient | 0.176 | 0.089 | 0.293 | 0.120 | 0.059 | 0.207 |
| El Salvador | \< 5 | Outpatient | 0.125 | 0.066 | 0.204 | 0.083 | 0.043 | 0.140 |
| El Salvador | \>= 5 | Outpatient | 0.139 | 0.071 | 0.235 | 0.093 | 0.046 | 0.165 |
| Equatorial Guinea | \< 5 | Inpatient | 0.160 | 0.083 | 0.259 | 0.112 | 0.054 | 0.193 |
| Equatorial Guinea | \>= 5 | Inpatient | 0.178 | 0.089 | 0.291 | 0.125 | 0.058 | 0.221 |
| Equatorial Guinea | \< 5 | Outpatient | 0.125 | 0.066 | 0.203 | 0.086 | 0.042 | 0.152 |
| Equatorial Guinea | \>= 5 | Outpatient | 0.140 | 0.072 | 0.232 | 0.098 | 0.046 | 0.175 |
| Eritrea | \< 5 | Inpatient | 0.168 | 0.099 | 0.257 | 0.147 | 0.078 | 0.242 |
| Eritrea | \>= 5 | Inpatient | 0.186 | 0.105 | 0.292 | 0.164 | 0.082 | 0.278 |
| Eritrea | \< 5 | Outpatient | 0.132 | 0.079 | 0.201 | 0.115 | 0.060 | 0.195 |
| Eritrea | \>= 5 | Outpatient | 0.147 | 0.083 | 0.232 | 0.129 | 0.064 | 0.223 |
| Estonia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Estonia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Estonia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Estonia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Eswatini | \< 5 | Inpatient | 0.158 | 0.093 | 0.240 | 0.127 | 0.075 | 0.197 |
| Eswatini | \>= 5 | Inpatient | 0.175 | 0.100 | 0.272 | 0.142 | 0.079 | 0.227 |
| Eswatini | \< 5 | Outpatient | 0.124 | 0.075 | 0.185 | 0.098 | 0.059 | 0.153 |
| Eswatini | \>= 5 | Outpatient | 0.138 | 0.080 | 0.214 | 0.111 | 0.063 | 0.178 |
| Ethiopia | \< 5 | Inpatient | 0.168 | 0.084 | 0.273 | 0.138 | 0.067 | 0.238 |
| Ethiopia | \>= 5 | Inpatient | 0.186 | 0.090 | 0.307 | 0.154 | 0.072 | 0.272 |
| Ethiopia | \< 5 | Outpatient | 0.132 | 0.067 | 0.219 | 0.108 | 0.052 | 0.191 |
| Ethiopia | \>= 5 | Outpatient | 0.147 | 0.071 | 0.248 | 0.121 | 0.056 | 0.219 |
| Fiji | \< 5 | Inpatient | 0.061 | 0.025 | 0.126 | 0.054 | 0.019 | 0.123 |
| Fiji | \>= 5 | Inpatient | 0.069 | 0.029 | 0.143 | 0.061 | 0.021 | 0.140 |
| Fiji | \< 5 | Outpatient | 0.047 | 0.020 | 0.092 | 0.041 | 0.015 | 0.095 |
| Fiji | \>= 5 | Outpatient | 0.053 | 0.023 | 0.106 | 0.046 | 0.016 | 0.107 |
| Finland | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Finland | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Finland | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Finland | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| France | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| France | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| France | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| France | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Gabon | \< 5 | Inpatient | 0.167 | 0.079 | 0.289 | 0.118 | 0.052 | 0.215 |
| Gabon | \>= 5 | Inpatient | 0.185 | 0.086 | 0.323 | 0.132 | 0.056 | 0.245 |
| Gabon | \< 5 | Outpatient | 0.132 | 0.064 | 0.227 | 0.091 | 0.041 | 0.169 |
| Gabon | \>= 5 | Outpatient | 0.147 | 0.069 | 0.257 | 0.103 | 0.045 | 0.193 |
| Gambia, The | \< 5 | Inpatient | 0.184 | 0.101 | 0.301 | 0.158 | 0.079 | 0.279 |
| Gambia, The | \>= 5 | Inpatient | 0.203 | 0.107 | 0.337 | 0.175 | 0.085 | 0.312 |
| Gambia, The | \< 5 | Outpatient | 0.145 | 0.081 | 0.245 | 0.124 | 0.062 | 0.226 |
| Gambia, The | \>= 5 | Outpatient | 0.161 | 0.086 | 0.275 | 0.138 | 0.066 | 0.258 |
| Georgia | \< 5 | Inpatient | 0.092 | 0.033 | 0.191 | 0.076 | 0.024 | 0.162 |
| Georgia | \>= 5 | Inpatient | 0.104 | 0.037 | 0.211 | 0.085 | 0.027 | 0.183 |
| Georgia | \< 5 | Outpatient | 0.071 | 0.026 | 0.144 | 0.058 | 0.019 | 0.125 |
| Georgia | \>= 5 | Outpatient | 0.080 | 0.029 | 0.162 | 0.066 | 0.022 | 0.142 |
| Germany | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Germany | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Germany | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Germany | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Ghana | \< 5 | Inpatient | 0.154 | 0.079 | 0.245 | 0.131 | 0.065 | 0.223 |
| Ghana | \>= 5 | Inpatient | 0.171 | 0.087 | 0.278 | 0.147 | 0.070 | 0.256 |
| Ghana | \< 5 | Outpatient | 0.120 | 0.062 | 0.190 | 0.102 | 0.051 | 0.176 |
| Ghana | \>= 5 | Outpatient | 0.134 | 0.068 | 0.220 | 0.115 | 0.055 | 0.201 |
| Greece | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Greece | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Greece | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Greece | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Grenada | \< 5 | Inpatient | 0.130 | 0.066 | 0.218 | 0.090 | 0.040 | 0.167 |
| Grenada | \>= 5 | Inpatient | 0.145 | 0.072 | 0.247 | 0.101 | 0.044 | 0.194 |
| Grenada | \< 5 | Outpatient | 0.101 | 0.052 | 0.170 | 0.069 | 0.031 | 0.128 |
| Grenada | \>= 5 | Outpatient | 0.114 | 0.056 | 0.195 | 0.078 | 0.034 | 0.148 |
| Guatemala | \< 5 | Inpatient | 0.144 | 0.071 | 0.253 | 0.112 | 0.052 | 0.208 |
| Guatemala | \>= 5 | Inpatient | 0.161 | 0.076 | 0.284 | 0.126 | 0.056 | 0.238 |
| Guatemala | \< 5 | Outpatient | 0.113 | 0.056 | 0.198 | 0.087 | 0.040 | 0.163 |
| Guatemala | \>= 5 | Outpatient | 0.126 | 0.060 | 0.226 | 0.098 | 0.043 | 0.188 |
| Guinea | \< 5 | Inpatient | 0.158 | 0.093 | 0.240 | 0.110 | 0.060 | 0.177 |
| Guinea | \>= 5 | Inpatient | 0.175 | 0.100 | 0.272 | 0.123 | 0.066 | 0.204 |
| Guinea | \< 5 | Outpatient | 0.124 | 0.075 | 0.185 | 0.085 | 0.048 | 0.136 |
| Guinea | \>= 5 | Outpatient | 0.138 | 0.080 | 0.214 | 0.096 | 0.052 | 0.158 |
| Guinea-Bissau | \< 5 | Inpatient | 0.174 | 0.090 | 0.289 | 0.147 | 0.073 | 0.260 |
| Guinea-Bissau | \>= 5 | Inpatient | 0.193 | 0.096 | 0.326 | 0.164 | 0.078 | 0.294 |
| Guinea-Bissau | \< 5 | Outpatient | 0.137 | 0.072 | 0.229 | 0.115 | 0.057 | 0.206 |
| Guinea-Bissau | \>= 5 | Outpatient | 0.153 | 0.076 | 0.260 | 0.129 | 0.060 | 0.238 |
| Guyana | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Guyana | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Guyana | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Guyana | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Haiti | \< 5 | Inpatient | 0.133 | 0.058 | 0.244 | 0.103 | 0.044 | 0.193 |
| Haiti | \>= 5 | Inpatient | 0.149 | 0.064 | 0.277 | 0.115 | 0.048 | 0.223 |
| Haiti | \< 5 | Outpatient | 0.104 | 0.046 | 0.195 | 0.080 | 0.034 | 0.152 |
| Haiti | \>= 5 | Outpatient | 0.117 | 0.050 | 0.222 | 0.090 | 0.037 | 0.176 |
| Honduras | \< 5 | Inpatient | 0.165 | 0.074 | 0.297 | 0.112 | 0.048 | 0.210 |
| Honduras | \>= 5 | Inpatient | 0.183 | 0.078 | 0.333 | 0.125 | 0.051 | 0.244 |
| Honduras | \< 5 | Outpatient | 0.130 | 0.057 | 0.240 | 0.087 | 0.037 | 0.167 |
| Honduras | \>= 5 | Outpatient | 0.145 | 0.062 | 0.273 | 0.098 | 0.039 | 0.193 |
| Hungary | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Hungary | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Hungary | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Hungary | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Iceland | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Iceland | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Iceland | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Iceland | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| India | \< 5 | Inpatient | 0.139 | 0.079 | 0.221 | 0.118 | 0.057 | 0.213 |
| India | \>= 5 | Inpatient | 0.155 | 0.086 | 0.248 | 0.132 | 0.061 | 0.240 |
| India | \< 5 | Outpatient | 0.109 | 0.063 | 0.169 | 0.092 | 0.044 | 0.166 |
| India | \>= 5 | Outpatient | 0.122 | 0.069 | 0.195 | 0.103 | 0.048 | 0.192 |
| Indonesia | \< 5 | Inpatient | 0.128 | 0.058 | 0.229 | 0.088 | 0.036 | 0.171 |
| Indonesia | \>= 5 | Inpatient | 0.142 | 0.064 | 0.254 | 0.099 | 0.040 | 0.193 |
| Indonesia | \< 5 | Outpatient | 0.099 | 0.046 | 0.178 | 0.068 | 0.028 | 0.133 |
| Indonesia | \>= 5 | Outpatient | 0.111 | 0.051 | 0.202 | 0.077 | 0.031 | 0.152 |
| Iran, Islamic Rep. | \< 5 | Inpatient | 0.177 | 0.090 | 0.297 | 0.125 | 0.058 | 0.225 |
| Iran, Islamic Rep. | \>= 5 | Inpatient | 0.195 | 0.098 | 0.330 | 0.139 | 0.064 | 0.252 |
| Iran, Islamic Rep. | \< 5 | Outpatient | 0.139 | 0.072 | 0.232 | 0.097 | 0.046 | 0.175 |
| Iran, Islamic Rep. | \>= 5 | Outpatient | 0.155 | 0.079 | 0.262 | 0.109 | 0.050 | 0.200 |
| Iraq | \< 5 | Inpatient | 0.169 | 0.091 | 0.275 | 0.130 | 0.068 | 0.220 |
| Iraq | \>= 5 | Inpatient | 0.187 | 0.098 | 0.308 | 0.145 | 0.073 | 0.250 |
| Iraq | \< 5 | Outpatient | 0.132 | 0.072 | 0.213 | 0.101 | 0.053 | 0.169 |
| Iraq | \>= 5 | Outpatient | 0.147 | 0.079 | 0.243 | 0.113 | 0.057 | 0.196 |
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
| Jamaica | \< 5 | Inpatient | 0.130 | 0.066 | 0.218 | 0.090 | 0.040 | 0.167 |
| Jamaica | \>= 5 | Inpatient | 0.145 | 0.072 | 0.247 | 0.101 | 0.044 | 0.194 |
| Jamaica | \< 5 | Outpatient | 0.101 | 0.052 | 0.170 | 0.069 | 0.031 | 0.128 |
| Jamaica | \>= 5 | Outpatient | 0.114 | 0.056 | 0.195 | 0.078 | 0.034 | 0.148 |
| Japan | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Japan | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Japan | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Japan | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Jordan | \< 5 | Inpatient | 0.169 | 0.091 | 0.275 | 0.141 | 0.070 | 0.241 |
| Jordan | \>= 5 | Inpatient | 0.187 | 0.098 | 0.308 | 0.157 | 0.076 | 0.275 |
| Jordan | \< 5 | Outpatient | 0.132 | 0.072 | 0.213 | 0.110 | 0.055 | 0.189 |
| Jordan | \>= 5 | Outpatient | 0.147 | 0.079 | 0.243 | 0.123 | 0.059 | 0.218 |
| Kazakhstan | \< 5 | Inpatient | 0.092 | 0.033 | 0.191 | 0.062 | 0.020 | 0.134 |
| Kazakhstan | \>= 5 | Inpatient | 0.104 | 0.037 | 0.211 | 0.070 | 0.023 | 0.151 |
| Kazakhstan | \< 5 | Outpatient | 0.071 | 0.026 | 0.144 | 0.048 | 0.016 | 0.102 |
| Kazakhstan | \>= 5 | Outpatient | 0.080 | 0.029 | 0.162 | 0.054 | 0.018 | 0.116 |
| Kenya | \< 5 | Inpatient | 0.159 | 0.091 | 0.244 | 0.137 | 0.073 | 0.229 |
| Kenya | \>= 5 | Inpatient | 0.176 | 0.099 | 0.274 | 0.153 | 0.078 | 0.259 |
| Kenya | \< 5 | Outpatient | 0.124 | 0.073 | 0.189 | 0.107 | 0.057 | 0.180 |
| Kenya | \>= 5 | Outpatient | 0.139 | 0.079 | 0.216 | 0.120 | 0.061 | 0.204 |
| Kiribati | \< 5 | Inpatient | 0.077 | 0.034 | 0.152 | 0.066 | 0.024 | 0.149 |
| Kiribati | \>= 5 | Inpatient | 0.087 | 0.037 | 0.173 | 0.075 | 0.027 | 0.168 |
| Kiribati | \< 5 | Outpatient | 0.059 | 0.027 | 0.116 | 0.051 | 0.019 | 0.114 |
| Kiribati | \>= 5 | Outpatient | 0.067 | 0.030 | 0.133 | 0.058 | 0.021 | 0.130 |
| Korea, Dem. People’s Rep. | \< 5 | Inpatient | 0.133 | 0.075 | 0.208 | 0.091 | 0.046 | 0.156 |
| Korea, Dem. People’s Rep. | \>= 5 | Inpatient | 0.148 | 0.082 | 0.237 | 0.103 | 0.051 | 0.179 |
| Korea, Dem. People’s Rep. | \< 5 | Outpatient | 0.103 | 0.060 | 0.159 | 0.070 | 0.037 | 0.120 |
| Korea, Dem. People’s Rep. | \>= 5 | Outpatient | 0.116 | 0.064 | 0.184 | 0.079 | 0.040 | 0.138 |
| Korea, Rep. | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Korea, Rep. | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Korea, Rep. | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Korea, Rep. | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kuwait | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kuwait | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kuwait | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kuwait | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kyrgyz Republic | \< 5 | Inpatient | 0.100 | 0.033 | 0.200 | 0.077 | 0.025 | 0.159 |
| Kyrgyz Republic | \>= 5 | Inpatient | 0.112 | 0.038 | 0.226 | 0.087 | 0.027 | 0.180 |
| Kyrgyz Republic | \< 5 | Outpatient | 0.077 | 0.027 | 0.153 | 0.059 | 0.019 | 0.121 |
| Kyrgyz Republic | \>= 5 | Outpatient | 0.087 | 0.031 | 0.173 | 0.067 | 0.022 | 0.139 |
| Lao PDR | \< 5 | Inpatient | 0.077 | 0.034 | 0.152 | 0.052 | 0.021 | 0.110 |
| Lao PDR | \>= 5 | Inpatient | 0.087 | 0.037 | 0.173 | 0.059 | 0.023 | 0.124 |
| Lao PDR | \< 5 | Outpatient | 0.059 | 0.027 | 0.116 | 0.039 | 0.016 | 0.082 |
| Lao PDR | \>= 5 | Outpatient | 0.067 | 0.030 | 0.133 | 0.045 | 0.018 | 0.095 |
| Latvia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Latvia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Latvia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Latvia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Lebanon | \< 5 | Inpatient | 0.174 | 0.085 | 0.303 | 0.123 | 0.055 | 0.232 |
| Lebanon | \>= 5 | Inpatient | 0.192 | 0.094 | 0.335 | 0.137 | 0.060 | 0.261 |
| Lebanon | \< 5 | Outpatient | 0.137 | 0.067 | 0.236 | 0.095 | 0.043 | 0.178 |
| Lebanon | \>= 5 | Outpatient | 0.152 | 0.074 | 0.268 | 0.107 | 0.047 | 0.203 |
| Lesotho | \< 5 | Inpatient | 0.158 | 0.093 | 0.240 | 0.134 | 0.075 | 0.216 |
| Lesotho | \>= 5 | Inpatient | 0.175 | 0.100 | 0.272 | 0.150 | 0.079 | 0.248 |
| Lesotho | \< 5 | Outpatient | 0.124 | 0.075 | 0.185 | 0.105 | 0.059 | 0.170 |
| Lesotho | \>= 5 | Outpatient | 0.138 | 0.080 | 0.214 | 0.118 | 0.063 | 0.195 |
| Liberia | \< 5 | Inpatient | 0.168 | 0.099 | 0.257 | 0.136 | 0.078 | 0.215 |
| Liberia | \>= 5 | Inpatient | 0.186 | 0.105 | 0.292 | 0.152 | 0.082 | 0.249 |
| Liberia | \< 5 | Outpatient | 0.132 | 0.079 | 0.201 | 0.106 | 0.061 | 0.171 |
| Liberia | \>= 5 | Outpatient | 0.147 | 0.083 | 0.232 | 0.119 | 0.064 | 0.200 |
| Libya | \< 5 | Inpatient | 0.169 | 0.080 | 0.294 | 0.140 | 0.063 | 0.257 |
| Libya | \>= 5 | Inpatient | 0.187 | 0.087 | 0.326 | 0.156 | 0.067 | 0.290 |
| Libya | \< 5 | Outpatient | 0.133 | 0.062 | 0.231 | 0.110 | 0.049 | 0.203 |
| Libya | \>= 5 | Outpatient | 0.148 | 0.068 | 0.260 | 0.123 | 0.053 | 0.230 |
| Lithuania | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Lithuania | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Lithuania | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Lithuania | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Luxembourg | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Luxembourg | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Luxembourg | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Luxembourg | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Madagascar | \< 5 | Inpatient | 0.167 | 0.082 | 0.276 | 0.136 | 0.066 | 0.234 |
| Madagascar | \>= 5 | Inpatient | 0.185 | 0.087 | 0.310 | 0.152 | 0.070 | 0.267 |
| Madagascar | \< 5 | Outpatient | 0.131 | 0.064 | 0.218 | 0.106 | 0.051 | 0.184 |
| Madagascar | \>= 5 | Outpatient | 0.147 | 0.069 | 0.248 | 0.119 | 0.055 | 0.211 |
| Malawi | \< 5 | Inpatient | 0.167 | 0.083 | 0.273 | 0.145 | 0.069 | 0.251 |
| Malawi | \>= 5 | Inpatient | 0.186 | 0.089 | 0.307 | 0.161 | 0.073 | 0.291 |
| Malawi | \< 5 | Outpatient | 0.132 | 0.066 | 0.217 | 0.113 | 0.053 | 0.204 |
| Malawi | \>= 5 | Outpatient | 0.147 | 0.070 | 0.248 | 0.127 | 0.057 | 0.233 |
| Malaysia | \< 5 | Inpatient | 0.059 | 0.022 | 0.127 | 0.039 | 0.014 | 0.088 |
| Malaysia | \>= 5 | Inpatient | 0.067 | 0.025 | 0.144 | 0.044 | 0.016 | 0.099 |
| Malaysia | \< 5 | Outpatient | 0.045 | 0.017 | 0.095 | 0.030 | 0.011 | 0.065 |
| Malaysia | \>= 5 | Outpatient | 0.051 | 0.020 | 0.107 | 0.034 | 0.013 | 0.074 |
| Maldives | \< 5 | Inpatient | 0.128 | 0.062 | 0.218 | 0.088 | 0.039 | 0.162 |
| Maldives | \>= 5 | Inpatient | 0.143 | 0.070 | 0.243 | 0.099 | 0.044 | 0.184 |
| Maldives | \< 5 | Outpatient | 0.100 | 0.050 | 0.169 | 0.068 | 0.031 | 0.126 |
| Maldives | \>= 5 | Outpatient | 0.112 | 0.055 | 0.193 | 0.077 | 0.034 | 0.146 |
| Mali | \< 5 | Inpatient | 0.173 | 0.093 | 0.279 | 0.141 | 0.074 | 0.241 |
| Mali | \>= 5 | Inpatient | 0.191 | 0.101 | 0.314 | 0.158 | 0.079 | 0.275 |
| Mali | \< 5 | Outpatient | 0.136 | 0.074 | 0.221 | 0.111 | 0.058 | 0.192 |
| Mali | \>= 5 | Outpatient | 0.151 | 0.079 | 0.250 | 0.124 | 0.061 | 0.219 |
| Malta | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Malta | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Malta | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Malta | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Marshall Islands | \< 5 | Inpatient | 0.066 | 0.027 | 0.135 | 0.047 | 0.019 | 0.103 |
| Marshall Islands | \>= 5 | Inpatient | 0.074 | 0.030 | 0.153 | 0.054 | 0.021 | 0.117 |
| Marshall Islands | \< 5 | Outpatient | 0.050 | 0.021 | 0.100 | 0.036 | 0.015 | 0.078 |
| Marshall Islands | \>= 5 | Outpatient | 0.057 | 0.024 | 0.115 | 0.041 | 0.017 | 0.088 |
| Mauritania | \< 5 | Inpatient | 0.158 | 0.093 | 0.240 | 0.133 | 0.075 | 0.213 |
| Mauritania | \>= 5 | Inpatient | 0.175 | 0.100 | 0.272 | 0.149 | 0.079 | 0.245 |
| Mauritania | \< 5 | Outpatient | 0.124 | 0.075 | 0.185 | 0.104 | 0.059 | 0.168 |
| Mauritania | \>= 5 | Outpatient | 0.138 | 0.080 | 0.214 | 0.117 | 0.063 | 0.192 |
| Mauritius | \< 5 | Inpatient | 0.160 | 0.083 | 0.259 | 0.138 | 0.070 | 0.232 |
| Mauritius | \>= 5 | Inpatient | 0.178 | 0.089 | 0.291 | 0.154 | 0.073 | 0.267 |
| Mauritius | \< 5 | Outpatient | 0.125 | 0.066 | 0.203 | 0.107 | 0.054 | 0.184 |
| Mauritius | \>= 5 | Outpatient | 0.140 | 0.072 | 0.232 | 0.121 | 0.058 | 0.212 |
| Mexico | \< 5 | Inpatient | 0.160 | 0.077 | 0.276 | 0.106 | 0.050 | 0.189 |
| Mexico | \>= 5 | Inpatient | 0.178 | 0.083 | 0.311 | 0.119 | 0.054 | 0.220 |
| Mexico | \< 5 | Outpatient | 0.126 | 0.062 | 0.220 | 0.082 | 0.039 | 0.150 |
| Mexico | \>= 5 | Outpatient | 0.140 | 0.065 | 0.250 | 0.093 | 0.042 | 0.174 |
| Micronesia, Fed. Sts. | \< 5 | Inpatient | 0.093 | 0.036 | 0.193 | 0.057 | 0.023 | 0.120 |
| Micronesia, Fed. Sts. | \>= 5 | Inpatient | 0.105 | 0.040 | 0.217 | 0.065 | 0.026 | 0.137 |
| Micronesia, Fed. Sts. | \< 5 | Outpatient | 0.072 | 0.029 | 0.148 | 0.044 | 0.018 | 0.092 |
| Micronesia, Fed. Sts. | \>= 5 | Outpatient | 0.081 | 0.032 | 0.169 | 0.050 | 0.021 | 0.106 |
| Moldova | \< 5 | Inpatient | 0.092 | 0.033 | 0.191 | 0.072 | 0.024 | 0.151 |
| Moldova | \>= 5 | Inpatient | 0.104 | 0.037 | 0.211 | 0.081 | 0.027 | 0.171 |
| Moldova | \< 5 | Outpatient | 0.071 | 0.026 | 0.144 | 0.055 | 0.019 | 0.116 |
| Moldova | \>= 5 | Outpatient | 0.080 | 0.029 | 0.162 | 0.062 | 0.021 | 0.133 |
| Monaco | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Monaco | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Monaco | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Monaco | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Mongolia | \< 5 | Inpatient | 0.077 | 0.034 | 0.152 | 0.052 | 0.021 | 0.110 |
| Mongolia | \>= 5 | Inpatient | 0.087 | 0.037 | 0.173 | 0.059 | 0.023 | 0.124 |
| Mongolia | \< 5 | Outpatient | 0.059 | 0.027 | 0.116 | 0.039 | 0.016 | 0.082 |
| Mongolia | \>= 5 | Outpatient | 0.067 | 0.030 | 0.133 | 0.045 | 0.018 | 0.095 |
| Montenegro | \< 5 | Inpatient | 0.092 | 0.033 | 0.191 | 0.062 | 0.020 | 0.134 |
| Montenegro | \>= 5 | Inpatient | 0.104 | 0.037 | 0.211 | 0.070 | 0.023 | 0.151 |
| Montenegro | \< 5 | Outpatient | 0.071 | 0.026 | 0.144 | 0.048 | 0.016 | 0.102 |
| Montenegro | \>= 5 | Outpatient | 0.080 | 0.029 | 0.162 | 0.054 | 0.018 | 0.116 |
| Morocco | \< 5 | Inpatient | 0.172 | 0.082 | 0.300 | 0.152 | 0.064 | 0.287 |
| Morocco | \>= 5 | Inpatient | 0.191 | 0.090 | 0.331 | 0.169 | 0.069 | 0.317 |
| Morocco | \< 5 | Outpatient | 0.136 | 0.065 | 0.238 | 0.119 | 0.050 | 0.232 |
| Morocco | \>= 5 | Outpatient | 0.151 | 0.070 | 0.265 | 0.133 | 0.054 | 0.257 |
| Mozambique | \< 5 | Inpatient | 0.177 | 0.096 | 0.287 | 0.142 | 0.074 | 0.242 |
| Mozambique | \>= 5 | Inpatient | 0.195 | 0.102 | 0.324 | 0.158 | 0.079 | 0.277 |
| Mozambique | \< 5 | Outpatient | 0.139 | 0.076 | 0.230 | 0.111 | 0.057 | 0.195 |
| Mozambique | \>= 5 | Outpatient | 0.155 | 0.081 | 0.262 | 0.125 | 0.061 | 0.221 |
| Myanmar | \< 5 | Inpatient | 0.133 | 0.075 | 0.208 | 0.103 | 0.053 | 0.170 |
| Myanmar | \>= 5 | Inpatient | 0.148 | 0.082 | 0.237 | 0.115 | 0.058 | 0.195 |
| Myanmar | \< 5 | Outpatient | 0.103 | 0.060 | 0.159 | 0.079 | 0.042 | 0.131 |
| Myanmar | \>= 5 | Outpatient | 0.116 | 0.064 | 0.184 | 0.090 | 0.046 | 0.153 |
| Namibia | \< 5 | Inpatient | 0.160 | 0.083 | 0.259 | 0.133 | 0.070 | 0.220 |
| Namibia | \>= 5 | Inpatient | 0.178 | 0.089 | 0.291 | 0.148 | 0.073 | 0.253 |
| Namibia | \< 5 | Outpatient | 0.125 | 0.066 | 0.203 | 0.103 | 0.054 | 0.173 |
| Namibia | \>= 5 | Outpatient | 0.140 | 0.072 | 0.232 | 0.116 | 0.058 | 0.199 |
| Nauru | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nauru | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nauru | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nauru | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nepal | \< 5 | Inpatient | 0.123 | 0.057 | 0.208 | 0.084 | 0.037 | 0.151 |
| Nepal | \>= 5 | Inpatient | 0.138 | 0.063 | 0.234 | 0.095 | 0.041 | 0.174 |
| Nepal | \< 5 | Outpatient | 0.095 | 0.045 | 0.156 | 0.065 | 0.028 | 0.115 |
| Nepal | \>= 5 | Outpatient | 0.107 | 0.049 | 0.181 | 0.073 | 0.031 | 0.132 |
| Netherlands | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Netherlands | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Netherlands | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Netherlands | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| New Zealand | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| New Zealand | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| New Zealand | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| New Zealand | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nicaragua | \< 5 | Inpatient | 0.159 | 0.064 | 0.304 | 0.110 | 0.042 | 0.223 |
| Nicaragua | \>= 5 | Inpatient | 0.177 | 0.067 | 0.338 | 0.123 | 0.045 | 0.253 |
| Nicaragua | \< 5 | Outpatient | 0.126 | 0.050 | 0.245 | 0.085 | 0.032 | 0.174 |
| Nicaragua | \>= 5 | Outpatient | 0.140 | 0.052 | 0.275 | 0.096 | 0.035 | 0.200 |
| Niger | \< 5 | Inpatient | 0.168 | 0.099 | 0.257 | 0.142 | 0.078 | 0.230 |
| Niger | \>= 5 | Inpatient | 0.186 | 0.105 | 0.292 | 0.159 | 0.083 | 0.265 |
| Niger | \< 5 | Outpatient | 0.132 | 0.079 | 0.201 | 0.111 | 0.061 | 0.183 |
| Niger | \>= 5 | Outpatient | 0.147 | 0.083 | 0.232 | 0.125 | 0.065 | 0.213 |
| Nigeria | \< 5 | Inpatient | 0.160 | 0.081 | 0.263 | 0.111 | 0.054 | 0.198 |
| Nigeria | \>= 5 | Inpatient | 0.177 | 0.087 | 0.293 | 0.125 | 0.058 | 0.223 |
| Nigeria | \< 5 | Outpatient | 0.125 | 0.065 | 0.207 | 0.086 | 0.042 | 0.154 |
| Nigeria | \>= 5 | Outpatient | 0.140 | 0.070 | 0.235 | 0.097 | 0.046 | 0.175 |
| Niue | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Niue | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Niue | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Niue | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| North Macedonia | \< 5 | Inpatient | 0.092 | 0.033 | 0.191 | 0.073 | 0.024 | 0.155 |
| North Macedonia | \>= 5 | Inpatient | 0.104 | 0.037 | 0.211 | 0.082 | 0.027 | 0.174 |
| North Macedonia | \< 5 | Outpatient | 0.071 | 0.026 | 0.144 | 0.056 | 0.019 | 0.118 |
| North Macedonia | \>= 5 | Outpatient | 0.080 | 0.029 | 0.162 | 0.063 | 0.021 | 0.135 |
| Norway | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Norway | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Norway | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Norway | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Oman | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Oman | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Oman | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Oman | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Pakistan | \< 5 | Inpatient | 0.180 | 0.094 | 0.302 | 0.155 | 0.072 | 0.283 |
| Pakistan | \>= 5 | Inpatient | 0.199 | 0.102 | 0.334 | 0.172 | 0.078 | 0.316 |
| Pakistan | \< 5 | Outpatient | 0.142 | 0.075 | 0.237 | 0.121 | 0.056 | 0.223 |
| Pakistan | \>= 5 | Outpatient | 0.158 | 0.082 | 0.267 | 0.135 | 0.060 | 0.249 |
| Palau | \< 5 | Inpatient | 0.076 | 0.027 | 0.166 | 0.051 | 0.019 | 0.115 |
| Palau | \>= 5 | Inpatient | 0.086 | 0.031 | 0.184 | 0.058 | 0.021 | 0.130 |
| Palau | \< 5 | Outpatient | 0.058 | 0.022 | 0.126 | 0.039 | 0.015 | 0.089 |
| Palau | \>= 5 | Outpatient | 0.066 | 0.024 | 0.143 | 0.044 | 0.016 | 0.099 |
| Panama | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Panama | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Panama | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Panama | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Papua New Guinea | \< 5 | Inpatient | 0.081 | 0.032 | 0.178 | 0.055 | 0.020 | 0.127 |
| Papua New Guinea | \>= 5 | Inpatient | 0.092 | 0.036 | 0.200 | 0.062 | 0.022 | 0.145 |
| Papua New Guinea | \< 5 | Outpatient | 0.062 | 0.026 | 0.137 | 0.042 | 0.016 | 0.098 |
| Papua New Guinea | \>= 5 | Outpatient | 0.070 | 0.029 | 0.155 | 0.047 | 0.018 | 0.111 |
| Paraguay | \< 5 | Inpatient | 0.147 | 0.081 | 0.230 | 0.108 | 0.055 | 0.180 |
| Paraguay | \>= 5 | Inpatient | 0.163 | 0.088 | 0.265 | 0.121 | 0.059 | 0.210 |
| Paraguay | \< 5 | Outpatient | 0.115 | 0.064 | 0.181 | 0.083 | 0.043 | 0.142 |
| Paraguay | \>= 5 | Outpatient | 0.128 | 0.069 | 0.211 | 0.094 | 0.046 | 0.166 |
| Peru | \< 5 | Inpatient | 0.148 | 0.077 | 0.241 | 0.103 | 0.050 | 0.178 |
| Peru | \>= 5 | Inpatient | 0.165 | 0.082 | 0.277 | 0.116 | 0.054 | 0.207 |
| Peru | \< 5 | Outpatient | 0.116 | 0.060 | 0.192 | 0.080 | 0.039 | 0.140 |
| Peru | \>= 5 | Outpatient | 0.130 | 0.064 | 0.220 | 0.090 | 0.041 | 0.162 |
| Philippines | \< 5 | Inpatient | 0.077 | 0.034 | 0.152 | 0.052 | 0.021 | 0.110 |
| Philippines | \>= 5 | Inpatient | 0.087 | 0.037 | 0.173 | 0.059 | 0.023 | 0.124 |
| Philippines | \< 5 | Outpatient | 0.059 | 0.027 | 0.116 | 0.039 | 0.016 | 0.082 |
| Philippines | \>= 5 | Outpatient | 0.067 | 0.030 | 0.133 | 0.045 | 0.018 | 0.095 |
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
| Russian Federation | \< 5 | Inpatient | 0.092 | 0.033 | 0.191 | 0.062 | 0.020 | 0.134 |
| Russian Federation | \>= 5 | Inpatient | 0.104 | 0.037 | 0.211 | 0.070 | 0.023 | 0.151 |
| Russian Federation | \< 5 | Outpatient | 0.071 | 0.026 | 0.144 | 0.048 | 0.016 | 0.102 |
| Russian Federation | \>= 5 | Outpatient | 0.080 | 0.029 | 0.162 | 0.054 | 0.018 | 0.116 |
| Rwanda | \< 5 | Inpatient | 0.168 | 0.099 | 0.257 | 0.145 | 0.078 | 0.237 |
| Rwanda | \>= 5 | Inpatient | 0.186 | 0.105 | 0.292 | 0.162 | 0.082 | 0.272 |
| Rwanda | \< 5 | Outpatient | 0.132 | 0.079 | 0.201 | 0.114 | 0.060 | 0.190 |
| Rwanda | \>= 5 | Outpatient | 0.147 | 0.083 | 0.232 | 0.127 | 0.064 | 0.219 |
| Samoa | \< 5 | Inpatient | 0.077 | 0.034 | 0.152 | 0.052 | 0.021 | 0.110 |
| Samoa | \>= 5 | Inpatient | 0.087 | 0.037 | 0.173 | 0.059 | 0.023 | 0.124 |
| Samoa | \< 5 | Outpatient | 0.059 | 0.027 | 0.116 | 0.039 | 0.016 | 0.082 |
| Samoa | \>= 5 | Outpatient | 0.067 | 0.030 | 0.133 | 0.045 | 0.018 | 0.095 |
| San Marino | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| San Marino | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| San Marino | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| San Marino | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| São Tomé and Principe | \< 5 | Inpatient | 0.158 | 0.093 | 0.240 | 0.134 | 0.075 | 0.216 |
| São Tomé and Principe | \>= 5 | Inpatient | 0.175 | 0.100 | 0.272 | 0.150 | 0.079 | 0.248 |
| São Tomé and Principe | \< 5 | Outpatient | 0.124 | 0.075 | 0.185 | 0.105 | 0.059 | 0.170 |
| São Tomé and Principe | \>= 5 | Outpatient | 0.138 | 0.080 | 0.214 | 0.118 | 0.063 | 0.195 |
| Saudi Arabia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Saudi Arabia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Saudi Arabia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Saudi Arabia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Senegal | \< 5 | Inpatient | 0.166 | 0.087 | 0.280 | 0.143 | 0.070 | 0.258 |
| Senegal | \>= 5 | Inpatient | 0.183 | 0.095 | 0.311 | 0.160 | 0.076 | 0.291 |
| Senegal | \< 5 | Outpatient | 0.130 | 0.070 | 0.222 | 0.112 | 0.055 | 0.209 |
| Senegal | \>= 5 | Outpatient | 0.145 | 0.075 | 0.246 | 0.126 | 0.059 | 0.232 |
| Serbia | \< 5 | Inpatient | 0.092 | 0.033 | 0.191 | 0.062 | 0.020 | 0.134 |
| Serbia | \>= 5 | Inpatient | 0.104 | 0.037 | 0.211 | 0.070 | 0.023 | 0.151 |
| Serbia | \< 5 | Outpatient | 0.071 | 0.026 | 0.144 | 0.048 | 0.016 | 0.102 |
| Serbia | \>= 5 | Outpatient | 0.080 | 0.029 | 0.162 | 0.054 | 0.018 | 0.116 |
| Seychelles | \< 5 | Inpatient | 0.160 | 0.083 | 0.259 | 0.140 | 0.069 | 0.238 |
| Seychelles | \>= 5 | Inpatient | 0.178 | 0.089 | 0.291 | 0.156 | 0.073 | 0.272 |
| Seychelles | \< 5 | Outpatient | 0.125 | 0.066 | 0.203 | 0.109 | 0.054 | 0.190 |
| Seychelles | \>= 5 | Outpatient | 0.140 | 0.072 | 0.232 | 0.123 | 0.058 | 0.217 |
| Sierra Leone | \< 5 | Inpatient | 0.168 | 0.099 | 0.257 | 0.145 | 0.078 | 0.236 |
| Sierra Leone | \>= 5 | Inpatient | 0.186 | 0.105 | 0.292 | 0.161 | 0.082 | 0.272 |
| Sierra Leone | \< 5 | Outpatient | 0.132 | 0.079 | 0.201 | 0.113 | 0.060 | 0.189 |
| Sierra Leone | \>= 5 | Outpatient | 0.147 | 0.083 | 0.232 | 0.127 | 0.064 | 0.218 |
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
| Solomon Islands | \< 5 | Inpatient | 0.077 | 0.034 | 0.152 | 0.058 | 0.024 | 0.120 |
| Solomon Islands | \>= 5 | Inpatient | 0.087 | 0.037 | 0.173 | 0.065 | 0.026 | 0.137 |
| Solomon Islands | \< 5 | Outpatient | 0.059 | 0.027 | 0.116 | 0.044 | 0.018 | 0.092 |
| Solomon Islands | \>= 5 | Outpatient | 0.067 | 0.030 | 0.133 | 0.050 | 0.021 | 0.106 |
| Somalia | \< 5 | Inpatient | 0.169 | 0.070 | 0.314 | 0.119 | 0.046 | 0.241 |
| Somalia | \>= 5 | Inpatient | 0.186 | 0.078 | 0.351 | 0.133 | 0.051 | 0.273 |
| Somalia | \< 5 | Outpatient | 0.133 | 0.055 | 0.251 | 0.093 | 0.036 | 0.189 |
| Somalia | \>= 5 | Outpatient | 0.148 | 0.061 | 0.285 | 0.104 | 0.040 | 0.218 |
| South Africa | \< 5 | Inpatient | 0.178 | 0.090 | 0.288 | 0.132 | 0.064 | 0.225 |
| South Africa | \>= 5 | Inpatient | 0.197 | 0.096 | 0.325 | 0.147 | 0.068 | 0.259 |
| South Africa | \< 5 | Outpatient | 0.140 | 0.073 | 0.227 | 0.103 | 0.051 | 0.175 |
| South Africa | \>= 5 | Outpatient | 0.156 | 0.077 | 0.259 | 0.115 | 0.054 | 0.204 |
| South Sudan | \< 5 | Inpatient | 0.168 | 0.099 | 0.257 | 0.118 | 0.063 | 0.194 |
| South Sudan | \>= 5 | Inpatient | 0.186 | 0.105 | 0.292 | 0.132 | 0.068 | 0.223 |
| South Sudan | \< 5 | Outpatient | 0.132 | 0.079 | 0.201 | 0.091 | 0.050 | 0.150 |
| South Sudan | \>= 5 | Outpatient | 0.147 | 0.083 | 0.232 | 0.103 | 0.053 | 0.176 |
| Spain | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Spain | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Spain | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Spain | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Sri Lanka | \< 5 | Inpatient | 0.133 | 0.075 | 0.208 | 0.091 | 0.046 | 0.156 |
| Sri Lanka | \>= 5 | Inpatient | 0.148 | 0.082 | 0.237 | 0.103 | 0.051 | 0.179 |
| Sri Lanka | \< 5 | Outpatient | 0.103 | 0.060 | 0.159 | 0.070 | 0.037 | 0.120 |
| Sri Lanka | \>= 5 | Outpatient | 0.116 | 0.064 | 0.184 | 0.079 | 0.040 | 0.138 |
| St. Kitts and Nevis | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| St. Kitts and Nevis | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| St. Kitts and Nevis | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| St. Kitts and Nevis | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| St. Lucia | \< 5 | Inpatient | 0.130 | 0.066 | 0.218 | 0.090 | 0.040 | 0.167 |
| St. Lucia | \>= 5 | Inpatient | 0.145 | 0.072 | 0.247 | 0.101 | 0.044 | 0.194 |
| St. Lucia | \< 5 | Outpatient | 0.101 | 0.052 | 0.170 | 0.069 | 0.031 | 0.128 |
| St. Lucia | \>= 5 | Outpatient | 0.114 | 0.056 | 0.195 | 0.078 | 0.034 | 0.148 |
| St. Vincent and the Grenadines | \< 5 | Inpatient | 0.130 | 0.066 | 0.218 | 0.090 | 0.040 | 0.167 |
| St. Vincent and the Grenadines | \>= 5 | Inpatient | 0.145 | 0.072 | 0.247 | 0.101 | 0.044 | 0.194 |
| St. Vincent and the Grenadines | \< 5 | Outpatient | 0.101 | 0.052 | 0.170 | 0.069 | 0.031 | 0.128 |
| St. Vincent and the Grenadines | \>= 5 | Outpatient | 0.114 | 0.056 | 0.195 | 0.078 | 0.034 | 0.148 |
| Sudan | \< 5 | Inpatient | 0.173 | 0.067 | 0.333 | 0.151 | 0.054 | 0.306 |
| Sudan | \>= 5 | Inpatient | 0.191 | 0.073 | 0.369 | 0.167 | 0.060 | 0.338 |
| Sudan | \< 5 | Outpatient | 0.137 | 0.052 | 0.265 | 0.118 | 0.042 | 0.248 |
| Sudan | \>= 5 | Outpatient | 0.152 | 0.058 | 0.298 | 0.132 | 0.047 | 0.276 |
| Suriname | \< 5 | Inpatient | 0.130 | 0.066 | 0.218 | 0.090 | 0.040 | 0.167 |
| Suriname | \>= 5 | Inpatient | 0.145 | 0.072 | 0.247 | 0.101 | 0.044 | 0.194 |
| Suriname | \< 5 | Outpatient | 0.101 | 0.052 | 0.170 | 0.069 | 0.031 | 0.128 |
| Suriname | \>= 5 | Outpatient | 0.114 | 0.056 | 0.195 | 0.078 | 0.034 | 0.148 |
| Sweden | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Sweden | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Sweden | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Sweden | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Switzerland | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Switzerland | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Switzerland | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Switzerland | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Syrian Arab Republic | \< 5 | Inpatient | 0.169 | 0.070 | 0.314 | 0.119 | 0.046 | 0.241 |
| Syrian Arab Republic | \>= 5 | Inpatient | 0.186 | 0.078 | 0.351 | 0.133 | 0.051 | 0.273 |
| Syrian Arab Republic | \< 5 | Outpatient | 0.133 | 0.055 | 0.251 | 0.093 | 0.036 | 0.189 |
| Syrian Arab Republic | \>= 5 | Outpatient | 0.148 | 0.061 | 0.285 | 0.104 | 0.040 | 0.218 |
| Tajikistan | \< 5 | Inpatient | 0.100 | 0.033 | 0.200 | 0.087 | 0.025 | 0.186 |
| Tajikistan | \>= 5 | Inpatient | 0.112 | 0.038 | 0.226 | 0.098 | 0.028 | 0.212 |
| Tajikistan | \< 5 | Outpatient | 0.077 | 0.027 | 0.153 | 0.067 | 0.020 | 0.144 |
| Tajikistan | \>= 5 | Outpatient | 0.087 | 0.031 | 0.173 | 0.075 | 0.023 | 0.165 |
| Tanzania | \< 5 | Inpatient | 0.154 | 0.082 | 0.243 | 0.132 | 0.067 | 0.221 |
| Tanzania | \>= 5 | Inpatient | 0.171 | 0.088 | 0.277 | 0.147 | 0.071 | 0.256 |
| Tanzania | \< 5 | Outpatient | 0.121 | 0.066 | 0.191 | 0.103 | 0.052 | 0.175 |
| Tanzania | \>= 5 | Outpatient | 0.135 | 0.070 | 0.219 | 0.116 | 0.056 | 0.203 |
| Thailand | \< 5 | Inpatient | 0.129 | 0.058 | 0.227 | 0.089 | 0.036 | 0.171 |
| Thailand | \>= 5 | Inpatient | 0.144 | 0.065 | 0.253 | 0.100 | 0.041 | 0.193 |
| Thailand | \< 5 | Outpatient | 0.100 | 0.046 | 0.176 | 0.068 | 0.029 | 0.132 |
| Thailand | \>= 5 | Outpatient | 0.112 | 0.052 | 0.200 | 0.077 | 0.032 | 0.151 |
| Timor-Leste | \< 5 | Inpatient | 0.133 | 0.075 | 0.208 | 0.091 | 0.046 | 0.156 |
| Timor-Leste | \>= 5 | Inpatient | 0.148 | 0.082 | 0.237 | 0.103 | 0.051 | 0.179 |
| Timor-Leste | \< 5 | Outpatient | 0.103 | 0.060 | 0.159 | 0.070 | 0.037 | 0.120 |
| Timor-Leste | \>= 5 | Outpatient | 0.116 | 0.064 | 0.184 | 0.079 | 0.040 | 0.138 |
| Togo | \< 5 | Inpatient | 0.168 | 0.099 | 0.257 | 0.140 | 0.078 | 0.225 |
| Togo | \>= 5 | Inpatient | 0.186 | 0.105 | 0.292 | 0.157 | 0.082 | 0.260 |
| Togo | \< 5 | Outpatient | 0.132 | 0.079 | 0.201 | 0.110 | 0.061 | 0.179 |
| Togo | \>= 5 | Outpatient | 0.147 | 0.083 | 0.232 | 0.123 | 0.065 | 0.209 |
| Tonga | \< 5 | Inpatient | 0.061 | 0.025 | 0.126 | 0.041 | 0.016 | 0.087 |
| Tonga | \>= 5 | Inpatient | 0.069 | 0.029 | 0.143 | 0.046 | 0.019 | 0.098 |
| Tonga | \< 5 | Outpatient | 0.047 | 0.020 | 0.092 | 0.031 | 0.013 | 0.064 |
| Tonga | \>= 5 | Outpatient | 0.053 | 0.023 | 0.106 | 0.035 | 0.015 | 0.073 |
| Trinidad and Tobago | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Trinidad and Tobago | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Trinidad and Tobago | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Trinidad and Tobago | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Tunisia | \< 5 | Inpatient | 0.169 | 0.091 | 0.275 | 0.118 | 0.058 | 0.207 |
| Tunisia | \>= 5 | Inpatient | 0.187 | 0.098 | 0.308 | 0.132 | 0.064 | 0.237 |
| Tunisia | \< 5 | Outpatient | 0.132 | 0.072 | 0.213 | 0.092 | 0.046 | 0.158 |
| Tunisia | \>= 5 | Outpatient | 0.147 | 0.079 | 0.243 | 0.103 | 0.051 | 0.184 |
| Türkiye | \< 5 | Inpatient | 0.085 | 0.030 | 0.178 | 0.057 | 0.019 | 0.126 |
| Türkiye | \>= 5 | Inpatient | 0.096 | 0.033 | 0.199 | 0.065 | 0.021 | 0.140 |
| Türkiye | \< 5 | Outpatient | 0.065 | 0.023 | 0.136 | 0.043 | 0.015 | 0.095 |
| Türkiye | \>= 5 | Outpatient | 0.074 | 0.026 | 0.154 | 0.049 | 0.017 | 0.109 |
| Turkmenistan | \< 5 | Inpatient | 0.092 | 0.033 | 0.191 | 0.080 | 0.025 | 0.177 |
| Turkmenistan | \>= 5 | Inpatient | 0.104 | 0.037 | 0.211 | 0.090 | 0.027 | 0.197 |
| Turkmenistan | \< 5 | Outpatient | 0.071 | 0.026 | 0.144 | 0.062 | 0.019 | 0.138 |
| Turkmenistan | \>= 5 | Outpatient | 0.080 | 0.029 | 0.162 | 0.070 | 0.021 | 0.155 |
| Tuvalu | \< 5 | Inpatient | 0.061 | 0.025 | 0.126 | 0.041 | 0.016 | 0.087 |
| Tuvalu | \>= 5 | Inpatient | 0.069 | 0.029 | 0.143 | 0.046 | 0.019 | 0.098 |
| Tuvalu | \< 5 | Outpatient | 0.047 | 0.020 | 0.092 | 0.031 | 0.013 | 0.064 |
| Tuvalu | \>= 5 | Outpatient | 0.053 | 0.023 | 0.106 | 0.035 | 0.015 | 0.073 |
| Uganda | \< 5 | Inpatient | 0.168 | 0.099 | 0.257 | 0.144 | 0.078 | 0.234 |
| Uganda | \>= 5 | Inpatient | 0.186 | 0.105 | 0.292 | 0.161 | 0.082 | 0.270 |
| Uganda | \< 5 | Outpatient | 0.132 | 0.079 | 0.201 | 0.113 | 0.060 | 0.187 |
| Uganda | \>= 5 | Outpatient | 0.147 | 0.083 | 0.232 | 0.126 | 0.064 | 0.217 |
| Ukraine | \< 5 | Inpatient | 0.100 | 0.033 | 0.200 | 0.068 | 0.021 | 0.143 |
| Ukraine | \>= 5 | Inpatient | 0.112 | 0.038 | 0.226 | 0.076 | 0.024 | 0.161 |
| Ukraine | \< 5 | Outpatient | 0.077 | 0.027 | 0.153 | 0.052 | 0.017 | 0.108 |
| Ukraine | \>= 5 | Outpatient | 0.087 | 0.031 | 0.173 | 0.059 | 0.019 | 0.124 |
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
| Uzbekistan | \< 5 | Inpatient | 0.100 | 0.033 | 0.200 | 0.086 | 0.026 | 0.182 |
| Uzbekistan | \>= 5 | Inpatient | 0.112 | 0.038 | 0.226 | 0.096 | 0.028 | 0.207 |
| Uzbekistan | \< 5 | Outpatient | 0.077 | 0.027 | 0.153 | 0.066 | 0.020 | 0.141 |
| Uzbekistan | \>= 5 | Outpatient | 0.087 | 0.031 | 0.173 | 0.074 | 0.023 | 0.161 |
| Vanuatu | \< 5 | Inpatient | 0.077 | 0.034 | 0.152 | 0.052 | 0.021 | 0.110 |
| Vanuatu | \>= 5 | Inpatient | 0.087 | 0.037 | 0.173 | 0.059 | 0.023 | 0.124 |
| Vanuatu | \< 5 | Outpatient | 0.059 | 0.027 | 0.116 | 0.039 | 0.016 | 0.082 |
| Vanuatu | \>= 5 | Outpatient | 0.067 | 0.030 | 0.133 | 0.045 | 0.018 | 0.095 |
| Venezuela | \< 5 | Inpatient | 0.156 | 0.070 | 0.290 | 0.099 | 0.036 | 0.210 |
| Venezuela | \>= 5 | Inpatient | 0.173 | 0.076 | 0.326 | 0.111 | 0.040 | 0.238 |
| Venezuela | \< 5 | Outpatient | 0.123 | 0.054 | 0.236 | 0.077 | 0.028 | 0.168 |
| Venezuela | \>= 5 | Outpatient | 0.137 | 0.059 | 0.265 | 0.086 | 0.031 | 0.191 |
| Vietnam | \< 5 | Inpatient | 0.075 | 0.033 | 0.143 | 0.050 | 0.020 | 0.104 |
| Vietnam | \>= 5 | Inpatient | 0.084 | 0.036 | 0.166 | 0.057 | 0.022 | 0.120 |
| Vietnam | \< 5 | Outpatient | 0.057 | 0.026 | 0.109 | 0.038 | 0.016 | 0.078 |
| Vietnam | \>= 5 | Outpatient | 0.065 | 0.029 | 0.128 | 0.043 | 0.017 | 0.092 |
| Yemen, Rep. | \< 5 | Inpatient | 0.169 | 0.070 | 0.314 | 0.134 | 0.054 | 0.261 |
| Yemen, Rep. | \>= 5 | Inpatient | 0.186 | 0.078 | 0.351 | 0.150 | 0.059 | 0.293 |
| Yemen, Rep. | \< 5 | Outpatient | 0.133 | 0.055 | 0.251 | 0.105 | 0.042 | 0.206 |
| Yemen, Rep. | \>= 5 | Outpatient | 0.148 | 0.061 | 0.285 | 0.117 | 0.047 | 0.237 |
| Zambia | \< 5 | Inpatient | 0.165 | 0.085 | 0.277 | 0.141 | 0.070 | 0.258 |
| Zambia | \>= 5 | Inpatient | 0.183 | 0.093 | 0.308 | 0.157 | 0.074 | 0.288 |
| Zambia | \< 5 | Outpatient | 0.130 | 0.069 | 0.217 | 0.111 | 0.055 | 0.206 |
| Zambia | \>= 5 | Outpatient | 0.144 | 0.074 | 0.243 | 0.124 | 0.059 | 0.230 |
| Zimbabwe | \< 5 | Inpatient | 0.157 | 0.079 | 0.259 | 0.135 | 0.065 | 0.234 |
| Zimbabwe | \>= 5 | Inpatient | 0.175 | 0.085 | 0.290 | 0.150 | 0.070 | 0.268 |
| Zimbabwe | \< 5 | Outpatient | 0.123 | 0.063 | 0.203 | 0.105 | 0.052 | 0.186 |
| Zimbabwe | \>= 5 | Outpatient | 0.138 | 0.068 | 0.231 | 0.118 | 0.055 | 0.212 |

Estimated `r params$Pathogen` proportion by country, 2010 vs 2020

# Session info

``` r
sessioninfo::session_info()
```

    ## Warning in system2("quarto", "-V", stdout = TRUE, env = paste0("TMPDIR=", : running command '"quarto"
    ## TMPDIR=C:/Users/fbbu6966/AppData/Local/Temp/RtmpKiRM7N/file3b207c8028d7 -V' had status 1

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
    ##  quarto   ERROR: Unknown command "TMPDIR=C:/Users/fbbu6966/AppData/Local/Temp/RtmpKiRM7N/file3b207c8028d7". Did you mean command "update"? @ C:\\Users\\fbbu6966\\AppData\\Local\\Programs\\RStudio\\RESOUR~1\\app\\bin\\quarto\\bin\\quarto.exe
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
