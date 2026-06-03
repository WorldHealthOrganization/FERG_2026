Global proportion of Norovirus (Diarrheal) • ASYMPT+ROTA • Estimate
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

    ## Warning: There were 20 divergent transitions after warmup. Increasing adapt_delta above 0.95 may help. See
    ## http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup

    ##  Family: gaussian 
    ##   Links: mu = identity; sigma = identity 
    ## Formula: yi | se(sei) ~ 1 + YEAR + SYNDROMTYPE + AGE + REFERENCE + COVERAGE + STRAIN + (1 | REG2) + (1 | REG2:SUB2) + (1 | REG2:SUB2:COUNTRY) + (1 | REG2:SUB2:COUNTRY:ID) + (1 | REG2:SUB2:COUNTRY:ID:DTP_ID) 
    ##    Data: subset(es, as.integer(FLAG) == 1) (Number of observations: 503) 
    ##   Draws: 5 chains, each with iter = 5000; warmup = 3000; thin = 1;
    ##          total post-warmup draws = 10000
    ## 
    ## Multilevel Hyperparameters:
    ## ~REG2 (Number of levels: 6) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.17      0.15     0.01     0.56 1.00     5504     5992
    ## 
    ## ~REG2:SUB2 (Number of levels: 11) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.14      0.11     0.00     0.41 1.00     4657     6686
    ## 
    ## ~REG2:SUB2:COUNTRY (Number of levels: 41) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.17      0.11     0.01     0.40 1.00     1713     4179
    ## 
    ## ~REG2:SUB2:COUNTRY:ID (Number of levels: 134) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.62      0.06     0.51     0.75 1.00     3589     5011
    ## 
    ## ~REG2:SUB2:COUNTRY:ID:DTP_ID (Number of levels: 503) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.43      0.03     0.38     0.48 1.00     3252     6045
    ## 
    ## Regression Coefficients:
    ##                       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept               -40.40     32.08  -104.41    21.85 1.00     6510     7290
    ## YEAR                      0.02      0.02    -0.01     0.05 1.00     6503     7291
    ## SYNDROMTYPEInpatient      0.48      0.09     0.31     0.64 1.00     7398     7534
    ## SYNDROMTYPEOutpatient     0.36      0.07     0.22     0.49 1.00     6565     6960
    ## AGEAgebelow5              0.42      0.11     0.19     0.64 1.00     7480     7459
    ## AGEMixedages              0.15      0.16    -0.16     0.47 1.00     8376     7245
    ## REFERENCEOther           -0.48      0.34    -1.15     0.20 1.00     8689     7813
    ## COVERAGE                  0.00      0.00    -0.00     0.00 1.00     7056     7308
    ## STRAINGII                -0.50      0.36    -1.21     0.20 1.00     8115     6869
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

    ## [1] 5412

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

    ## [1] 11088

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

    ## [1] 1188

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
    ##  $ VAL_MEAN     : num  0.032 0.0465 0.0325 0.0472 0.033 ...
    ##  $ VAL_LWR      : num  0.017 0.0254 0.0175 0.0262 0.0179 ...
    ##  $ VAL_UPR      : num  0.0537 0.0758 0.0539 0.0759 0.0543 ...
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
    ##  $ VAL_MEAN     : num  0.032 0.0465 0.0325 0.0472 0.033 ...
    ##  $ VAL_LWR      : num  0.017 0.0254 0.0175 0.0262 0.0179 ...
    ##  $ VAL_UPR      : num  0.0537 0.0758 0.0539 0.0759 0.0543 ...
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

\[1\] 0.00 0.01 0.02 0.03 0.04 0.05 0.06

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

\[1\] 0.00 0.01 0.02 0.03 0.04

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

\[1\] 0.00 0.01 0.02 0.03 0.04 0.05 0.06

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

\[1\] 0.00 0.01 0.02 0.03 0.04 0.05 0.06

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

\[1\] 0.00 0.01 0.02 0.03 0.04 0.05

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

\[1\] 0.00 0.02 0.04 0.06 0.08 0.10

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

\[1\] 0.00 0.01 0.02 0.03 0.04 0.05 0.06

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
| Afghanistan | \< 5 | Inpatient | 0.054 | 0.032 | 0.083 | 0.064 | 0.036 | 0.101 |
| Afghanistan | \>= 5 | Inpatient | 0.038 | 0.022 | 0.059 | 0.045 | 0.024 | 0.073 |
| Afghanistan | \< 5 | Outpatient | 0.039 | 0.022 | 0.059 | 0.046 | 0.025 | 0.072 |
| Afghanistan | \>= 5 | Outpatient | 0.027 | 0.015 | 0.041 | 0.032 | 0.017 | 0.052 |
| Albania | \< 5 | Inpatient | 0.058 | 0.030 | 0.097 | 0.069 | 0.033 | 0.121 |
| Albania | \>= 5 | Inpatient | 0.040 | 0.020 | 0.070 | 0.049 | 0.022 | 0.089 |
| Albania | \< 5 | Outpatient | 0.041 | 0.021 | 0.068 | 0.049 | 0.023 | 0.086 |
| Albania | \>= 5 | Outpatient | 0.028 | 0.014 | 0.049 | 0.034 | 0.015 | 0.063 |
| Algeria | \< 5 | Inpatient | 0.050 | 0.030 | 0.073 | 0.059 | 0.033 | 0.092 |
| Algeria | \>= 5 | Inpatient | 0.035 | 0.020 | 0.053 | 0.041 | 0.022 | 0.066 |
| Algeria | \< 5 | Outpatient | 0.036 | 0.021 | 0.052 | 0.042 | 0.023 | 0.065 |
| Algeria | \>= 5 | Outpatient | 0.024 | 0.014 | 0.037 | 0.029 | 0.015 | 0.046 |
| Andorra | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Andorra | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Andorra | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Andorra | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Angola | \< 5 | Inpatient | 0.050 | 0.030 | 0.073 | 0.059 | 0.034 | 0.090 |
| Angola | \>= 5 | Inpatient | 0.035 | 0.020 | 0.053 | 0.041 | 0.023 | 0.065 |
| Angola | \< 5 | Outpatient | 0.036 | 0.021 | 0.052 | 0.042 | 0.024 | 0.063 |
| Angola | \>= 5 | Outpatient | 0.024 | 0.014 | 0.037 | 0.029 | 0.016 | 0.046 |
| Antigua and Barbuda | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Antigua and Barbuda | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Antigua and Barbuda | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Antigua and Barbuda | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Argentina | \< 5 | Inpatient | 0.058 | 0.033 | 0.091 | 0.068 | 0.038 | 0.110 |
| Argentina | \>= 5 | Inpatient | 0.040 | 0.022 | 0.065 | 0.048 | 0.026 | 0.080 |
| Argentina | \< 5 | Outpatient | 0.041 | 0.023 | 0.065 | 0.049 | 0.027 | 0.077 |
| Argentina | \>= 5 | Outpatient | 0.028 | 0.016 | 0.046 | 0.034 | 0.018 | 0.056 |
| Armenia | \< 5 | Inpatient | 0.057 | 0.033 | 0.091 | 0.069 | 0.036 | 0.114 |
| Armenia | \>= 5 | Inpatient | 0.040 | 0.022 | 0.065 | 0.048 | 0.024 | 0.083 |
| Armenia | \< 5 | Outpatient | 0.041 | 0.023 | 0.064 | 0.049 | 0.025 | 0.081 |
| Armenia | \>= 5 | Outpatient | 0.028 | 0.015 | 0.046 | 0.034 | 0.017 | 0.059 |
| Australia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Australia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Australia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Australia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Austria | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Austria | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Austria | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Austria | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Azerbaijan | \< 5 | Inpatient | 0.057 | 0.033 | 0.091 | 0.068 | 0.036 | 0.113 |
| Azerbaijan | \>= 5 | Inpatient | 0.040 | 0.022 | 0.065 | 0.047 | 0.024 | 0.082 |
| Azerbaijan | \< 5 | Outpatient | 0.041 | 0.023 | 0.064 | 0.048 | 0.025 | 0.079 |
| Azerbaijan | \>= 5 | Outpatient | 0.028 | 0.015 | 0.046 | 0.034 | 0.017 | 0.057 |
| Bahamas, The | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahamas, The | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahamas, The | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahamas, The | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahrain | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahrain | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahrain | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahrain | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bangladesh | \< 5 | Inpatient | 0.059 | 0.033 | 0.093 | 0.069 | 0.038 | 0.112 |
| Bangladesh | \>= 5 | Inpatient | 0.041 | 0.023 | 0.066 | 0.048 | 0.026 | 0.080 |
| Bangladesh | \< 5 | Outpatient | 0.042 | 0.023 | 0.066 | 0.049 | 0.026 | 0.080 |
| Bangladesh | \>= 5 | Outpatient | 0.029 | 0.016 | 0.047 | 0.034 | 0.018 | 0.057 |
| Barbados | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Barbados | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Barbados | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Barbados | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belarus | \< 5 | Inpatient | 0.057 | 0.033 | 0.091 | 0.068 | 0.036 | 0.113 |
| Belarus | \>= 5 | Inpatient | 0.040 | 0.022 | 0.065 | 0.047 | 0.024 | 0.082 |
| Belarus | \< 5 | Outpatient | 0.041 | 0.023 | 0.064 | 0.048 | 0.025 | 0.079 |
| Belarus | \>= 5 | Outpatient | 0.028 | 0.015 | 0.046 | 0.034 | 0.017 | 0.057 |
| Belgium | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belgium | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belgium | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belgium | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belize | \< 5 | Inpatient | 0.058 | 0.033 | 0.091 | 0.068 | 0.035 | 0.117 |
| Belize | \>= 5 | Inpatient | 0.040 | 0.022 | 0.065 | 0.048 | 0.024 | 0.085 |
| Belize | \< 5 | Outpatient | 0.041 | 0.023 | 0.065 | 0.049 | 0.025 | 0.084 |
| Belize | \>= 5 | Outpatient | 0.028 | 0.016 | 0.046 | 0.034 | 0.017 | 0.060 |
| Benin | \< 5 | Inpatient | 0.050 | 0.030 | 0.073 | 0.060 | 0.034 | 0.091 |
| Benin | \>= 5 | Inpatient | 0.035 | 0.020 | 0.053 | 0.041 | 0.023 | 0.066 |
| Benin | \< 5 | Outpatient | 0.036 | 0.021 | 0.052 | 0.042 | 0.024 | 0.064 |
| Benin | \>= 5 | Outpatient | 0.024 | 0.014 | 0.037 | 0.029 | 0.016 | 0.046 |
| Bhutan | \< 5 | Inpatient | 0.052 | 0.031 | 0.075 | 0.061 | 0.034 | 0.096 |
| Bhutan | \>= 5 | Inpatient | 0.036 | 0.021 | 0.054 | 0.043 | 0.023 | 0.069 |
| Bhutan | \< 5 | Outpatient | 0.037 | 0.022 | 0.055 | 0.044 | 0.024 | 0.069 |
| Bhutan | \>= 5 | Outpatient | 0.025 | 0.014 | 0.039 | 0.030 | 0.016 | 0.049 |
| Bolivia | \< 5 | Inpatient | 0.070 | 0.037 | 0.124 | 0.081 | 0.041 | 0.145 |
| Bolivia | \>= 5 | Inpatient | 0.049 | 0.025 | 0.092 | 0.058 | 0.028 | 0.108 |
| Bolivia | \< 5 | Outpatient | 0.050 | 0.026 | 0.088 | 0.058 | 0.029 | 0.103 |
| Bolivia | \>= 5 | Outpatient | 0.035 | 0.017 | 0.065 | 0.041 | 0.020 | 0.076 |
| Bosnia and Herzegovina | \< 5 | Inpatient | 0.057 | 0.033 | 0.091 | 0.068 | 0.036 | 0.113 |
| Bosnia and Herzegovina | \>= 5 | Inpatient | 0.040 | 0.022 | 0.065 | 0.047 | 0.024 | 0.082 |
| Bosnia and Herzegovina | \< 5 | Outpatient | 0.041 | 0.023 | 0.064 | 0.048 | 0.025 | 0.079 |
| Bosnia and Herzegovina | \>= 5 | Outpatient | 0.028 | 0.015 | 0.046 | 0.034 | 0.017 | 0.057 |
| Botswana | \< 5 | Inpatient | 0.051 | 0.030 | 0.077 | 0.060 | 0.034 | 0.095 |
| Botswana | \>= 5 | Inpatient | 0.035 | 0.020 | 0.056 | 0.042 | 0.022 | 0.070 |
| Botswana | \< 5 | Outpatient | 0.036 | 0.021 | 0.054 | 0.043 | 0.024 | 0.067 |
| Botswana | \>= 5 | Outpatient | 0.025 | 0.014 | 0.039 | 0.030 | 0.016 | 0.049 |
| Brazil | \< 5 | Inpatient | 0.055 | 0.028 | 0.091 | 0.065 | 0.032 | 0.109 |
| Brazil | \>= 5 | Inpatient | 0.038 | 0.019 | 0.065 | 0.045 | 0.021 | 0.080 |
| Brazil | \< 5 | Outpatient | 0.039 | 0.020 | 0.064 | 0.046 | 0.023 | 0.077 |
| Brazil | \>= 5 | Outpatient | 0.027 | 0.013 | 0.046 | 0.032 | 0.015 | 0.055 |
| Brunei Darussalam | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Brunei Darussalam | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Brunei Darussalam | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Brunei Darussalam | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bulgaria | \< 5 | Inpatient | 0.057 | 0.033 | 0.091 | 0.068 | 0.037 | 0.111 |
| Bulgaria | \>= 5 | Inpatient | 0.040 | 0.022 | 0.065 | 0.048 | 0.025 | 0.081 |
| Bulgaria | \< 5 | Outpatient | 0.041 | 0.023 | 0.064 | 0.048 | 0.026 | 0.078 |
| Bulgaria | \>= 5 | Outpatient | 0.028 | 0.015 | 0.046 | 0.034 | 0.018 | 0.056 |
| Burkina Faso | \< 5 | Inpatient | 0.048 | 0.029 | 0.071 | 0.058 | 0.031 | 0.092 |
| Burkina Faso | \>= 5 | Inpatient | 0.033 | 0.019 | 0.051 | 0.040 | 0.021 | 0.066 |
| Burkina Faso | \< 5 | Outpatient | 0.034 | 0.020 | 0.051 | 0.041 | 0.022 | 0.065 |
| Burkina Faso | \>= 5 | Outpatient | 0.023 | 0.013 | 0.036 | 0.028 | 0.014 | 0.047 |
| Burundi | \< 5 | Inpatient | 0.048 | 0.029 | 0.071 | 0.058 | 0.031 | 0.092 |
| Burundi | \>= 5 | Inpatient | 0.033 | 0.019 | 0.051 | 0.040 | 0.020 | 0.067 |
| Burundi | \< 5 | Outpatient | 0.034 | 0.020 | 0.051 | 0.041 | 0.022 | 0.066 |
| Burundi | \>= 5 | Outpatient | 0.023 | 0.013 | 0.036 | 0.029 | 0.014 | 0.047 |
| Cabo Verde | \< 5 | Inpatient | 0.049 | 0.026 | 0.079 | 0.058 | 0.029 | 0.097 |
| Cabo Verde | \>= 5 | Inpatient | 0.034 | 0.018 | 0.057 | 0.041 | 0.020 | 0.070 |
| Cabo Verde | \< 5 | Outpatient | 0.035 | 0.018 | 0.056 | 0.042 | 0.020 | 0.069 |
| Cabo Verde | \>= 5 | Outpatient | 0.024 | 0.012 | 0.040 | 0.029 | 0.013 | 0.049 |
| Cambodia | \< 5 | Inpatient | 0.050 | 0.025 | 0.082 | 0.059 | 0.028 | 0.101 |
| Cambodia | \>= 5 | Inpatient | 0.035 | 0.017 | 0.058 | 0.041 | 0.019 | 0.073 |
| Cambodia | \< 5 | Outpatient | 0.035 | 0.018 | 0.058 | 0.042 | 0.020 | 0.072 |
| Cambodia | \>= 5 | Outpatient | 0.024 | 0.012 | 0.041 | 0.029 | 0.013 | 0.051 |
| Cameroon | \< 5 | Inpatient | 0.048 | 0.025 | 0.077 | 0.057 | 0.029 | 0.094 |
| Cameroon | \>= 5 | Inpatient | 0.033 | 0.017 | 0.055 | 0.040 | 0.020 | 0.068 |
| Cameroon | \< 5 | Outpatient | 0.034 | 0.017 | 0.055 | 0.041 | 0.020 | 0.067 |
| Cameroon | \>= 5 | Outpatient | 0.023 | 0.011 | 0.039 | 0.028 | 0.013 | 0.048 |
| Canada | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Canada | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Canada | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Canada | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Central African Republic | \< 5 | Inpatient | 0.049 | 0.026 | 0.080 | 0.058 | 0.028 | 0.097 |
| Central African Republic | \>= 5 | Inpatient | 0.034 | 0.017 | 0.057 | 0.040 | 0.019 | 0.070 |
| Central African Republic | \< 5 | Outpatient | 0.034 | 0.018 | 0.056 | 0.041 | 0.020 | 0.069 |
| Central African Republic | \>= 5 | Outpatient | 0.024 | 0.012 | 0.040 | 0.028 | 0.013 | 0.050 |
| Chad | \< 5 | Inpatient | 0.048 | 0.029 | 0.071 | 0.057 | 0.031 | 0.091 |
| Chad | \>= 5 | Inpatient | 0.033 | 0.019 | 0.051 | 0.040 | 0.021 | 0.065 |
| Chad | \< 5 | Outpatient | 0.034 | 0.020 | 0.051 | 0.041 | 0.022 | 0.064 |
| Chad | \>= 5 | Outpatient | 0.023 | 0.013 | 0.036 | 0.028 | 0.014 | 0.046 |
| Chile | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Chile | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Chile | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Chile | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| China | \< 5 | Inpatient | 0.061 | 0.036 | 0.090 | 0.071 | 0.041 | 0.109 |
| China | \>= 5 | Inpatient | 0.042 | 0.025 | 0.064 | 0.050 | 0.028 | 0.079 |
| China | \< 5 | Outpatient | 0.043 | 0.026 | 0.063 | 0.051 | 0.030 | 0.076 |
| China | \>= 5 | Outpatient | 0.030 | 0.018 | 0.044 | 0.035 | 0.020 | 0.054 |
| Colombia | \< 5 | Inpatient | 0.060 | 0.032 | 0.096 | 0.070 | 0.036 | 0.116 |
| Colombia | \>= 5 | Inpatient | 0.042 | 0.022 | 0.070 | 0.049 | 0.025 | 0.085 |
| Colombia | \< 5 | Outpatient | 0.042 | 0.023 | 0.068 | 0.050 | 0.026 | 0.082 |
| Colombia | \>= 5 | Outpatient | 0.029 | 0.015 | 0.050 | 0.035 | 0.017 | 0.060 |
| Comoros | \< 5 | Inpatient | 0.050 | 0.030 | 0.073 | 0.059 | 0.033 | 0.092 |
| Comoros | \>= 5 | Inpatient | 0.035 | 0.020 | 0.053 | 0.041 | 0.022 | 0.066 |
| Comoros | \< 5 | Outpatient | 0.036 | 0.021 | 0.052 | 0.042 | 0.023 | 0.065 |
| Comoros | \>= 5 | Outpatient | 0.024 | 0.014 | 0.037 | 0.029 | 0.015 | 0.046 |
| Congo, Dem. Rep. | \< 5 | Inpatient | 0.048 | 0.029 | 0.071 | 0.057 | 0.032 | 0.089 |
| Congo, Dem. Rep. | \>= 5 | Inpatient | 0.033 | 0.019 | 0.051 | 0.040 | 0.021 | 0.064 |
| Congo, Dem. Rep. | \< 5 | Outpatient | 0.034 | 0.020 | 0.051 | 0.041 | 0.022 | 0.063 |
| Congo, Dem. Rep. | \>= 5 | Outpatient | 0.023 | 0.013 | 0.036 | 0.028 | 0.015 | 0.045 |
| Congo, Rep. | \< 5 | Inpatient | 0.050 | 0.030 | 0.073 | 0.060 | 0.034 | 0.091 |
| Congo, Rep. | \>= 5 | Inpatient | 0.035 | 0.020 | 0.053 | 0.041 | 0.023 | 0.066 |
| Congo, Rep. | \< 5 | Outpatient | 0.036 | 0.021 | 0.052 | 0.042 | 0.024 | 0.064 |
| Congo, Rep. | \>= 5 | Outpatient | 0.024 | 0.014 | 0.037 | 0.029 | 0.016 | 0.046 |
| Cook Islands | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Cook Islands | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Cook Islands | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Cook Islands | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Costa Rica | \< 5 | Inpatient | 0.058 | 0.033 | 0.091 | 0.069 | 0.038 | 0.109 |
| Costa Rica | \>= 5 | Inpatient | 0.040 | 0.022 | 0.065 | 0.048 | 0.025 | 0.080 |
| Costa Rica | \< 5 | Outpatient | 0.041 | 0.023 | 0.065 | 0.049 | 0.027 | 0.078 |
| Costa Rica | \>= 5 | Outpatient | 0.028 | 0.016 | 0.046 | 0.034 | 0.018 | 0.057 |
| Côte d’Ivoire | \< 5 | Inpatient | 0.050 | 0.030 | 0.073 | 0.060 | 0.034 | 0.091 |
| Côte d’Ivoire | \>= 5 | Inpatient | 0.035 | 0.020 | 0.053 | 0.042 | 0.023 | 0.066 |
| Côte d’Ivoire | \< 5 | Outpatient | 0.036 | 0.021 | 0.052 | 0.042 | 0.024 | 0.064 |
| Côte d’Ivoire | \>= 5 | Outpatient | 0.024 | 0.014 | 0.037 | 0.029 | 0.016 | 0.047 |
| Croatia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Croatia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Croatia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Croatia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Cuba | \< 5 | Inpatient | 0.056 | 0.027 | 0.093 | 0.066 | 0.029 | 0.120 |
| Cuba | \>= 5 | Inpatient | 0.039 | 0.018 | 0.067 | 0.046 | 0.020 | 0.086 |
| Cuba | \< 5 | Outpatient | 0.039 | 0.019 | 0.066 | 0.047 | 0.020 | 0.085 |
| Cuba | \>= 5 | Outpatient | 0.027 | 0.013 | 0.047 | 0.033 | 0.014 | 0.061 |
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
| Djibouti | \< 5 | Inpatient | 0.054 | 0.032 | 0.084 | 0.065 | 0.036 | 0.103 |
| Djibouti | \>= 5 | Inpatient | 0.038 | 0.021 | 0.060 | 0.045 | 0.024 | 0.075 |
| Djibouti | \< 5 | Outpatient | 0.039 | 0.022 | 0.059 | 0.046 | 0.025 | 0.074 |
| Djibouti | \>= 5 | Outpatient | 0.027 | 0.015 | 0.042 | 0.032 | 0.017 | 0.053 |
| Dominica | \< 5 | Inpatient | 0.058 | 0.033 | 0.091 | 0.068 | 0.035 | 0.117 |
| Dominica | \>= 5 | Inpatient | 0.040 | 0.022 | 0.065 | 0.048 | 0.024 | 0.085 |
| Dominica | \< 5 | Outpatient | 0.041 | 0.023 | 0.065 | 0.049 | 0.025 | 0.084 |
| Dominica | \>= 5 | Outpatient | 0.028 | 0.016 | 0.046 | 0.034 | 0.017 | 0.060 |
| Dominican Republic | \< 5 | Inpatient | 0.058 | 0.033 | 0.091 | 0.068 | 0.038 | 0.109 |
| Dominican Republic | \>= 5 | Inpatient | 0.040 | 0.022 | 0.065 | 0.048 | 0.026 | 0.080 |
| Dominican Republic | \< 5 | Outpatient | 0.041 | 0.023 | 0.065 | 0.049 | 0.027 | 0.077 |
| Dominican Republic | \>= 5 | Outpatient | 0.028 | 0.016 | 0.046 | 0.034 | 0.018 | 0.056 |
| Ecuador | \< 5 | Inpatient | 0.059 | 0.033 | 0.091 | 0.068 | 0.038 | 0.109 |
| Ecuador | \>= 5 | Inpatient | 0.041 | 0.022 | 0.065 | 0.048 | 0.026 | 0.080 |
| Ecuador | \< 5 | Outpatient | 0.042 | 0.023 | 0.065 | 0.049 | 0.027 | 0.077 |
| Ecuador | \>= 5 | Outpatient | 0.029 | 0.015 | 0.047 | 0.034 | 0.018 | 0.056 |
| Egypt, Arab Rep. | \< 5 | Inpatient | 0.055 | 0.030 | 0.090 | 0.065 | 0.032 | 0.111 |
| Egypt, Arab Rep. | \>= 5 | Inpatient | 0.038 | 0.020 | 0.064 | 0.045 | 0.022 | 0.081 |
| Egypt, Arab Rep. | \< 5 | Outpatient | 0.039 | 0.020 | 0.063 | 0.046 | 0.022 | 0.078 |
| Egypt, Arab Rep. | \>= 5 | Outpatient | 0.027 | 0.014 | 0.045 | 0.032 | 0.015 | 0.057 |
| El Salvador | \< 5 | Inpatient | 0.058 | 0.033 | 0.090 | 0.068 | 0.038 | 0.109 |
| El Salvador | \>= 5 | Inpatient | 0.041 | 0.022 | 0.065 | 0.048 | 0.026 | 0.080 |
| El Salvador | \< 5 | Outpatient | 0.042 | 0.023 | 0.064 | 0.049 | 0.027 | 0.077 |
| El Salvador | \>= 5 | Outpatient | 0.029 | 0.015 | 0.046 | 0.034 | 0.018 | 0.056 |
| Equatorial Guinea | \< 5 | Inpatient | 0.051 | 0.030 | 0.077 | 0.060 | 0.033 | 0.097 |
| Equatorial Guinea | \>= 5 | Inpatient | 0.035 | 0.020 | 0.056 | 0.042 | 0.022 | 0.069 |
| Equatorial Guinea | \< 5 | Outpatient | 0.036 | 0.021 | 0.054 | 0.043 | 0.023 | 0.068 |
| Equatorial Guinea | \>= 5 | Outpatient | 0.025 | 0.014 | 0.039 | 0.029 | 0.015 | 0.048 |
| Eritrea | \< 5 | Inpatient | 0.048 | 0.029 | 0.071 | 0.058 | 0.031 | 0.092 |
| Eritrea | \>= 5 | Inpatient | 0.033 | 0.019 | 0.051 | 0.040 | 0.020 | 0.067 |
| Eritrea | \< 5 | Outpatient | 0.034 | 0.020 | 0.051 | 0.041 | 0.022 | 0.066 |
| Eritrea | \>= 5 | Outpatient | 0.023 | 0.013 | 0.036 | 0.029 | 0.014 | 0.047 |
| Estonia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Estonia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Estonia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Estonia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Eswatini | \< 5 | Inpatient | 0.050 | 0.030 | 0.073 | 0.060 | 0.034 | 0.091 |
| Eswatini | \>= 5 | Inpatient | 0.035 | 0.020 | 0.053 | 0.041 | 0.023 | 0.066 |
| Eswatini | \< 5 | Outpatient | 0.036 | 0.021 | 0.052 | 0.042 | 0.024 | 0.064 |
| Eswatini | \>= 5 | Outpatient | 0.024 | 0.014 | 0.037 | 0.029 | 0.016 | 0.046 |
| Ethiopia | \< 5 | Inpatient | 0.048 | 0.029 | 0.071 | 0.057 | 0.032 | 0.090 |
| Ethiopia | \>= 5 | Inpatient | 0.033 | 0.019 | 0.051 | 0.040 | 0.021 | 0.065 |
| Ethiopia | \< 5 | Outpatient | 0.034 | 0.020 | 0.051 | 0.041 | 0.023 | 0.063 |
| Ethiopia | \>= 5 | Outpatient | 0.023 | 0.013 | 0.036 | 0.028 | 0.015 | 0.045 |
| Fiji | \< 5 | Inpatient | 0.057 | 0.033 | 0.086 | 0.068 | 0.035 | 0.110 |
| Fiji | \>= 5 | Inpatient | 0.039 | 0.022 | 0.061 | 0.048 | 0.024 | 0.080 |
| Fiji | \< 5 | Outpatient | 0.040 | 0.023 | 0.061 | 0.048 | 0.025 | 0.078 |
| Fiji | \>= 5 | Outpatient | 0.028 | 0.016 | 0.043 | 0.034 | 0.017 | 0.056 |
| Finland | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Finland | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Finland | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Finland | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| France | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| France | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| France | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| France | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Gabon | \< 5 | Inpatient | 0.053 | 0.029 | 0.085 | 0.062 | 0.033 | 0.103 |
| Gabon | \>= 5 | Inpatient | 0.036 | 0.019 | 0.061 | 0.043 | 0.022 | 0.074 |
| Gabon | \< 5 | Outpatient | 0.037 | 0.020 | 0.060 | 0.044 | 0.023 | 0.072 |
| Gabon | \>= 5 | Outpatient | 0.026 | 0.013 | 0.043 | 0.031 | 0.015 | 0.052 |
| Gambia, The | \< 5 | Inpatient | 0.048 | 0.026 | 0.075 | 0.057 | 0.029 | 0.093 |
| Gambia, The | \>= 5 | Inpatient | 0.033 | 0.018 | 0.053 | 0.040 | 0.019 | 0.068 |
| Gambia, The | \< 5 | Outpatient | 0.034 | 0.018 | 0.053 | 0.040 | 0.021 | 0.067 |
| Gambia, The | \>= 5 | Outpatient | 0.023 | 0.012 | 0.038 | 0.028 | 0.014 | 0.048 |
| Georgia | \< 5 | Inpatient | 0.057 | 0.033 | 0.091 | 0.068 | 0.037 | 0.112 |
| Georgia | \>= 5 | Inpatient | 0.040 | 0.022 | 0.065 | 0.048 | 0.024 | 0.082 |
| Georgia | \< 5 | Outpatient | 0.041 | 0.023 | 0.064 | 0.049 | 0.026 | 0.080 |
| Georgia | \>= 5 | Outpatient | 0.028 | 0.015 | 0.046 | 0.034 | 0.017 | 0.057 |
| Germany | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Germany | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Germany | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Germany | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Ghana | \< 5 | Inpatient | 0.047 | 0.026 | 0.074 | 0.057 | 0.029 | 0.093 |
| Ghana | \>= 5 | Inpatient | 0.033 | 0.017 | 0.053 | 0.040 | 0.019 | 0.068 |
| Ghana | \< 5 | Outpatient | 0.034 | 0.018 | 0.053 | 0.040 | 0.020 | 0.066 |
| Ghana | \>= 5 | Outpatient | 0.023 | 0.012 | 0.037 | 0.028 | 0.013 | 0.048 |
| Greece | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Greece | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Greece | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Greece | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Grenada | \< 5 | Inpatient | 0.058 | 0.033 | 0.091 | 0.068 | 0.035 | 0.117 |
| Grenada | \>= 5 | Inpatient | 0.040 | 0.022 | 0.065 | 0.048 | 0.024 | 0.085 |
| Grenada | \< 5 | Outpatient | 0.041 | 0.023 | 0.065 | 0.049 | 0.025 | 0.084 |
| Grenada | \>= 5 | Outpatient | 0.028 | 0.016 | 0.046 | 0.034 | 0.017 | 0.060 |
| Guatemala | \< 5 | Inpatient | 0.060 | 0.033 | 0.095 | 0.071 | 0.036 | 0.116 |
| Guatemala | \>= 5 | Inpatient | 0.042 | 0.022 | 0.068 | 0.050 | 0.025 | 0.085 |
| Guatemala | \< 5 | Outpatient | 0.042 | 0.024 | 0.067 | 0.050 | 0.027 | 0.082 |
| Guatemala | \>= 5 | Outpatient | 0.029 | 0.016 | 0.048 | 0.035 | 0.018 | 0.059 |
| Guinea | \< 5 | Inpatient | 0.050 | 0.030 | 0.073 | 0.059 | 0.033 | 0.092 |
| Guinea | \>= 5 | Inpatient | 0.035 | 0.020 | 0.053 | 0.041 | 0.022 | 0.066 |
| Guinea | \< 5 | Outpatient | 0.036 | 0.021 | 0.052 | 0.042 | 0.023 | 0.065 |
| Guinea | \>= 5 | Outpatient | 0.024 | 0.014 | 0.037 | 0.029 | 0.015 | 0.046 |
| Guinea-Bissau | \< 5 | Inpatient | 0.052 | 0.028 | 0.088 | 0.062 | 0.031 | 0.106 |
| Guinea-Bissau | \>= 5 | Inpatient | 0.036 | 0.019 | 0.063 | 0.043 | 0.021 | 0.077 |
| Guinea-Bissau | \< 5 | Outpatient | 0.037 | 0.020 | 0.061 | 0.044 | 0.023 | 0.075 |
| Guinea-Bissau | \>= 5 | Outpatient | 0.025 | 0.013 | 0.044 | 0.031 | 0.015 | 0.055 |
| Guyana | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Guyana | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Guyana | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Guyana | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Haiti | \< 5 | Inpatient | 0.064 | 0.036 | 0.111 | 0.076 | 0.040 | 0.130 |
| Haiti | \>= 5 | Inpatient | 0.045 | 0.024 | 0.080 | 0.054 | 0.027 | 0.096 |
| Haiti | \< 5 | Outpatient | 0.046 | 0.025 | 0.078 | 0.054 | 0.029 | 0.092 |
| Haiti | \>= 5 | Outpatient | 0.032 | 0.017 | 0.056 | 0.038 | 0.019 | 0.067 |
| Honduras | \< 5 | Inpatient | 0.065 | 0.036 | 0.109 | 0.076 | 0.041 | 0.128 |
| Honduras | \>= 5 | Inpatient | 0.046 | 0.024 | 0.080 | 0.054 | 0.028 | 0.095 |
| Honduras | \< 5 | Outpatient | 0.047 | 0.025 | 0.077 | 0.054 | 0.029 | 0.091 |
| Honduras | \>= 5 | Outpatient | 0.032 | 0.016 | 0.056 | 0.038 | 0.019 | 0.066 |
| Hungary | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Hungary | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Hungary | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Hungary | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Iceland | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Iceland | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Iceland | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Iceland | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| India | \< 5 | Inpatient | 0.049 | 0.029 | 0.073 | 0.059 | 0.031 | 0.095 |
| India | \>= 5 | Inpatient | 0.034 | 0.019 | 0.053 | 0.041 | 0.021 | 0.070 |
| India | \< 5 | Outpatient | 0.035 | 0.020 | 0.053 | 0.042 | 0.022 | 0.069 |
| India | \>= 5 | Outpatient | 0.024 | 0.013 | 0.037 | 0.029 | 0.014 | 0.050 |
| Indonesia | \< 5 | Inpatient | 0.056 | 0.031 | 0.091 | 0.066 | 0.034 | 0.111 |
| Indonesia | \>= 5 | Inpatient | 0.039 | 0.021 | 0.065 | 0.046 | 0.023 | 0.081 |
| Indonesia | \< 5 | Outpatient | 0.040 | 0.021 | 0.064 | 0.047 | 0.024 | 0.079 |
| Indonesia | \>= 5 | Outpatient | 0.027 | 0.014 | 0.046 | 0.033 | 0.016 | 0.057 |
| Iran, Islamic Rep. | \< 5 | Inpatient | 0.052 | 0.028 | 0.084 | 0.061 | 0.030 | 0.103 |
| Iran, Islamic Rep. | \>= 5 | Inpatient | 0.036 | 0.018 | 0.060 | 0.043 | 0.020 | 0.074 |
| Iran, Islamic Rep. | \< 5 | Outpatient | 0.037 | 0.019 | 0.060 | 0.044 | 0.021 | 0.074 |
| Iran, Islamic Rep. | \>= 5 | Outpatient | 0.025 | 0.012 | 0.042 | 0.030 | 0.014 | 0.053 |
| Iraq | \< 5 | Inpatient | 0.054 | 0.032 | 0.084 | 0.064 | 0.036 | 0.102 |
| Iraq | \>= 5 | Inpatient | 0.038 | 0.021 | 0.060 | 0.045 | 0.024 | 0.073 |
| Iraq | \< 5 | Outpatient | 0.039 | 0.022 | 0.059 | 0.046 | 0.025 | 0.073 |
| Iraq | \>= 5 | Outpatient | 0.027 | 0.015 | 0.042 | 0.032 | 0.017 | 0.052 |
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
| Jamaica | \< 5 | Inpatient | 0.058 | 0.033 | 0.091 | 0.068 | 0.035 | 0.117 |
| Jamaica | \>= 5 | Inpatient | 0.040 | 0.022 | 0.065 | 0.048 | 0.024 | 0.085 |
| Jamaica | \< 5 | Outpatient | 0.041 | 0.023 | 0.065 | 0.049 | 0.025 | 0.084 |
| Jamaica | \>= 5 | Outpatient | 0.028 | 0.016 | 0.046 | 0.034 | 0.017 | 0.060 |
| Japan | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Japan | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Japan | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Japan | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Jordan | \< 5 | Inpatient | 0.054 | 0.032 | 0.084 | 0.065 | 0.035 | 0.104 |
| Jordan | \>= 5 | Inpatient | 0.038 | 0.021 | 0.060 | 0.045 | 0.024 | 0.075 |
| Jordan | \< 5 | Outpatient | 0.039 | 0.022 | 0.059 | 0.046 | 0.025 | 0.075 |
| Jordan | \>= 5 | Outpatient | 0.027 | 0.015 | 0.042 | 0.032 | 0.017 | 0.053 |
| Kazakhstan | \< 5 | Inpatient | 0.057 | 0.033 | 0.091 | 0.068 | 0.036 | 0.113 |
| Kazakhstan | \>= 5 | Inpatient | 0.040 | 0.022 | 0.065 | 0.047 | 0.024 | 0.082 |
| Kazakhstan | \< 5 | Outpatient | 0.041 | 0.023 | 0.064 | 0.048 | 0.025 | 0.079 |
| Kazakhstan | \>= 5 | Outpatient | 0.028 | 0.015 | 0.046 | 0.034 | 0.017 | 0.057 |
| Kenya | \< 5 | Inpatient | 0.051 | 0.029 | 0.077 | 0.061 | 0.033 | 0.098 |
| Kenya | \>= 5 | Inpatient | 0.035 | 0.020 | 0.055 | 0.043 | 0.022 | 0.071 |
| Kenya | \< 5 | Outpatient | 0.036 | 0.021 | 0.055 | 0.043 | 0.023 | 0.070 |
| Kenya | \>= 5 | Outpatient | 0.025 | 0.014 | 0.039 | 0.030 | 0.015 | 0.050 |
| Kiribati | \< 5 | Inpatient | 0.052 | 0.029 | 0.080 | 0.062 | 0.032 | 0.102 |
| Kiribati | \>= 5 | Inpatient | 0.036 | 0.019 | 0.057 | 0.043 | 0.021 | 0.073 |
| Kiribati | \< 5 | Outpatient | 0.037 | 0.020 | 0.057 | 0.044 | 0.022 | 0.073 |
| Kiribati | \>= 5 | Outpatient | 0.025 | 0.013 | 0.040 | 0.031 | 0.015 | 0.052 |
| Korea, Dem. People’s Rep. | \< 5 | Inpatient | 0.052 | 0.031 | 0.075 | 0.061 | 0.034 | 0.096 |
| Korea, Dem. People’s Rep. | \>= 5 | Inpatient | 0.036 | 0.021 | 0.054 | 0.043 | 0.023 | 0.069 |
| Korea, Dem. People’s Rep. | \< 5 | Outpatient | 0.037 | 0.022 | 0.055 | 0.044 | 0.024 | 0.069 |
| Korea, Dem. People’s Rep. | \>= 5 | Outpatient | 0.025 | 0.014 | 0.039 | 0.030 | 0.016 | 0.049 |
| Korea, Rep. | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Korea, Rep. | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Korea, Rep. | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Korea, Rep. | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kuwait | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kuwait | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kuwait | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kuwait | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kyrgyz Republic | \< 5 | Inpatient | 0.056 | 0.032 | 0.088 | 0.066 | 0.037 | 0.108 |
| Kyrgyz Republic | \>= 5 | Inpatient | 0.039 | 0.022 | 0.063 | 0.047 | 0.025 | 0.078 |
| Kyrgyz Republic | \< 5 | Outpatient | 0.040 | 0.023 | 0.063 | 0.047 | 0.026 | 0.076 |
| Kyrgyz Republic | \>= 5 | Outpatient | 0.027 | 0.015 | 0.044 | 0.033 | 0.017 | 0.055 |
| Lao PDR | \< 5 | Inpatient | 0.052 | 0.029 | 0.080 | 0.061 | 0.032 | 0.100 |
| Lao PDR | \>= 5 | Inpatient | 0.036 | 0.019 | 0.057 | 0.043 | 0.022 | 0.071 |
| Lao PDR | \< 5 | Outpatient | 0.037 | 0.020 | 0.057 | 0.043 | 0.023 | 0.070 |
| Lao PDR | \>= 5 | Outpatient | 0.025 | 0.013 | 0.040 | 0.030 | 0.015 | 0.050 |
| Latvia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Latvia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Latvia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Latvia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Lebanon | \< 5 | Inpatient | 0.054 | 0.028 | 0.088 | 0.063 | 0.031 | 0.107 |
| Lebanon | \>= 5 | Inpatient | 0.037 | 0.019 | 0.063 | 0.044 | 0.021 | 0.078 |
| Lebanon | \< 5 | Outpatient | 0.038 | 0.020 | 0.062 | 0.045 | 0.022 | 0.076 |
| Lebanon | \>= 5 | Outpatient | 0.026 | 0.013 | 0.044 | 0.031 | 0.015 | 0.054 |
| Lesotho | \< 5 | Inpatient | 0.050 | 0.030 | 0.073 | 0.060 | 0.034 | 0.093 |
| Lesotho | \>= 5 | Inpatient | 0.035 | 0.020 | 0.053 | 0.042 | 0.022 | 0.068 |
| Lesotho | \< 5 | Outpatient | 0.036 | 0.021 | 0.052 | 0.043 | 0.023 | 0.066 |
| Lesotho | \>= 5 | Outpatient | 0.024 | 0.014 | 0.037 | 0.030 | 0.016 | 0.048 |
| Liberia | \< 5 | Inpatient | 0.048 | 0.029 | 0.071 | 0.057 | 0.032 | 0.089 |
| Liberia | \>= 5 | Inpatient | 0.033 | 0.019 | 0.051 | 0.040 | 0.021 | 0.064 |
| Liberia | \< 5 | Outpatient | 0.034 | 0.020 | 0.051 | 0.041 | 0.023 | 0.063 |
| Liberia | \>= 5 | Outpatient | 0.023 | 0.013 | 0.036 | 0.028 | 0.015 | 0.045 |
| Libya | \< 5 | Inpatient | 0.054 | 0.032 | 0.084 | 0.065 | 0.036 | 0.104 |
| Libya | \>= 5 | Inpatient | 0.038 | 0.021 | 0.060 | 0.045 | 0.024 | 0.075 |
| Libya | \< 5 | Outpatient | 0.039 | 0.022 | 0.059 | 0.046 | 0.025 | 0.074 |
| Libya | \>= 5 | Outpatient | 0.027 | 0.015 | 0.042 | 0.032 | 0.017 | 0.053 |
| Lithuania | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Lithuania | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Lithuania | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Lithuania | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Luxembourg | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Luxembourg | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Luxembourg | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Luxembourg | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Madagascar | \< 5 | Inpatient | 0.048 | 0.029 | 0.071 | 0.057 | 0.032 | 0.089 |
| Madagascar | \>= 5 | Inpatient | 0.033 | 0.019 | 0.051 | 0.040 | 0.021 | 0.064 |
| Madagascar | \< 5 | Outpatient | 0.034 | 0.020 | 0.051 | 0.041 | 0.023 | 0.063 |
| Madagascar | \>= 5 | Outpatient | 0.023 | 0.013 | 0.036 | 0.028 | 0.015 | 0.045 |
| Malawi | \< 5 | Inpatient | 0.048 | 0.026 | 0.075 | 0.058 | 0.029 | 0.095 |
| Malawi | \>= 5 | Inpatient | 0.033 | 0.018 | 0.054 | 0.040 | 0.019 | 0.069 |
| Malawi | \< 5 | Outpatient | 0.034 | 0.018 | 0.054 | 0.041 | 0.020 | 0.068 |
| Malawi | \>= 5 | Outpatient | 0.023 | 0.012 | 0.038 | 0.028 | 0.013 | 0.049 |
| Malaysia | \< 5 | Inpatient | 0.057 | 0.033 | 0.086 | 0.067 | 0.036 | 0.107 |
| Malaysia | \>= 5 | Inpatient | 0.039 | 0.022 | 0.061 | 0.047 | 0.025 | 0.076 |
| Malaysia | \< 5 | Outpatient | 0.040 | 0.023 | 0.061 | 0.047 | 0.026 | 0.075 |
| Malaysia | \>= 5 | Outpatient | 0.028 | 0.016 | 0.043 | 0.033 | 0.017 | 0.053 |
| Maldives | \< 5 | Inpatient | 0.052 | 0.031 | 0.076 | 0.061 | 0.033 | 0.096 |
| Maldives | \>= 5 | Inpatient | 0.036 | 0.021 | 0.054 | 0.042 | 0.022 | 0.069 |
| Maldives | \< 5 | Outpatient | 0.037 | 0.021 | 0.055 | 0.043 | 0.023 | 0.070 |
| Maldives | \>= 5 | Outpatient | 0.025 | 0.014 | 0.039 | 0.030 | 0.016 | 0.050 |
| Mali | \< 5 | Inpatient | 0.045 | 0.025 | 0.071 | 0.054 | 0.028 | 0.088 |
| Mali | \>= 5 | Inpatient | 0.031 | 0.016 | 0.050 | 0.037 | 0.018 | 0.064 |
| Mali | \< 5 | Outpatient | 0.032 | 0.017 | 0.051 | 0.038 | 0.019 | 0.063 |
| Mali | \>= 5 | Outpatient | 0.022 | 0.011 | 0.036 | 0.026 | 0.013 | 0.045 |
| Malta | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Malta | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Malta | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Malta | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Marshall Islands | \< 5 | Inpatient | 0.057 | 0.033 | 0.087 | 0.067 | 0.037 | 0.106 |
| Marshall Islands | \>= 5 | Inpatient | 0.039 | 0.022 | 0.062 | 0.047 | 0.025 | 0.077 |
| Marshall Islands | \< 5 | Outpatient | 0.040 | 0.023 | 0.061 | 0.048 | 0.027 | 0.074 |
| Marshall Islands | \>= 5 | Outpatient | 0.028 | 0.016 | 0.043 | 0.033 | 0.018 | 0.053 |
| Mauritania | \< 5 | Inpatient | 0.050 | 0.030 | 0.073 | 0.060 | 0.034 | 0.093 |
| Mauritania | \>= 5 | Inpatient | 0.035 | 0.020 | 0.053 | 0.042 | 0.022 | 0.068 |
| Mauritania | \< 5 | Outpatient | 0.036 | 0.021 | 0.052 | 0.043 | 0.024 | 0.066 |
| Mauritania | \>= 5 | Outpatient | 0.024 | 0.014 | 0.037 | 0.030 | 0.016 | 0.047 |
| Mauritius | \< 5 | Inpatient | 0.051 | 0.030 | 0.077 | 0.061 | 0.033 | 0.096 |
| Mauritius | \>= 5 | Inpatient | 0.035 | 0.020 | 0.056 | 0.042 | 0.022 | 0.070 |
| Mauritius | \< 5 | Outpatient | 0.036 | 0.021 | 0.054 | 0.043 | 0.023 | 0.068 |
| Mauritius | \>= 5 | Outpatient | 0.025 | 0.014 | 0.039 | 0.030 | 0.016 | 0.049 |
| Mexico | \< 5 | Inpatient | 0.058 | 0.033 | 0.090 | 0.068 | 0.038 | 0.109 |
| Mexico | \>= 5 | Inpatient | 0.041 | 0.022 | 0.065 | 0.048 | 0.026 | 0.080 |
| Mexico | \< 5 | Outpatient | 0.042 | 0.023 | 0.064 | 0.049 | 0.027 | 0.077 |
| Mexico | \>= 5 | Outpatient | 0.029 | 0.015 | 0.046 | 0.034 | 0.018 | 0.056 |
| Micronesia, Fed. Sts. | \< 5 | Inpatient | 0.053 | 0.027 | 0.085 | 0.061 | 0.033 | 0.098 |
| Micronesia, Fed. Sts. | \>= 5 | Inpatient | 0.036 | 0.018 | 0.061 | 0.043 | 0.022 | 0.070 |
| Micronesia, Fed. Sts. | \< 5 | Outpatient | 0.037 | 0.019 | 0.060 | 0.044 | 0.023 | 0.069 |
| Micronesia, Fed. Sts. | \>= 5 | Outpatient | 0.026 | 0.013 | 0.043 | 0.030 | 0.015 | 0.050 |
| Moldova | \< 5 | Inpatient | 0.057 | 0.033 | 0.091 | 0.068 | 0.037 | 0.111 |
| Moldova | \>= 5 | Inpatient | 0.040 | 0.022 | 0.065 | 0.048 | 0.025 | 0.081 |
| Moldova | \< 5 | Outpatient | 0.041 | 0.023 | 0.064 | 0.048 | 0.026 | 0.078 |
| Moldova | \>= 5 | Outpatient | 0.028 | 0.015 | 0.046 | 0.034 | 0.017 | 0.056 |
| Monaco | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Monaco | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Monaco | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Monaco | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Mongolia | \< 5 | Inpatient | 0.052 | 0.029 | 0.080 | 0.061 | 0.032 | 0.100 |
| Mongolia | \>= 5 | Inpatient | 0.036 | 0.019 | 0.057 | 0.043 | 0.022 | 0.071 |
| Mongolia | \< 5 | Outpatient | 0.037 | 0.020 | 0.057 | 0.043 | 0.023 | 0.070 |
| Mongolia | \>= 5 | Outpatient | 0.025 | 0.013 | 0.040 | 0.030 | 0.015 | 0.050 |
| Montenegro | \< 5 | Inpatient | 0.057 | 0.033 | 0.091 | 0.068 | 0.036 | 0.113 |
| Montenegro | \>= 5 | Inpatient | 0.040 | 0.022 | 0.065 | 0.047 | 0.024 | 0.082 |
| Montenegro | \< 5 | Outpatient | 0.041 | 0.023 | 0.064 | 0.048 | 0.025 | 0.079 |
| Montenegro | \>= 5 | Outpatient | 0.028 | 0.015 | 0.046 | 0.034 | 0.017 | 0.057 |
| Morocco | \< 5 | Inpatient | 0.054 | 0.028 | 0.087 | 0.064 | 0.031 | 0.110 |
| Morocco | \>= 5 | Inpatient | 0.037 | 0.019 | 0.062 | 0.045 | 0.021 | 0.080 |
| Morocco | \< 5 | Outpatient | 0.038 | 0.020 | 0.062 | 0.046 | 0.022 | 0.080 |
| Morocco | \>= 5 | Outpatient | 0.026 | 0.013 | 0.044 | 0.032 | 0.015 | 0.057 |
| Mozambique | \< 5 | Inpatient | 0.045 | 0.024 | 0.072 | 0.054 | 0.027 | 0.089 |
| Mozambique | \>= 5 | Inpatient | 0.031 | 0.016 | 0.051 | 0.038 | 0.018 | 0.065 |
| Mozambique | \< 5 | Outpatient | 0.032 | 0.017 | 0.051 | 0.039 | 0.019 | 0.063 |
| Mozambique | \>= 5 | Outpatient | 0.022 | 0.011 | 0.036 | 0.027 | 0.013 | 0.046 |
| Myanmar | \< 5 | Inpatient | 0.052 | 0.031 | 0.075 | 0.062 | 0.035 | 0.095 |
| Myanmar | \>= 5 | Inpatient | 0.036 | 0.021 | 0.054 | 0.043 | 0.024 | 0.069 |
| Myanmar | \< 5 | Outpatient | 0.037 | 0.022 | 0.055 | 0.044 | 0.025 | 0.068 |
| Myanmar | \>= 5 | Outpatient | 0.025 | 0.014 | 0.039 | 0.031 | 0.016 | 0.049 |
| Namibia | \< 5 | Inpatient | 0.051 | 0.030 | 0.077 | 0.060 | 0.034 | 0.094 |
| Namibia | \>= 5 | Inpatient | 0.035 | 0.020 | 0.056 | 0.042 | 0.023 | 0.069 |
| Namibia | \< 5 | Outpatient | 0.036 | 0.021 | 0.054 | 0.043 | 0.024 | 0.067 |
| Namibia | \>= 5 | Outpatient | 0.025 | 0.014 | 0.039 | 0.030 | 0.016 | 0.048 |
| Nauru | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nauru | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nauru | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nauru | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nepal | \< 5 | Inpatient | 0.049 | 0.027 | 0.076 | 0.058 | 0.029 | 0.095 |
| Nepal | \>= 5 | Inpatient | 0.034 | 0.018 | 0.054 | 0.040 | 0.019 | 0.068 |
| Nepal | \< 5 | Outpatient | 0.035 | 0.018 | 0.055 | 0.041 | 0.020 | 0.068 |
| Nepal | \>= 5 | Outpatient | 0.024 | 0.012 | 0.039 | 0.028 | 0.013 | 0.049 |
| Netherlands | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Netherlands | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Netherlands | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Netherlands | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| New Zealand | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| New Zealand | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| New Zealand | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| New Zealand | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nicaragua | \< 5 | Inpatient | 0.072 | 0.038 | 0.126 | 0.083 | 0.043 | 0.143 |
| Nicaragua | \>= 5 | Inpatient | 0.051 | 0.025 | 0.093 | 0.060 | 0.029 | 0.108 |
| Nicaragua | \< 5 | Outpatient | 0.051 | 0.026 | 0.088 | 0.060 | 0.030 | 0.102 |
| Nicaragua | \>= 5 | Outpatient | 0.036 | 0.017 | 0.065 | 0.042 | 0.021 | 0.076 |
| Niger | \< 5 | Inpatient | 0.048 | 0.029 | 0.071 | 0.058 | 0.032 | 0.091 |
| Niger | \>= 5 | Inpatient | 0.033 | 0.019 | 0.051 | 0.040 | 0.021 | 0.066 |
| Niger | \< 5 | Outpatient | 0.034 | 0.020 | 0.051 | 0.041 | 0.022 | 0.064 |
| Niger | \>= 5 | Outpatient | 0.023 | 0.013 | 0.036 | 0.028 | 0.015 | 0.046 |
| Nigeria | \< 5 | Inpatient | 0.047 | 0.025 | 0.074 | 0.056 | 0.027 | 0.091 |
| Nigeria | \>= 5 | Inpatient | 0.032 | 0.017 | 0.053 | 0.039 | 0.018 | 0.066 |
| Nigeria | \< 5 | Outpatient | 0.033 | 0.017 | 0.052 | 0.040 | 0.019 | 0.065 |
| Nigeria | \>= 5 | Outpatient | 0.023 | 0.011 | 0.038 | 0.027 | 0.013 | 0.046 |
| Niue | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Niue | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Niue | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Niue | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| North Macedonia | \< 5 | Inpatient | 0.057 | 0.033 | 0.091 | 0.068 | 0.037 | 0.111 |
| North Macedonia | \>= 5 | Inpatient | 0.040 | 0.022 | 0.065 | 0.048 | 0.025 | 0.081 |
| North Macedonia | \< 5 | Outpatient | 0.041 | 0.023 | 0.064 | 0.048 | 0.026 | 0.078 |
| North Macedonia | \>= 5 | Outpatient | 0.028 | 0.015 | 0.046 | 0.034 | 0.017 | 0.056 |
| Norway | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Norway | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Norway | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Norway | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Oman | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Oman | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Oman | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Oman | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Pakistan | \< 5 | Inpatient | 0.061 | 0.033 | 0.099 | 0.072 | 0.037 | 0.123 |
| Pakistan | \>= 5 | Inpatient | 0.042 | 0.023 | 0.072 | 0.051 | 0.025 | 0.090 |
| Pakistan | \< 5 | Outpatient | 0.043 | 0.023 | 0.071 | 0.052 | 0.026 | 0.088 |
| Pakistan | \>= 5 | Outpatient | 0.030 | 0.016 | 0.051 | 0.036 | 0.018 | 0.065 |
| Palau | \< 5 | Inpatient | 0.058 | 0.031 | 0.093 | 0.067 | 0.036 | 0.108 |
| Palau | \>= 5 | Inpatient | 0.040 | 0.021 | 0.067 | 0.047 | 0.024 | 0.079 |
| Palau | \< 5 | Outpatient | 0.041 | 0.021 | 0.066 | 0.048 | 0.026 | 0.076 |
| Palau | \>= 5 | Outpatient | 0.028 | 0.014 | 0.047 | 0.034 | 0.017 | 0.055 |
| Panama | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Panama | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Panama | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Panama | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Papua New Guinea | \< 5 | Inpatient | 0.051 | 0.026 | 0.085 | 0.061 | 0.029 | 0.105 |
| Papua New Guinea | \>= 5 | Inpatient | 0.036 | 0.017 | 0.061 | 0.042 | 0.019 | 0.076 |
| Papua New Guinea | \< 5 | Outpatient | 0.037 | 0.018 | 0.060 | 0.043 | 0.020 | 0.074 |
| Papua New Guinea | \>= 5 | Outpatient | 0.025 | 0.012 | 0.043 | 0.030 | 0.013 | 0.053 |
| Paraguay | \< 5 | Inpatient | 0.058 | 0.034 | 0.088 | 0.068 | 0.038 | 0.109 |
| Paraguay | \>= 5 | Inpatient | 0.040 | 0.023 | 0.063 | 0.048 | 0.026 | 0.080 |
| Paraguay | \< 5 | Outpatient | 0.041 | 0.024 | 0.062 | 0.049 | 0.027 | 0.077 |
| Paraguay | \>= 5 | Outpatient | 0.028 | 0.016 | 0.045 | 0.034 | 0.018 | 0.056 |
| Peru | \< 5 | Inpatient | 0.061 | 0.033 | 0.099 | 0.072 | 0.038 | 0.119 |
| Peru | \>= 5 | Inpatient | 0.043 | 0.022 | 0.072 | 0.051 | 0.025 | 0.087 |
| Peru | \< 5 | Outpatient | 0.043 | 0.023 | 0.071 | 0.051 | 0.027 | 0.085 |
| Peru | \>= 5 | Outpatient | 0.030 | 0.016 | 0.051 | 0.036 | 0.018 | 0.062 |
| Philippines | \< 5 | Inpatient | 0.052 | 0.029 | 0.080 | 0.061 | 0.032 | 0.100 |
| Philippines | \>= 5 | Inpatient | 0.036 | 0.019 | 0.057 | 0.043 | 0.022 | 0.071 |
| Philippines | \< 5 | Outpatient | 0.037 | 0.020 | 0.057 | 0.043 | 0.023 | 0.070 |
| Philippines | \>= 5 | Outpatient | 0.025 | 0.013 | 0.040 | 0.030 | 0.015 | 0.050 |
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
| Russian Federation | \< 5 | Inpatient | 0.057 | 0.030 | 0.092 | 0.067 | 0.033 | 0.114 |
| Russian Federation | \>= 5 | Inpatient | 0.039 | 0.020 | 0.066 | 0.047 | 0.022 | 0.083 |
| Russian Federation | \< 5 | Outpatient | 0.040 | 0.020 | 0.066 | 0.047 | 0.023 | 0.082 |
| Russian Federation | \>= 5 | Outpatient | 0.028 | 0.014 | 0.047 | 0.033 | 0.015 | 0.058 |
| Rwanda | \< 5 | Inpatient | 0.048 | 0.029 | 0.071 | 0.058 | 0.031 | 0.092 |
| Rwanda | \>= 5 | Inpatient | 0.033 | 0.019 | 0.051 | 0.040 | 0.020 | 0.066 |
| Rwanda | \< 5 | Outpatient | 0.034 | 0.020 | 0.051 | 0.041 | 0.022 | 0.065 |
| Rwanda | \>= 5 | Outpatient | 0.023 | 0.013 | 0.036 | 0.028 | 0.014 | 0.047 |
| Samoa | \< 5 | Inpatient | 0.052 | 0.029 | 0.080 | 0.061 | 0.032 | 0.100 |
| Samoa | \>= 5 | Inpatient | 0.036 | 0.019 | 0.057 | 0.043 | 0.022 | 0.071 |
| Samoa | \< 5 | Outpatient | 0.037 | 0.020 | 0.057 | 0.043 | 0.023 | 0.070 |
| Samoa | \>= 5 | Outpatient | 0.025 | 0.013 | 0.040 | 0.030 | 0.015 | 0.050 |
| San Marino | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| San Marino | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| San Marino | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| San Marino | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| São Tomé and Principe | \< 5 | Inpatient | 0.050 | 0.030 | 0.073 | 0.060 | 0.034 | 0.093 |
| São Tomé and Principe | \>= 5 | Inpatient | 0.035 | 0.020 | 0.053 | 0.042 | 0.022 | 0.068 |
| São Tomé and Principe | \< 5 | Outpatient | 0.036 | 0.021 | 0.052 | 0.043 | 0.023 | 0.066 |
| São Tomé and Principe | \>= 5 | Outpatient | 0.024 | 0.014 | 0.037 | 0.030 | 0.016 | 0.048 |
| Saudi Arabia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Saudi Arabia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Saudi Arabia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Saudi Arabia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Senegal | \< 5 | Inpatient | 0.050 | 0.030 | 0.073 | 0.060 | 0.033 | 0.094 |
| Senegal | \>= 5 | Inpatient | 0.035 | 0.020 | 0.053 | 0.042 | 0.022 | 0.069 |
| Senegal | \< 5 | Outpatient | 0.036 | 0.021 | 0.052 | 0.043 | 0.023 | 0.067 |
| Senegal | \>= 5 | Outpatient | 0.024 | 0.014 | 0.037 | 0.030 | 0.015 | 0.048 |
| Serbia | \< 5 | Inpatient | 0.057 | 0.033 | 0.091 | 0.068 | 0.036 | 0.113 |
| Serbia | \>= 5 | Inpatient | 0.040 | 0.022 | 0.065 | 0.047 | 0.024 | 0.082 |
| Serbia | \< 5 | Outpatient | 0.041 | 0.023 | 0.064 | 0.048 | 0.025 | 0.079 |
| Serbia | \>= 5 | Outpatient | 0.028 | 0.015 | 0.046 | 0.034 | 0.017 | 0.057 |
| Seychelles | \< 5 | Inpatient | 0.051 | 0.030 | 0.077 | 0.061 | 0.033 | 0.097 |
| Seychelles | \>= 5 | Inpatient | 0.035 | 0.020 | 0.056 | 0.042 | 0.022 | 0.070 |
| Seychelles | \< 5 | Outpatient | 0.036 | 0.021 | 0.054 | 0.043 | 0.023 | 0.069 |
| Seychelles | \>= 5 | Outpatient | 0.025 | 0.014 | 0.039 | 0.030 | 0.015 | 0.050 |
| Sierra Leone | \< 5 | Inpatient | 0.048 | 0.029 | 0.071 | 0.058 | 0.031 | 0.092 |
| Sierra Leone | \>= 5 | Inpatient | 0.033 | 0.019 | 0.051 | 0.040 | 0.020 | 0.066 |
| Sierra Leone | \< 5 | Outpatient | 0.034 | 0.020 | 0.051 | 0.041 | 0.022 | 0.065 |
| Sierra Leone | \>= 5 | Outpatient | 0.023 | 0.013 | 0.036 | 0.028 | 0.014 | 0.047 |
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
| Solomon Islands | \< 5 | Inpatient | 0.052 | 0.029 | 0.080 | 0.061 | 0.033 | 0.098 |
| Solomon Islands | \>= 5 | Inpatient | 0.036 | 0.019 | 0.057 | 0.043 | 0.022 | 0.070 |
| Solomon Islands | \< 5 | Outpatient | 0.037 | 0.020 | 0.057 | 0.044 | 0.023 | 0.069 |
| Solomon Islands | \>= 5 | Outpatient | 0.025 | 0.013 | 0.040 | 0.030 | 0.015 | 0.050 |
| Somalia | \< 5 | Inpatient | 0.054 | 0.032 | 0.083 | 0.064 | 0.035 | 0.103 |
| Somalia | \>= 5 | Inpatient | 0.038 | 0.022 | 0.059 | 0.045 | 0.024 | 0.074 |
| Somalia | \< 5 | Outpatient | 0.039 | 0.022 | 0.059 | 0.046 | 0.025 | 0.074 |
| Somalia | \>= 5 | Outpatient | 0.027 | 0.015 | 0.041 | 0.032 | 0.017 | 0.052 |
| South Africa | \< 5 | Inpatient | 0.050 | 0.028 | 0.078 | 0.059 | 0.032 | 0.093 |
| South Africa | \>= 5 | Inpatient | 0.034 | 0.019 | 0.056 | 0.041 | 0.022 | 0.068 |
| South Africa | \< 5 | Outpatient | 0.035 | 0.020 | 0.055 | 0.042 | 0.023 | 0.067 |
| South Africa | \>= 5 | Outpatient | 0.024 | 0.013 | 0.040 | 0.029 | 0.015 | 0.048 |
| South Sudan | \< 5 | Inpatient | 0.048 | 0.029 | 0.071 | 0.057 | 0.031 | 0.091 |
| South Sudan | \>= 5 | Inpatient | 0.033 | 0.019 | 0.051 | 0.040 | 0.021 | 0.065 |
| South Sudan | \< 5 | Outpatient | 0.034 | 0.020 | 0.051 | 0.041 | 0.022 | 0.064 |
| South Sudan | \>= 5 | Outpatient | 0.023 | 0.013 | 0.036 | 0.028 | 0.014 | 0.046 |
| Spain | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Spain | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Spain | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Spain | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Sri Lanka | \< 5 | Inpatient | 0.052 | 0.031 | 0.075 | 0.061 | 0.034 | 0.096 |
| Sri Lanka | \>= 5 | Inpatient | 0.036 | 0.021 | 0.054 | 0.043 | 0.023 | 0.069 |
| Sri Lanka | \< 5 | Outpatient | 0.037 | 0.022 | 0.055 | 0.044 | 0.024 | 0.069 |
| Sri Lanka | \>= 5 | Outpatient | 0.025 | 0.014 | 0.039 | 0.030 | 0.016 | 0.049 |
| St. Kitts and Nevis | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| St. Kitts and Nevis | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| St. Kitts and Nevis | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| St. Kitts and Nevis | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| St. Lucia | \< 5 | Inpatient | 0.058 | 0.033 | 0.091 | 0.068 | 0.035 | 0.117 |
| St. Lucia | \>= 5 | Inpatient | 0.040 | 0.022 | 0.065 | 0.048 | 0.024 | 0.085 |
| St. Lucia | \< 5 | Outpatient | 0.041 | 0.023 | 0.065 | 0.049 | 0.025 | 0.084 |
| St. Lucia | \>= 5 | Outpatient | 0.028 | 0.016 | 0.046 | 0.034 | 0.017 | 0.060 |
| St. Vincent and the Grenadines | \< 5 | Inpatient | 0.058 | 0.033 | 0.091 | 0.068 | 0.035 | 0.117 |
| St. Vincent and the Grenadines | \>= 5 | Inpatient | 0.040 | 0.022 | 0.065 | 0.048 | 0.024 | 0.085 |
| St. Vincent and the Grenadines | \< 5 | Outpatient | 0.041 | 0.023 | 0.065 | 0.049 | 0.025 | 0.084 |
| St. Vincent and the Grenadines | \>= 5 | Outpatient | 0.028 | 0.016 | 0.046 | 0.034 | 0.017 | 0.060 |
| Sudan | \< 5 | Inpatient | 0.054 | 0.032 | 0.083 | 0.065 | 0.035 | 0.104 |
| Sudan | \>= 5 | Inpatient | 0.038 | 0.022 | 0.059 | 0.046 | 0.024 | 0.076 |
| Sudan | \< 5 | Outpatient | 0.039 | 0.022 | 0.059 | 0.046 | 0.025 | 0.076 |
| Sudan | \>= 5 | Outpatient | 0.027 | 0.015 | 0.041 | 0.032 | 0.016 | 0.054 |
| Suriname | \< 5 | Inpatient | 0.058 | 0.033 | 0.091 | 0.068 | 0.035 | 0.117 |
| Suriname | \>= 5 | Inpatient | 0.040 | 0.022 | 0.065 | 0.048 | 0.024 | 0.085 |
| Suriname | \< 5 | Outpatient | 0.041 | 0.023 | 0.065 | 0.049 | 0.025 | 0.084 |
| Suriname | \>= 5 | Outpatient | 0.028 | 0.016 | 0.046 | 0.034 | 0.017 | 0.060 |
| Sweden | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Sweden | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Sweden | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Sweden | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Switzerland | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Switzerland | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Switzerland | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Switzerland | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Syrian Arab Republic | \< 5 | Inpatient | 0.054 | 0.032 | 0.083 | 0.064 | 0.035 | 0.103 |
| Syrian Arab Republic | \>= 5 | Inpatient | 0.038 | 0.022 | 0.059 | 0.045 | 0.024 | 0.074 |
| Syrian Arab Republic | \< 5 | Outpatient | 0.039 | 0.022 | 0.059 | 0.046 | 0.025 | 0.074 |
| Syrian Arab Republic | \>= 5 | Outpatient | 0.027 | 0.015 | 0.041 | 0.032 | 0.017 | 0.052 |
| Tajikistan | \< 5 | Inpatient | 0.056 | 0.032 | 0.088 | 0.067 | 0.036 | 0.110 |
| Tajikistan | \>= 5 | Inpatient | 0.039 | 0.022 | 0.063 | 0.047 | 0.024 | 0.080 |
| Tajikistan | \< 5 | Outpatient | 0.040 | 0.023 | 0.063 | 0.048 | 0.025 | 0.079 |
| Tajikistan | \>= 5 | Outpatient | 0.027 | 0.015 | 0.044 | 0.033 | 0.017 | 0.057 |
| Tanzania | \< 5 | Inpatient | 0.055 | 0.030 | 0.089 | 0.065 | 0.034 | 0.109 |
| Tanzania | \>= 5 | Inpatient | 0.038 | 0.020 | 0.064 | 0.046 | 0.023 | 0.080 |
| Tanzania | \< 5 | Outpatient | 0.039 | 0.021 | 0.063 | 0.046 | 0.024 | 0.077 |
| Tanzania | \>= 5 | Outpatient | 0.027 | 0.014 | 0.045 | 0.032 | 0.016 | 0.056 |
| Thailand | \< 5 | Inpatient | 0.047 | 0.027 | 0.072 | 0.056 | 0.029 | 0.091 |
| Thailand | \>= 5 | Inpatient | 0.033 | 0.018 | 0.051 | 0.039 | 0.020 | 0.065 |
| Thailand | \< 5 | Outpatient | 0.034 | 0.018 | 0.053 | 0.040 | 0.020 | 0.067 |
| Thailand | \>= 5 | Outpatient | 0.023 | 0.012 | 0.037 | 0.028 | 0.014 | 0.047 |
| Timor-Leste | \< 5 | Inpatient | 0.052 | 0.031 | 0.075 | 0.061 | 0.034 | 0.096 |
| Timor-Leste | \>= 5 | Inpatient | 0.036 | 0.021 | 0.054 | 0.043 | 0.023 | 0.069 |
| Timor-Leste | \< 5 | Outpatient | 0.037 | 0.022 | 0.055 | 0.044 | 0.024 | 0.069 |
| Timor-Leste | \>= 5 | Outpatient | 0.025 | 0.014 | 0.039 | 0.030 | 0.016 | 0.049 |
| Togo | \< 5 | Inpatient | 0.048 | 0.029 | 0.071 | 0.058 | 0.032 | 0.090 |
| Togo | \>= 5 | Inpatient | 0.033 | 0.019 | 0.051 | 0.040 | 0.021 | 0.065 |
| Togo | \< 5 | Outpatient | 0.034 | 0.020 | 0.051 | 0.041 | 0.022 | 0.064 |
| Togo | \>= 5 | Outpatient | 0.023 | 0.013 | 0.036 | 0.028 | 0.015 | 0.046 |
| Tonga | \< 5 | Inpatient | 0.057 | 0.033 | 0.086 | 0.067 | 0.036 | 0.107 |
| Tonga | \>= 5 | Inpatient | 0.039 | 0.022 | 0.061 | 0.047 | 0.025 | 0.076 |
| Tonga | \< 5 | Outpatient | 0.040 | 0.023 | 0.061 | 0.047 | 0.026 | 0.075 |
| Tonga | \>= 5 | Outpatient | 0.028 | 0.016 | 0.043 | 0.033 | 0.017 | 0.053 |
| Trinidad and Tobago | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Trinidad and Tobago | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Trinidad and Tobago | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Trinidad and Tobago | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Tunisia | \< 5 | Inpatient | 0.054 | 0.032 | 0.084 | 0.064 | 0.034 | 0.104 |
| Tunisia | \>= 5 | Inpatient | 0.038 | 0.021 | 0.060 | 0.045 | 0.024 | 0.075 |
| Tunisia | \< 5 | Outpatient | 0.039 | 0.022 | 0.059 | 0.046 | 0.025 | 0.074 |
| Tunisia | \>= 5 | Outpatient | 0.027 | 0.015 | 0.042 | 0.032 | 0.017 | 0.053 |
| Türkiye | \< 5 | Inpatient | 0.062 | 0.034 | 0.101 | 0.072 | 0.038 | 0.123 |
| Türkiye | \>= 5 | Inpatient | 0.043 | 0.023 | 0.073 | 0.051 | 0.026 | 0.089 |
| Türkiye | \< 5 | Outpatient | 0.044 | 0.024 | 0.071 | 0.051 | 0.027 | 0.086 |
| Türkiye | \>= 5 | Outpatient | 0.030 | 0.016 | 0.051 | 0.036 | 0.018 | 0.062 |
| Turkmenistan | \< 5 | Inpatient | 0.057 | 0.033 | 0.091 | 0.069 | 0.036 | 0.115 |
| Turkmenistan | \>= 5 | Inpatient | 0.040 | 0.022 | 0.065 | 0.048 | 0.024 | 0.084 |
| Turkmenistan | \< 5 | Outpatient | 0.041 | 0.023 | 0.064 | 0.049 | 0.025 | 0.082 |
| Turkmenistan | \>= 5 | Outpatient | 0.028 | 0.015 | 0.046 | 0.034 | 0.017 | 0.059 |
| Tuvalu | \< 5 | Inpatient | 0.057 | 0.033 | 0.086 | 0.067 | 0.036 | 0.107 |
| Tuvalu | \>= 5 | Inpatient | 0.039 | 0.022 | 0.061 | 0.047 | 0.025 | 0.076 |
| Tuvalu | \< 5 | Outpatient | 0.040 | 0.023 | 0.061 | 0.047 | 0.026 | 0.075 |
| Tuvalu | \>= 5 | Outpatient | 0.028 | 0.016 | 0.043 | 0.033 | 0.017 | 0.053 |
| Uganda | \< 5 | Inpatient | 0.048 | 0.029 | 0.071 | 0.058 | 0.031 | 0.092 |
| Uganda | \>= 5 | Inpatient | 0.033 | 0.019 | 0.051 | 0.040 | 0.021 | 0.066 |
| Uganda | \< 5 | Outpatient | 0.034 | 0.020 | 0.051 | 0.041 | 0.022 | 0.065 |
| Uganda | \>= 5 | Outpatient | 0.023 | 0.013 | 0.036 | 0.028 | 0.014 | 0.047 |
| Ukraine | \< 5 | Inpatient | 0.056 | 0.032 | 0.088 | 0.066 | 0.036 | 0.110 |
| Ukraine | \>= 5 | Inpatient | 0.039 | 0.022 | 0.063 | 0.046 | 0.024 | 0.079 |
| Ukraine | \< 5 | Outpatient | 0.040 | 0.023 | 0.063 | 0.047 | 0.025 | 0.078 |
| Ukraine | \>= 5 | Outpatient | 0.027 | 0.015 | 0.044 | 0.033 | 0.017 | 0.055 |
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
| Uzbekistan | \< 5 | Inpatient | 0.056 | 0.032 | 0.088 | 0.067 | 0.036 | 0.110 |
| Uzbekistan | \>= 5 | Inpatient | 0.039 | 0.022 | 0.063 | 0.047 | 0.024 | 0.080 |
| Uzbekistan | \< 5 | Outpatient | 0.040 | 0.023 | 0.063 | 0.048 | 0.025 | 0.079 |
| Uzbekistan | \>= 5 | Outpatient | 0.027 | 0.015 | 0.044 | 0.033 | 0.017 | 0.057 |
| Vanuatu | \< 5 | Inpatient | 0.052 | 0.029 | 0.080 | 0.061 | 0.032 | 0.100 |
| Vanuatu | \>= 5 | Inpatient | 0.036 | 0.019 | 0.057 | 0.043 | 0.022 | 0.071 |
| Vanuatu | \< 5 | Outpatient | 0.037 | 0.020 | 0.057 | 0.043 | 0.023 | 0.070 |
| Vanuatu | \>= 5 | Outpatient | 0.025 | 0.013 | 0.040 | 0.030 | 0.015 | 0.050 |
| Venezuela | \< 5 | Inpatient | 0.065 | 0.033 | 0.111 | 0.076 | 0.035 | 0.141 |
| Venezuela | \>= 5 | Inpatient | 0.045 | 0.022 | 0.082 | 0.054 | 0.024 | 0.105 |
| Venezuela | \< 5 | Outpatient | 0.046 | 0.023 | 0.079 | 0.054 | 0.024 | 0.101 |
| Venezuela | \>= 5 | Outpatient | 0.032 | 0.015 | 0.058 | 0.038 | 0.016 | 0.074 |
| Vietnam | \< 5 | Inpatient | 0.051 | 0.027 | 0.083 | 0.061 | 0.030 | 0.103 |
| Vietnam | \>= 5 | Inpatient | 0.035 | 0.018 | 0.059 | 0.042 | 0.020 | 0.073 |
| Vietnam | \< 5 | Outpatient | 0.036 | 0.019 | 0.059 | 0.043 | 0.021 | 0.073 |
| Vietnam | \>= 5 | Outpatient | 0.025 | 0.012 | 0.042 | 0.030 | 0.014 | 0.052 |
| Yemen, Rep. | \< 5 | Inpatient | 0.054 | 0.032 | 0.083 | 0.064 | 0.036 | 0.101 |
| Yemen, Rep. | \>= 5 | Inpatient | 0.038 | 0.022 | 0.059 | 0.045 | 0.024 | 0.073 |
| Yemen, Rep. | \< 5 | Outpatient | 0.039 | 0.022 | 0.059 | 0.046 | 0.026 | 0.072 |
| Yemen, Rep. | \>= 5 | Outpatient | 0.027 | 0.015 | 0.041 | 0.032 | 0.017 | 0.052 |
| Zambia | \< 5 | Inpatient | 0.053 | 0.029 | 0.086 | 0.063 | 0.033 | 0.106 |
| Zambia | \>= 5 | Inpatient | 0.037 | 0.019 | 0.062 | 0.044 | 0.022 | 0.078 |
| Zambia | \< 5 | Outpatient | 0.038 | 0.020 | 0.061 | 0.045 | 0.023 | 0.076 |
| Zambia | \>= 5 | Outpatient | 0.026 | 0.013 | 0.044 | 0.031 | 0.015 | 0.055 |
| Zimbabwe | \< 5 | Inpatient | 0.053 | 0.029 | 0.088 | 0.063 | 0.032 | 0.108 |
| Zimbabwe | \>= 5 | Inpatient | 0.037 | 0.019 | 0.064 | 0.044 | 0.022 | 0.080 |
| Zimbabwe | \< 5 | Outpatient | 0.038 | 0.020 | 0.062 | 0.045 | 0.023 | 0.077 |
| Zimbabwe | \>= 5 | Outpatient | 0.026 | 0.013 | 0.044 | 0.031 | 0.015 | 0.056 |

Estimated `r params$Pathogen` proportion by country, 2010 vs 2020

# Session info

``` r
sessioninfo::session_info()
```

    ## Warning in system2("quarto", "-V", stdout = TRUE, env = paste0("TMPDIR=", : running command '"quarto"
    ## TMPDIR=C:/Users/fbbu6966/AppData/Local/Temp/RtmpMtu9zo/file40b815737a4a -V' had status 1

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
    ##  quarto   ERROR: Unknown command "TMPDIR=C:/Users/fbbu6966/AppData/Local/Temp/RtmpMtu9zo/file40b815737a4a". Did you mean command "update"? @ C:\\Users\\fbbu6966\\AppData\\Local\\Programs\\RStudio\\RESOUR~1\\app\\bin\\quarto\\bin\\quarto.exe
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
