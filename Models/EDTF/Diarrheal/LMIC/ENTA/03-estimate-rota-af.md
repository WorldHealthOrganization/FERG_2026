Global proportion of Entamoeba (Diarrheal) • ASYMPT+ROTA • Estimate
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

    ## Warning: There were 9 divergent transitions after warmup. Increasing adapt_delta above 0.95 may help. See
    ## http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup

    ##  Family: gaussian 
    ##   Links: mu = identity; sigma = identity 
    ## Formula: yi | se(sei) ~ 1 + YEAR + SYNDROMTYPE + AGE + REFERENCE + COVERAGE + STRAIN + (1 | REG2) + (1 | REG2:SUB2) + (1 | REG2:SUB2:COUNTRY) + (1 | REG2:SUB2:COUNTRY:ID) + (1 | REG2:SUB2:COUNTRY:ID:DTP_ID) 
    ##    Data: subset(es, as.integer(FLAG) == 1) (Number of observations: 351) 
    ##   Draws: 5 chains, each with iter = 5000; warmup = 3000; thin = 1;
    ##          total post-warmup draws = 10000
    ## 
    ## Multilevel Hyperparameters:
    ## ~REG2 (Number of levels: 6) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.37      0.31     0.01     1.15 1.00     2885     5301
    ## 
    ## ~REG2:SUB2 (Number of levels: 12) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.27      0.23     0.01     0.88 1.00     2714     4200
    ## 
    ## ~REG2:SUB2:COUNTRY (Number of levels: 39) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.23      0.18     0.01     0.66 1.00     1807     2741
    ## 
    ## ~REG2:SUB2:COUNTRY:ID (Number of levels: 83) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     1.35      0.13     1.11     1.62 1.00     3331     5538
    ## 
    ## ~REG2:SUB2:COUNTRY:ID:DTP_ID (Number of levels: 351) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.48      0.05     0.39     0.59 1.00     3891     6275
    ## 
    ## Regression Coefficients:
    ##                          Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept                   26.27     65.18  -100.70   157.08 1.00     2173     3703
    ## YEAR                        -0.02      0.03    -0.08     0.05 1.00     2176     3705
    ## SYNDROMTYPEInpatient         1.14      0.20     0.76     1.54 1.00     6082     7384
    ## SYNDROMTYPEOutpatient        0.38      0.13     0.12     0.64 1.00     5361     6590
    ## AGEAgebelow5                -0.64      0.25    -1.13    -0.13 1.00     3709     5428
    ## AGEMixedages                 0.13      0.31    -0.48     0.75 1.00     3783     5469
    ## REFERENCEOther               0.93      0.38     0.18     1.69 1.00     1780     3879
    ## COVERAGE                     0.01      0.01    -0.01     0.02 1.00     2488     4702
    ## STRAINhistolytica_dispar     0.70      0.55    -0.40     1.78 1.00     2915     4706
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

    ## [1] 5148

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
    ##  $ VAL_MEAN     : num  0.0261 0.014 0.0255 0.0136 0.0249 ...
    ##  $ VAL_LWR      : num  0.00507 0.00287 0.00514 0.00292 0.00523 ...
    ##  $ VAL_UPR      : num  0.083 0.0443 0.0794 0.0426 0.0755 ...
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
    ##  $ VAL_MEAN     : num  0.0261 0.014 0.0255 0.0136 0.0249 ...
    ##  $ VAL_LWR      : num  0.00507 0.00287 0.00514 0.00292 0.00523 ...
    ##  $ VAL_UPR      : num  0.083 0.0443 0.0794 0.0426 0.0755 ...
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

\[1\] 0.000 0.005 0.010 0.015 0.020 0.025 0.030 0.035

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

\[1\] 0.000 0.001 0.002 0.003 0.004 0.005 0.006 0.007

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

\[1\] 0.000 0.005 0.010 0.015 0.020

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

\[1\] 0.000 0.001 0.002 0.003 0.004

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

\[1\] 0.00 0.01 0.02 0.03 0.04

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

\[1\] 0.000 0.002 0.004 0.006 0.008 0.010

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

\[1\] 0.000 0.005 0.010 0.015 0.020 0.025

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

\[1\] 0.000 0.001 0.002 0.003 0.004 0.005

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
| Afghanistan | \< 5 | Inpatient | 0.011 | 0.003 | 0.032 | 0.014 | 0.004 | 0.040 |
| Afghanistan | \>= 5 | Inpatient | 0.021 | 0.005 | 0.060 | 0.026 | 0.006 | 0.075 |
| Afghanistan | \< 5 | Outpatient | 0.002 | 0.000 | 0.007 | 0.003 | 0.001 | 0.009 |
| Afghanistan | \>= 5 | Outpatient | 0.005 | 0.001 | 0.014 | 0.006 | 0.001 | 0.017 |
| Albania | \< 5 | Inpatient | 0.012 | 0.003 | 0.037 | 0.020 | 0.004 | 0.069 |
| Albania | \>= 5 | Inpatient | 0.023 | 0.005 | 0.072 | 0.038 | 0.006 | 0.127 |
| Albania | \< 5 | Outpatient | 0.003 | 0.000 | 0.008 | 0.005 | 0.001 | 0.017 |
| Albania | \>= 5 | Outpatient | 0.005 | 0.001 | 0.016 | 0.009 | 0.001 | 0.032 |
| Algeria | \< 5 | Inpatient | 0.008 | 0.003 | 0.017 | 0.007 | 0.002 | 0.016 |
| Algeria | \>= 5 | Inpatient | 0.015 | 0.005 | 0.034 | 0.014 | 0.004 | 0.033 |
| Algeria | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.002 | 0.000 | 0.004 |
| Algeria | \>= 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.003 | 0.001 | 0.008 |
| Andorra | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Andorra | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Andorra | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Andorra | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Angola | \< 5 | Inpatient | 0.008 | 0.002 | 0.019 | 0.009 | 0.003 | 0.021 |
| Angola | \>= 5 | Inpatient | 0.015 | 0.004 | 0.038 | 0.017 | 0.005 | 0.043 |
| Angola | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.002 | 0.000 | 0.005 |
| Angola | \>= 5 | Outpatient | 0.003 | 0.001 | 0.009 | 0.004 | 0.001 | 0.010 |
| Antigua and Barbuda | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Antigua and Barbuda | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Antigua and Barbuda | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Antigua and Barbuda | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Argentina | \< 5 | Inpatient | 0.009 | 0.002 | 0.022 | 0.011 | 0.003 | 0.029 |
| Argentina | \>= 5 | Inpatient | 0.016 | 0.004 | 0.043 | 0.022 | 0.006 | 0.057 |
| Argentina | \< 5 | Outpatient | 0.002 | 0.000 | 0.005 | 0.002 | 0.000 | 0.007 |
| Argentina | \>= 5 | Outpatient | 0.004 | 0.001 | 0.010 | 0.005 | 0.001 | 0.013 |
| Armenia | \< 5 | Inpatient | 0.013 | 0.003 | 0.037 | 0.020 | 0.004 | 0.064 |
| Armenia | \>= 5 | Inpatient | 0.023 | 0.006 | 0.071 | 0.037 | 0.007 | 0.119 |
| Armenia | \< 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.004 | 0.001 | 0.016 |
| Armenia | \>= 5 | Outpatient | 0.005 | 0.001 | 0.016 | 0.008 | 0.001 | 0.030 |
| Australia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Australia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Australia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Australia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Austria | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Austria | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Austria | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Austria | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Azerbaijan | \< 5 | Inpatient | 0.013 | 0.003 | 0.037 | 0.011 | 0.003 | 0.034 |
| Azerbaijan | \>= 5 | Inpatient | 0.023 | 0.006 | 0.071 | 0.021 | 0.005 | 0.065 |
| Azerbaijan | \< 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.002 | 0.000 | 0.008 |
| Azerbaijan | \>= 5 | Outpatient | 0.005 | 0.001 | 0.016 | 0.005 | 0.001 | 0.015 |
| Bahamas, The | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahamas, The | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahamas, The | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahamas, The | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahrain | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahrain | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahrain | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahrain | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bangladesh | \< 5 | Inpatient | 0.011 | 0.004 | 0.024 | 0.010 | 0.003 | 0.025 |
| Bangladesh | \>= 5 | Inpatient | 0.021 | 0.007 | 0.047 | 0.019 | 0.005 | 0.049 |
| Bangladesh | \< 5 | Outpatient | 0.002 | 0.001 | 0.005 | 0.002 | 0.000 | 0.006 |
| Bangladesh | \>= 5 | Outpatient | 0.004 | 0.001 | 0.011 | 0.004 | 0.001 | 0.011 |
| Barbados | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Barbados | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Barbados | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Barbados | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belarus | \< 5 | Inpatient | 0.013 | 0.003 | 0.037 | 0.011 | 0.003 | 0.034 |
| Belarus | \>= 5 | Inpatient | 0.023 | 0.006 | 0.071 | 0.021 | 0.005 | 0.065 |
| Belarus | \< 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.002 | 0.000 | 0.008 |
| Belarus | \>= 5 | Outpatient | 0.005 | 0.001 | 0.016 | 0.005 | 0.001 | 0.015 |
| Belgium | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belgium | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belgium | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belgium | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belize | \< 5 | Inpatient | 0.009 | 0.002 | 0.022 | 0.008 | 0.002 | 0.022 |
| Belize | \>= 5 | Inpatient | 0.016 | 0.004 | 0.043 | 0.015 | 0.003 | 0.043 |
| Belize | \< 5 | Outpatient | 0.002 | 0.000 | 0.005 | 0.002 | 0.000 | 0.005 |
| Belize | \>= 5 | Outpatient | 0.004 | 0.001 | 0.010 | 0.003 | 0.000 | 0.010 |
| Benin | \< 5 | Inpatient | 0.008 | 0.003 | 0.017 | 0.010 | 0.004 | 0.023 |
| Benin | \>= 5 | Inpatient | 0.015 | 0.005 | 0.034 | 0.020 | 0.006 | 0.046 |
| Benin | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.002 | 0.001 | 0.005 |
| Benin | \>= 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.004 | 0.001 | 0.011 |
| Bhutan | \< 5 | Inpatient | 0.011 | 0.004 | 0.023 | 0.010 | 0.003 | 0.024 |
| Bhutan | \>= 5 | Inpatient | 0.020 | 0.007 | 0.046 | 0.018 | 0.005 | 0.047 |
| Bhutan | \< 5 | Outpatient | 0.002 | 0.001 | 0.005 | 0.002 | 0.000 | 0.006 |
| Bhutan | \>= 5 | Outpatient | 0.004 | 0.001 | 0.010 | 0.004 | 0.001 | 0.011 |
| Bolivia | \< 5 | Inpatient | 0.014 | 0.003 | 0.037 | 0.011 | 0.002 | 0.028 |
| Bolivia | \>= 5 | Inpatient | 0.026 | 0.005 | 0.074 | 0.021 | 0.004 | 0.056 |
| Bolivia | \< 5 | Outpatient | 0.003 | 0.000 | 0.009 | 0.002 | 0.000 | 0.007 |
| Bolivia | \>= 5 | Outpatient | 0.006 | 0.001 | 0.018 | 0.004 | 0.001 | 0.013 |
| Bosnia and Herzegovina | \< 5 | Inpatient | 0.013 | 0.003 | 0.037 | 0.011 | 0.003 | 0.034 |
| Bosnia and Herzegovina | \>= 5 | Inpatient | 0.023 | 0.006 | 0.071 | 0.021 | 0.005 | 0.065 |
| Bosnia and Herzegovina | \< 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.002 | 0.000 | 0.008 |
| Bosnia and Herzegovina | \>= 5 | Outpatient | 0.005 | 0.001 | 0.016 | 0.005 | 0.001 | 0.015 |
| Botswana | \< 5 | Inpatient | 0.008 | 0.002 | 0.018 | 0.012 | 0.003 | 0.030 |
| Botswana | \>= 5 | Inpatient | 0.015 | 0.004 | 0.036 | 0.022 | 0.005 | 0.058 |
| Botswana | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.003 | 0.000 | 0.007 |
| Botswana | \>= 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.005 | 0.001 | 0.014 |
| Brazil | \< 5 | Inpatient | 0.014 | 0.004 | 0.037 | 0.012 | 0.003 | 0.030 |
| Brazil | \>= 5 | Inpatient | 0.027 | 0.006 | 0.074 | 0.022 | 0.005 | 0.061 |
| Brazil | \< 5 | Outpatient | 0.003 | 0.001 | 0.009 | 0.003 | 0.000 | 0.007 |
| Brazil | \>= 5 | Outpatient | 0.006 | 0.001 | 0.018 | 0.005 | 0.001 | 0.014 |
| Brunei Darussalam | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Brunei Darussalam | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Brunei Darussalam | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Brunei Darussalam | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bulgaria | \< 5 | Inpatient | 0.013 | 0.003 | 0.037 | 0.014 | 0.004 | 0.042 |
| Bulgaria | \>= 5 | Inpatient | 0.023 | 0.006 | 0.071 | 0.027 | 0.006 | 0.082 |
| Bulgaria | \< 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.003 | 0.001 | 0.010 |
| Bulgaria | \>= 5 | Outpatient | 0.005 | 0.001 | 0.016 | 0.006 | 0.001 | 0.019 |
| Burkina Faso | \< 5 | Inpatient | 0.009 | 0.003 | 0.021 | 0.014 | 0.004 | 0.035 |
| Burkina Faso | \>= 5 | Inpatient | 0.017 | 0.005 | 0.043 | 0.026 | 0.007 | 0.070 |
| Burkina Faso | \< 5 | Outpatient | 0.002 | 0.000 | 0.005 | 0.003 | 0.001 | 0.008 |
| Burkina Faso | \>= 5 | Outpatient | 0.004 | 0.001 | 0.010 | 0.006 | 0.001 | 0.017 |
| Burundi | \< 5 | Inpatient | 0.009 | 0.003 | 0.019 | 0.014 | 0.004 | 0.035 |
| Burundi | \>= 5 | Inpatient | 0.017 | 0.006 | 0.039 | 0.027 | 0.007 | 0.070 |
| Burundi | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.003 | 0.001 | 0.009 |
| Burundi | \>= 5 | Outpatient | 0.004 | 0.001 | 0.009 | 0.006 | 0.001 | 0.017 |
| Cabo Verde | \< 5 | Inpatient | 0.008 | 0.003 | 0.017 | 0.007 | 0.002 | 0.016 |
| Cabo Verde | \>= 5 | Inpatient | 0.015 | 0.005 | 0.034 | 0.014 | 0.004 | 0.033 |
| Cabo Verde | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.002 | 0.000 | 0.004 |
| Cabo Verde | \>= 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.003 | 0.001 | 0.008 |
| Cambodia | \< 5 | Inpatient | 0.007 | 0.001 | 0.020 | 0.006 | 0.001 | 0.018 |
| Cambodia | \>= 5 | Inpatient | 0.013 | 0.002 | 0.038 | 0.012 | 0.002 | 0.036 |
| Cambodia | \< 5 | Outpatient | 0.001 | 0.000 | 0.004 | 0.001 | 0.000 | 0.004 |
| Cambodia | \>= 5 | Outpatient | 0.003 | 0.000 | 0.009 | 0.003 | 0.000 | 0.008 |
| Cameroon | \< 5 | Inpatient | 0.009 | 0.003 | 0.024 | 0.012 | 0.004 | 0.031 |
| Cameroon | \>= 5 | Inpatient | 0.017 | 0.005 | 0.046 | 0.023 | 0.006 | 0.061 |
| Cameroon | \< 5 | Outpatient | 0.002 | 0.000 | 0.005 | 0.003 | 0.001 | 0.007 |
| Cameroon | \>= 5 | Outpatient | 0.004 | 0.001 | 0.010 | 0.005 | 0.001 | 0.014 |
| Canada | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Canada | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Canada | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Canada | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Central African Republic | \< 5 | Inpatient | 0.009 | 0.003 | 0.023 | 0.008 | 0.002 | 0.021 |
| Central African Republic | \>= 5 | Inpatient | 0.018 | 0.005 | 0.045 | 0.016 | 0.004 | 0.043 |
| Central African Republic | \< 5 | Outpatient | 0.002 | 0.000 | 0.006 | 0.002 | 0.000 | 0.005 |
| Central African Republic | \>= 5 | Outpatient | 0.004 | 0.001 | 0.011 | 0.003 | 0.001 | 0.010 |
| Chad | \< 5 | Inpatient | 0.009 | 0.003 | 0.019 | 0.008 | 0.003 | 0.018 |
| Chad | \>= 5 | Inpatient | 0.017 | 0.006 | 0.039 | 0.015 | 0.005 | 0.037 |
| Chad | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.002 | 0.000 | 0.004 |
| Chad | \>= 5 | Outpatient | 0.004 | 0.001 | 0.009 | 0.003 | 0.001 | 0.009 |
| Chile | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Chile | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Chile | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Chile | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| China | \< 5 | Inpatient | 0.008 | 0.002 | 0.022 | 0.007 | 0.001 | 0.020 |
| China | \>= 5 | Inpatient | 0.015 | 0.003 | 0.043 | 0.013 | 0.002 | 0.039 |
| China | \< 5 | Outpatient | 0.002 | 0.000 | 0.005 | 0.001 | 0.000 | 0.005 |
| China | \>= 5 | Outpatient | 0.003 | 0.000 | 0.010 | 0.003 | 0.000 | 0.009 |
| Colombia | \< 5 | Inpatient | 0.015 | 0.004 | 0.041 | 0.014 | 0.003 | 0.039 |
| Colombia | \>= 5 | Inpatient | 0.028 | 0.007 | 0.079 | 0.026 | 0.006 | 0.076 |
| Colombia | \< 5 | Outpatient | 0.003 | 0.001 | 0.010 | 0.003 | 0.001 | 0.009 |
| Colombia | \>= 5 | Outpatient | 0.006 | 0.001 | 0.019 | 0.006 | 0.001 | 0.018 |
| Comoros | \< 5 | Inpatient | 0.008 | 0.003 | 0.017 | 0.007 | 0.002 | 0.016 |
| Comoros | \>= 5 | Inpatient | 0.015 | 0.005 | 0.034 | 0.014 | 0.004 | 0.033 |
| Comoros | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.002 | 0.000 | 0.004 |
| Comoros | \>= 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.003 | 0.001 | 0.008 |
| Congo, Dem. Rep. | \< 5 | Inpatient | 0.009 | 0.003 | 0.019 | 0.009 | 0.003 | 0.020 |
| Congo, Dem. Rep. | \>= 5 | Inpatient | 0.017 | 0.006 | 0.039 | 0.017 | 0.006 | 0.040 |
| Congo, Dem. Rep. | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.002 | 0.000 | 0.005 |
| Congo, Dem. Rep. | \>= 5 | Outpatient | 0.004 | 0.001 | 0.009 | 0.004 | 0.001 | 0.010 |
| Congo, Rep. | \< 5 | Inpatient | 0.008 | 0.003 | 0.017 | 0.010 | 0.004 | 0.022 |
| Congo, Rep. | \>= 5 | Inpatient | 0.015 | 0.005 | 0.034 | 0.019 | 0.006 | 0.045 |
| Congo, Rep. | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.002 | 0.000 | 0.005 |
| Congo, Rep. | \>= 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.004 | 0.001 | 0.011 |
| Cook Islands | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Cook Islands | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Cook Islands | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Cook Islands | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Costa Rica | \< 5 | Inpatient | 0.009 | 0.002 | 0.022 | 0.014 | 0.004 | 0.036 |
| Costa Rica | \>= 5 | Inpatient | 0.016 | 0.004 | 0.043 | 0.026 | 0.006 | 0.069 |
| Costa Rica | \< 5 | Outpatient | 0.002 | 0.000 | 0.005 | 0.003 | 0.001 | 0.008 |
| Costa Rica | \>= 5 | Outpatient | 0.004 | 0.001 | 0.010 | 0.006 | 0.001 | 0.017 |
| Côte d’Ivoire | \< 5 | Inpatient | 0.008 | 0.003 | 0.017 | 0.011 | 0.004 | 0.024 |
| Côte d’Ivoire | \>= 5 | Inpatient | 0.015 | 0.005 | 0.034 | 0.020 | 0.006 | 0.048 |
| Côte d’Ivoire | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.002 | 0.001 | 0.006 |
| Côte d’Ivoire | \>= 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.004 | 0.001 | 0.012 |
| Croatia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Croatia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Croatia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Croatia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Cuba | \< 5 | Inpatient | 0.009 | 0.002 | 0.022 | 0.008 | 0.002 | 0.022 |
| Cuba | \>= 5 | Inpatient | 0.016 | 0.004 | 0.043 | 0.015 | 0.003 | 0.043 |
| Cuba | \< 5 | Outpatient | 0.002 | 0.000 | 0.005 | 0.002 | 0.000 | 0.005 |
| Cuba | \>= 5 | Outpatient | 0.004 | 0.001 | 0.010 | 0.003 | 0.000 | 0.010 |
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
| Djibouti | \< 5 | Inpatient | 0.011 | 0.004 | 0.026 | 0.015 | 0.005 | 0.038 |
| Djibouti | \>= 5 | Inpatient | 0.021 | 0.007 | 0.052 | 0.028 | 0.008 | 0.076 |
| Djibouti | \< 5 | Outpatient | 0.002 | 0.001 | 0.006 | 0.003 | 0.001 | 0.009 |
| Djibouti | \>= 5 | Outpatient | 0.005 | 0.001 | 0.012 | 0.006 | 0.001 | 0.018 |
| Dominica | \< 5 | Inpatient | 0.009 | 0.002 | 0.022 | 0.008 | 0.002 | 0.022 |
| Dominica | \>= 5 | Inpatient | 0.016 | 0.004 | 0.043 | 0.015 | 0.003 | 0.043 |
| Dominica | \< 5 | Outpatient | 0.002 | 0.000 | 0.005 | 0.002 | 0.000 | 0.005 |
| Dominica | \>= 5 | Outpatient | 0.004 | 0.001 | 0.010 | 0.003 | 0.000 | 0.010 |
| Dominican Republic | \< 5 | Inpatient | 0.009 | 0.002 | 0.022 | 0.012 | 0.003 | 0.030 |
| Dominican Republic | \>= 5 | Inpatient | 0.016 | 0.004 | 0.043 | 0.022 | 0.006 | 0.059 |
| Dominican Republic | \< 5 | Outpatient | 0.002 | 0.000 | 0.005 | 0.003 | 0.001 | 0.007 |
| Dominican Republic | \>= 5 | Outpatient | 0.004 | 0.001 | 0.010 | 0.005 | 0.001 | 0.014 |
| Ecuador | \< 5 | Inpatient | 0.016 | 0.004 | 0.041 | 0.012 | 0.003 | 0.029 |
| Ecuador | \>= 5 | Inpatient | 0.029 | 0.007 | 0.079 | 0.022 | 0.006 | 0.058 |
| Ecuador | \< 5 | Outpatient | 0.003 | 0.001 | 0.010 | 0.003 | 0.001 | 0.007 |
| Ecuador | \>= 5 | Outpatient | 0.006 | 0.001 | 0.019 | 0.005 | 0.001 | 0.014 |
| Egypt, Arab Rep. | \< 5 | Inpatient | 0.011 | 0.003 | 0.028 | 0.010 | 0.002 | 0.027 |
| Egypt, Arab Rep. | \>= 5 | Inpatient | 0.021 | 0.005 | 0.057 | 0.019 | 0.004 | 0.053 |
| Egypt, Arab Rep. | \< 5 | Outpatient | 0.002 | 0.000 | 0.007 | 0.002 | 0.000 | 0.007 |
| Egypt, Arab Rep. | \>= 5 | Outpatient | 0.005 | 0.001 | 0.013 | 0.004 | 0.001 | 0.012 |
| El Salvador | \< 5 | Inpatient | 0.015 | 0.004 | 0.039 | 0.012 | 0.003 | 0.030 |
| El Salvador | \>= 5 | Inpatient | 0.029 | 0.007 | 0.076 | 0.022 | 0.006 | 0.059 |
| El Salvador | \< 5 | Outpatient | 0.003 | 0.001 | 0.009 | 0.003 | 0.001 | 0.007 |
| El Salvador | \>= 5 | Outpatient | 0.006 | 0.001 | 0.018 | 0.005 | 0.001 | 0.014 |
| Equatorial Guinea | \< 5 | Inpatient | 0.008 | 0.002 | 0.018 | 0.007 | 0.002 | 0.018 |
| Equatorial Guinea | \>= 5 | Inpatient | 0.015 | 0.004 | 0.036 | 0.013 | 0.003 | 0.035 |
| Equatorial Guinea | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.002 | 0.000 | 0.004 |
| Equatorial Guinea | \>= 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.003 | 0.000 | 0.008 |
| Eritrea | \< 5 | Inpatient | 0.009 | 0.003 | 0.019 | 0.014 | 0.004 | 0.035 |
| Eritrea | \>= 5 | Inpatient | 0.017 | 0.006 | 0.039 | 0.027 | 0.007 | 0.070 |
| Eritrea | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.003 | 0.001 | 0.009 |
| Eritrea | \>= 5 | Outpatient | 0.004 | 0.001 | 0.009 | 0.006 | 0.001 | 0.017 |
| Estonia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Estonia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Estonia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Estonia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Eswatini | \< 5 | Inpatient | 0.008 | 0.003 | 0.017 | 0.010 | 0.004 | 0.023 |
| Eswatini | \>= 5 | Inpatient | 0.015 | 0.005 | 0.034 | 0.020 | 0.006 | 0.046 |
| Eswatini | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.002 | 0.001 | 0.005 |
| Eswatini | \>= 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.004 | 0.001 | 0.011 |
| Ethiopia | \< 5 | Inpatient | 0.009 | 0.003 | 0.022 | 0.012 | 0.004 | 0.031 |
| Ethiopia | \>= 5 | Inpatient | 0.018 | 0.005 | 0.044 | 0.023 | 0.006 | 0.061 |
| Ethiopia | \< 5 | Outpatient | 0.002 | 0.000 | 0.005 | 0.003 | 0.001 | 0.007 |
| Ethiopia | \>= 5 | Outpatient | 0.004 | 0.001 | 0.010 | 0.005 | 0.001 | 0.015 |
| Fiji | \< 5 | Inpatient | 0.008 | 0.002 | 0.021 | 0.013 | 0.002 | 0.038 |
| Fiji | \>= 5 | Inpatient | 0.015 | 0.003 | 0.040 | 0.024 | 0.004 | 0.074 |
| Fiji | \< 5 | Outpatient | 0.002 | 0.000 | 0.005 | 0.003 | 0.000 | 0.009 |
| Fiji | \>= 5 | Outpatient | 0.003 | 0.001 | 0.009 | 0.005 | 0.001 | 0.018 |
| Finland | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Finland | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Finland | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Finland | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| France | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| France | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| France | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| France | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Gabon | \< 5 | Inpatient | 0.008 | 0.002 | 0.020 | 0.007 | 0.001 | 0.019 |
| Gabon | \>= 5 | Inpatient | 0.015 | 0.003 | 0.040 | 0.013 | 0.003 | 0.037 |
| Gabon | \< 5 | Outpatient | 0.002 | 0.000 | 0.005 | 0.002 | 0.000 | 0.004 |
| Gabon | \>= 5 | Outpatient | 0.003 | 0.001 | 0.009 | 0.003 | 0.000 | 0.009 |
| Gambia, The | \< 5 | Inpatient | 0.009 | 0.003 | 0.020 | 0.013 | 0.004 | 0.033 |
| Gambia, The | \>= 5 | Inpatient | 0.016 | 0.005 | 0.040 | 0.024 | 0.006 | 0.064 |
| Gambia, The | \< 5 | Outpatient | 0.002 | 0.000 | 0.005 | 0.003 | 0.001 | 0.008 |
| Gambia, The | \>= 5 | Outpatient | 0.004 | 0.001 | 0.010 | 0.005 | 0.001 | 0.016 |
| Georgia | \< 5 | Inpatient | 0.013 | 0.003 | 0.037 | 0.018 | 0.004 | 0.054 |
| Georgia | \>= 5 | Inpatient | 0.023 | 0.006 | 0.071 | 0.033 | 0.007 | 0.103 |
| Georgia | \< 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.004 | 0.001 | 0.013 |
| Georgia | \>= 5 | Outpatient | 0.005 | 0.001 | 0.016 | 0.007 | 0.001 | 0.025 |
| Germany | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Germany | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Germany | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Germany | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Ghana | \< 5 | Inpatient | 0.008 | 0.003 | 0.020 | 0.013 | 0.003 | 0.033 |
| Ghana | \>= 5 | Inpatient | 0.016 | 0.005 | 0.038 | 0.024 | 0.006 | 0.064 |
| Ghana | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.003 | 0.001 | 0.008 |
| Ghana | \>= 5 | Outpatient | 0.003 | 0.001 | 0.009 | 0.005 | 0.001 | 0.016 |
| Greece | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Greece | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Greece | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Greece | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Grenada | \< 5 | Inpatient | 0.009 | 0.002 | 0.022 | 0.008 | 0.002 | 0.022 |
| Grenada | \>= 5 | Inpatient | 0.016 | 0.004 | 0.043 | 0.015 | 0.003 | 0.043 |
| Grenada | \< 5 | Outpatient | 0.002 | 0.000 | 0.005 | 0.002 | 0.000 | 0.005 |
| Grenada | \>= 5 | Outpatient | 0.004 | 0.001 | 0.010 | 0.003 | 0.000 | 0.010 |
| Guatemala | \< 5 | Inpatient | 0.011 | 0.003 | 0.025 | 0.013 | 0.003 | 0.032 |
| Guatemala | \>= 5 | Inpatient | 0.020 | 0.006 | 0.050 | 0.024 | 0.006 | 0.063 |
| Guatemala | \< 5 | Outpatient | 0.002 | 0.000 | 0.006 | 0.003 | 0.001 | 0.008 |
| Guatemala | \>= 5 | Outpatient | 0.004 | 0.001 | 0.012 | 0.005 | 0.001 | 0.015 |
| Guinea | \< 5 | Inpatient | 0.008 | 0.003 | 0.017 | 0.007 | 0.002 | 0.016 |
| Guinea | \>= 5 | Inpatient | 0.015 | 0.005 | 0.034 | 0.014 | 0.004 | 0.033 |
| Guinea | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.002 | 0.000 | 0.004 |
| Guinea | \>= 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.003 | 0.001 | 0.008 |
| Guinea-Bissau | \< 5 | Inpatient | 0.010 | 0.003 | 0.024 | 0.014 | 0.004 | 0.035 |
| Guinea-Bissau | \>= 5 | Inpatient | 0.018 | 0.005 | 0.047 | 0.026 | 0.007 | 0.069 |
| Guinea-Bissau | \< 5 | Outpatient | 0.002 | 0.000 | 0.005 | 0.003 | 0.001 | 0.008 |
| Guinea-Bissau | \>= 5 | Outpatient | 0.004 | 0.001 | 0.011 | 0.006 | 0.001 | 0.017 |
| Guyana | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Guyana | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Guyana | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Guyana | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Haiti | \< 5 | Inpatient | 0.008 | 0.002 | 0.020 | 0.009 | 0.002 | 0.023 |
| Haiti | \>= 5 | Inpatient | 0.015 | 0.003 | 0.040 | 0.018 | 0.004 | 0.048 |
| Haiti | \< 5 | Outpatient | 0.002 | 0.000 | 0.005 | 0.002 | 0.000 | 0.005 |
| Haiti | \>= 5 | Outpatient | 0.003 | 0.001 | 0.009 | 0.004 | 0.001 | 0.011 |
| Honduras | \< 5 | Inpatient | 0.015 | 0.003 | 0.041 | 0.011 | 0.003 | 0.030 |
| Honduras | \>= 5 | Inpatient | 0.028 | 0.005 | 0.082 | 0.022 | 0.004 | 0.060 |
| Honduras | \< 5 | Outpatient | 0.003 | 0.000 | 0.010 | 0.002 | 0.000 | 0.007 |
| Honduras | \>= 5 | Outpatient | 0.006 | 0.001 | 0.019 | 0.005 | 0.001 | 0.014 |
| Hungary | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Hungary | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Hungary | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Hungary | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Iceland | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Iceland | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Iceland | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Iceland | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| India | \< 5 | Inpatient | 0.012 | 0.004 | 0.026 | 0.018 | 0.005 | 0.047 |
| India | \>= 5 | Inpatient | 0.022 | 0.008 | 0.050 | 0.033 | 0.008 | 0.090 |
| India | \< 5 | Outpatient | 0.003 | 0.001 | 0.006 | 0.004 | 0.001 | 0.011 |
| India | \>= 5 | Outpatient | 0.005 | 0.001 | 0.012 | 0.007 | 0.001 | 0.022 |
| Indonesia | \< 5 | Inpatient | 0.011 | 0.003 | 0.027 | 0.010 | 0.002 | 0.026 |
| Indonesia | \>= 5 | Inpatient | 0.020 | 0.006 | 0.052 | 0.018 | 0.004 | 0.051 |
| Indonesia | \< 5 | Outpatient | 0.002 | 0.001 | 0.006 | 0.002 | 0.000 | 0.006 |
| Indonesia | \>= 5 | Outpatient | 0.004 | 0.001 | 0.012 | 0.004 | 0.001 | 0.012 |
| Iran, Islamic Rep. | \< 5 | Inpatient | 0.011 | 0.004 | 0.026 | 0.010 | 0.003 | 0.026 |
| Iran, Islamic Rep. | \>= 5 | Inpatient | 0.021 | 0.007 | 0.052 | 0.019 | 0.005 | 0.050 |
| Iran, Islamic Rep. | \< 5 | Outpatient | 0.002 | 0.001 | 0.006 | 0.002 | 0.000 | 0.006 |
| Iran, Islamic Rep. | \>= 5 | Outpatient | 0.005 | 0.001 | 0.012 | 0.004 | 0.001 | 0.012 |
| Iraq | \< 5 | Inpatient | 0.012 | 0.004 | 0.034 | 0.014 | 0.004 | 0.038 |
| Iraq | \>= 5 | Inpatient | 0.023 | 0.006 | 0.066 | 0.026 | 0.007 | 0.074 |
| Iraq | \< 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.003 | 0.001 | 0.009 |
| Iraq | \>= 5 | Outpatient | 0.005 | 0.001 | 0.016 | 0.006 | 0.001 | 0.017 |
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
| Jamaica | \< 5 | Inpatient | 0.009 | 0.002 | 0.022 | 0.008 | 0.002 | 0.022 |
| Jamaica | \>= 5 | Inpatient | 0.016 | 0.004 | 0.043 | 0.015 | 0.003 | 0.043 |
| Jamaica | \< 5 | Outpatient | 0.002 | 0.000 | 0.005 | 0.002 | 0.000 | 0.005 |
| Jamaica | \>= 5 | Outpatient | 0.004 | 0.001 | 0.010 | 0.003 | 0.000 | 0.010 |
| Japan | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Japan | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Japan | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Japan | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Jordan | \< 5 | Inpatient | 0.012 | 0.003 | 0.031 | 0.016 | 0.004 | 0.045 |
| Jordan | \>= 5 | Inpatient | 0.022 | 0.006 | 0.059 | 0.031 | 0.007 | 0.086 |
| Jordan | \< 5 | Outpatient | 0.003 | 0.001 | 0.007 | 0.004 | 0.001 | 0.011 |
| Jordan | \>= 5 | Outpatient | 0.005 | 0.001 | 0.014 | 0.007 | 0.001 | 0.021 |
| Kazakhstan | \< 5 | Inpatient | 0.013 | 0.003 | 0.037 | 0.011 | 0.003 | 0.034 |
| Kazakhstan | \>= 5 | Inpatient | 0.023 | 0.006 | 0.071 | 0.021 | 0.005 | 0.065 |
| Kazakhstan | \< 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.002 | 0.000 | 0.008 |
| Kazakhstan | \>= 5 | Outpatient | 0.005 | 0.001 | 0.016 | 0.005 | 0.001 | 0.015 |
| Kenya | \< 5 | Inpatient | 0.008 | 0.003 | 0.018 | 0.013 | 0.004 | 0.033 |
| Kenya | \>= 5 | Inpatient | 0.015 | 0.005 | 0.036 | 0.024 | 0.006 | 0.064 |
| Kenya | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.003 | 0.001 | 0.008 |
| Kenya | \>= 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.005 | 0.001 | 0.015 |
| Kiribati | \< 5 | Inpatient | 0.007 | 0.001 | 0.019 | 0.011 | 0.002 | 0.032 |
| Kiribati | \>= 5 | Inpatient | 0.014 | 0.002 | 0.037 | 0.021 | 0.003 | 0.063 |
| Kiribati | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.002 | 0.000 | 0.008 |
| Kiribati | \>= 5 | Outpatient | 0.003 | 0.000 | 0.009 | 0.005 | 0.001 | 0.015 |
| Korea, Dem. People’s Rep. | \< 5 | Inpatient | 0.011 | 0.004 | 0.023 | 0.010 | 0.003 | 0.024 |
| Korea, Dem. People’s Rep. | \>= 5 | Inpatient | 0.020 | 0.007 | 0.046 | 0.018 | 0.005 | 0.047 |
| Korea, Dem. People’s Rep. | \< 5 | Outpatient | 0.002 | 0.001 | 0.005 | 0.002 | 0.000 | 0.006 |
| Korea, Dem. People’s Rep. | \>= 5 | Outpatient | 0.004 | 0.001 | 0.010 | 0.004 | 0.001 | 0.011 |
| Korea, Rep. | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Korea, Rep. | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Korea, Rep. | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Korea, Rep. | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kuwait | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kuwait | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kuwait | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kuwait | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kyrgyz Republic | \< 5 | Inpatient | 0.011 | 0.003 | 0.032 | 0.013 | 0.004 | 0.039 |
| Kyrgyz Republic | \>= 5 | Inpatient | 0.021 | 0.006 | 0.060 | 0.025 | 0.006 | 0.074 |
| Kyrgyz Republic | \< 5 | Outpatient | 0.002 | 0.000 | 0.007 | 0.003 | 0.001 | 0.009 |
| Kyrgyz Republic | \>= 5 | Outpatient | 0.005 | 0.001 | 0.014 | 0.006 | 0.001 | 0.017 |
| Lao PDR | \< 5 | Inpatient | 0.007 | 0.001 | 0.019 | 0.006 | 0.001 | 0.018 |
| Lao PDR | \>= 5 | Inpatient | 0.014 | 0.002 | 0.037 | 0.012 | 0.002 | 0.036 |
| Lao PDR | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.001 | 0.000 | 0.004 |
| Lao PDR | \>= 5 | Outpatient | 0.003 | 0.000 | 0.009 | 0.003 | 0.000 | 0.008 |
| Latvia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Latvia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Latvia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Latvia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Lebanon | \< 5 | Inpatient | 0.012 | 0.004 | 0.031 | 0.011 | 0.003 | 0.030 |
| Lebanon | \>= 5 | Inpatient | 0.023 | 0.006 | 0.061 | 0.020 | 0.005 | 0.057 |
| Lebanon | \< 5 | Outpatient | 0.003 | 0.001 | 0.007 | 0.002 | 0.000 | 0.007 |
| Lebanon | \>= 5 | Outpatient | 0.005 | 0.001 | 0.014 | 0.004 | 0.001 | 0.014 |
| Lesotho | \< 5 | Inpatient | 0.008 | 0.003 | 0.017 | 0.012 | 0.004 | 0.029 |
| Lesotho | \>= 5 | Inpatient | 0.015 | 0.005 | 0.034 | 0.023 | 0.006 | 0.057 |
| Lesotho | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.003 | 0.001 | 0.007 |
| Lesotho | \>= 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.005 | 0.001 | 0.014 |
| Liberia | \< 5 | Inpatient | 0.009 | 0.003 | 0.019 | 0.011 | 0.004 | 0.026 |
| Liberia | \>= 5 | Inpatient | 0.017 | 0.006 | 0.039 | 0.022 | 0.007 | 0.052 |
| Liberia | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.003 | 0.001 | 0.006 |
| Liberia | \>= 5 | Outpatient | 0.004 | 0.001 | 0.009 | 0.005 | 0.001 | 0.013 |
| Libya | \< 5 | Inpatient | 0.011 | 0.003 | 0.029 | 0.016 | 0.004 | 0.044 |
| Libya | \>= 5 | Inpatient | 0.022 | 0.006 | 0.057 | 0.030 | 0.007 | 0.085 |
| Libya | \< 5 | Outpatient | 0.002 | 0.001 | 0.007 | 0.003 | 0.001 | 0.010 |
| Libya | \>= 5 | Outpatient | 0.005 | 0.001 | 0.013 | 0.007 | 0.001 | 0.020 |
| Lithuania | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Lithuania | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Lithuania | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Lithuania | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Luxembourg | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Luxembourg | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Luxembourg | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Luxembourg | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Madagascar | \< 5 | Inpatient | 0.009 | 0.003 | 0.019 | 0.011 | 0.004 | 0.026 |
| Madagascar | \>= 5 | Inpatient | 0.017 | 0.006 | 0.039 | 0.022 | 0.007 | 0.052 |
| Madagascar | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.003 | 0.001 | 0.006 |
| Madagascar | \>= 5 | Outpatient | 0.004 | 0.001 | 0.009 | 0.005 | 0.001 | 0.013 |
| Malawi | \< 5 | Inpatient | 0.010 | 0.003 | 0.024 | 0.015 | 0.004 | 0.040 |
| Malawi | \>= 5 | Inpatient | 0.019 | 0.005 | 0.049 | 0.028 | 0.007 | 0.076 |
| Malawi | \< 5 | Outpatient | 0.002 | 0.000 | 0.006 | 0.003 | 0.001 | 0.010 |
| Malawi | \>= 5 | Outpatient | 0.004 | 0.001 | 0.012 | 0.006 | 0.001 | 0.019 |
| Malaysia | \< 5 | Inpatient | 0.008 | 0.002 | 0.021 | 0.007 | 0.001 | 0.019 |
| Malaysia | \>= 5 | Inpatient | 0.015 | 0.003 | 0.040 | 0.013 | 0.003 | 0.038 |
| Malaysia | \< 5 | Outpatient | 0.002 | 0.000 | 0.005 | 0.001 | 0.000 | 0.004 |
| Malaysia | \>= 5 | Outpatient | 0.003 | 0.001 | 0.009 | 0.003 | 0.000 | 0.008 |
| Maldives | \< 5 | Inpatient | 0.011 | 0.003 | 0.027 | 0.010 | 0.002 | 0.026 |
| Maldives | \>= 5 | Inpatient | 0.020 | 0.006 | 0.052 | 0.018 | 0.004 | 0.051 |
| Maldives | \< 5 | Outpatient | 0.002 | 0.001 | 0.006 | 0.002 | 0.000 | 0.006 |
| Maldives | \>= 5 | Outpatient | 0.004 | 0.001 | 0.012 | 0.004 | 0.001 | 0.012 |
| Mali | \< 5 | Inpatient | 0.009 | 0.003 | 0.022 | 0.012 | 0.004 | 0.030 |
| Mali | \>= 5 | Inpatient | 0.018 | 0.005 | 0.045 | 0.023 | 0.006 | 0.060 |
| Mali | \< 5 | Outpatient | 0.002 | 0.000 | 0.005 | 0.003 | 0.001 | 0.007 |
| Mali | \>= 5 | Outpatient | 0.004 | 0.001 | 0.010 | 0.005 | 0.001 | 0.014 |
| Malta | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Malta | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Malta | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Malta | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Marshall Islands | \< 5 | Inpatient | 0.009 | 0.002 | 0.024 | 0.010 | 0.002 | 0.026 |
| Marshall Islands | \>= 5 | Inpatient | 0.017 | 0.003 | 0.046 | 0.018 | 0.004 | 0.050 |
| Marshall Islands | \< 5 | Outpatient | 0.002 | 0.000 | 0.005 | 0.002 | 0.000 | 0.006 |
| Marshall Islands | \>= 5 | Outpatient | 0.004 | 0.001 | 0.011 | 0.004 | 0.001 | 0.012 |
| Mauritania | \< 5 | Inpatient | 0.008 | 0.003 | 0.017 | 0.012 | 0.004 | 0.028 |
| Mauritania | \>= 5 | Inpatient | 0.015 | 0.005 | 0.034 | 0.023 | 0.006 | 0.055 |
| Mauritania | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.003 | 0.001 | 0.007 |
| Mauritania | \>= 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.005 | 0.001 | 0.013 |
| Mauritius | \< 5 | Inpatient | 0.008 | 0.002 | 0.018 | 0.012 | 0.003 | 0.032 |
| Mauritius | \>= 5 | Inpatient | 0.015 | 0.004 | 0.036 | 0.023 | 0.005 | 0.062 |
| Mauritius | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.003 | 0.000 | 0.007 |
| Mauritius | \>= 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.005 | 0.001 | 0.015 |
| Mexico | \< 5 | Inpatient | 0.016 | 0.004 | 0.046 | 0.012 | 0.003 | 0.033 |
| Mexico | \>= 5 | Inpatient | 0.031 | 0.007 | 0.087 | 0.023 | 0.005 | 0.064 |
| Mexico | \< 5 | Outpatient | 0.004 | 0.001 | 0.011 | 0.003 | 0.000 | 0.008 |
| Mexico | \>= 5 | Outpatient | 0.007 | 0.001 | 0.021 | 0.005 | 0.001 | 0.015 |
| Micronesia, Fed. Sts. | \< 5 | Inpatient | 0.012 | 0.002 | 0.033 | 0.008 | 0.001 | 0.021 |
| Micronesia, Fed. Sts. | \>= 5 | Inpatient | 0.022 | 0.003 | 0.065 | 0.015 | 0.003 | 0.042 |
| Micronesia, Fed. Sts. | \< 5 | Outpatient | 0.003 | 0.000 | 0.008 | 0.002 | 0.000 | 0.005 |
| Micronesia, Fed. Sts. | \>= 5 | Outpatient | 0.005 | 0.001 | 0.015 | 0.003 | 0.000 | 0.010 |
| Moldova | \< 5 | Inpatient | 0.013 | 0.003 | 0.037 | 0.015 | 0.004 | 0.046 |
| Moldova | \>= 5 | Inpatient | 0.023 | 0.006 | 0.071 | 0.029 | 0.007 | 0.089 |
| Moldova | \< 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.003 | 0.001 | 0.011 |
| Moldova | \>= 5 | Outpatient | 0.005 | 0.001 | 0.016 | 0.006 | 0.001 | 0.021 |
| Monaco | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Monaco | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Monaco | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Monaco | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Mongolia | \< 5 | Inpatient | 0.007 | 0.001 | 0.019 | 0.006 | 0.001 | 0.018 |
| Mongolia | \>= 5 | Inpatient | 0.014 | 0.002 | 0.037 | 0.012 | 0.002 | 0.036 |
| Mongolia | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.001 | 0.000 | 0.004 |
| Mongolia | \>= 5 | Outpatient | 0.003 | 0.000 | 0.009 | 0.003 | 0.000 | 0.008 |
| Montenegro | \< 5 | Inpatient | 0.013 | 0.003 | 0.037 | 0.011 | 0.003 | 0.034 |
| Montenegro | \>= 5 | Inpatient | 0.023 | 0.006 | 0.071 | 0.021 | 0.005 | 0.065 |
| Montenegro | \< 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.002 | 0.000 | 0.008 |
| Montenegro | \>= 5 | Outpatient | 0.005 | 0.001 | 0.016 | 0.005 | 0.001 | 0.015 |
| Morocco | \< 5 | Inpatient | 0.011 | 0.003 | 0.028 | 0.018 | 0.004 | 0.052 |
| Morocco | \>= 5 | Inpatient | 0.021 | 0.005 | 0.055 | 0.033 | 0.007 | 0.099 |
| Morocco | \< 5 | Outpatient | 0.002 | 0.000 | 0.007 | 0.004 | 0.001 | 0.012 |
| Morocco | \>= 5 | Outpatient | 0.005 | 0.001 | 0.013 | 0.008 | 0.001 | 0.024 |
| Mozambique | \< 5 | Inpatient | 0.009 | 0.003 | 0.020 | 0.011 | 0.004 | 0.025 |
| Mozambique | \>= 5 | Inpatient | 0.017 | 0.005 | 0.041 | 0.021 | 0.006 | 0.051 |
| Mozambique | \< 5 | Outpatient | 0.002 | 0.000 | 0.005 | 0.002 | 0.000 | 0.006 |
| Mozambique | \>= 5 | Outpatient | 0.004 | 0.001 | 0.010 | 0.005 | 0.001 | 0.012 |
| Myanmar | \< 5 | Inpatient | 0.011 | 0.004 | 0.023 | 0.013 | 0.004 | 0.030 |
| Myanmar | \>= 5 | Inpatient | 0.020 | 0.007 | 0.046 | 0.024 | 0.007 | 0.060 |
| Myanmar | \< 5 | Outpatient | 0.002 | 0.001 | 0.005 | 0.003 | 0.001 | 0.007 |
| Myanmar | \>= 5 | Outpatient | 0.004 | 0.001 | 0.010 | 0.005 | 0.001 | 0.014 |
| Namibia | \< 5 | Inpatient | 0.008 | 0.002 | 0.018 | 0.011 | 0.003 | 0.027 |
| Namibia | \>= 5 | Inpatient | 0.015 | 0.004 | 0.036 | 0.021 | 0.005 | 0.053 |
| Namibia | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.002 | 0.000 | 0.006 |
| Namibia | \>= 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.005 | 0.001 | 0.013 |
| Nauru | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nauru | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nauru | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nauru | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nepal | \< 5 | Inpatient | 0.010 | 0.004 | 0.023 | 0.009 | 0.002 | 0.024 |
| Nepal | \>= 5 | Inpatient | 0.020 | 0.006 | 0.047 | 0.018 | 0.004 | 0.048 |
| Nepal | \< 5 | Outpatient | 0.002 | 0.001 | 0.005 | 0.002 | 0.000 | 0.006 |
| Nepal | \>= 5 | Outpatient | 0.004 | 0.001 | 0.011 | 0.004 | 0.001 | 0.011 |
| Netherlands | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Netherlands | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Netherlands | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Netherlands | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| New Zealand | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| New Zealand | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| New Zealand | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| New Zealand | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nicaragua | \< 5 | Inpatient | 0.016 | 0.003 | 0.046 | 0.013 | 0.003 | 0.037 |
| Nicaragua | \>= 5 | Inpatient | 0.030 | 0.005 | 0.093 | 0.024 | 0.005 | 0.073 |
| Nicaragua | \< 5 | Outpatient | 0.003 | 0.000 | 0.011 | 0.003 | 0.000 | 0.009 |
| Nicaragua | \>= 5 | Outpatient | 0.007 | 0.001 | 0.022 | 0.005 | 0.001 | 0.017 |
| Niger | \< 5 | Inpatient | 0.009 | 0.003 | 0.019 | 0.013 | 0.004 | 0.031 |
| Niger | \>= 5 | Inpatient | 0.017 | 0.006 | 0.039 | 0.025 | 0.007 | 0.061 |
| Niger | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.003 | 0.001 | 0.007 |
| Niger | \>= 5 | Outpatient | 0.004 | 0.001 | 0.009 | 0.005 | 0.001 | 0.015 |
| Nigeria | \< 5 | Inpatient | 0.008 | 0.003 | 0.019 | 0.007 | 0.002 | 0.019 |
| Nigeria | \>= 5 | Inpatient | 0.016 | 0.005 | 0.038 | 0.014 | 0.003 | 0.037 |
| Nigeria | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.002 | 0.000 | 0.004 |
| Nigeria | \>= 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.003 | 0.001 | 0.008 |
| Niue | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Niue | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Niue | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Niue | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| North Macedonia | \< 5 | Inpatient | 0.013 | 0.003 | 0.037 | 0.016 | 0.004 | 0.048 |
| North Macedonia | \>= 5 | Inpatient | 0.023 | 0.006 | 0.071 | 0.030 | 0.007 | 0.091 |
| North Macedonia | \< 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.003 | 0.001 | 0.011 |
| North Macedonia | \>= 5 | Outpatient | 0.005 | 0.001 | 0.016 | 0.007 | 0.001 | 0.022 |
| Norway | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Norway | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Norway | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Norway | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Oman | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Oman | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Oman | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Oman | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Pakistan | \< 5 | Inpatient | 0.011 | 0.003 | 0.027 | 0.016 | 0.004 | 0.045 |
| Pakistan | \>= 5 | Inpatient | 0.021 | 0.006 | 0.054 | 0.031 | 0.007 | 0.088 |
| Pakistan | \< 5 | Outpatient | 0.002 | 0.001 | 0.006 | 0.004 | 0.001 | 0.011 |
| Pakistan | \>= 5 | Outpatient | 0.005 | 0.001 | 0.013 | 0.007 | 0.001 | 0.021 |
| Palau | \< 5 | Inpatient | 0.013 | 0.002 | 0.039 | 0.012 | 0.002 | 0.033 |
| Palau | \>= 5 | Inpatient | 0.025 | 0.004 | 0.074 | 0.022 | 0.004 | 0.063 |
| Palau | \< 5 | Outpatient | 0.003 | 0.000 | 0.009 | 0.003 | 0.000 | 0.008 |
| Palau | \>= 5 | Outpatient | 0.005 | 0.001 | 0.018 | 0.005 | 0.001 | 0.015 |
| Panama | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Panama | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Panama | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Panama | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Papua New Guinea | \< 5 | Inpatient | 0.007 | 0.001 | 0.019 | 0.006 | 0.001 | 0.018 |
| Papua New Guinea | \>= 5 | Inpatient | 0.014 | 0.002 | 0.037 | 0.012 | 0.002 | 0.036 |
| Papua New Guinea | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.001 | 0.000 | 0.004 |
| Papua New Guinea | \>= 5 | Outpatient | 0.003 | 0.000 | 0.009 | 0.003 | 0.000 | 0.008 |
| Paraguay | \< 5 | Inpatient | 0.012 | 0.004 | 0.028 | 0.012 | 0.003 | 0.031 |
| Paraguay | \>= 5 | Inpatient | 0.023 | 0.006 | 0.057 | 0.023 | 0.006 | 0.060 |
| Paraguay | \< 5 | Outpatient | 0.003 | 0.001 | 0.007 | 0.003 | 0.001 | 0.007 |
| Paraguay | \>= 5 | Outpatient | 0.005 | 0.001 | 0.014 | 0.005 | 0.001 | 0.014 |
| Peru | \< 5 | Inpatient | 0.013 | 0.003 | 0.034 | 0.011 | 0.003 | 0.030 |
| Peru | \>= 5 | Inpatient | 0.025 | 0.006 | 0.068 | 0.021 | 0.005 | 0.059 |
| Peru | \< 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.002 | 0.000 | 0.007 |
| Peru | \>= 5 | Outpatient | 0.005 | 0.001 | 0.016 | 0.005 | 0.001 | 0.014 |
| Philippines | \< 5 | Inpatient | 0.007 | 0.001 | 0.019 | 0.006 | 0.001 | 0.018 |
| Philippines | \>= 5 | Inpatient | 0.014 | 0.002 | 0.037 | 0.012 | 0.002 | 0.036 |
| Philippines | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.001 | 0.000 | 0.004 |
| Philippines | \>= 5 | Outpatient | 0.003 | 0.000 | 0.009 | 0.003 | 0.000 | 0.008 |
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
| Russian Federation | \< 5 | Inpatient | 0.013 | 0.003 | 0.037 | 0.011 | 0.003 | 0.034 |
| Russian Federation | \>= 5 | Inpatient | 0.023 | 0.006 | 0.071 | 0.021 | 0.005 | 0.065 |
| Russian Federation | \< 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.002 | 0.000 | 0.008 |
| Russian Federation | \>= 5 | Outpatient | 0.005 | 0.001 | 0.016 | 0.005 | 0.001 | 0.015 |
| Rwanda | \< 5 | Inpatient | 0.009 | 0.003 | 0.019 | 0.014 | 0.004 | 0.034 |
| Rwanda | \>= 5 | Inpatient | 0.017 | 0.006 | 0.039 | 0.026 | 0.007 | 0.066 |
| Rwanda | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.003 | 0.001 | 0.008 |
| Rwanda | \>= 5 | Outpatient | 0.004 | 0.001 | 0.009 | 0.006 | 0.001 | 0.017 |
| Samoa | \< 5 | Inpatient | 0.007 | 0.001 | 0.019 | 0.006 | 0.001 | 0.018 |
| Samoa | \>= 5 | Inpatient | 0.014 | 0.002 | 0.037 | 0.012 | 0.002 | 0.036 |
| Samoa | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.001 | 0.000 | 0.004 |
| Samoa | \>= 5 | Outpatient | 0.003 | 0.000 | 0.009 | 0.003 | 0.000 | 0.008 |
| San Marino | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| San Marino | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| San Marino | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| San Marino | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| São Tomé and Principe | \< 5 | Inpatient | 0.008 | 0.003 | 0.017 | 0.012 | 0.004 | 0.029 |
| São Tomé and Principe | \>= 5 | Inpatient | 0.015 | 0.005 | 0.034 | 0.023 | 0.006 | 0.057 |
| São Tomé and Principe | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.003 | 0.001 | 0.007 |
| São Tomé and Principe | \>= 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.005 | 0.001 | 0.014 |
| Saudi Arabia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Saudi Arabia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Saudi Arabia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Saudi Arabia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Senegal | \< 5 | Inpatient | 0.008 | 0.003 | 0.017 | 0.013 | 0.004 | 0.031 |
| Senegal | \>= 5 | Inpatient | 0.015 | 0.005 | 0.034 | 0.024 | 0.006 | 0.061 |
| Senegal | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.003 | 0.001 | 0.007 |
| Senegal | \>= 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.005 | 0.001 | 0.015 |
| Serbia | \< 5 | Inpatient | 0.013 | 0.003 | 0.037 | 0.011 | 0.003 | 0.034 |
| Serbia | \>= 5 | Inpatient | 0.023 | 0.006 | 0.071 | 0.021 | 0.005 | 0.065 |
| Serbia | \< 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.002 | 0.000 | 0.008 |
| Serbia | \>= 5 | Outpatient | 0.005 | 0.001 | 0.016 | 0.005 | 0.001 | 0.015 |
| Seychelles | \< 5 | Inpatient | 0.008 | 0.002 | 0.018 | 0.013 | 0.003 | 0.034 |
| Seychelles | \>= 5 | Inpatient | 0.015 | 0.004 | 0.036 | 0.024 | 0.005 | 0.067 |
| Seychelles | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.003 | 0.000 | 0.008 |
| Seychelles | \>= 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.005 | 0.001 | 0.016 |
| Sierra Leone | \< 5 | Inpatient | 0.009 | 0.003 | 0.019 | 0.014 | 0.004 | 0.033 |
| Sierra Leone | \>= 5 | Inpatient | 0.017 | 0.006 | 0.039 | 0.026 | 0.007 | 0.066 |
| Sierra Leone | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.003 | 0.001 | 0.008 |
| Sierra Leone | \>= 5 | Outpatient | 0.004 | 0.001 | 0.009 | 0.006 | 0.001 | 0.016 |
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
| Solomon Islands | \< 5 | Inpatient | 0.007 | 0.001 | 0.019 | 0.008 | 0.001 | 0.022 |
| Solomon Islands | \>= 5 | Inpatient | 0.014 | 0.002 | 0.037 | 0.015 | 0.003 | 0.042 |
| Solomon Islands | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.002 | 0.000 | 0.005 |
| Solomon Islands | \>= 5 | Outpatient | 0.003 | 0.000 | 0.009 | 0.003 | 0.000 | 0.010 |
| Somalia | \< 5 | Inpatient | 0.011 | 0.003 | 0.032 | 0.010 | 0.002 | 0.030 |
| Somalia | \>= 5 | Inpatient | 0.021 | 0.005 | 0.060 | 0.019 | 0.004 | 0.056 |
| Somalia | \< 5 | Outpatient | 0.002 | 0.000 | 0.007 | 0.002 | 0.000 | 0.007 |
| Somalia | \>= 5 | Outpatient | 0.005 | 0.001 | 0.014 | 0.004 | 0.001 | 0.013 |
| South Africa | \< 5 | Inpatient | 0.011 | 0.003 | 0.029 | 0.011 | 0.003 | 0.030 |
| South Africa | \>= 5 | Inpatient | 0.021 | 0.005 | 0.058 | 0.021 | 0.005 | 0.059 |
| South Africa | \< 5 | Outpatient | 0.002 | 0.000 | 0.007 | 0.002 | 0.000 | 0.007 |
| South Africa | \>= 5 | Outpatient | 0.005 | 0.001 | 0.014 | 0.005 | 0.001 | 0.014 |
| South Sudan | \< 5 | Inpatient | 0.009 | 0.003 | 0.019 | 0.008 | 0.003 | 0.018 |
| South Sudan | \>= 5 | Inpatient | 0.017 | 0.006 | 0.039 | 0.015 | 0.005 | 0.037 |
| South Sudan | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.002 | 0.000 | 0.004 |
| South Sudan | \>= 5 | Outpatient | 0.004 | 0.001 | 0.009 | 0.003 | 0.001 | 0.009 |
| Spain | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Spain | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Spain | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Spain | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Sri Lanka | \< 5 | Inpatient | 0.011 | 0.004 | 0.023 | 0.010 | 0.003 | 0.024 |
| Sri Lanka | \>= 5 | Inpatient | 0.020 | 0.007 | 0.046 | 0.018 | 0.005 | 0.047 |
| Sri Lanka | \< 5 | Outpatient | 0.002 | 0.001 | 0.005 | 0.002 | 0.000 | 0.006 |
| Sri Lanka | \>= 5 | Outpatient | 0.004 | 0.001 | 0.010 | 0.004 | 0.001 | 0.011 |
| St. Kitts and Nevis | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| St. Kitts and Nevis | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| St. Kitts and Nevis | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| St. Kitts and Nevis | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| St. Lucia | \< 5 | Inpatient | 0.009 | 0.002 | 0.022 | 0.008 | 0.002 | 0.022 |
| St. Lucia | \>= 5 | Inpatient | 0.016 | 0.004 | 0.043 | 0.015 | 0.003 | 0.043 |
| St. Lucia | \< 5 | Outpatient | 0.002 | 0.000 | 0.005 | 0.002 | 0.000 | 0.005 |
| St. Lucia | \>= 5 | Outpatient | 0.004 | 0.001 | 0.010 | 0.003 | 0.000 | 0.010 |
| St. Vincent and the Grenadines | \< 5 | Inpatient | 0.009 | 0.002 | 0.022 | 0.008 | 0.002 | 0.022 |
| St. Vincent and the Grenadines | \>= 5 | Inpatient | 0.016 | 0.004 | 0.043 | 0.015 | 0.003 | 0.043 |
| St. Vincent and the Grenadines | \< 5 | Outpatient | 0.002 | 0.000 | 0.005 | 0.002 | 0.000 | 0.005 |
| St. Vincent and the Grenadines | \>= 5 | Outpatient | 0.004 | 0.001 | 0.010 | 0.003 | 0.000 | 0.010 |
| Sudan | \< 5 | Inpatient | 0.012 | 0.003 | 0.036 | 0.018 | 0.004 | 0.057 |
| Sudan | \>= 5 | Inpatient | 0.022 | 0.005 | 0.067 | 0.034 | 0.007 | 0.106 |
| Sudan | \< 5 | Outpatient | 0.003 | 0.000 | 0.008 | 0.004 | 0.001 | 0.014 |
| Sudan | \>= 5 | Outpatient | 0.005 | 0.001 | 0.016 | 0.008 | 0.001 | 0.025 |
| Suriname | \< 5 | Inpatient | 0.009 | 0.002 | 0.022 | 0.008 | 0.002 | 0.022 |
| Suriname | \>= 5 | Inpatient | 0.016 | 0.004 | 0.043 | 0.015 | 0.003 | 0.043 |
| Suriname | \< 5 | Outpatient | 0.002 | 0.000 | 0.005 | 0.002 | 0.000 | 0.005 |
| Suriname | \>= 5 | Outpatient | 0.004 | 0.001 | 0.010 | 0.003 | 0.000 | 0.010 |
| Sweden | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Sweden | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Sweden | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Sweden | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Switzerland | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Switzerland | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Switzerland | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Switzerland | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Syrian Arab Republic | \< 5 | Inpatient | 0.011 | 0.003 | 0.032 | 0.010 | 0.002 | 0.030 |
| Syrian Arab Republic | \>= 5 | Inpatient | 0.021 | 0.005 | 0.060 | 0.019 | 0.004 | 0.056 |
| Syrian Arab Republic | \< 5 | Outpatient | 0.002 | 0.000 | 0.007 | 0.002 | 0.000 | 0.007 |
| Syrian Arab Republic | \>= 5 | Outpatient | 0.005 | 0.001 | 0.014 | 0.004 | 0.001 | 0.013 |
| Tajikistan | \< 5 | Inpatient | 0.011 | 0.003 | 0.032 | 0.018 | 0.004 | 0.058 |
| Tajikistan | \>= 5 | Inpatient | 0.021 | 0.006 | 0.060 | 0.034 | 0.007 | 0.106 |
| Tajikistan | \< 5 | Outpatient | 0.002 | 0.000 | 0.007 | 0.004 | 0.001 | 0.013 |
| Tajikistan | \>= 5 | Outpatient | 0.005 | 0.001 | 0.014 | 0.008 | 0.001 | 0.026 |
| Tanzania | \< 5 | Inpatient | 0.008 | 0.003 | 0.017 | 0.012 | 0.003 | 0.030 |
| Tanzania | \>= 5 | Inpatient | 0.015 | 0.004 | 0.035 | 0.022 | 0.006 | 0.059 |
| Tanzania | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.003 | 0.000 | 0.007 |
| Tanzania | \>= 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.005 | 0.001 | 0.015 |
| Thailand | \< 5 | Inpatient | 0.011 | 0.003 | 0.029 | 0.010 | 0.002 | 0.029 |
| Thailand | \>= 5 | Inpatient | 0.021 | 0.006 | 0.057 | 0.019 | 0.004 | 0.055 |
| Thailand | \< 5 | Outpatient | 0.002 | 0.000 | 0.007 | 0.002 | 0.000 | 0.007 |
| Thailand | \>= 5 | Outpatient | 0.005 | 0.001 | 0.013 | 0.004 | 0.001 | 0.013 |
| Timor-Leste | \< 5 | Inpatient | 0.011 | 0.004 | 0.023 | 0.010 | 0.003 | 0.024 |
| Timor-Leste | \>= 5 | Inpatient | 0.020 | 0.007 | 0.046 | 0.018 | 0.005 | 0.047 |
| Timor-Leste | \< 5 | Outpatient | 0.002 | 0.001 | 0.005 | 0.002 | 0.000 | 0.006 |
| Timor-Leste | \>= 5 | Outpatient | 0.004 | 0.001 | 0.010 | 0.004 | 0.001 | 0.011 |
| Togo | \< 5 | Inpatient | 0.009 | 0.003 | 0.019 | 0.013 | 0.004 | 0.029 |
| Togo | \>= 5 | Inpatient | 0.017 | 0.006 | 0.039 | 0.024 | 0.007 | 0.058 |
| Togo | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.003 | 0.001 | 0.007 |
| Togo | \>= 5 | Outpatient | 0.004 | 0.001 | 0.009 | 0.005 | 0.001 | 0.014 |
| Tonga | \< 5 | Inpatient | 0.008 | 0.002 | 0.021 | 0.007 | 0.001 | 0.019 |
| Tonga | \>= 5 | Inpatient | 0.015 | 0.003 | 0.040 | 0.013 | 0.003 | 0.038 |
| Tonga | \< 5 | Outpatient | 0.002 | 0.000 | 0.005 | 0.001 | 0.000 | 0.004 |
| Tonga | \>= 5 | Outpatient | 0.003 | 0.001 | 0.009 | 0.003 | 0.000 | 0.008 |
| Trinidad and Tobago | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Trinidad and Tobago | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Trinidad and Tobago | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Trinidad and Tobago | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Tunisia | \< 5 | Inpatient | 0.011 | 0.004 | 0.026 | 0.010 | 0.003 | 0.026 |
| Tunisia | \>= 5 | Inpatient | 0.021 | 0.007 | 0.052 | 0.019 | 0.005 | 0.050 |
| Tunisia | \< 5 | Outpatient | 0.002 | 0.001 | 0.006 | 0.002 | 0.000 | 0.006 |
| Tunisia | \>= 5 | Outpatient | 0.005 | 0.001 | 0.012 | 0.004 | 0.001 | 0.012 |
| Türkiye | \< 5 | Inpatient | 0.014 | 0.004 | 0.043 | 0.012 | 0.003 | 0.039 |
| Türkiye | \>= 5 | Inpatient | 0.026 | 0.006 | 0.081 | 0.023 | 0.005 | 0.073 |
| Türkiye | \< 5 | Outpatient | 0.003 | 0.001 | 0.010 | 0.003 | 0.000 | 0.009 |
| Türkiye | \>= 5 | Outpatient | 0.006 | 0.001 | 0.019 | 0.005 | 0.001 | 0.017 |
| Turkmenistan | \< 5 | Inpatient | 0.013 | 0.003 | 0.037 | 0.020 | 0.004 | 0.066 |
| Turkmenistan | \>= 5 | Inpatient | 0.023 | 0.006 | 0.071 | 0.038 | 0.007 | 0.123 |
| Turkmenistan | \< 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.004 | 0.001 | 0.016 |
| Turkmenistan | \>= 5 | Outpatient | 0.005 | 0.001 | 0.016 | 0.008 | 0.001 | 0.031 |
| Tuvalu | \< 5 | Inpatient | 0.008 | 0.002 | 0.021 | 0.007 | 0.001 | 0.019 |
| Tuvalu | \>= 5 | Inpatient | 0.015 | 0.003 | 0.040 | 0.013 | 0.003 | 0.038 |
| Tuvalu | \< 5 | Outpatient | 0.002 | 0.000 | 0.005 | 0.001 | 0.000 | 0.004 |
| Tuvalu | \>= 5 | Outpatient | 0.003 | 0.001 | 0.009 | 0.003 | 0.000 | 0.008 |
| Uganda | \< 5 | Inpatient | 0.009 | 0.003 | 0.019 | 0.013 | 0.004 | 0.033 |
| Uganda | \>= 5 | Inpatient | 0.017 | 0.006 | 0.039 | 0.026 | 0.007 | 0.064 |
| Uganda | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.003 | 0.001 | 0.008 |
| Uganda | \>= 5 | Outpatient | 0.004 | 0.001 | 0.009 | 0.006 | 0.001 | 0.016 |
| Ukraine | \< 5 | Inpatient | 0.011 | 0.003 | 0.032 | 0.010 | 0.002 | 0.030 |
| Ukraine | \>= 5 | Inpatient | 0.021 | 0.006 | 0.060 | 0.019 | 0.004 | 0.057 |
| Ukraine | \< 5 | Outpatient | 0.002 | 0.000 | 0.007 | 0.002 | 0.000 | 0.007 |
| Ukraine | \>= 5 | Outpatient | 0.005 | 0.001 | 0.014 | 0.004 | 0.001 | 0.013 |
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
| Uzbekistan | \< 5 | Inpatient | 0.011 | 0.003 | 0.032 | 0.017 | 0.004 | 0.054 |
| Uzbekistan | \>= 5 | Inpatient | 0.021 | 0.006 | 0.060 | 0.032 | 0.007 | 0.102 |
| Uzbekistan | \< 5 | Outpatient | 0.002 | 0.000 | 0.007 | 0.004 | 0.001 | 0.013 |
| Uzbekistan | \>= 5 | Outpatient | 0.005 | 0.001 | 0.014 | 0.007 | 0.001 | 0.024 |
| Vanuatu | \< 5 | Inpatient | 0.007 | 0.001 | 0.019 | 0.006 | 0.001 | 0.018 |
| Vanuatu | \>= 5 | Inpatient | 0.014 | 0.002 | 0.037 | 0.012 | 0.002 | 0.036 |
| Vanuatu | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.001 | 0.000 | 0.004 |
| Vanuatu | \>= 5 | Outpatient | 0.003 | 0.000 | 0.009 | 0.003 | 0.000 | 0.008 |
| Venezuela | \< 5 | Inpatient | 0.010 | 0.002 | 0.026 | 0.007 | 0.001 | 0.021 |
| Venezuela | \>= 5 | Inpatient | 0.019 | 0.004 | 0.052 | 0.013 | 0.002 | 0.042 |
| Venezuela | \< 5 | Outpatient | 0.002 | 0.000 | 0.006 | 0.001 | 0.000 | 0.005 |
| Venezuela | \>= 5 | Outpatient | 0.004 | 0.001 | 0.012 | 0.003 | 0.000 | 0.009 |
| Vietnam | \< 5 | Inpatient | 0.007 | 0.001 | 0.019 | 0.006 | 0.001 | 0.018 |
| Vietnam | \>= 5 | Inpatient | 0.014 | 0.002 | 0.037 | 0.012 | 0.002 | 0.036 |
| Vietnam | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.001 | 0.000 | 0.004 |
| Vietnam | \>= 5 | Outpatient | 0.003 | 0.000 | 0.009 | 0.003 | 0.000 | 0.008 |
| Yemen, Rep. | \< 5 | Inpatient | 0.011 | 0.003 | 0.032 | 0.014 | 0.004 | 0.039 |
| Yemen, Rep. | \>= 5 | Inpatient | 0.021 | 0.005 | 0.060 | 0.026 | 0.006 | 0.073 |
| Yemen, Rep. | \< 5 | Outpatient | 0.002 | 0.000 | 0.007 | 0.003 | 0.001 | 0.009 |
| Yemen, Rep. | \>= 5 | Outpatient | 0.005 | 0.001 | 0.014 | 0.006 | 0.001 | 0.017 |
| Zambia | \< 5 | Inpatient | 0.008 | 0.003 | 0.019 | 0.012 | 0.003 | 0.032 |
| Zambia | \>= 5 | Inpatient | 0.016 | 0.005 | 0.039 | 0.024 | 0.006 | 0.063 |
| Zambia | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.003 | 0.001 | 0.007 |
| Zambia | \>= 5 | Outpatient | 0.003 | 0.001 | 0.009 | 0.005 | 0.001 | 0.015 |
| Zimbabwe | \< 5 | Inpatient | 0.008 | 0.003 | 0.017 | 0.012 | 0.004 | 0.030 |
| Zimbabwe | \>= 5 | Inpatient | 0.015 | 0.005 | 0.034 | 0.023 | 0.006 | 0.058 |
| Zimbabwe | \< 5 | Outpatient | 0.002 | 0.000 | 0.004 | 0.003 | 0.001 | 0.007 |
| Zimbabwe | \>= 5 | Outpatient | 0.003 | 0.001 | 0.008 | 0.005 | 0.001 | 0.014 |

Estimated `r params$Pathogen` proportion by country, 2010 vs 2020

# Session info

``` r
sessioninfo::session_info()
```

    ## Warning in system2("quarto", "-V", stdout = TRUE, env = paste0("TMPDIR=", : running command '"quarto"
    ## TMPDIR=C:/Users/fbbu6966/AppData/Local/Temp/Rtmpgxvndi/file2db4f5730 -V' had status 1

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
    ##  quarto   ERROR: Unknown command "TMPDIR=C:/Users/fbbu6966/AppData/Local/Temp/Rtmpgxvndi/file2db4f5730". Did you mean command "render"? @ C:\\Users\\fbbu6966\\AppData\\Local\\Programs\\RStudio\\RESOUR~1\\app\\bin\\quarto\\bin\\quarto.exe
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
