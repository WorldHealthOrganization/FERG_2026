Global proportion of Campylobacter (Diarrheal) • SYMPT+ROTA • Estimate
proportion
================
fbbu6966
2025-08-31

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
image <- paste0("03-estimate-rota_files/figure-gfm/imputation_map.png")
cat("![](", image ,")")
```

![](03-estimate-rota_files/figure-gfm/imputation_map.png)

``` r
fit_brms_reg_s <- readRDS(paste0(params$Dir, "/fit_brms_reg_s.rds"))
print(fit_brms_reg_s)
```

    ## Warning: There were 13 divergent transitions after warmup. Increasing adapt_delta above 0.95 may help. See
    ## http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup

    ##  Family: gaussian 
    ##   Links: mu = identity; sigma = identity 
    ## Formula: yi | se(sei) ~ 1 + YEAR + SYNDROMTYPE + AGE + REFERENCE + COVERAGE + STRAIN + (1 | REG2) + (1 | REG2:SUB2) + (1 | REG2:SUB2:COUNTRY) + (1 | REG2:SUB2:COUNTRY:ID) + (1 | REG2:SUB2:COUNTRY:ID:DTP_ID) 
    ##    Data: subset(es, as.integer(FLAG) == 1) (Number of observations: 328) 
    ##   Draws: 5 chains, each with iter = 5000; warmup = 3000; thin = 1;
    ##          total post-warmup draws = 10000
    ## 
    ## Multilevel Hyperparameters:
    ## ~REG2 (Number of levels: 6) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.22      0.19     0.01     0.72 1.00     6020     6486
    ## 
    ## ~REG2:SUB2 (Number of levels: 12) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.17      0.15     0.01     0.55 1.00     4254     3793
    ## 
    ## ~REG2:SUB2:COUNTRY (Number of levels: 47) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.43      0.18     0.07     0.77 1.01      995     1404
    ## 
    ## ~REG2:SUB2:COUNTRY:ID (Number of levels: 122) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.85      0.09     0.67     1.04 1.00     2309     3016
    ## 
    ## ~REG2:SUB2:COUNTRY:ID:DTP_ID (Number of levels: 328) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.48      0.04     0.41     0.56 1.00     3451     5877
    ## 
    ## Regression Coefficients:
    ##                        Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept                101.04     38.88    23.79   178.32 1.00     4210     5963
    ## YEAR                      -0.05      0.02    -0.09    -0.01 1.00     4205     6019
    ## SYNDROMTYPEOutpatient      0.14      0.12    -0.11     0.38 1.00     5357     5959
    ## AGEAgebelow5               0.47      0.15     0.18     0.77 1.00     5380     6382
    ## AGEMixedages               0.60      0.21     0.18     1.01 1.00     6608     7317
    ## REFERENCEOther            -1.95      0.23    -2.40    -1.49 1.00     3963     5757
    ## COVERAGE                   0.00      0.00    -0.01     0.01 1.00     5838     5891
    ## STRAINcoli                -1.18      1.06    -3.24     0.90 1.00     7222     7161
    ## STRAINjejuni              -0.53      0.32    -1.15     0.09 1.00     4685     5597
    ## STRAINUndifferentiated     0.12      0.24    -0.36     0.59 1.00     4575     5268
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
read_xlsx("endemic_countries.xlsx")%>%
  select(REG2, SUB2, ISO3, Country, edtf_diarrheal) %>% 
  rename(COUNTRY=ISO3, COUNTRY_LABEL = Country) %>%
  mutate(DISEASEFREE = edtf_diarrheal)
```

    ## New names:
    ## • `` -> `...33`

``` r
kable(
  caption = "High income countries, estimates not made with this approach",
  row.names = FALSE,
  subset(zero_cases, edtf_diarrheal==0)[, 4])
```

| COUNTRY_LABEL        |
|:---------------------|
| Andorra              |
| Antigua and Barbuda  |
| Australia            |
| Austria              |
| Bahamas, The         |
| Bahrain              |
| Barbados             |
| Belgium              |
| Brunei Darussalam    |
| Canada               |
| Chile                |
| Cook Islands         |
| Croatia              |
| Cyprus               |
| Czech Republic       |
| Denmark              |
| Estonia              |
| Finland              |
| France               |
| Germany              |
| Greece               |
| Guyana               |
| Hungary              |
| Iceland              |
| Ireland              |
| Israel               |
| Italy                |
| Japan                |
| Korea, Rep.          |
| Kuwait               |
| Latvia               |
| Lithuania            |
| Luxembourg           |
| Malta                |
| Monaco               |
| Nauru                |
| Netherlands          |
| New Zealand          |
| Niue                 |
| Norway               |
| Oman                 |
| Panama               |
| Poland               |
| Portugal             |
| Qatar                |
| Romania              |
| San Marino           |
| Saudi Arabia         |
| Singapore            |
| Slovak Republic      |
| Slovenia             |
| Spain                |
| St. Kitts and Nevis  |
| Sweden               |
| Switzerland          |
| Trinidad and Tobago  |
| United Arab Emirates |
| United Kingdom       |
| United States        |
| Uruguay              |

High income countries, estimates not made with this approach

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
    SYNDROMTYPE = 1:2,
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
         levels = 1:2, 
         labels = c("Inpatient", "Outpatient"))
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
b0_ad_si <- draws_fit$b_Intercept
b0_ad_so <- b0_ad_si + draws_fit$b_SYNDROMTYPEOutpatient
b0_ch_si <- b0_ad_si + draws_fit$b_AGEAgebelow5
b0_ch_so <- b0_ch_si + draws_fit$b_SYNDROMTYPEOutpatient
b0 <-
  function(age, type) {
    b0 <-
      paste(
        "b0",
        c("ad", "ch")[age],
        c("si", "so")[type],
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

    ## [1] 4136

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

    ## [1] 5280

``` r
for (x in id) {
  sim_all$SIM[x] <- -Inf
}

## data in subregion
id <- which(as.integer(sim_all$ESTIMATES) == 3); length(id)
```

    ## [1] 7304

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

    ## [1] 352

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

# sim_all_as <- subset(sim_all, SYNDROMTYPE == "Asymptomatic")
# sim_all_si <- subset(sim_all, SYNDROMTYPE == "Inpatient")
# sim_all_so <- subset(sim_all, SYNDROMTYPE == "Outpatient")
# 
# sim_all_si$PROP <- mapply(`-`, sim_all_si$PROP, sim_all_as$PROP)
# sim_all_so$PROP <- mapply(`-`, sim_all_so$PROP, sim_all_as$PROP)
# 
# sim_all <- rbind(sim_all_si, sim_all_so)

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
    ##  $ YEAR         : int  2000 2000 2000 2000 2001 2001 2001 2001 2002 2002 ...
    ##  $ REG2         : chr  "EMR" "EMR" "EMR" "EMR" ...
    ##  $ SUB2         : chr  "EMRD" "EMRD" "EMRD" "EMRD" ...
    ##  $ COUNTRY      : chr  "AFG" "AFG" "AFG" "AFG" ...
    ##  $ SYNDROMTYPE  : Factor w/ 2 levels "Inpatient","Outpatient": 1 2 2 1 2 1 2 1 1 1 ...
    ##  $ AGE          : Factor w/ 2 levels "Age above or equal 5",..: 1 1 2 2 2 1 1 2 2 1 ...
    ##  $ VAL_MEAN     : num  0.198 0.22 0.308 0.28 0.297 ...
    ##  $ VAL_LWR      : num  0.0816 0.0938 0.1472 0.1298 0.1425 ...
    ##  $ VAL_UPR      : num  0.367 0.398 0.505 0.474 0.488 ...
    ##  $ LOCATION     : chr  "Country" "Country" "Country" "Country" ...
    ##  $ LOCATION_NAME: chr  "Afghanistan" "Afghanistan" "Afghanistan" "Afghanistan" ...
    ##  $ METRIC       : chr  "Proportion" "Proportion" "Proportion" "Proportion" ...

``` r
## compile all
all_est <- all_cnt_prop
str(all_est)
```

    ## 'data.frame':    17072 obs. of  12 variables:
    ##  $ YEAR         : int  2000 2000 2000 2000 2001 2001 2001 2001 2002 2002 ...
    ##  $ REG2         : chr  "EMR" "EMR" "EMR" "EMR" ...
    ##  $ SUB2         : chr  "EMRD" "EMRD" "EMRD" "EMRD" ...
    ##  $ COUNTRY      : chr  "AFG" "AFG" "AFG" "AFG" ...
    ##  $ SYNDROMTYPE  : Factor w/ 2 levels "Inpatient","Outpatient": 1 2 2 1 2 1 2 1 1 1 ...
    ##  $ AGE          : Factor w/ 2 levels "Age above or equal 5",..: 1 1 2 2 2 1 1 2 2 1 ...
    ##  $ VAL_MEAN     : num  0.198 0.22 0.308 0.28 0.297 ...
    ##  $ VAL_LWR      : num  0.0816 0.0938 0.1472 0.1298 0.1425 ...
    ##  $ VAL_UPR      : num  0.367 0.398 0.505 0.474 0.488 ...
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
image <- paste0("03-estimate-rota_files/figure-gfm/r_world1.png")
cat("![](",image,")")
```

![](03-estimate-rota_files/figure-gfm/r_world1.png)

``` r
png(paste0(params$PlotDir, "/r_world2.png"), width=960, height=480)
plot_world(subset(all_cnt_prop, YEAR == 2010 & AGE == "Age above or equal 5" & SYNDROMTYPE == "Outpatient"),  
           "COUNTRY", "VAL_MEAN", legend.title = "Proportion", diseasefree = zero_cases, text.width = 23)
```

\[1\] 0.00 0.05 0.10 0.15 0.20 0.25

``` r
dev.off()
```

png 2

``` r
setwd(params$Dir)
image <- paste0("03-estimate-rota_files/figure-gfm/r_world2.png")
cat("![](",image,")")
```

![](03-estimate-rota_files/figure-gfm/r_world2.png)

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
image <- paste0("03-estimate-rota_files/figure-gfm/r_world3.png")
cat("![](",image,")")
```

![](03-estimate-rota_files/figure-gfm/r_world3.png)

``` r
png(paste0(params$PlotDir, "/r_world4.png"), width=960, height=480)
plot_world(subset(all_cnt_prop, YEAR == 2010 & AGE == "Age below 5" & SYNDROMTYPE == "Outpatient"), 
           "COUNTRY", "VAL_MEAN", legend.title = "Proportion", diseasefree = zero_cases, text.width = 23)
```

\[1\] 0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35

``` r
dev.off()
```

png 2

``` r
setwd(params$Dir)
image <- paste0("03-estimate-rota_files/figure-gfm/r_world4.png")
cat("![](",image,")")
```

![](03-estimate-rota_files/figure-gfm/r_world4.png)

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
image <- paste0("03-estimate-rota_files/figure-gfm/r_world5.png")
cat("![](",image,")")
```

![](03-estimate-rota_files/figure-gfm/r_world5.png)

``` r
png(paste0(params$PlotDir, "/r_world6.png"), width=960, height=480)
plot_world(subset(all_cnt_prop, YEAR == 2020 & AGE == "Age above or equal 5" & SYNDROMTYPE == "Outpatient"),  
           "COUNTRY", "VAL_MEAN", legend.title = "Proportion", diseasefree = zero_cases, text.width = 23)
```

\[1\] 0.00 0.05 0.10 0.15 0.20

``` r
dev.off()
```

png 2

``` r
setwd(params$Dir)
image <- paste0("03-estimate-rota_files/figure-gfm/r_world6.png")
cat("![](",image,")")
```

![](03-estimate-rota_files/figure-gfm/r_world6.png)

``` r
png(paste0(params$PlotDir, "/r_world7.png"), width=960, height=480)
plot_world(subset(all_cnt_prop, YEAR == 2020 & AGE == "Age below 5" & SYNDROMTYPE == "Inpatient"),
           "COUNTRY", "VAL_MEAN", legend.title = "Proportion", diseasefree = zero_cases, text.width = 23)
```

\[1\] 0.00 0.05 0.10 0.15 0.20 0.25

``` r
dev.off()
```

png 2

``` r
setwd(params$Dir)
image <- paste0("03-estimate-rota_files/figure-gfm/r_world7.png")
cat("![](",image,")")
```

![](03-estimate-rota_files/figure-gfm/r_world7.png)

``` r
png(paste0(params$PlotDir, "/r_world8.png"), width=960, height=480)
plot_world(subset(all_cnt_prop, YEAR == 2020 & AGE == "Age below 5" & SYNDROMTYPE == "Outpatient"),  
           "COUNTRY", "VAL_MEAN", legend.title = "Proportion", diseasefree = zero_cases, text.width = 23)
```

\[1\] 0.00 0.05 0.10 0.15 0.20 0.25 0.30

``` r
dev.off()
```

png 2

``` r
setwd(params$Dir)
image <- paste0("03-estimate-rota_files/figure-gfm/r_world8.png")
cat("![](",image,")")
```

![](03-estimate-rota_files/figure-gfm/r_world8.png)

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
image <- paste0("03-estimate-rota_files/figure-gfm/r_CI1.png")
cat("![](",image,")")
```

![](03-estimate-rota_files/figure-gfm/r_CI1.png)

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
image <- paste0("03-estimate-rota_files/figure-gfm/r_CI2.png")
cat("![](",image,")")
```

![](03-estimate-rota_files/figure-gfm/r_CI2.png)

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
image <- paste0("03-estimate-rota_files/figure-gfm/r_CI3.png")
cat("![](",image,")")
```

![](03-estimate-rota_files/figure-gfm/r_CI3.png)

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
image <- paste0("03-estimate-rota_files/figure-gfm/r_CI4.png")
cat("![](",image,")")
```

![](03-estimate-rota_files/figure-gfm/r_CI4.png)

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
image <- paste0("03-estimate-rota_files/figure-gfm/r_CI5.png")
cat("![](",image,")")
```

![](03-estimate-rota_files/figure-gfm/r_CI5.png)

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
image <- paste0("03-estimate-rota_files/figure-gfm/r_CI6.png")
cat("![](",image,")")
```

![](03-estimate-rota_files/figure-gfm/r_CI6.png)

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
| Afghanistan | \< 5 | Inpatient | 0.189 | 0.090 | 0.326 | 0.135 | 0.061 | 0.243 |
| Afghanistan | \>= 5 | Inpatient | 0.129 | 0.057 | 0.237 | 0.090 | 0.038 | 0.172 |
| Afghanistan | \< 5 | Outpatient | 0.211 | 0.103 | 0.353 | 0.152 | 0.071 | 0.267 |
| Afghanistan | \>= 5 | Outpatient | 0.144 | 0.066 | 0.262 | 0.101 | 0.044 | 0.190 |
| Albania | \< 5 | Inpatient | 0.197 | 0.073 | 0.392 | 0.152 | 0.049 | 0.335 |
| Albania | \>= 5 | Inpatient | 0.135 | 0.047 | 0.288 | 0.102 | 0.030 | 0.241 |
| Albania | \< 5 | Outpatient | 0.219 | 0.085 | 0.427 | 0.170 | 0.056 | 0.364 |
| Albania | \>= 5 | Outpatient | 0.151 | 0.054 | 0.315 | 0.115 | 0.035 | 0.266 |
| Algeria | \< 5 | Inpatient | 0.198 | 0.118 | 0.301 | 0.130 | 0.069 | 0.217 |
| Algeria | \>= 5 | Inpatient | 0.134 | 0.074 | 0.223 | 0.086 | 0.043 | 0.155 |
| Algeria | \< 5 | Outpatient | 0.220 | 0.134 | 0.328 | 0.147 | 0.079 | 0.238 |
| Algeria | \>= 5 | Outpatient | 0.151 | 0.085 | 0.244 | 0.097 | 0.050 | 0.170 |
| Andorra | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Andorra | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Andorra | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Andorra | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Angola | \< 5 | Inpatient | 0.183 | 0.057 | 0.354 | 0.128 | 0.037 | 0.267 |
| Angola | \>= 5 | Inpatient | 0.125 | 0.035 | 0.266 | 0.085 | 0.023 | 0.189 |
| Angola | \< 5 | Outpatient | 0.204 | 0.065 | 0.390 | 0.144 | 0.042 | 0.294 |
| Angola | \>= 5 | Outpatient | 0.140 | 0.040 | 0.291 | 0.096 | 0.026 | 0.208 |
| Antigua and Barbuda | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Antigua and Barbuda | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Antigua and Barbuda | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Antigua and Barbuda | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Argentina | \< 5 | Inpatient | 0.262 | 0.110 | 0.502 | 0.197 | 0.076 | 0.418 |
| Argentina | \>= 5 | Inpatient | 0.185 | 0.070 | 0.390 | 0.136 | 0.047 | 0.312 |
| Argentina | \< 5 | Outpatient | 0.288 | 0.126 | 0.537 | 0.219 | 0.088 | 0.447 |
| Argentina | \>= 5 | Outpatient | 0.206 | 0.081 | 0.419 | 0.152 | 0.055 | 0.335 |
| Armenia | \< 5 | Inpatient | 0.184 | 0.092 | 0.312 | 0.140 | 0.060 | 0.262 |
| Armenia | \>= 5 | Inpatient | 0.125 | 0.058 | 0.224 | 0.093 | 0.037 | 0.187 |
| Armenia | \< 5 | Outpatient | 0.206 | 0.106 | 0.338 | 0.156 | 0.069 | 0.288 |
| Armenia | \>= 5 | Outpatient | 0.140 | 0.067 | 0.246 | 0.105 | 0.043 | 0.207 |
| Australia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Australia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Australia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Australia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Austria | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Austria | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Austria | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Austria | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Azerbaijan | \< 5 | Inpatient | 0.184 | 0.092 | 0.312 | 0.121 | 0.055 | 0.223 |
| Azerbaijan | \>= 5 | Inpatient | 0.125 | 0.058 | 0.224 | 0.080 | 0.034 | 0.155 |
| Azerbaijan | \< 5 | Outpatient | 0.206 | 0.106 | 0.338 | 0.136 | 0.063 | 0.244 |
| Azerbaijan | \>= 5 | Outpatient | 0.140 | 0.067 | 0.246 | 0.090 | 0.040 | 0.172 |
| Bahamas, The | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahamas, The | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahamas, The | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahamas, The | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahrain | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahrain | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahrain | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahrain | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bangladesh | \< 5 | Inpatient | 0.225 | 0.125 | 0.363 | 0.150 | 0.074 | 0.269 |
| Bangladesh | \>= 5 | Inpatient | 0.155 | 0.078 | 0.272 | 0.101 | 0.046 | 0.193 |
| Bangladesh | \< 5 | Outpatient | 0.249 | 0.141 | 0.391 | 0.169 | 0.084 | 0.294 |
| Bangladesh | \>= 5 | Outpatient | 0.174 | 0.089 | 0.294 | 0.114 | 0.052 | 0.211 |
| Barbados | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Barbados | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Barbados | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Barbados | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belarus | \< 5 | Inpatient | 0.184 | 0.092 | 0.312 | 0.121 | 0.055 | 0.223 |
| Belarus | \>= 5 | Inpatient | 0.125 | 0.058 | 0.224 | 0.080 | 0.034 | 0.155 |
| Belarus | \< 5 | Outpatient | 0.206 | 0.106 | 0.338 | 0.136 | 0.063 | 0.244 |
| Belarus | \>= 5 | Outpatient | 0.140 | 0.067 | 0.246 | 0.090 | 0.040 | 0.172 |
| Belgium | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belgium | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belgium | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belgium | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belize | \< 5 | Inpatient | 0.197 | 0.106 | 0.323 | 0.130 | 0.061 | 0.240 |
| Belize | \>= 5 | Inpatient | 0.134 | 0.066 | 0.237 | 0.086 | 0.038 | 0.170 |
| Belize | \< 5 | Outpatient | 0.219 | 0.121 | 0.351 | 0.147 | 0.070 | 0.264 |
| Belize | \>= 5 | Outpatient | 0.150 | 0.077 | 0.260 | 0.098 | 0.045 | 0.185 |
| Benin | \< 5 | Inpatient | 0.198 | 0.118 | 0.301 | 0.142 | 0.080 | 0.230 |
| Benin | \>= 5 | Inpatient | 0.134 | 0.074 | 0.223 | 0.095 | 0.049 | 0.163 |
| Benin | \< 5 | Outpatient | 0.220 | 0.134 | 0.328 | 0.160 | 0.092 | 0.251 |
| Benin | \>= 5 | Outpatient | 0.151 | 0.085 | 0.244 | 0.107 | 0.056 | 0.182 |
| Bhutan | \< 5 | Inpatient | 0.194 | 0.108 | 0.311 | 0.128 | 0.064 | 0.224 |
| Bhutan | \>= 5 | Inpatient | 0.132 | 0.067 | 0.226 | 0.085 | 0.039 | 0.158 |
| Bhutan | \< 5 | Outpatient | 0.216 | 0.123 | 0.338 | 0.144 | 0.073 | 0.249 |
| Bhutan | \>= 5 | Outpatient | 0.148 | 0.079 | 0.249 | 0.096 | 0.046 | 0.175 |
| Bolivia | \< 5 | Inpatient | 0.207 | 0.094 | 0.359 | 0.134 | 0.058 | 0.240 |
| Bolivia | \>= 5 | Inpatient | 0.142 | 0.058 | 0.266 | 0.089 | 0.036 | 0.170 |
| Bolivia | \< 5 | Outpatient | 0.230 | 0.107 | 0.388 | 0.150 | 0.068 | 0.262 |
| Bolivia | \>= 5 | Outpatient | 0.159 | 0.066 | 0.290 | 0.100 | 0.042 | 0.188 |
| Bosnia and Herzegovina | \< 5 | Inpatient | 0.184 | 0.092 | 0.312 | 0.121 | 0.055 | 0.223 |
| Bosnia and Herzegovina | \>= 5 | Inpatient | 0.125 | 0.058 | 0.224 | 0.080 | 0.034 | 0.155 |
| Bosnia and Herzegovina | \< 5 | Outpatient | 0.206 | 0.106 | 0.338 | 0.136 | 0.063 | 0.244 |
| Bosnia and Herzegovina | \>= 5 | Outpatient | 0.140 | 0.067 | 0.246 | 0.090 | 0.040 | 0.172 |
| Botswana | \< 5 | Inpatient | 0.189 | 0.103 | 0.303 | 0.141 | 0.071 | 0.241 |
| Botswana | \>= 5 | Inpatient | 0.128 | 0.064 | 0.221 | 0.094 | 0.043 | 0.174 |
| Botswana | \< 5 | Outpatient | 0.211 | 0.117 | 0.331 | 0.158 | 0.081 | 0.265 |
| Botswana | \>= 5 | Outpatient | 0.144 | 0.075 | 0.241 | 0.106 | 0.050 | 0.192 |
| Brazil | \< 5 | Inpatient | 0.181 | 0.073 | 0.337 | 0.117 | 0.044 | 0.235 |
| Brazil | \>= 5 | Inpatient | 0.123 | 0.044 | 0.250 | 0.078 | 0.027 | 0.165 |
| Brazil | \< 5 | Outpatient | 0.202 | 0.083 | 0.370 | 0.132 | 0.050 | 0.258 |
| Brazil | \>= 5 | Outpatient | 0.139 | 0.051 | 0.275 | 0.088 | 0.031 | 0.184 |
| Brunei Darussalam | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Brunei Darussalam | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Brunei Darussalam | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Brunei Darussalam | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bulgaria | \< 5 | Inpatient | 0.185 | 0.067 | 0.366 | 0.130 | 0.044 | 0.277 |
| Bulgaria | \>= 5 | Inpatient | 0.126 | 0.042 | 0.262 | 0.086 | 0.027 | 0.191 |
| Bulgaria | \< 5 | Outpatient | 0.206 | 0.078 | 0.392 | 0.146 | 0.051 | 0.302 |
| Bulgaria | \>= 5 | Outpatient | 0.141 | 0.049 | 0.288 | 0.097 | 0.032 | 0.213 |
| Burkina Faso | \< 5 | Inpatient | 0.189 | 0.067 | 0.364 | 0.142 | 0.045 | 0.300 |
| Burkina Faso | \>= 5 | Inpatient | 0.129 | 0.042 | 0.271 | 0.096 | 0.028 | 0.217 |
| Burkina Faso | \< 5 | Outpatient | 0.210 | 0.075 | 0.398 | 0.160 | 0.051 | 0.331 |
| Burkina Faso | \>= 5 | Outpatient | 0.145 | 0.047 | 0.299 | 0.108 | 0.031 | 0.242 |
| Burundi | \< 5 | Inpatient | 0.200 | 0.118 | 0.312 | 0.153 | 0.075 | 0.268 |
| Burundi | \>= 5 | Inpatient | 0.136 | 0.074 | 0.228 | 0.102 | 0.046 | 0.193 |
| Burundi | \< 5 | Outpatient | 0.223 | 0.134 | 0.339 | 0.171 | 0.086 | 0.295 |
| Burundi | \>= 5 | Outpatient | 0.153 | 0.085 | 0.250 | 0.115 | 0.053 | 0.212 |
| Cabo Verde | \< 5 | Inpatient | 0.189 | 0.071 | 0.356 | 0.125 | 0.043 | 0.256 |
| Cabo Verde | \>= 5 | Inpatient | 0.129 | 0.045 | 0.262 | 0.083 | 0.026 | 0.181 |
| Cabo Verde | \< 5 | Outpatient | 0.210 | 0.081 | 0.388 | 0.140 | 0.049 | 0.281 |
| Cabo Verde | \>= 5 | Outpatient | 0.144 | 0.051 | 0.291 | 0.093 | 0.030 | 0.201 |
| Cambodia | \< 5 | Inpatient | 0.167 | 0.067 | 0.318 | 0.109 | 0.040 | 0.227 |
| Cambodia | \>= 5 | Inpatient | 0.112 | 0.042 | 0.231 | 0.072 | 0.025 | 0.158 |
| Cambodia | \< 5 | Outpatient | 0.186 | 0.077 | 0.349 | 0.123 | 0.046 | 0.249 |
| Cambodia | \>= 5 | Outpatient | 0.126 | 0.049 | 0.255 | 0.081 | 0.029 | 0.175 |
| Cameroon | \< 5 | Inpatient | 0.239 | 0.107 | 0.443 | 0.176 | 0.076 | 0.350 |
| Cameroon | \>= 5 | Inpatient | 0.167 | 0.068 | 0.339 | 0.119 | 0.048 | 0.258 |
| Cameroon | \< 5 | Outpatient | 0.264 | 0.123 | 0.479 | 0.196 | 0.087 | 0.382 |
| Cameroon | \>= 5 | Outpatient | 0.186 | 0.079 | 0.370 | 0.134 | 0.055 | 0.283 |
| Canada | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Canada | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Canada | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Canada | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Central African Republic | \< 5 | Inpatient | 0.189 | 0.070 | 0.360 | 0.125 | 0.043 | 0.262 |
| Central African Republic | \>= 5 | Inpatient | 0.129 | 0.044 | 0.268 | 0.083 | 0.027 | 0.184 |
| Central African Republic | \< 5 | Outpatient | 0.210 | 0.080 | 0.396 | 0.141 | 0.049 | 0.291 |
| Central African Republic | \>= 5 | Outpatient | 0.145 | 0.050 | 0.295 | 0.094 | 0.030 | 0.207 |
| Chad | \< 5 | Inpatient | 0.200 | 0.118 | 0.312 | 0.132 | 0.069 | 0.223 |
| Chad | \>= 5 | Inpatient | 0.136 | 0.074 | 0.228 | 0.087 | 0.043 | 0.158 |
| Chad | \< 5 | Outpatient | 0.223 | 0.134 | 0.339 | 0.149 | 0.079 | 0.247 |
| Chad | \>= 5 | Outpatient | 0.153 | 0.085 | 0.250 | 0.099 | 0.050 | 0.175 |
| Chile | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Chile | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Chile | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Chile | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| China | \< 5 | Inpatient | 0.152 | 0.077 | 0.259 | 0.097 | 0.050 | 0.168 |
| China | \>= 5 | Inpatient | 0.102 | 0.047 | 0.184 | 0.063 | 0.030 | 0.115 |
| China | \< 5 | Outpatient | 0.171 | 0.087 | 0.286 | 0.110 | 0.057 | 0.186 |
| China | \>= 5 | Outpatient | 0.115 | 0.055 | 0.203 | 0.072 | 0.035 | 0.128 |
| Colombia | \< 5 | Inpatient | 0.213 | 0.091 | 0.384 | 0.144 | 0.057 | 0.276 |
| Colombia | \>= 5 | Inpatient | 0.147 | 0.056 | 0.289 | 0.096 | 0.036 | 0.198 |
| Colombia | \< 5 | Outpatient | 0.236 | 0.104 | 0.413 | 0.161 | 0.066 | 0.303 |
| Colombia | \>= 5 | Outpatient | 0.164 | 0.064 | 0.315 | 0.109 | 0.041 | 0.220 |
| Comoros | \< 5 | Inpatient | 0.198 | 0.118 | 0.301 | 0.130 | 0.069 | 0.217 |
| Comoros | \>= 5 | Inpatient | 0.134 | 0.074 | 0.223 | 0.086 | 0.043 | 0.155 |
| Comoros | \< 5 | Outpatient | 0.220 | 0.134 | 0.328 | 0.147 | 0.079 | 0.238 |
| Comoros | \>= 5 | Outpatient | 0.151 | 0.085 | 0.244 | 0.097 | 0.050 | 0.170 |
| Congo, Dem. Rep. | \< 5 | Inpatient | 0.200 | 0.118 | 0.312 | 0.137 | 0.076 | 0.224 |
| Congo, Dem. Rep. | \>= 5 | Inpatient | 0.136 | 0.074 | 0.228 | 0.091 | 0.046 | 0.158 |
| Congo, Dem. Rep. | \< 5 | Outpatient | 0.223 | 0.134 | 0.339 | 0.154 | 0.086 | 0.249 |
| Congo, Dem. Rep. | \>= 5 | Outpatient | 0.153 | 0.085 | 0.250 | 0.102 | 0.054 | 0.175 |
| Congo, Rep. | \< 5 | Inpatient | 0.198 | 0.118 | 0.301 | 0.141 | 0.080 | 0.228 |
| Congo, Rep. | \>= 5 | Inpatient | 0.134 | 0.074 | 0.223 | 0.094 | 0.049 | 0.161 |
| Congo, Rep. | \< 5 | Outpatient | 0.220 | 0.134 | 0.328 | 0.159 | 0.092 | 0.250 |
| Congo, Rep. | \>= 5 | Outpatient | 0.151 | 0.085 | 0.244 | 0.106 | 0.056 | 0.180 |
| Cook Islands | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Cook Islands | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Cook Islands | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Cook Islands | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Costa Rica | \< 5 | Inpatient | 0.197 | 0.106 | 0.323 | 0.149 | 0.074 | 0.263 |
| Costa Rica | \>= 5 | Inpatient | 0.134 | 0.066 | 0.237 | 0.100 | 0.045 | 0.188 |
| Costa Rica | \< 5 | Outpatient | 0.219 | 0.121 | 0.351 | 0.167 | 0.085 | 0.291 |
| Costa Rica | \>= 5 | Outpatient | 0.150 | 0.077 | 0.260 | 0.113 | 0.052 | 0.208 |
| Côte d’Ivoire | \< 5 | Inpatient | 0.198 | 0.118 | 0.301 | 0.144 | 0.080 | 0.234 |
| Côte d’Ivoire | \>= 5 | Inpatient | 0.134 | 0.074 | 0.223 | 0.096 | 0.049 | 0.167 |
| Côte d’Ivoire | \< 5 | Outpatient | 0.220 | 0.134 | 0.328 | 0.161 | 0.092 | 0.255 |
| Côte d’Ivoire | \>= 5 | Outpatient | 0.151 | 0.085 | 0.244 | 0.108 | 0.056 | 0.186 |
| Croatia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Croatia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Croatia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Croatia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Cuba | \< 5 | Inpatient | 0.159 | 0.046 | 0.323 | 0.105 | 0.026 | 0.239 |
| Cuba | \>= 5 | Inpatient | 0.107 | 0.028 | 0.236 | 0.069 | 0.016 | 0.170 |
| Cuba | \< 5 | Outpatient | 0.177 | 0.052 | 0.353 | 0.118 | 0.030 | 0.262 |
| Cuba | \>= 5 | Outpatient | 0.121 | 0.033 | 0.257 | 0.078 | 0.018 | 0.187 |
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
| Djibouti | \< 5 | Inpatient | 0.192 | 0.104 | 0.315 | 0.140 | 0.069 | 0.243 |
| Djibouti | \>= 5 | Inpatient | 0.130 | 0.065 | 0.228 | 0.093 | 0.042 | 0.174 |
| Djibouti | \< 5 | Outpatient | 0.213 | 0.119 | 0.341 | 0.157 | 0.079 | 0.268 |
| Djibouti | \>= 5 | Outpatient | 0.146 | 0.076 | 0.249 | 0.105 | 0.049 | 0.192 |
| Dominica | \< 5 | Inpatient | 0.197 | 0.106 | 0.323 | 0.130 | 0.061 | 0.240 |
| Dominica | \>= 5 | Inpatient | 0.134 | 0.066 | 0.237 | 0.086 | 0.038 | 0.170 |
| Dominica | \< 5 | Outpatient | 0.219 | 0.121 | 0.351 | 0.147 | 0.070 | 0.264 |
| Dominica | \>= 5 | Outpatient | 0.150 | 0.077 | 0.260 | 0.098 | 0.045 | 0.185 |
| Dominican Republic | \< 5 | Inpatient | 0.197 | 0.106 | 0.323 | 0.144 | 0.074 | 0.248 |
| Dominican Republic | \>= 5 | Inpatient | 0.134 | 0.066 | 0.237 | 0.096 | 0.045 | 0.177 |
| Dominican Republic | \< 5 | Outpatient | 0.219 | 0.121 | 0.351 | 0.162 | 0.085 | 0.273 |
| Dominican Republic | \>= 5 | Outpatient | 0.150 | 0.077 | 0.260 | 0.109 | 0.053 | 0.194 |
| Ecuador | \< 5 | Inpatient | 0.224 | 0.114 | 0.376 | 0.144 | 0.074 | 0.246 |
| Ecuador | \>= 5 | Inpatient | 0.155 | 0.070 | 0.280 | 0.096 | 0.045 | 0.176 |
| Ecuador | \< 5 | Outpatient | 0.248 | 0.130 | 0.406 | 0.161 | 0.085 | 0.272 |
| Ecuador | \>= 5 | Outpatient | 0.173 | 0.082 | 0.306 | 0.108 | 0.053 | 0.194 |
| Egypt, Arab Rep. | \< 5 | Inpatient | 0.152 | 0.063 | 0.284 | 0.099 | 0.038 | 0.200 |
| Egypt, Arab Rep. | \>= 5 | Inpatient | 0.102 | 0.038 | 0.206 | 0.065 | 0.023 | 0.141 |
| Egypt, Arab Rep. | \< 5 | Outpatient | 0.171 | 0.072 | 0.311 | 0.112 | 0.042 | 0.220 |
| Egypt, Arab Rep. | \>= 5 | Outpatient | 0.115 | 0.044 | 0.227 | 0.074 | 0.026 | 0.156 |
| El Salvador | \< 5 | Inpatient | 0.223 | 0.115 | 0.369 | 0.144 | 0.074 | 0.248 |
| El Salvador | \>= 5 | Inpatient | 0.154 | 0.071 | 0.276 | 0.096 | 0.045 | 0.177 |
| El Salvador | \< 5 | Outpatient | 0.247 | 0.132 | 0.401 | 0.162 | 0.085 | 0.273 |
| El Salvador | \>= 5 | Outpatient | 0.172 | 0.083 | 0.301 | 0.109 | 0.053 | 0.194 |
| Equatorial Guinea | \< 5 | Inpatient | 0.189 | 0.103 | 0.303 | 0.124 | 0.060 | 0.217 |
| Equatorial Guinea | \>= 5 | Inpatient | 0.128 | 0.064 | 0.221 | 0.082 | 0.038 | 0.154 |
| Equatorial Guinea | \< 5 | Outpatient | 0.211 | 0.117 | 0.331 | 0.140 | 0.069 | 0.239 |
| Equatorial Guinea | \>= 5 | Outpatient | 0.144 | 0.075 | 0.241 | 0.093 | 0.044 | 0.169 |
| Eritrea | \< 5 | Inpatient | 0.200 | 0.118 | 0.312 | 0.153 | 0.075 | 0.268 |
| Eritrea | \>= 5 | Inpatient | 0.136 | 0.074 | 0.228 | 0.102 | 0.046 | 0.193 |
| Eritrea | \< 5 | Outpatient | 0.223 | 0.134 | 0.339 | 0.171 | 0.086 | 0.295 |
| Eritrea | \>= 5 | Outpatient | 0.153 | 0.085 | 0.250 | 0.115 | 0.053 | 0.212 |
| Estonia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Estonia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Estonia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Estonia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Eswatini | \< 5 | Inpatient | 0.198 | 0.118 | 0.301 | 0.142 | 0.080 | 0.230 |
| Eswatini | \>= 5 | Inpatient | 0.134 | 0.074 | 0.223 | 0.095 | 0.049 | 0.163 |
| Eswatini | \< 5 | Outpatient | 0.220 | 0.134 | 0.328 | 0.160 | 0.092 | 0.251 |
| Eswatini | \>= 5 | Outpatient | 0.151 | 0.085 | 0.244 | 0.107 | 0.056 | 0.182 |
| Ethiopia | \< 5 | Inpatient | 0.252 | 0.111 | 0.478 | 0.190 | 0.074 | 0.404 |
| Ethiopia | \>= 5 | Inpatient | 0.177 | 0.070 | 0.366 | 0.130 | 0.046 | 0.301 |
| Ethiopia | \< 5 | Outpatient | 0.278 | 0.126 | 0.511 | 0.211 | 0.084 | 0.430 |
| Ethiopia | \>= 5 | Outpatient | 0.197 | 0.080 | 0.398 | 0.146 | 0.053 | 0.326 |
| Fiji | \< 5 | Inpatient | 0.159 | 0.073 | 0.265 | 0.121 | 0.047 | 0.230 |
| Fiji | \>= 5 | Inpatient | 0.107 | 0.046 | 0.190 | 0.080 | 0.029 | 0.163 |
| Fiji | \< 5 | Outpatient | 0.178 | 0.084 | 0.292 | 0.136 | 0.054 | 0.256 |
| Fiji | \>= 5 | Outpatient | 0.120 | 0.053 | 0.211 | 0.090 | 0.034 | 0.181 |
| Finland | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Finland | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Finland | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Finland | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| France | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| France | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| France | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| France | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Gabon | \< 5 | Inpatient | 0.166 | 0.056 | 0.318 | 0.108 | 0.034 | 0.226 |
| Gabon | \>= 5 | Inpatient | 0.112 | 0.036 | 0.233 | 0.071 | 0.021 | 0.159 |
| Gabon | \< 5 | Outpatient | 0.185 | 0.064 | 0.347 | 0.122 | 0.039 | 0.249 |
| Gabon | \>= 5 | Outpatient | 0.126 | 0.041 | 0.254 | 0.081 | 0.024 | 0.176 |
| Gambia, The | \< 5 | Inpatient | 0.221 | 0.103 | 0.396 | 0.168 | 0.068 | 0.332 |
| Gambia, The | \>= 5 | Inpatient | 0.153 | 0.065 | 0.295 | 0.113 | 0.042 | 0.239 |
| Gambia, The | \< 5 | Outpatient | 0.245 | 0.116 | 0.426 | 0.187 | 0.077 | 0.361 |
| Gambia, The | \>= 5 | Outpatient | 0.171 | 0.074 | 0.322 | 0.128 | 0.048 | 0.264 |
| Georgia | \< 5 | Inpatient | 0.184 | 0.092 | 0.312 | 0.135 | 0.061 | 0.248 |
| Georgia | \>= 5 | Inpatient | 0.125 | 0.058 | 0.224 | 0.090 | 0.038 | 0.175 |
| Georgia | \< 5 | Outpatient | 0.206 | 0.106 | 0.338 | 0.152 | 0.070 | 0.272 |
| Georgia | \>= 5 | Outpatient | 0.140 | 0.067 | 0.246 | 0.102 | 0.044 | 0.194 |
| Germany | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Germany | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Germany | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Germany | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Ghana | \< 5 | Inpatient | 0.210 | 0.094 | 0.377 | 0.159 | 0.065 | 0.312 |
| Ghana | \>= 5 | Inpatient | 0.144 | 0.059 | 0.281 | 0.107 | 0.040 | 0.230 |
| Ghana | \< 5 | Outpatient | 0.233 | 0.107 | 0.408 | 0.178 | 0.076 | 0.337 |
| Ghana | \>= 5 | Outpatient | 0.162 | 0.068 | 0.306 | 0.120 | 0.047 | 0.247 |
| Greece | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Greece | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Greece | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Greece | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Grenada | \< 5 | Inpatient | 0.197 | 0.106 | 0.323 | 0.130 | 0.061 | 0.240 |
| Grenada | \>= 5 | Inpatient | 0.134 | 0.066 | 0.237 | 0.086 | 0.038 | 0.170 |
| Grenada | \< 5 | Outpatient | 0.219 | 0.121 | 0.351 | 0.147 | 0.070 | 0.264 |
| Grenada | \>= 5 | Outpatient | 0.150 | 0.077 | 0.260 | 0.098 | 0.045 | 0.185 |
| Guatemala | \< 5 | Inpatient | 0.203 | 0.078 | 0.388 | 0.145 | 0.051 | 0.296 |
| Guatemala | \>= 5 | Inpatient | 0.140 | 0.049 | 0.288 | 0.097 | 0.032 | 0.212 |
| Guatemala | \< 5 | Outpatient | 0.226 | 0.089 | 0.414 | 0.162 | 0.059 | 0.323 |
| Guatemala | \>= 5 | Outpatient | 0.156 | 0.056 | 0.316 | 0.109 | 0.036 | 0.235 |
| Guinea | \< 5 | Inpatient | 0.198 | 0.118 | 0.301 | 0.130 | 0.069 | 0.217 |
| Guinea | \>= 5 | Inpatient | 0.134 | 0.074 | 0.223 | 0.086 | 0.043 | 0.155 |
| Guinea | \< 5 | Outpatient | 0.220 | 0.134 | 0.328 | 0.147 | 0.079 | 0.238 |
| Guinea | \>= 5 | Outpatient | 0.151 | 0.085 | 0.244 | 0.097 | 0.050 | 0.170 |
| Guinea-Bissau | \< 5 | Inpatient | 0.255 | 0.117 | 0.482 | 0.194 | 0.076 | 0.402 |
| Guinea-Bissau | \>= 5 | Inpatient | 0.179 | 0.074 | 0.371 | 0.133 | 0.047 | 0.301 |
| Guinea-Bissau | \< 5 | Outpatient | 0.281 | 0.132 | 0.509 | 0.216 | 0.088 | 0.434 |
| Guinea-Bissau | \>= 5 | Outpatient | 0.200 | 0.084 | 0.399 | 0.149 | 0.054 | 0.327 |
| Guyana | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Guyana | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Guyana | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Guyana | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Haiti | \< 5 | Inpatient | 0.183 | 0.085 | 0.310 | 0.128 | 0.057 | 0.228 |
| Haiti | \>= 5 | Inpatient | 0.124 | 0.053 | 0.227 | 0.085 | 0.035 | 0.162 |
| Haiti | \< 5 | Outpatient | 0.204 | 0.098 | 0.340 | 0.144 | 0.066 | 0.250 |
| Haiti | \>= 5 | Outpatient | 0.140 | 0.062 | 0.250 | 0.096 | 0.041 | 0.179 |
| Honduras | \< 5 | Inpatient | 0.211 | 0.092 | 0.372 | 0.136 | 0.059 | 0.245 |
| Honduras | \>= 5 | Inpatient | 0.145 | 0.057 | 0.276 | 0.090 | 0.036 | 0.175 |
| Honduras | \< 5 | Outpatient | 0.234 | 0.105 | 0.401 | 0.152 | 0.068 | 0.268 |
| Honduras | \>= 5 | Outpatient | 0.162 | 0.065 | 0.301 | 0.102 | 0.041 | 0.194 |
| Hungary | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Hungary | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Hungary | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Hungary | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Iceland | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Iceland | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Iceland | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Iceland | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| India | \< 5 | Inpatient | 0.216 | 0.117 | 0.349 | 0.163 | 0.074 | 0.305 |
| India | \>= 5 | Inpatient | 0.148 | 0.075 | 0.258 | 0.110 | 0.045 | 0.216 |
| India | \< 5 | Outpatient | 0.240 | 0.132 | 0.384 | 0.182 | 0.083 | 0.334 |
| India | \>= 5 | Outpatient | 0.166 | 0.085 | 0.280 | 0.124 | 0.052 | 0.240 |
| Indonesia | \< 5 | Inpatient | 0.139 | 0.052 | 0.266 | 0.090 | 0.030 | 0.187 |
| Indonesia | \>= 5 | Inpatient | 0.093 | 0.032 | 0.190 | 0.059 | 0.018 | 0.128 |
| Indonesia | \< 5 | Outpatient | 0.156 | 0.060 | 0.293 | 0.102 | 0.035 | 0.207 |
| Indonesia | \>= 5 | Outpatient | 0.105 | 0.038 | 0.210 | 0.067 | 0.021 | 0.143 |
| Iran, Islamic Rep. | \< 5 | Inpatient | 0.158 | 0.055 | 0.310 | 0.103 | 0.032 | 0.222 |
| Iran, Islamic Rep. | \>= 5 | Inpatient | 0.106 | 0.035 | 0.225 | 0.068 | 0.020 | 0.154 |
| Iran, Islamic Rep. | \< 5 | Outpatient | 0.176 | 0.064 | 0.337 | 0.116 | 0.038 | 0.246 |
| Iran, Islamic Rep. | \>= 5 | Outpatient | 0.119 | 0.040 | 0.247 | 0.077 | 0.024 | 0.171 |
| Iraq | \< 5 | Inpatient | 0.192 | 0.104 | 0.315 | 0.133 | 0.069 | 0.227 |
| Iraq | \>= 5 | Inpatient | 0.130 | 0.065 | 0.228 | 0.088 | 0.042 | 0.161 |
| Iraq | \< 5 | Outpatient | 0.213 | 0.119 | 0.341 | 0.150 | 0.079 | 0.250 |
| Iraq | \>= 5 | Outpatient | 0.146 | 0.076 | 0.249 | 0.100 | 0.049 | 0.177 |
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
| Jamaica | \< 5 | Inpatient | 0.197 | 0.106 | 0.323 | 0.130 | 0.061 | 0.240 |
| Jamaica | \>= 5 | Inpatient | 0.134 | 0.066 | 0.237 | 0.086 | 0.038 | 0.170 |
| Jamaica | \< 5 | Outpatient | 0.219 | 0.121 | 0.351 | 0.147 | 0.070 | 0.264 |
| Jamaica | \>= 5 | Outpatient | 0.150 | 0.077 | 0.260 | 0.098 | 0.045 | 0.185 |
| Japan | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Japan | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Japan | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Japan | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Jordan | \< 5 | Inpatient | 0.192 | 0.104 | 0.315 | 0.141 | 0.069 | 0.245 |
| Jordan | \>= 5 | Inpatient | 0.130 | 0.065 | 0.228 | 0.094 | 0.042 | 0.175 |
| Jordan | \< 5 | Outpatient | 0.213 | 0.119 | 0.341 | 0.158 | 0.079 | 0.269 |
| Jordan | \>= 5 | Outpatient | 0.146 | 0.076 | 0.249 | 0.106 | 0.049 | 0.194 |
| Kazakhstan | \< 5 | Inpatient | 0.184 | 0.092 | 0.312 | 0.121 | 0.055 | 0.223 |
| Kazakhstan | \>= 5 | Inpatient | 0.125 | 0.058 | 0.224 | 0.080 | 0.034 | 0.155 |
| Kazakhstan | \< 5 | Outpatient | 0.206 | 0.106 | 0.338 | 0.136 | 0.063 | 0.244 |
| Kazakhstan | \>= 5 | Outpatient | 0.140 | 0.067 | 0.246 | 0.090 | 0.040 | 0.172 |
| Kenya | \< 5 | Inpatient | 0.250 | 0.143 | 0.393 | 0.193 | 0.094 | 0.342 |
| Kenya | \>= 5 | Inpatient | 0.174 | 0.093 | 0.291 | 0.131 | 0.059 | 0.251 |
| Kenya | \< 5 | Outpatient | 0.277 | 0.162 | 0.423 | 0.215 | 0.107 | 0.371 |
| Kenya | \>= 5 | Outpatient | 0.194 | 0.106 | 0.315 | 0.147 | 0.068 | 0.272 |
| Kiribati | \< 5 | Inpatient | 0.161 | 0.077 | 0.269 | 0.121 | 0.048 | 0.226 |
| Kiribati | \>= 5 | Inpatient | 0.108 | 0.048 | 0.194 | 0.080 | 0.030 | 0.162 |
| Kiribati | \< 5 | Outpatient | 0.181 | 0.088 | 0.296 | 0.136 | 0.056 | 0.248 |
| Kiribati | \>= 5 | Outpatient | 0.122 | 0.055 | 0.213 | 0.090 | 0.034 | 0.178 |
| Korea, Dem. People’s Rep. | \< 5 | Inpatient | 0.194 | 0.108 | 0.311 | 0.128 | 0.064 | 0.224 |
| Korea, Dem. People’s Rep. | \>= 5 | Inpatient | 0.132 | 0.067 | 0.226 | 0.085 | 0.039 | 0.158 |
| Korea, Dem. People’s Rep. | \< 5 | Outpatient | 0.216 | 0.123 | 0.338 | 0.144 | 0.073 | 0.249 |
| Korea, Dem. People’s Rep. | \>= 5 | Outpatient | 0.148 | 0.079 | 0.249 | 0.096 | 0.046 | 0.175 |
| Korea, Rep. | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Korea, Rep. | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Korea, Rep. | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Korea, Rep. | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kuwait | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kuwait | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kuwait | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kuwait | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kyrgyz Republic | \< 5 | Inpatient | 0.184 | 0.096 | 0.300 | 0.130 | 0.063 | 0.226 |
| Kyrgyz Republic | \>= 5 | Inpatient | 0.125 | 0.060 | 0.219 | 0.086 | 0.039 | 0.160 |
| Kyrgyz Republic | \< 5 | Outpatient | 0.205 | 0.109 | 0.329 | 0.146 | 0.072 | 0.248 |
| Kyrgyz Republic | \>= 5 | Outpatient | 0.140 | 0.070 | 0.240 | 0.097 | 0.045 | 0.176 |
| Lao PDR | \< 5 | Inpatient | 0.161 | 0.077 | 0.269 | 0.105 | 0.046 | 0.190 |
| Lao PDR | \>= 5 | Inpatient | 0.108 | 0.048 | 0.194 | 0.069 | 0.028 | 0.132 |
| Lao PDR | \< 5 | Outpatient | 0.181 | 0.088 | 0.296 | 0.118 | 0.052 | 0.210 |
| Lao PDR | \>= 5 | Outpatient | 0.122 | 0.055 | 0.213 | 0.078 | 0.032 | 0.145 |
| Latvia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Latvia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Latvia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Latvia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Lebanon | \< 5 | Inpatient | 0.207 | 0.084 | 0.391 | 0.137 | 0.052 | 0.285 |
| Lebanon | \>= 5 | Inpatient | 0.142 | 0.053 | 0.292 | 0.092 | 0.033 | 0.204 |
| Lebanon | \< 5 | Outpatient | 0.229 | 0.097 | 0.422 | 0.154 | 0.059 | 0.311 |
| Lebanon | \>= 5 | Outpatient | 0.159 | 0.063 | 0.316 | 0.103 | 0.038 | 0.222 |
| Lesotho | \< 5 | Inpatient | 0.198 | 0.118 | 0.301 | 0.148 | 0.077 | 0.250 |
| Lesotho | \>= 5 | Inpatient | 0.134 | 0.074 | 0.223 | 0.099 | 0.047 | 0.179 |
| Lesotho | \< 5 | Outpatient | 0.220 | 0.134 | 0.328 | 0.166 | 0.090 | 0.273 |
| Lesotho | \>= 5 | Outpatient | 0.151 | 0.085 | 0.244 | 0.112 | 0.055 | 0.198 |
| Liberia | \< 5 | Inpatient | 0.200 | 0.118 | 0.312 | 0.145 | 0.079 | 0.240 |
| Liberia | \>= 5 | Inpatient | 0.136 | 0.074 | 0.228 | 0.096 | 0.048 | 0.171 |
| Liberia | \< 5 | Outpatient | 0.223 | 0.134 | 0.339 | 0.162 | 0.090 | 0.266 |
| Liberia | \>= 5 | Outpatient | 0.153 | 0.085 | 0.250 | 0.109 | 0.055 | 0.189 |
| Libya | \< 5 | Inpatient | 0.192 | 0.104 | 0.315 | 0.140 | 0.069 | 0.243 |
| Libya | \>= 5 | Inpatient | 0.130 | 0.065 | 0.228 | 0.093 | 0.042 | 0.174 |
| Libya | \< 5 | Outpatient | 0.213 | 0.119 | 0.341 | 0.157 | 0.079 | 0.268 |
| Libya | \>= 5 | Outpatient | 0.146 | 0.076 | 0.249 | 0.105 | 0.049 | 0.192 |
| Lithuania | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Lithuania | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Lithuania | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Lithuania | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Luxembourg | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Luxembourg | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Luxembourg | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Luxembourg | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Madagascar | \< 5 | Inpatient | 0.185 | 0.063 | 0.360 | 0.134 | 0.042 | 0.279 |
| Madagascar | \>= 5 | Inpatient | 0.126 | 0.039 | 0.264 | 0.089 | 0.026 | 0.198 |
| Madagascar | \< 5 | Outpatient | 0.206 | 0.072 | 0.392 | 0.150 | 0.048 | 0.308 |
| Madagascar | \>= 5 | Outpatient | 0.142 | 0.046 | 0.289 | 0.101 | 0.030 | 0.222 |
| Malawi | \< 5 | Inpatient | 0.190 | 0.084 | 0.335 | 0.144 | 0.058 | 0.282 |
| Malawi | \>= 5 | Inpatient | 0.130 | 0.053 | 0.249 | 0.096 | 0.036 | 0.203 |
| Malawi | \< 5 | Outpatient | 0.212 | 0.094 | 0.370 | 0.161 | 0.065 | 0.310 |
| Malawi | \>= 5 | Outpatient | 0.146 | 0.060 | 0.274 | 0.109 | 0.041 | 0.228 |
| Malaysia | \< 5 | Inpatient | 0.132 | 0.040 | 0.273 | 0.085 | 0.023 | 0.187 |
| Malaysia | \>= 5 | Inpatient | 0.088 | 0.025 | 0.194 | 0.056 | 0.014 | 0.130 |
| Malaysia | \< 5 | Outpatient | 0.148 | 0.045 | 0.301 | 0.096 | 0.027 | 0.211 |
| Malaysia | \>= 5 | Outpatient | 0.100 | 0.028 | 0.214 | 0.063 | 0.016 | 0.144 |
| Maldives | \< 5 | Inpatient | 0.183 | 0.094 | 0.296 | 0.120 | 0.056 | 0.216 |
| Maldives | \>= 5 | Inpatient | 0.123 | 0.059 | 0.216 | 0.079 | 0.035 | 0.150 |
| Maldives | \< 5 | Outpatient | 0.204 | 0.108 | 0.325 | 0.135 | 0.064 | 0.237 |
| Maldives | \>= 5 | Outpatient | 0.139 | 0.068 | 0.237 | 0.089 | 0.040 | 0.167 |
| Mali | \< 5 | Inpatient | 0.201 | 0.087 | 0.364 | 0.147 | 0.059 | 0.287 |
| Mali | \>= 5 | Inpatient | 0.138 | 0.055 | 0.268 | 0.098 | 0.037 | 0.206 |
| Mali | \< 5 | Outpatient | 0.223 | 0.099 | 0.392 | 0.165 | 0.067 | 0.314 |
| Mali | \>= 5 | Outpatient | 0.154 | 0.063 | 0.297 | 0.111 | 0.042 | 0.229 |
| Malta | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Malta | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Malta | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Malta | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Marshall Islands | \< 5 | Inpatient | 0.165 | 0.076 | 0.271 | 0.112 | 0.048 | 0.197 |
| Marshall Islands | \>= 5 | Inpatient | 0.111 | 0.048 | 0.197 | 0.074 | 0.030 | 0.138 |
| Marshall Islands | \< 5 | Outpatient | 0.184 | 0.087 | 0.301 | 0.126 | 0.056 | 0.220 |
| Marshall Islands | \>= 5 | Outpatient | 0.125 | 0.055 | 0.218 | 0.084 | 0.035 | 0.154 |
| Mauritania | \< 5 | Inpatient | 0.198 | 0.118 | 0.301 | 0.147 | 0.078 | 0.247 |
| Mauritania | \>= 5 | Inpatient | 0.134 | 0.074 | 0.223 | 0.098 | 0.048 | 0.177 |
| Mauritania | \< 5 | Outpatient | 0.220 | 0.134 | 0.328 | 0.165 | 0.090 | 0.270 |
| Mauritania | \>= 5 | Outpatient | 0.151 | 0.085 | 0.244 | 0.111 | 0.055 | 0.196 |
| Mauritius | \< 5 | Inpatient | 0.189 | 0.103 | 0.303 | 0.142 | 0.070 | 0.248 |
| Mauritius | \>= 5 | Inpatient | 0.128 | 0.064 | 0.221 | 0.095 | 0.043 | 0.180 |
| Mauritius | \< 5 | Outpatient | 0.211 | 0.117 | 0.331 | 0.160 | 0.080 | 0.274 |
| Mauritius | \>= 5 | Outpatient | 0.144 | 0.075 | 0.241 | 0.107 | 0.050 | 0.198 |
| Mexico | \< 5 | Inpatient | 0.265 | 0.113 | 0.480 | 0.174 | 0.074 | 0.339 |
| Mexico | \>= 5 | Inpatient | 0.187 | 0.072 | 0.373 | 0.118 | 0.046 | 0.249 |
| Mexico | \< 5 | Outpatient | 0.291 | 0.128 | 0.516 | 0.194 | 0.084 | 0.367 |
| Mexico | \>= 5 | Outpatient | 0.208 | 0.082 | 0.405 | 0.133 | 0.052 | 0.273 |
| Micronesia, Fed. Sts. | \< 5 | Inpatient | 0.180 | 0.078 | 0.313 | 0.111 | 0.050 | 0.194 |
| Micronesia, Fed. Sts. | \>= 5 | Inpatient | 0.122 | 0.048 | 0.230 | 0.073 | 0.031 | 0.137 |
| Micronesia, Fed. Sts. | \< 5 | Outpatient | 0.200 | 0.089 | 0.343 | 0.125 | 0.057 | 0.216 |
| Micronesia, Fed. Sts. | \>= 5 | Outpatient | 0.137 | 0.057 | 0.251 | 0.083 | 0.035 | 0.153 |
| Moldova | \< 5 | Inpatient | 0.184 | 0.092 | 0.312 | 0.131 | 0.061 | 0.237 |
| Moldova | \>= 5 | Inpatient | 0.125 | 0.058 | 0.224 | 0.087 | 0.038 | 0.165 |
| Moldova | \< 5 | Outpatient | 0.206 | 0.106 | 0.338 | 0.147 | 0.070 | 0.260 |
| Moldova | \>= 5 | Outpatient | 0.140 | 0.067 | 0.246 | 0.098 | 0.044 | 0.181 |
| Monaco | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Monaco | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Monaco | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Monaco | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Mongolia | \< 5 | Inpatient | 0.161 | 0.077 | 0.269 | 0.105 | 0.046 | 0.190 |
| Mongolia | \>= 5 | Inpatient | 0.108 | 0.048 | 0.194 | 0.069 | 0.028 | 0.132 |
| Mongolia | \< 5 | Outpatient | 0.181 | 0.088 | 0.296 | 0.118 | 0.052 | 0.210 |
| Mongolia | \>= 5 | Outpatient | 0.122 | 0.055 | 0.213 | 0.078 | 0.032 | 0.145 |
| Montenegro | \< 5 | Inpatient | 0.184 | 0.092 | 0.312 | 0.121 | 0.055 | 0.223 |
| Montenegro | \>= 5 | Inpatient | 0.125 | 0.058 | 0.224 | 0.080 | 0.034 | 0.155 |
| Montenegro | \< 5 | Outpatient | 0.206 | 0.106 | 0.338 | 0.136 | 0.063 | 0.244 |
| Montenegro | \>= 5 | Outpatient | 0.140 | 0.067 | 0.246 | 0.090 | 0.040 | 0.172 |
| Morocco | \< 5 | Inpatient | 0.202 | 0.075 | 0.396 | 0.154 | 0.052 | 0.328 |
| Morocco | \>= 5 | Inpatient | 0.139 | 0.047 | 0.297 | 0.104 | 0.032 | 0.241 |
| Morocco | \< 5 | Outpatient | 0.224 | 0.086 | 0.431 | 0.172 | 0.059 | 0.359 |
| Morocco | \>= 5 | Outpatient | 0.155 | 0.053 | 0.328 | 0.117 | 0.036 | 0.266 |
| Mozambique | \< 5 | Inpatient | 0.227 | 0.106 | 0.405 | 0.165 | 0.071 | 0.318 |
| Mozambique | \>= 5 | Inpatient | 0.157 | 0.066 | 0.302 | 0.112 | 0.044 | 0.230 |
| Mozambique | \< 5 | Outpatient | 0.252 | 0.121 | 0.437 | 0.185 | 0.081 | 0.349 |
| Mozambique | \>= 5 | Outpatient | 0.176 | 0.076 | 0.331 | 0.126 | 0.050 | 0.254 |
| Myanmar | \< 5 | Inpatient | 0.194 | 0.108 | 0.311 | 0.137 | 0.071 | 0.235 |
| Myanmar | \>= 5 | Inpatient | 0.132 | 0.067 | 0.226 | 0.091 | 0.044 | 0.168 |
| Myanmar | \< 5 | Outpatient | 0.216 | 0.123 | 0.338 | 0.154 | 0.081 | 0.260 |
| Myanmar | \>= 5 | Outpatient | 0.148 | 0.079 | 0.249 | 0.103 | 0.050 | 0.185 |
| Namibia | \< 5 | Inpatient | 0.189 | 0.103 | 0.303 | 0.139 | 0.071 | 0.235 |
| Namibia | \>= 5 | Inpatient | 0.128 | 0.064 | 0.221 | 0.092 | 0.044 | 0.169 |
| Namibia | \< 5 | Outpatient | 0.211 | 0.117 | 0.331 | 0.156 | 0.080 | 0.258 |
| Namibia | \>= 5 | Outpatient | 0.144 | 0.075 | 0.241 | 0.104 | 0.050 | 0.186 |
| Nauru | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nauru | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nauru | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nauru | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nepal | \< 5 | Inpatient | 0.185 | 0.084 | 0.322 | 0.121 | 0.051 | 0.229 |
| Nepal | \>= 5 | Inpatient | 0.125 | 0.053 | 0.234 | 0.080 | 0.031 | 0.161 |
| Nepal | \< 5 | Outpatient | 0.206 | 0.095 | 0.351 | 0.137 | 0.057 | 0.255 |
| Nepal | \>= 5 | Outpatient | 0.141 | 0.060 | 0.258 | 0.091 | 0.036 | 0.181 |
| Netherlands | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Netherlands | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Netherlands | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Netherlands | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| New Zealand | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| New Zealand | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| New Zealand | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| New Zealand | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nicaragua | \< 5 | Inpatient | 0.181 | 0.049 | 0.377 | 0.117 | 0.030 | 0.257 |
| Nicaragua | \>= 5 | Inpatient | 0.124 | 0.030 | 0.281 | 0.078 | 0.018 | 0.184 |
| Nicaragua | \< 5 | Outpatient | 0.201 | 0.057 | 0.411 | 0.132 | 0.035 | 0.282 |
| Nicaragua | \>= 5 | Outpatient | 0.139 | 0.035 | 0.310 | 0.088 | 0.021 | 0.203 |
| Niger | \< 5 | Inpatient | 0.200 | 0.118 | 0.312 | 0.149 | 0.077 | 0.255 |
| Niger | \>= 5 | Inpatient | 0.136 | 0.074 | 0.228 | 0.100 | 0.047 | 0.183 |
| Niger | \< 5 | Outpatient | 0.223 | 0.134 | 0.339 | 0.167 | 0.088 | 0.280 |
| Niger | \>= 5 | Outpatient | 0.153 | 0.085 | 0.250 | 0.112 | 0.055 | 0.201 |
| Nigeria | \< 5 | Inpatient | 0.198 | 0.118 | 0.301 | 0.130 | 0.069 | 0.217 |
| Nigeria | \>= 5 | Inpatient | 0.134 | 0.074 | 0.223 | 0.086 | 0.043 | 0.155 |
| Nigeria | \< 5 | Outpatient | 0.220 | 0.134 | 0.328 | 0.147 | 0.079 | 0.238 |
| Nigeria | \>= 5 | Outpatient | 0.151 | 0.085 | 0.244 | 0.097 | 0.050 | 0.170 |
| Niue | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Niue | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Niue | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Niue | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| North Macedonia | \< 5 | Inpatient | 0.184 | 0.092 | 0.312 | 0.132 | 0.061 | 0.240 |
| North Macedonia | \>= 5 | Inpatient | 0.125 | 0.058 | 0.224 | 0.088 | 0.038 | 0.167 |
| North Macedonia | \< 5 | Outpatient | 0.206 | 0.106 | 0.338 | 0.148 | 0.070 | 0.262 |
| North Macedonia | \>= 5 | Outpatient | 0.140 | 0.067 | 0.246 | 0.099 | 0.044 | 0.184 |
| Norway | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Norway | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Norway | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Norway | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Oman | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Oman | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Oman | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Oman | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Pakistan | \< 5 | Inpatient | 0.309 | 0.144 | 0.533 | 0.241 | 0.095 | 0.469 |
| Pakistan | \>= 5 | Inpatient | 0.222 | 0.094 | 0.423 | 0.169 | 0.060 | 0.362 |
| Pakistan | \< 5 | Outpatient | 0.337 | 0.166 | 0.563 | 0.266 | 0.110 | 0.501 |
| Pakistan | \>= 5 | Outpatient | 0.245 | 0.107 | 0.452 | 0.188 | 0.069 | 0.388 |
| Palau | \< 5 | Inpatient | 0.179 | 0.074 | 0.318 | 0.117 | 0.048 | 0.215 |
| Palau | \>= 5 | Inpatient | 0.121 | 0.046 | 0.234 | 0.077 | 0.029 | 0.153 |
| Palau | \< 5 | Outpatient | 0.200 | 0.085 | 0.348 | 0.132 | 0.055 | 0.239 |
| Palau | \>= 5 | Outpatient | 0.137 | 0.054 | 0.258 | 0.088 | 0.034 | 0.169 |
| Panama | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Panama | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Panama | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Panama | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Papua New Guinea | \< 5 | Inpatient | 0.161 | 0.077 | 0.269 | 0.105 | 0.046 | 0.190 |
| Papua New Guinea | \>= 5 | Inpatient | 0.108 | 0.048 | 0.194 | 0.069 | 0.028 | 0.132 |
| Papua New Guinea | \< 5 | Outpatient | 0.181 | 0.088 | 0.296 | 0.118 | 0.052 | 0.210 |
| Papua New Guinea | \>= 5 | Outpatient | 0.122 | 0.055 | 0.213 | 0.078 | 0.032 | 0.145 |
| Paraguay | \< 5 | Inpatient | 0.212 | 0.118 | 0.337 | 0.145 | 0.074 | 0.249 |
| Paraguay | \>= 5 | Inpatient | 0.145 | 0.073 | 0.247 | 0.097 | 0.045 | 0.179 |
| Paraguay | \< 5 | Outpatient | 0.235 | 0.135 | 0.366 | 0.163 | 0.085 | 0.275 |
| Paraguay | \>= 5 | Outpatient | 0.162 | 0.085 | 0.271 | 0.109 | 0.053 | 0.196 |
| Peru | \< 5 | Inpatient | 0.307 | 0.142 | 0.540 | 0.214 | 0.092 | 0.409 |
| Peru | \>= 5 | Inpatient | 0.221 | 0.092 | 0.431 | 0.148 | 0.058 | 0.308 |
| Peru | \< 5 | Outpatient | 0.336 | 0.165 | 0.569 | 0.237 | 0.105 | 0.439 |
| Peru | \>= 5 | Outpatient | 0.244 | 0.105 | 0.459 | 0.165 | 0.067 | 0.334 |
| Philippines | \< 5 | Inpatient | 0.161 | 0.077 | 0.269 | 0.105 | 0.046 | 0.190 |
| Philippines | \>= 5 | Inpatient | 0.108 | 0.048 | 0.194 | 0.069 | 0.028 | 0.132 |
| Philippines | \< 5 | Outpatient | 0.181 | 0.088 | 0.296 | 0.118 | 0.052 | 0.210 |
| Philippines | \>= 5 | Outpatient | 0.122 | 0.055 | 0.213 | 0.078 | 0.032 | 0.145 |
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
| Russian Federation | \< 5 | Inpatient | 0.184 | 0.092 | 0.312 | 0.121 | 0.055 | 0.223 |
| Russian Federation | \>= 5 | Inpatient | 0.125 | 0.058 | 0.224 | 0.080 | 0.034 | 0.155 |
| Russian Federation | \< 5 | Outpatient | 0.206 | 0.106 | 0.338 | 0.136 | 0.063 | 0.244 |
| Russian Federation | \>= 5 | Outpatient | 0.140 | 0.067 | 0.246 | 0.090 | 0.040 | 0.172 |
| Rwanda | \< 5 | Inpatient | 0.200 | 0.118 | 0.312 | 0.151 | 0.076 | 0.263 |
| Rwanda | \>= 5 | Inpatient | 0.136 | 0.074 | 0.228 | 0.101 | 0.046 | 0.189 |
| Rwanda | \< 5 | Outpatient | 0.223 | 0.134 | 0.339 | 0.170 | 0.087 | 0.288 |
| Rwanda | \>= 5 | Outpatient | 0.153 | 0.085 | 0.250 | 0.114 | 0.054 | 0.207 |
| Samoa | \< 5 | Inpatient | 0.161 | 0.077 | 0.269 | 0.105 | 0.046 | 0.190 |
| Samoa | \>= 5 | Inpatient | 0.108 | 0.048 | 0.194 | 0.069 | 0.028 | 0.132 |
| Samoa | \< 5 | Outpatient | 0.181 | 0.088 | 0.296 | 0.118 | 0.052 | 0.210 |
| Samoa | \>= 5 | Outpatient | 0.122 | 0.055 | 0.213 | 0.078 | 0.032 | 0.145 |
| San Marino | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| San Marino | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| San Marino | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| San Marino | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| São Tomé and Principe | \< 5 | Inpatient | 0.198 | 0.118 | 0.301 | 0.148 | 0.077 | 0.250 |
| São Tomé and Principe | \>= 5 | Inpatient | 0.134 | 0.074 | 0.223 | 0.099 | 0.047 | 0.179 |
| São Tomé and Principe | \< 5 | Outpatient | 0.220 | 0.134 | 0.328 | 0.166 | 0.090 | 0.273 |
| São Tomé and Principe | \>= 5 | Outpatient | 0.151 | 0.085 | 0.244 | 0.112 | 0.055 | 0.198 |
| Saudi Arabia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Saudi Arabia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Saudi Arabia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Saudi Arabia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Senegal | \< 5 | Inpatient | 0.177 | 0.059 | 0.339 | 0.133 | 0.040 | 0.279 |
| Senegal | \>= 5 | Inpatient | 0.120 | 0.036 | 0.246 | 0.089 | 0.024 | 0.199 |
| Senegal | \< 5 | Outpatient | 0.197 | 0.067 | 0.365 | 0.149 | 0.045 | 0.304 |
| Senegal | \>= 5 | Outpatient | 0.135 | 0.042 | 0.271 | 0.101 | 0.028 | 0.218 |
| Serbia | \< 5 | Inpatient | 0.184 | 0.092 | 0.312 | 0.121 | 0.055 | 0.223 |
| Serbia | \>= 5 | Inpatient | 0.125 | 0.058 | 0.224 | 0.080 | 0.034 | 0.155 |
| Serbia | \< 5 | Outpatient | 0.206 | 0.106 | 0.338 | 0.136 | 0.063 | 0.244 |
| Serbia | \>= 5 | Outpatient | 0.140 | 0.067 | 0.246 | 0.090 | 0.040 | 0.172 |
| Seychelles | \< 5 | Inpatient | 0.189 | 0.103 | 0.303 | 0.144 | 0.070 | 0.255 |
| Seychelles | \>= 5 | Inpatient | 0.128 | 0.064 | 0.221 | 0.096 | 0.042 | 0.185 |
| Seychelles | \< 5 | Outpatient | 0.211 | 0.117 | 0.331 | 0.162 | 0.080 | 0.279 |
| Seychelles | \>= 5 | Outpatient | 0.144 | 0.075 | 0.241 | 0.108 | 0.049 | 0.204 |
| Sierra Leone | \< 5 | Inpatient | 0.200 | 0.118 | 0.312 | 0.151 | 0.076 | 0.262 |
| Sierra Leone | \>= 5 | Inpatient | 0.136 | 0.074 | 0.228 | 0.101 | 0.046 | 0.189 |
| Sierra Leone | \< 5 | Outpatient | 0.223 | 0.134 | 0.339 | 0.169 | 0.087 | 0.287 |
| Sierra Leone | \>= 5 | Outpatient | 0.153 | 0.085 | 0.250 | 0.114 | 0.054 | 0.206 |
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
| Solomon Islands | \< 5 | Inpatient | 0.161 | 0.077 | 0.269 | 0.111 | 0.050 | 0.194 |
| Solomon Islands | \>= 5 | Inpatient | 0.108 | 0.048 | 0.194 | 0.073 | 0.031 | 0.138 |
| Solomon Islands | \< 5 | Outpatient | 0.181 | 0.088 | 0.296 | 0.125 | 0.057 | 0.216 |
| Solomon Islands | \>= 5 | Outpatient | 0.122 | 0.055 | 0.213 | 0.083 | 0.035 | 0.153 |
| Somalia | \< 5 | Inpatient | 0.189 | 0.090 | 0.326 | 0.125 | 0.054 | 0.238 |
| Somalia | \>= 5 | Inpatient | 0.129 | 0.057 | 0.237 | 0.082 | 0.033 | 0.164 |
| Somalia | \< 5 | Outpatient | 0.211 | 0.103 | 0.353 | 0.140 | 0.061 | 0.257 |
| Somalia | \>= 5 | Outpatient | 0.144 | 0.066 | 0.262 | 0.093 | 0.039 | 0.182 |
| South Africa | \< 5 | Inpatient | 0.213 | 0.107 | 0.357 | 0.146 | 0.069 | 0.260 |
| South Africa | \>= 5 | Inpatient | 0.147 | 0.066 | 0.267 | 0.098 | 0.043 | 0.189 |
| South Africa | \< 5 | Outpatient | 0.237 | 0.121 | 0.388 | 0.164 | 0.078 | 0.285 |
| South Africa | \>= 5 | Outpatient | 0.164 | 0.077 | 0.294 | 0.110 | 0.049 | 0.206 |
| South Sudan | \< 5 | Inpatient | 0.200 | 0.118 | 0.312 | 0.132 | 0.069 | 0.223 |
| South Sudan | \>= 5 | Inpatient | 0.136 | 0.074 | 0.228 | 0.087 | 0.043 | 0.158 |
| South Sudan | \< 5 | Outpatient | 0.223 | 0.134 | 0.339 | 0.149 | 0.079 | 0.247 |
| South Sudan | \>= 5 | Outpatient | 0.153 | 0.085 | 0.250 | 0.099 | 0.050 | 0.175 |
| Spain | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Spain | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Spain | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Spain | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Sri Lanka | \< 5 | Inpatient | 0.194 | 0.108 | 0.311 | 0.128 | 0.064 | 0.224 |
| Sri Lanka | \>= 5 | Inpatient | 0.132 | 0.067 | 0.226 | 0.085 | 0.039 | 0.158 |
| Sri Lanka | \< 5 | Outpatient | 0.216 | 0.123 | 0.338 | 0.144 | 0.073 | 0.249 |
| Sri Lanka | \>= 5 | Outpatient | 0.148 | 0.079 | 0.249 | 0.096 | 0.046 | 0.175 |
| St. Kitts and Nevis | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| St. Kitts and Nevis | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| St. Kitts and Nevis | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| St. Kitts and Nevis | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| St. Lucia | \< 5 | Inpatient | 0.197 | 0.106 | 0.323 | 0.130 | 0.061 | 0.240 |
| St. Lucia | \>= 5 | Inpatient | 0.134 | 0.066 | 0.237 | 0.086 | 0.038 | 0.170 |
| St. Lucia | \< 5 | Outpatient | 0.219 | 0.121 | 0.351 | 0.147 | 0.070 | 0.264 |
| St. Lucia | \>= 5 | Outpatient | 0.150 | 0.077 | 0.260 | 0.098 | 0.045 | 0.185 |
| St. Vincent and the Grenadines | \< 5 | Inpatient | 0.197 | 0.106 | 0.323 | 0.130 | 0.061 | 0.240 |
| St. Vincent and the Grenadines | \>= 5 | Inpatient | 0.134 | 0.066 | 0.237 | 0.086 | 0.038 | 0.170 |
| St. Vincent and the Grenadines | \< 5 | Outpatient | 0.219 | 0.121 | 0.351 | 0.147 | 0.070 | 0.264 |
| St. Vincent and the Grenadines | \>= 5 | Outpatient | 0.150 | 0.077 | 0.260 | 0.098 | 0.045 | 0.185 |
| Sudan | \< 5 | Inpatient | 0.192 | 0.064 | 0.392 | 0.145 | 0.047 | 0.312 |
| Sudan | \>= 5 | Inpatient | 0.132 | 0.040 | 0.294 | 0.097 | 0.029 | 0.224 |
| Sudan | \< 5 | Outpatient | 0.214 | 0.075 | 0.421 | 0.162 | 0.055 | 0.336 |
| Sudan | \>= 5 | Outpatient | 0.148 | 0.047 | 0.320 | 0.110 | 0.034 | 0.246 |
| Suriname | \< 5 | Inpatient | 0.197 | 0.106 | 0.323 | 0.130 | 0.061 | 0.240 |
| Suriname | \>= 5 | Inpatient | 0.134 | 0.066 | 0.237 | 0.086 | 0.038 | 0.170 |
| Suriname | \< 5 | Outpatient | 0.219 | 0.121 | 0.351 | 0.147 | 0.070 | 0.264 |
| Suriname | \>= 5 | Outpatient | 0.150 | 0.077 | 0.260 | 0.098 | 0.045 | 0.185 |
| Sweden | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Sweden | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Sweden | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Sweden | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Switzerland | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Switzerland | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Switzerland | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Switzerland | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Syrian Arab Republic | \< 5 | Inpatient | 0.189 | 0.090 | 0.326 | 0.125 | 0.054 | 0.238 |
| Syrian Arab Republic | \>= 5 | Inpatient | 0.129 | 0.057 | 0.237 | 0.082 | 0.033 | 0.164 |
| Syrian Arab Republic | \< 5 | Outpatient | 0.211 | 0.103 | 0.353 | 0.140 | 0.061 | 0.257 |
| Syrian Arab Republic | \>= 5 | Outpatient | 0.144 | 0.066 | 0.262 | 0.093 | 0.039 | 0.182 |
| Tajikistan | \< 5 | Inpatient | 0.184 | 0.096 | 0.300 | 0.140 | 0.062 | 0.260 |
| Tajikistan | \>= 5 | Inpatient | 0.125 | 0.060 | 0.219 | 0.093 | 0.038 | 0.186 |
| Tajikistan | \< 5 | Outpatient | 0.205 | 0.109 | 0.329 | 0.157 | 0.071 | 0.285 |
| Tajikistan | \>= 5 | Outpatient | 0.140 | 0.070 | 0.240 | 0.105 | 0.044 | 0.204 |
| Tanzania | \< 5 | Inpatient | 0.173 | 0.079 | 0.296 | 0.129 | 0.055 | 0.239 |
| Tanzania | \>= 5 | Inpatient | 0.117 | 0.049 | 0.216 | 0.086 | 0.033 | 0.169 |
| Tanzania | \< 5 | Outpatient | 0.193 | 0.090 | 0.326 | 0.145 | 0.064 | 0.262 |
| Tanzania | \>= 5 | Outpatient | 0.131 | 0.057 | 0.239 | 0.097 | 0.039 | 0.190 |
| Thailand | \< 5 | Inpatient | 0.216 | 0.104 | 0.376 | 0.144 | 0.061 | 0.284 |
| Thailand | \>= 5 | Inpatient | 0.148 | 0.064 | 0.283 | 0.096 | 0.038 | 0.201 |
| Thailand | \< 5 | Outpatient | 0.239 | 0.117 | 0.412 | 0.162 | 0.070 | 0.309 |
| Thailand | \>= 5 | Outpatient | 0.166 | 0.074 | 0.309 | 0.109 | 0.044 | 0.221 |
| Timor-Leste | \< 5 | Inpatient | 0.194 | 0.108 | 0.311 | 0.128 | 0.064 | 0.224 |
| Timor-Leste | \>= 5 | Inpatient | 0.132 | 0.067 | 0.226 | 0.085 | 0.039 | 0.158 |
| Timor-Leste | \< 5 | Outpatient | 0.216 | 0.123 | 0.338 | 0.144 | 0.073 | 0.249 |
| Timor-Leste | \>= 5 | Outpatient | 0.148 | 0.079 | 0.249 | 0.096 | 0.046 | 0.175 |
| Togo | \< 5 | Inpatient | 0.200 | 0.118 | 0.312 | 0.148 | 0.078 | 0.250 |
| Togo | \>= 5 | Inpatient | 0.136 | 0.074 | 0.228 | 0.099 | 0.048 | 0.179 |
| Togo | \< 5 | Outpatient | 0.223 | 0.134 | 0.339 | 0.166 | 0.088 | 0.276 |
| Togo | \>= 5 | Outpatient | 0.153 | 0.085 | 0.250 | 0.111 | 0.055 | 0.198 |
| Tonga | \< 5 | Inpatient | 0.159 | 0.073 | 0.265 | 0.103 | 0.045 | 0.183 |
| Tonga | \>= 5 | Inpatient | 0.107 | 0.046 | 0.190 | 0.067 | 0.027 | 0.127 |
| Tonga | \< 5 | Outpatient | 0.178 | 0.084 | 0.292 | 0.116 | 0.051 | 0.203 |
| Tonga | \>= 5 | Outpatient | 0.120 | 0.053 | 0.211 | 0.076 | 0.031 | 0.143 |
| Trinidad and Tobago | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Trinidad and Tobago | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Trinidad and Tobago | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Trinidad and Tobago | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Tunisia | \< 5 | Inpatient | 0.192 | 0.104 | 0.315 | 0.126 | 0.063 | 0.227 |
| Tunisia | \>= 5 | Inpatient | 0.130 | 0.065 | 0.228 | 0.083 | 0.038 | 0.159 |
| Tunisia | \< 5 | Outpatient | 0.213 | 0.119 | 0.341 | 0.142 | 0.072 | 0.247 |
| Tunisia | \>= 5 | Outpatient | 0.146 | 0.076 | 0.249 | 0.094 | 0.045 | 0.174 |
| Türkiye | \< 5 | Inpatient | 0.183 | 0.079 | 0.333 | 0.120 | 0.048 | 0.234 |
| Türkiye | \>= 5 | Inpatient | 0.124 | 0.049 | 0.241 | 0.079 | 0.030 | 0.164 |
| Türkiye | \< 5 | Outpatient | 0.204 | 0.091 | 0.364 | 0.135 | 0.055 | 0.261 |
| Türkiye | \>= 5 | Outpatient | 0.139 | 0.057 | 0.264 | 0.089 | 0.034 | 0.182 |
| Turkmenistan | \< 5 | Inpatient | 0.184 | 0.092 | 0.312 | 0.140 | 0.060 | 0.265 |
| Turkmenistan | \>= 5 | Inpatient | 0.125 | 0.058 | 0.224 | 0.094 | 0.037 | 0.189 |
| Turkmenistan | \< 5 | Outpatient | 0.206 | 0.106 | 0.338 | 0.157 | 0.069 | 0.292 |
| Turkmenistan | \>= 5 | Outpatient | 0.140 | 0.067 | 0.246 | 0.106 | 0.043 | 0.209 |
| Tuvalu | \< 5 | Inpatient | 0.159 | 0.073 | 0.265 | 0.103 | 0.045 | 0.183 |
| Tuvalu | \>= 5 | Inpatient | 0.107 | 0.046 | 0.190 | 0.067 | 0.027 | 0.127 |
| Tuvalu | \< 5 | Outpatient | 0.178 | 0.084 | 0.292 | 0.116 | 0.051 | 0.203 |
| Tuvalu | \>= 5 | Outpatient | 0.120 | 0.053 | 0.211 | 0.076 | 0.031 | 0.143 |
| Uganda | \< 5 | Inpatient | 0.200 | 0.118 | 0.312 | 0.150 | 0.077 | 0.261 |
| Uganda | \>= 5 | Inpatient | 0.136 | 0.074 | 0.228 | 0.101 | 0.047 | 0.187 |
| Uganda | \< 5 | Outpatient | 0.223 | 0.134 | 0.339 | 0.169 | 0.088 | 0.285 |
| Uganda | \>= 5 | Outpatient | 0.153 | 0.085 | 0.250 | 0.114 | 0.054 | 0.205 |
| Ukraine | \< 5 | Inpatient | 0.184 | 0.096 | 0.300 | 0.121 | 0.056 | 0.218 |
| Ukraine | \>= 5 | Inpatient | 0.125 | 0.060 | 0.219 | 0.080 | 0.035 | 0.150 |
| Ukraine | \< 5 | Outpatient | 0.205 | 0.109 | 0.329 | 0.136 | 0.065 | 0.238 |
| Ukraine | \>= 5 | Outpatient | 0.140 | 0.070 | 0.240 | 0.090 | 0.041 | 0.166 |
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
| Uzbekistan | \< 5 | Inpatient | 0.184 | 0.096 | 0.300 | 0.139 | 0.062 | 0.255 |
| Uzbekistan | \>= 5 | Inpatient | 0.125 | 0.060 | 0.219 | 0.092 | 0.038 | 0.182 |
| Uzbekistan | \< 5 | Outpatient | 0.205 | 0.109 | 0.329 | 0.156 | 0.071 | 0.279 |
| Uzbekistan | \>= 5 | Outpatient | 0.140 | 0.070 | 0.240 | 0.104 | 0.045 | 0.200 |
| Vanuatu | \< 5 | Inpatient | 0.161 | 0.077 | 0.269 | 0.105 | 0.046 | 0.190 |
| Vanuatu | \>= 5 | Inpatient | 0.108 | 0.048 | 0.194 | 0.069 | 0.028 | 0.132 |
| Vanuatu | \< 5 | Outpatient | 0.181 | 0.088 | 0.296 | 0.118 | 0.052 | 0.210 |
| Vanuatu | \>= 5 | Outpatient | 0.122 | 0.055 | 0.213 | 0.078 | 0.032 | 0.145 |
| Venezuela | \< 5 | Inpatient | 0.205 | 0.078 | 0.392 | 0.129 | 0.041 | 0.281 |
| Venezuela | \>= 5 | Inpatient | 0.141 | 0.049 | 0.294 | 0.086 | 0.025 | 0.201 |
| Venezuela | \< 5 | Outpatient | 0.228 | 0.091 | 0.420 | 0.145 | 0.047 | 0.313 |
| Venezuela | \>= 5 | Outpatient | 0.158 | 0.056 | 0.319 | 0.097 | 0.029 | 0.223 |
| Vietnam | \< 5 | Inpatient | 0.137 | 0.058 | 0.253 | 0.088 | 0.034 | 0.177 |
| Vietnam | \>= 5 | Inpatient | 0.091 | 0.036 | 0.179 | 0.058 | 0.021 | 0.123 |
| Vietnam | \< 5 | Outpatient | 0.153 | 0.067 | 0.278 | 0.100 | 0.039 | 0.195 |
| Vietnam | \>= 5 | Outpatient | 0.103 | 0.042 | 0.198 | 0.065 | 0.024 | 0.135 |
| Yemen, Rep. | \< 5 | Inpatient | 0.189 | 0.090 | 0.326 | 0.134 | 0.061 | 0.241 |
| Yemen, Rep. | \>= 5 | Inpatient | 0.129 | 0.057 | 0.237 | 0.089 | 0.038 | 0.171 |
| Yemen, Rep. | \< 5 | Outpatient | 0.211 | 0.103 | 0.353 | 0.151 | 0.071 | 0.264 |
| Yemen, Rep. | \>= 5 | Outpatient | 0.144 | 0.066 | 0.262 | 0.101 | 0.044 | 0.189 |
| Zambia | \< 5 | Inpatient | 0.224 | 0.101 | 0.408 | 0.170 | 0.068 | 0.347 |
| Zambia | \>= 5 | Inpatient | 0.155 | 0.063 | 0.308 | 0.116 | 0.042 | 0.254 |
| Zambia | \< 5 | Outpatient | 0.248 | 0.116 | 0.439 | 0.190 | 0.077 | 0.380 |
| Zambia | \>= 5 | Outpatient | 0.173 | 0.072 | 0.337 | 0.130 | 0.048 | 0.279 |
| Zimbabwe | \< 5 | Inpatient | 0.232 | 0.101 | 0.442 | 0.177 | 0.069 | 0.365 |
| Zimbabwe | \>= 5 | Inpatient | 0.162 | 0.064 | 0.338 | 0.121 | 0.043 | 0.270 |
| Zimbabwe | \< 5 | Outpatient | 0.257 | 0.115 | 0.470 | 0.197 | 0.080 | 0.394 |
| Zimbabwe | \>= 5 | Outpatient | 0.180 | 0.074 | 0.365 | 0.135 | 0.050 | 0.291 |

Estimated `r params$Pathogen` proportion by country, 2010 vs 2020

# Session info

``` r
sessioninfo::session_info()
```

    ## Warning in system2("quarto", "-V", stdout = TRUE, env = paste0("TMPDIR=", : running command '"quarto"
    ## TMPDIR=C:/Users/fbbu6966/AppData/Local/Temp/RtmpYbH8Nj/file366463dc2435 -V' had status 1

    ## ─ Session info ──────────────────────────────────────────────────────────────────────────────────────────────────────────
    ##  setting  value
    ##  version  R version 4.5.0 (2025-04-11 ucrt)
    ##  os       Windows 10 x64 (build 19045)
    ##  system   x86_64, mingw32
    ##  ui       RStudio
    ##  language (EN)
    ##  collate  English_United States.utf8
    ##  ctype    English_United States.utf8
    ##  tz       Europe/Brussels
    ##  date     2025-08-31
    ##  rstudio  2025.05.0+496 Mariposa Orchid (desktop)
    ##  pandoc   3.4 @ C:/Users/fbbu6966/AppData/Local/Programs/RStudio/resources/app/bin/quarto/bin/tools/ (via rmarkdown)
    ##  quarto   ERROR: Unknown command "TMPDIR=C:/Users/fbbu6966/AppData/Local/Temp/RtmpYbH8Nj/file366463dc2435". Did you mean command "update"? @ C:\\Users\\fbbu6966\\AppData\\Local\\Programs\\RStudio\\RESOUR~1\\app\\bin\\quarto\\bin\\quarto.exe
    ## 
    ## ─ Packages ──────────────────────────────────────────────────────────────────────────────────────────────────────────────
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
    ## ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

``` r
# Save dataset for report created for expert to receive feedback
# save(all_cnt_rt, file="./00-Report_FB/all_cnt_rt.Rdata")
# save(all_glb_prop, file="./00-Report_FB/all_glb_prop.Rdata")
# save(all_reg_prop, file="./00-Report_FB/all_reg_prop.Rdata")
# save(all_reg_rt, file="./00-Report_FB/all_reg_rt.Rdata")
# save(all_sub_nr, file="./00-Report_FB/all_sub_nr.Rdata")
# save(all_sub_rt, file="./00-Report_FB/all_sub_rt.Rdata")
```
