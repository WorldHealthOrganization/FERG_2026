Global proportion of STEC (Diarrheal) • SYMPT+ROTA • Estimate proportion
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
image <- paste0("03-estimate-rota_files/figure-gfm/imputation_map.png")
cat("![](", image ,")")
```

![](03-estimate-rota_files/figure-gfm/imputation_map.png)

``` r
fit_brms_reg_s <- readRDS(paste0(params$Dir, "/fit_brms_reg_s.rds"))
print(fit_brms_reg_s)
```

    ## Warning: There were 7 divergent transitions after warmup. Increasing adapt_delta above 0.95 may help. See
    ## http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup

    ##  Family: gaussian 
    ##   Links: mu = identity; sigma = identity 
    ## Formula: yi | se(sei) ~ 1 + YEAR + SYNDROMTYPE + AGE + REFERENCE + COVERAGE + STRAIN + (1 | REG2) + (1 | REG2:SUB2) + (1 | REG2:SUB2:COUNTRY) + (1 | REG2:SUB2:COUNTRY:ID) + (1 | REG2:SUB2:COUNTRY:ID:DTP_ID) 
    ##    Data: subset(es, as.integer(FLAG) == 1) (Number of observations: 237) 
    ##   Draws: 5 chains, each with iter = 5000; warmup = 3000; thin = 1;
    ##          total post-warmup draws = 10000
    ## 
    ## Multilevel Hyperparameters:
    ## ~REG2 (Number of levels: 6) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.62      0.40     0.03     1.53 1.00     4393     5428
    ## 
    ## ~REG2:SUB2 (Number of levels: 10) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.36      0.28     0.01     1.04 1.00     4947     5762
    ## 
    ## ~REG2:SUB2:COUNTRY (Number of levels: 36) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.43      0.26     0.02     0.95 1.00     1776     4227
    ## 
    ## ~REG2:SUB2:COUNTRY:ID (Number of levels: 79) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     1.04      0.14     0.79     1.32 1.00     4566     6120
    ## 
    ## ~REG2:SUB2:COUNTRY:ID:DTP_ID (Number of levels: 237) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.44      0.06     0.33     0.56 1.00     4650     6924
    ## 
    ## Regression Coefficients:
    ##                       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept               -64.31     72.31  -207.08    77.38 1.00     7647     7134
    ## YEAR                      0.03      0.04    -0.04     0.10 1.00     7646     7177
    ## SYNDROMTYPEOutpatient     0.01      0.19    -0.37     0.38 1.00    11844     7951
    ## AGEAgebelow5             -0.01      0.20    -0.40     0.36 1.00     9387     7363
    ## AGEMixedages              0.05      0.26    -0.46     0.56 1.00     9801     8123
    ## REFERENCEOther            0.71      0.97    -1.19     2.63 1.00    10295     7910
    ## COVERAGE                 -0.00      0.01    -0.01     0.01 1.00     9326     7757
    ## STRAINother              -1.93      0.58    -3.08    -0.78 1.00     7983     7757
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

kable(
  caption = "High income countries, estimates not made with this approach",
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

    ## [1] 3168

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

    ## [1] 6864

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

    ## [1] 1760

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
    ##  $ VAL_MEAN     : num  0.0306 0.0305 0.03 0.03 0.0304 ...
    ##  $ VAL_LWR      : num  0.00676 0.00705 0.00729 0.00696 0.00788 ...
    ##  $ VAL_UPR      : num  0.093 0.0922 0.0889 0.089 0.0878 ...
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
    ##  $ VAL_MEAN     : num  0.0306 0.0305 0.03 0.03 0.0304 ...
    ##  $ VAL_LWR      : num  0.00676 0.00705 0.00729 0.00696 0.00788 ...
    ##  $ VAL_UPR      : num  0.093 0.0922 0.0889 0.089 0.0878 ...
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

\[1\] 0.00 0.01 0.02 0.03 0.04 0.05 0.06 0.07

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

\[1\] 0.00 0.01 0.02 0.03 0.04 0.05 0.06 0.07

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

\[1\] 0.00 0.01 0.02 0.03 0.04 0.05 0.06 0.07

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

\[1\] 0.00 0.01 0.02 0.03 0.04 0.05 0.06 0.07

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

\[1\] 0.00 0.02 0.04 0.06 0.08 0.10

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

\[1\] 0.00 0.02 0.04 0.06 0.08 0.10

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

\[1\] 0.00 0.02 0.04 0.06 0.08 0.10

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

\[1\] 0.00 0.02 0.04 0.06 0.08 0.10

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
| Afghanistan | \< 5 | Inpatient | 0.037 | 0.013 | 0.090 | 0.048 | 0.014 | 0.123 |
| Afghanistan | \>= 5 | Inpatient | 0.038 | 0.012 | 0.095 | 0.049 | 0.014 | 0.131 |
| Afghanistan | \< 5 | Outpatient | 0.037 | 0.013 | 0.089 | 0.049 | 0.015 | 0.123 |
| Afghanistan | \>= 5 | Outpatient | 0.037 | 0.013 | 0.091 | 0.050 | 0.014 | 0.129 |
| Albania | \< 5 | Inpatient | 0.043 | 0.010 | 0.138 | 0.055 | 0.011 | 0.180 |
| Albania | \>= 5 | Inpatient | 0.044 | 0.009 | 0.143 | 0.056 | 0.011 | 0.184 |
| Albania | \< 5 | Outpatient | 0.043 | 0.010 | 0.138 | 0.056 | 0.011 | 0.181 |
| Albania | \>= 5 | Outpatient | 0.044 | 0.010 | 0.141 | 0.057 | 0.011 | 0.188 |
| Algeria | \< 5 | Inpatient | 0.036 | 0.017 | 0.066 | 0.050 | 0.019 | 0.105 |
| Algeria | \>= 5 | Inpatient | 0.037 | 0.016 | 0.072 | 0.051 | 0.018 | 0.112 |
| Algeria | \< 5 | Outpatient | 0.036 | 0.018 | 0.062 | 0.050 | 0.019 | 0.103 |
| Algeria | \>= 5 | Outpatient | 0.037 | 0.017 | 0.067 | 0.051 | 0.019 | 0.110 |
| Andorra | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Andorra | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Andorra | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Andorra | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Angola | \< 5 | Inpatient | 0.036 | 0.017 | 0.066 | 0.047 | 0.020 | 0.092 |
| Angola | \>= 5 | Inpatient | 0.037 | 0.016 | 0.072 | 0.048 | 0.019 | 0.098 |
| Angola | \< 5 | Outpatient | 0.036 | 0.018 | 0.062 | 0.047 | 0.021 | 0.090 |
| Angola | \>= 5 | Outpatient | 0.037 | 0.017 | 0.067 | 0.048 | 0.020 | 0.096 |
| Antigua and Barbuda | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Antigua and Barbuda | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Antigua and Barbuda | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Antigua and Barbuda | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Argentina | \< 5 | Inpatient | 0.024 | 0.005 | 0.067 | 0.031 | 0.007 | 0.093 |
| Argentina | \>= 5 | Inpatient | 0.024 | 0.005 | 0.072 | 0.032 | 0.006 | 0.099 |
| Argentina | \< 5 | Outpatient | 0.024 | 0.006 | 0.064 | 0.031 | 0.007 | 0.093 |
| Argentina | \>= 5 | Outpatient | 0.024 | 0.005 | 0.068 | 0.032 | 0.007 | 0.096 |
| Armenia | \< 5 | Inpatient | 0.043 | 0.010 | 0.138 | 0.055 | 0.011 | 0.177 |
| Armenia | \>= 5 | Inpatient | 0.044 | 0.009 | 0.143 | 0.056 | 0.011 | 0.182 |
| Armenia | \< 5 | Outpatient | 0.043 | 0.010 | 0.138 | 0.055 | 0.011 | 0.178 |
| Armenia | \>= 5 | Outpatient | 0.044 | 0.010 | 0.141 | 0.057 | 0.011 | 0.185 |
| Australia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Australia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Australia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Australia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Austria | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Austria | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Austria | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Austria | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Azerbaijan | \< 5 | Inpatient | 0.043 | 0.010 | 0.138 | 0.057 | 0.013 | 0.179 |
| Azerbaijan | \>= 5 | Inpatient | 0.044 | 0.009 | 0.143 | 0.058 | 0.012 | 0.183 |
| Azerbaijan | \< 5 | Outpatient | 0.043 | 0.010 | 0.138 | 0.057 | 0.013 | 0.176 |
| Azerbaijan | \>= 5 | Outpatient | 0.044 | 0.010 | 0.141 | 0.058 | 0.013 | 0.180 |
| Bahamas, The | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahamas, The | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahamas, The | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahamas, The | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahrain | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahrain | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahrain | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bahrain | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bangladesh | \< 5 | Inpatient | 0.024 | 0.009 | 0.054 | 0.035 | 0.009 | 0.092 |
| Bangladesh | \>= 5 | Inpatient | 0.025 | 0.009 | 0.057 | 0.036 | 0.009 | 0.097 |
| Bangladesh | \< 5 | Outpatient | 0.024 | 0.010 | 0.051 | 0.035 | 0.010 | 0.089 |
| Bangladesh | \>= 5 | Outpatient | 0.025 | 0.009 | 0.054 | 0.036 | 0.009 | 0.094 |
| Barbados | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Barbados | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Barbados | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Barbados | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belarus | \< 5 | Inpatient | 0.043 | 0.010 | 0.138 | 0.057 | 0.013 | 0.179 |
| Belarus | \>= 5 | Inpatient | 0.044 | 0.009 | 0.143 | 0.058 | 0.012 | 0.183 |
| Belarus | \< 5 | Outpatient | 0.043 | 0.010 | 0.138 | 0.057 | 0.013 | 0.176 |
| Belarus | \>= 5 | Outpatient | 0.044 | 0.010 | 0.141 | 0.058 | 0.013 | 0.180 |
| Belgium | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belgium | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belgium | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belgium | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Belize | \< 5 | Inpatient | 0.023 | 0.008 | 0.052 | 0.034 | 0.008 | 0.094 |
| Belize | \>= 5 | Inpatient | 0.024 | 0.008 | 0.055 | 0.035 | 0.008 | 0.098 |
| Belize | \< 5 | Outpatient | 0.023 | 0.009 | 0.048 | 0.034 | 0.008 | 0.090 |
| Belize | \>= 5 | Outpatient | 0.024 | 0.009 | 0.051 | 0.035 | 0.008 | 0.095 |
| Benin | \< 5 | Inpatient | 0.036 | 0.017 | 0.066 | 0.047 | 0.020 | 0.092 |
| Benin | \>= 5 | Inpatient | 0.037 | 0.016 | 0.072 | 0.048 | 0.018 | 0.100 |
| Benin | \< 5 | Outpatient | 0.036 | 0.018 | 0.062 | 0.047 | 0.021 | 0.092 |
| Benin | \>= 5 | Outpatient | 0.037 | 0.017 | 0.067 | 0.048 | 0.019 | 0.098 |
| Bhutan | \< 5 | Inpatient | 0.022 | 0.009 | 0.045 | 0.032 | 0.009 | 0.078 |
| Bhutan | \>= 5 | Inpatient | 0.023 | 0.009 | 0.048 | 0.032 | 0.009 | 0.083 |
| Bhutan | \< 5 | Outpatient | 0.022 | 0.010 | 0.042 | 0.032 | 0.010 | 0.075 |
| Bhutan | \>= 5 | Outpatient | 0.023 | 0.009 | 0.045 | 0.032 | 0.009 | 0.079 |
| Bolivia | \< 5 | Inpatient | 0.024 | 0.007 | 0.061 | 0.033 | 0.009 | 0.087 |
| Bolivia | \>= 5 | Inpatient | 0.025 | 0.007 | 0.066 | 0.034 | 0.008 | 0.093 |
| Bolivia | \< 5 | Outpatient | 0.024 | 0.007 | 0.059 | 0.033 | 0.009 | 0.085 |
| Bolivia | \>= 5 | Outpatient | 0.025 | 0.007 | 0.064 | 0.034 | 0.009 | 0.091 |
| Bosnia and Herzegovina | \< 5 | Inpatient | 0.043 | 0.010 | 0.138 | 0.057 | 0.013 | 0.179 |
| Bosnia and Herzegovina | \>= 5 | Inpatient | 0.044 | 0.009 | 0.143 | 0.058 | 0.012 | 0.183 |
| Bosnia and Herzegovina | \< 5 | Outpatient | 0.043 | 0.010 | 0.138 | 0.057 | 0.013 | 0.176 |
| Bosnia and Herzegovina | \>= 5 | Outpatient | 0.044 | 0.010 | 0.141 | 0.058 | 0.013 | 0.180 |
| Botswana | \< 5 | Inpatient | 0.036 | 0.013 | 0.082 | 0.047 | 0.015 | 0.110 |
| Botswana | \>= 5 | Inpatient | 0.037 | 0.013 | 0.088 | 0.048 | 0.014 | 0.119 |
| Botswana | \< 5 | Outpatient | 0.036 | 0.014 | 0.079 | 0.047 | 0.015 | 0.110 |
| Botswana | \>= 5 | Outpatient | 0.037 | 0.013 | 0.085 | 0.048 | 0.015 | 0.118 |
| Brazil | \< 5 | Inpatient | 0.026 | 0.009 | 0.059 | 0.036 | 0.011 | 0.092 |
| Brazil | \>= 5 | Inpatient | 0.026 | 0.009 | 0.063 | 0.037 | 0.010 | 0.098 |
| Brazil | \< 5 | Outpatient | 0.026 | 0.010 | 0.057 | 0.036 | 0.011 | 0.090 |
| Brazil | \>= 5 | Outpatient | 0.026 | 0.009 | 0.060 | 0.037 | 0.010 | 0.096 |
| Brunei Darussalam | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Brunei Darussalam | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Brunei Darussalam | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Brunei Darussalam | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Bulgaria | \< 5 | Inpatient | 0.043 | 0.007 | 0.148 | 0.054 | 0.010 | 0.173 |
| Bulgaria | \>= 5 | Inpatient | 0.044 | 0.007 | 0.152 | 0.055 | 0.009 | 0.181 |
| Bulgaria | \< 5 | Outpatient | 0.043 | 0.007 | 0.147 | 0.054 | 0.010 | 0.176 |
| Bulgaria | \>= 5 | Outpatient | 0.044 | 0.007 | 0.151 | 0.055 | 0.010 | 0.182 |
| Burkina Faso | \< 5 | Inpatient | 0.029 | 0.006 | 0.068 | 0.038 | 0.007 | 0.103 |
| Burkina Faso | \>= 5 | Inpatient | 0.029 | 0.006 | 0.074 | 0.039 | 0.007 | 0.111 |
| Burkina Faso | \< 5 | Outpatient | 0.029 | 0.006 | 0.066 | 0.039 | 0.007 | 0.104 |
| Burkina Faso | \>= 5 | Outpatient | 0.029 | 0.006 | 0.072 | 0.040 | 0.007 | 0.111 |
| Burundi | \< 5 | Inpatient | 0.032 | 0.014 | 0.059 | 0.042 | 0.013 | 0.098 |
| Burundi | \>= 5 | Inpatient | 0.032 | 0.013 | 0.065 | 0.043 | 0.012 | 0.105 |
| Burundi | \< 5 | Outpatient | 0.032 | 0.015 | 0.057 | 0.043 | 0.013 | 0.100 |
| Burundi | \>= 5 | Outpatient | 0.032 | 0.014 | 0.063 | 0.044 | 0.013 | 0.106 |
| Cabo Verde | \< 5 | Inpatient | 0.037 | 0.010 | 0.086 | 0.050 | 0.013 | 0.130 |
| Cabo Verde | \>= 5 | Inpatient | 0.037 | 0.009 | 0.089 | 0.051 | 0.013 | 0.135 |
| Cabo Verde | \< 5 | Outpatient | 0.036 | 0.010 | 0.084 | 0.050 | 0.013 | 0.126 |
| Cabo Verde | \>= 5 | Outpatient | 0.037 | 0.010 | 0.085 | 0.051 | 0.013 | 0.132 |
| Cambodia | \< 5 | Inpatient | 0.016 | 0.003 | 0.041 | 0.022 | 0.003 | 0.067 |
| Cambodia | \>= 5 | Inpatient | 0.016 | 0.003 | 0.043 | 0.023 | 0.003 | 0.070 |
| Cambodia | \< 5 | Outpatient | 0.016 | 0.003 | 0.039 | 0.022 | 0.004 | 0.065 |
| Cambodia | \>= 5 | Outpatient | 0.016 | 0.003 | 0.040 | 0.022 | 0.004 | 0.067 |
| Cameroon | \< 5 | Inpatient | 0.037 | 0.011 | 0.088 | 0.048 | 0.014 | 0.116 |
| Cameroon | \>= 5 | Inpatient | 0.038 | 0.010 | 0.092 | 0.049 | 0.014 | 0.123 |
| Cameroon | \< 5 | Outpatient | 0.037 | 0.011 | 0.084 | 0.048 | 0.015 | 0.116 |
| Cameroon | \>= 5 | Outpatient | 0.038 | 0.011 | 0.088 | 0.049 | 0.014 | 0.122 |
| Canada | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Canada | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Canada | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Canada | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Central African Republic | \< 5 | Inpatient | 0.032 | 0.014 | 0.059 | 0.044 | 0.015 | 0.100 |
| Central African Republic | \>= 5 | Inpatient | 0.032 | 0.013 | 0.065 | 0.045 | 0.014 | 0.107 |
| Central African Republic | \< 5 | Outpatient | 0.032 | 0.015 | 0.057 | 0.044 | 0.016 | 0.096 |
| Central African Republic | \>= 5 | Outpatient | 0.032 | 0.014 | 0.063 | 0.045 | 0.015 | 0.103 |
| Chad | \< 5 | Inpatient | 0.032 | 0.014 | 0.059 | 0.044 | 0.015 | 0.100 |
| Chad | \>= 5 | Inpatient | 0.032 | 0.013 | 0.065 | 0.045 | 0.014 | 0.107 |
| Chad | \< 5 | Outpatient | 0.032 | 0.015 | 0.057 | 0.044 | 0.016 | 0.096 |
| Chad | \>= 5 | Outpatient | 0.032 | 0.014 | 0.063 | 0.045 | 0.015 | 0.103 |
| Chile | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Chile | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Chile | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Chile | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| China | \< 5 | Inpatient | 0.007 | 0.002 | 0.017 | 0.010 | 0.003 | 0.025 |
| China | \>= 5 | Inpatient | 0.007 | 0.002 | 0.017 | 0.010 | 0.003 | 0.025 |
| China | \< 5 | Outpatient | 0.007 | 0.002 | 0.015 | 0.010 | 0.003 | 0.023 |
| China | \>= 5 | Outpatient | 0.007 | 0.002 | 0.016 | 0.010 | 0.003 | 0.024 |
| Colombia | \< 5 | Inpatient | 0.022 | 0.007 | 0.052 | 0.030 | 0.008 | 0.076 |
| Colombia | \>= 5 | Inpatient | 0.022 | 0.006 | 0.054 | 0.031 | 0.008 | 0.083 |
| Colombia | \< 5 | Outpatient | 0.022 | 0.007 | 0.050 | 0.030 | 0.009 | 0.074 |
| Colombia | \>= 5 | Outpatient | 0.022 | 0.007 | 0.053 | 0.031 | 0.008 | 0.080 |
| Comoros | \< 5 | Inpatient | 0.036 | 0.017 | 0.066 | 0.050 | 0.019 | 0.105 |
| Comoros | \>= 5 | Inpatient | 0.037 | 0.016 | 0.072 | 0.051 | 0.018 | 0.112 |
| Comoros | \< 5 | Outpatient | 0.036 | 0.018 | 0.062 | 0.050 | 0.019 | 0.103 |
| Comoros | \>= 5 | Outpatient | 0.037 | 0.017 | 0.067 | 0.051 | 0.019 | 0.110 |
| Congo, Dem. Rep. | \< 5 | Inpatient | 0.032 | 0.014 | 0.059 | 0.043 | 0.016 | 0.090 |
| Congo, Dem. Rep. | \>= 5 | Inpatient | 0.032 | 0.013 | 0.065 | 0.044 | 0.015 | 0.096 |
| Congo, Dem. Rep. | \< 5 | Outpatient | 0.032 | 0.015 | 0.057 | 0.043 | 0.016 | 0.087 |
| Congo, Dem. Rep. | \>= 5 | Outpatient | 0.032 | 0.014 | 0.063 | 0.044 | 0.016 | 0.095 |
| Congo, Rep. | \< 5 | Inpatient | 0.036 | 0.017 | 0.066 | 0.047 | 0.020 | 0.092 |
| Congo, Rep. | \>= 5 | Inpatient | 0.037 | 0.016 | 0.072 | 0.048 | 0.019 | 0.099 |
| Congo, Rep. | \< 5 | Outpatient | 0.036 | 0.018 | 0.062 | 0.047 | 0.021 | 0.091 |
| Congo, Rep. | \>= 5 | Outpatient | 0.037 | 0.017 | 0.067 | 0.048 | 0.019 | 0.098 |
| Cook Islands | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Cook Islands | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Cook Islands | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Cook Islands | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Costa Rica | \< 5 | Inpatient | 0.023 | 0.008 | 0.052 | 0.030 | 0.010 | 0.072 |
| Costa Rica | \>= 5 | Inpatient | 0.024 | 0.008 | 0.055 | 0.031 | 0.009 | 0.077 |
| Costa Rica | \< 5 | Outpatient | 0.023 | 0.009 | 0.048 | 0.030 | 0.010 | 0.070 |
| Costa Rica | \>= 5 | Outpatient | 0.024 | 0.009 | 0.051 | 0.031 | 0.009 | 0.075 |
| Côte d’Ivoire | \< 5 | Inpatient | 0.036 | 0.017 | 0.066 | 0.047 | 0.020 | 0.093 |
| Côte d’Ivoire | \>= 5 | Inpatient | 0.037 | 0.016 | 0.072 | 0.048 | 0.018 | 0.101 |
| Côte d’Ivoire | \< 5 | Outpatient | 0.036 | 0.018 | 0.062 | 0.047 | 0.020 | 0.093 |
| Côte d’Ivoire | \>= 5 | Outpatient | 0.037 | 0.017 | 0.067 | 0.048 | 0.019 | 0.100 |
| Croatia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Croatia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Croatia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Croatia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Cuba | \< 5 | Inpatient | 0.023 | 0.008 | 0.052 | 0.034 | 0.008 | 0.094 |
| Cuba | \>= 5 | Inpatient | 0.024 | 0.008 | 0.055 | 0.035 | 0.008 | 0.098 |
| Cuba | \< 5 | Outpatient | 0.023 | 0.009 | 0.048 | 0.034 | 0.008 | 0.090 |
| Cuba | \>= 5 | Outpatient | 0.024 | 0.009 | 0.051 | 0.035 | 0.008 | 0.095 |
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
| Djibouti | \< 5 | Inpatient | 0.042 | 0.015 | 0.093 | 0.055 | 0.017 | 0.132 |
| Djibouti | \>= 5 | Inpatient | 0.043 | 0.015 | 0.099 | 0.056 | 0.016 | 0.142 |
| Djibouti | \< 5 | Outpatient | 0.042 | 0.016 | 0.090 | 0.055 | 0.017 | 0.131 |
| Djibouti | \>= 5 | Outpatient | 0.042 | 0.016 | 0.094 | 0.056 | 0.017 | 0.140 |
| Dominica | \< 5 | Inpatient | 0.023 | 0.008 | 0.052 | 0.034 | 0.008 | 0.094 |
| Dominica | \>= 5 | Inpatient | 0.024 | 0.008 | 0.055 | 0.035 | 0.008 | 0.098 |
| Dominica | \< 5 | Outpatient | 0.023 | 0.009 | 0.048 | 0.034 | 0.008 | 0.090 |
| Dominica | \>= 5 | Outpatient | 0.024 | 0.009 | 0.051 | 0.035 | 0.008 | 0.095 |
| Dominican Republic | \< 5 | Inpatient | 0.023 | 0.008 | 0.052 | 0.030 | 0.010 | 0.071 |
| Dominican Republic | \>= 5 | Inpatient | 0.024 | 0.008 | 0.055 | 0.031 | 0.009 | 0.076 |
| Dominican Republic | \< 5 | Outpatient | 0.023 | 0.009 | 0.048 | 0.030 | 0.011 | 0.068 |
| Dominican Republic | \>= 5 | Outpatient | 0.024 | 0.009 | 0.051 | 0.031 | 0.010 | 0.074 |
| Ecuador | \< 5 | Inpatient | 0.022 | 0.008 | 0.048 | 0.030 | 0.010 | 0.070 |
| Ecuador | \>= 5 | Inpatient | 0.023 | 0.008 | 0.051 | 0.031 | 0.009 | 0.076 |
| Ecuador | \< 5 | Outpatient | 0.022 | 0.008 | 0.046 | 0.030 | 0.010 | 0.068 |
| Ecuador | \>= 5 | Outpatient | 0.022 | 0.008 | 0.049 | 0.031 | 0.010 | 0.074 |
| Egypt, Arab Rep. | \< 5 | Inpatient | 0.042 | 0.015 | 0.093 | 0.057 | 0.018 | 0.135 |
| Egypt, Arab Rep. | \>= 5 | Inpatient | 0.043 | 0.015 | 0.099 | 0.058 | 0.017 | 0.142 |
| Egypt, Arab Rep. | \< 5 | Outpatient | 0.042 | 0.016 | 0.090 | 0.057 | 0.019 | 0.130 |
| Egypt, Arab Rep. | \>= 5 | Outpatient | 0.042 | 0.016 | 0.094 | 0.058 | 0.018 | 0.137 |
| El Salvador | \< 5 | Inpatient | 0.022 | 0.008 | 0.047 | 0.030 | 0.010 | 0.071 |
| El Salvador | \>= 5 | Inpatient | 0.022 | 0.008 | 0.050 | 0.031 | 0.009 | 0.076 |
| El Salvador | \< 5 | Outpatient | 0.022 | 0.009 | 0.045 | 0.030 | 0.011 | 0.068 |
| El Salvador | \>= 5 | Outpatient | 0.022 | 0.008 | 0.048 | 0.031 | 0.010 | 0.074 |
| Equatorial Guinea | \< 5 | Inpatient | 0.036 | 0.013 | 0.082 | 0.050 | 0.015 | 0.119 |
| Equatorial Guinea | \>= 5 | Inpatient | 0.037 | 0.013 | 0.088 | 0.051 | 0.015 | 0.126 |
| Equatorial Guinea | \< 5 | Outpatient | 0.036 | 0.014 | 0.079 | 0.050 | 0.016 | 0.118 |
| Equatorial Guinea | \>= 5 | Outpatient | 0.037 | 0.013 | 0.085 | 0.051 | 0.015 | 0.123 |
| Eritrea | \< 5 | Inpatient | 0.032 | 0.014 | 0.059 | 0.042 | 0.013 | 0.098 |
| Eritrea | \>= 5 | Inpatient | 0.032 | 0.013 | 0.065 | 0.043 | 0.012 | 0.105 |
| Eritrea | \< 5 | Outpatient | 0.032 | 0.015 | 0.057 | 0.043 | 0.013 | 0.100 |
| Eritrea | \>= 5 | Outpatient | 0.032 | 0.014 | 0.063 | 0.044 | 0.013 | 0.106 |
| Estonia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Estonia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Estonia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Estonia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Eswatini | \< 5 | Inpatient | 0.036 | 0.017 | 0.066 | 0.047 | 0.020 | 0.092 |
| Eswatini | \>= 5 | Inpatient | 0.037 | 0.016 | 0.072 | 0.048 | 0.018 | 0.100 |
| Eswatini | \< 5 | Outpatient | 0.036 | 0.018 | 0.062 | 0.047 | 0.021 | 0.092 |
| Eswatini | \>= 5 | Outpatient | 0.037 | 0.017 | 0.067 | 0.048 | 0.019 | 0.098 |
| Ethiopia | \< 5 | Inpatient | 0.032 | 0.014 | 0.059 | 0.042 | 0.015 | 0.090 |
| Ethiopia | \>= 5 | Inpatient | 0.032 | 0.013 | 0.065 | 0.043 | 0.014 | 0.098 |
| Ethiopia | \< 5 | Outpatient | 0.032 | 0.015 | 0.057 | 0.042 | 0.015 | 0.090 |
| Ethiopia | \>= 5 | Outpatient | 0.032 | 0.014 | 0.063 | 0.043 | 0.014 | 0.097 |
| Fiji | \< 5 | Inpatient | 0.012 | 0.002 | 0.034 | 0.016 | 0.003 | 0.056 |
| Fiji | \>= 5 | Inpatient | 0.012 | 0.002 | 0.034 | 0.016 | 0.003 | 0.057 |
| Fiji | \< 5 | Outpatient | 0.012 | 0.003 | 0.032 | 0.016 | 0.003 | 0.053 |
| Fiji | \>= 5 | Outpatient | 0.012 | 0.003 | 0.033 | 0.016 | 0.003 | 0.055 |
| Finland | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Finland | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Finland | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Finland | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| France | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| France | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| France | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| France | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Gabon | \< 5 | Inpatient | 0.048 | 0.013 | 0.136 | 0.064 | 0.016 | 0.184 |
| Gabon | \>= 5 | Inpatient | 0.048 | 0.012 | 0.139 | 0.066 | 0.016 | 0.191 |
| Gabon | \< 5 | Outpatient | 0.047 | 0.013 | 0.134 | 0.064 | 0.016 | 0.183 |
| Gabon | \>= 5 | Outpatient | 0.048 | 0.013 | 0.135 | 0.065 | 0.016 | 0.186 |
| Gambia, The | \< 5 | Inpatient | 0.038 | 0.014 | 0.087 | 0.051 | 0.014 | 0.137 |
| Gambia, The | \>= 5 | Inpatient | 0.039 | 0.013 | 0.094 | 0.053 | 0.013 | 0.146 |
| Gambia, The | \< 5 | Outpatient | 0.038 | 0.014 | 0.087 | 0.052 | 0.014 | 0.139 |
| Gambia, The | \>= 5 | Outpatient | 0.039 | 0.013 | 0.090 | 0.053 | 0.013 | 0.147 |
| Georgia | \< 5 | Inpatient | 0.043 | 0.010 | 0.138 | 0.055 | 0.012 | 0.171 |
| Georgia | \>= 5 | Inpatient | 0.044 | 0.009 | 0.143 | 0.056 | 0.012 | 0.176 |
| Georgia | \< 5 | Outpatient | 0.043 | 0.010 | 0.138 | 0.055 | 0.012 | 0.172 |
| Georgia | \>= 5 | Outpatient | 0.044 | 0.010 | 0.141 | 0.056 | 0.012 | 0.177 |
| Germany | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Germany | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Germany | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Germany | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Ghana | \< 5 | Inpatient | 0.042 | 0.015 | 0.098 | 0.054 | 0.017 | 0.139 |
| Ghana | \>= 5 | Inpatient | 0.043 | 0.014 | 0.104 | 0.056 | 0.015 | 0.145 |
| Ghana | \< 5 | Outpatient | 0.041 | 0.015 | 0.093 | 0.055 | 0.017 | 0.137 |
| Ghana | \>= 5 | Outpatient | 0.042 | 0.015 | 0.096 | 0.056 | 0.016 | 0.144 |
| Greece | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Greece | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Greece | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Greece | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Grenada | \< 5 | Inpatient | 0.023 | 0.008 | 0.052 | 0.034 | 0.008 | 0.094 |
| Grenada | \>= 5 | Inpatient | 0.024 | 0.008 | 0.055 | 0.035 | 0.008 | 0.098 |
| Grenada | \< 5 | Outpatient | 0.023 | 0.009 | 0.048 | 0.034 | 0.008 | 0.090 |
| Grenada | \>= 5 | Outpatient | 0.024 | 0.009 | 0.051 | 0.035 | 0.008 | 0.095 |
| Guatemala | \< 5 | Inpatient | 0.022 | 0.006 | 0.053 | 0.030 | 0.007 | 0.082 |
| Guatemala | \>= 5 | Inpatient | 0.022 | 0.006 | 0.056 | 0.030 | 0.007 | 0.087 |
| Guatemala | \< 5 | Outpatient | 0.022 | 0.006 | 0.051 | 0.030 | 0.007 | 0.079 |
| Guatemala | \>= 5 | Outpatient | 0.022 | 0.006 | 0.054 | 0.030 | 0.007 | 0.085 |
| Guinea | \< 5 | Inpatient | 0.036 | 0.017 | 0.066 | 0.050 | 0.019 | 0.105 |
| Guinea | \>= 5 | Inpatient | 0.037 | 0.016 | 0.072 | 0.051 | 0.018 | 0.112 |
| Guinea | \< 5 | Outpatient | 0.036 | 0.018 | 0.062 | 0.050 | 0.019 | 0.103 |
| Guinea | \>= 5 | Outpatient | 0.037 | 0.017 | 0.067 | 0.051 | 0.019 | 0.110 |
| Guinea-Bissau | \< 5 | Inpatient | 0.033 | 0.009 | 0.078 | 0.043 | 0.011 | 0.115 |
| Guinea-Bissau | \>= 5 | Inpatient | 0.033 | 0.009 | 0.083 | 0.044 | 0.010 | 0.125 |
| Guinea-Bissau | \< 5 | Outpatient | 0.033 | 0.009 | 0.075 | 0.043 | 0.011 | 0.116 |
| Guinea-Bissau | \>= 5 | Outpatient | 0.033 | 0.009 | 0.080 | 0.045 | 0.010 | 0.124 |
| Guyana | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Guyana | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Guyana | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Guyana | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Haiti | \< 5 | Inpatient | 0.025 | 0.007 | 0.064 | 0.034 | 0.009 | 0.089 |
| Haiti | \>= 5 | Inpatient | 0.026 | 0.007 | 0.068 | 0.035 | 0.008 | 0.097 |
| Haiti | \< 5 | Outpatient | 0.025 | 0.008 | 0.061 | 0.034 | 0.009 | 0.088 |
| Haiti | \>= 5 | Outpatient | 0.026 | 0.008 | 0.064 | 0.035 | 0.009 | 0.093 |
| Honduras | \< 5 | Inpatient | 0.024 | 0.007 | 0.063 | 0.033 | 0.009 | 0.088 |
| Honduras | \>= 5 | Inpatient | 0.025 | 0.006 | 0.068 | 0.034 | 0.008 | 0.094 |
| Honduras | \< 5 | Outpatient | 0.024 | 0.007 | 0.062 | 0.033 | 0.009 | 0.085 |
| Honduras | \>= 5 | Outpatient | 0.025 | 0.007 | 0.066 | 0.034 | 0.009 | 0.092 |
| Hungary | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Hungary | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Hungary | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Hungary | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Iceland | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Iceland | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Iceland | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Iceland | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| India | \< 5 | Inpatient | 0.025 | 0.011 | 0.050 | 0.035 | 0.010 | 0.088 |
| India | \>= 5 | Inpatient | 0.026 | 0.010 | 0.054 | 0.035 | 0.009 | 0.093 |
| India | \< 5 | Outpatient | 0.025 | 0.011 | 0.050 | 0.035 | 0.010 | 0.089 |
| India | \>= 5 | Outpatient | 0.026 | 0.011 | 0.052 | 0.036 | 0.009 | 0.095 |
| Indonesia | \< 5 | Inpatient | 0.025 | 0.005 | 0.080 | 0.037 | 0.006 | 0.130 |
| Indonesia | \>= 5 | Inpatient | 0.026 | 0.005 | 0.082 | 0.037 | 0.006 | 0.136 |
| Indonesia | \< 5 | Outpatient | 0.025 | 0.005 | 0.078 | 0.037 | 0.006 | 0.128 |
| Indonesia | \>= 5 | Outpatient | 0.026 | 0.005 | 0.081 | 0.037 | 0.006 | 0.131 |
| Iran, Islamic Rep. | \< 5 | Inpatient | 0.063 | 0.022 | 0.147 | 0.086 | 0.024 | 0.219 |
| Iran, Islamic Rep. | \>= 5 | Inpatient | 0.064 | 0.021 | 0.152 | 0.087 | 0.024 | 0.226 |
| Iran, Islamic Rep. | \< 5 | Outpatient | 0.063 | 0.023 | 0.142 | 0.086 | 0.025 | 0.216 |
| Iran, Islamic Rep. | \>= 5 | Outpatient | 0.063 | 0.022 | 0.145 | 0.087 | 0.025 | 0.222 |
| Iraq | \< 5 | Inpatient | 0.042 | 0.015 | 0.093 | 0.055 | 0.019 | 0.127 |
| Iraq | \>= 5 | Inpatient | 0.043 | 0.015 | 0.099 | 0.056 | 0.018 | 0.136 |
| Iraq | \< 5 | Outpatient | 0.042 | 0.016 | 0.090 | 0.055 | 0.019 | 0.124 |
| Iraq | \>= 5 | Outpatient | 0.042 | 0.016 | 0.094 | 0.056 | 0.018 | 0.133 |
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
| Jamaica | \< 5 | Inpatient | 0.023 | 0.008 | 0.052 | 0.034 | 0.008 | 0.094 |
| Jamaica | \>= 5 | Inpatient | 0.024 | 0.008 | 0.055 | 0.035 | 0.008 | 0.098 |
| Jamaica | \< 5 | Outpatient | 0.023 | 0.009 | 0.048 | 0.034 | 0.008 | 0.090 |
| Jamaica | \>= 5 | Outpatient | 0.024 | 0.009 | 0.051 | 0.035 | 0.008 | 0.095 |
| Japan | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Japan | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Japan | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Japan | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Jordan | \< 5 | Inpatient | 0.042 | 0.015 | 0.093 | 0.055 | 0.017 | 0.133 |
| Jordan | \>= 5 | Inpatient | 0.043 | 0.015 | 0.099 | 0.056 | 0.016 | 0.143 |
| Jordan | \< 5 | Outpatient | 0.042 | 0.016 | 0.090 | 0.055 | 0.017 | 0.134 |
| Jordan | \>= 5 | Outpatient | 0.042 | 0.016 | 0.094 | 0.056 | 0.016 | 0.140 |
| Kazakhstan | \< 5 | Inpatient | 0.043 | 0.010 | 0.138 | 0.057 | 0.013 | 0.179 |
| Kazakhstan | \>= 5 | Inpatient | 0.044 | 0.009 | 0.143 | 0.058 | 0.012 | 0.183 |
| Kazakhstan | \< 5 | Outpatient | 0.043 | 0.010 | 0.138 | 0.057 | 0.013 | 0.176 |
| Kazakhstan | \>= 5 | Outpatient | 0.044 | 0.010 | 0.141 | 0.058 | 0.013 | 0.180 |
| Kenya | \< 5 | Inpatient | 0.037 | 0.016 | 0.073 | 0.050 | 0.017 | 0.117 |
| Kenya | \>= 5 | Inpatient | 0.038 | 0.015 | 0.077 | 0.051 | 0.016 | 0.124 |
| Kenya | \< 5 | Outpatient | 0.037 | 0.017 | 0.070 | 0.050 | 0.016 | 0.118 |
| Kenya | \>= 5 | Outpatient | 0.038 | 0.016 | 0.073 | 0.051 | 0.016 | 0.124 |
| Kiribati | \< 5 | Inpatient | 0.016 | 0.003 | 0.041 | 0.021 | 0.003 | 0.064 |
| Kiribati | \>= 5 | Inpatient | 0.016 | 0.003 | 0.043 | 0.022 | 0.003 | 0.066 |
| Kiribati | \< 5 | Outpatient | 0.016 | 0.003 | 0.039 | 0.021 | 0.003 | 0.063 |
| Kiribati | \>= 5 | Outpatient | 0.016 | 0.003 | 0.040 | 0.022 | 0.003 | 0.066 |
| Korea, Dem. People’s Rep. | \< 5 | Inpatient | 0.022 | 0.009 | 0.045 | 0.032 | 0.009 | 0.078 |
| Korea, Dem. People’s Rep. | \>= 5 | Inpatient | 0.023 | 0.009 | 0.048 | 0.032 | 0.009 | 0.083 |
| Korea, Dem. People’s Rep. | \< 5 | Outpatient | 0.022 | 0.010 | 0.042 | 0.032 | 0.010 | 0.075 |
| Korea, Dem. People’s Rep. | \>= 5 | Outpatient | 0.023 | 0.009 | 0.045 | 0.032 | 0.009 | 0.079 |
| Korea, Rep. | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Korea, Rep. | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Korea, Rep. | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Korea, Rep. | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kuwait | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kuwait | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kuwait | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kuwait | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Kyrgyz Republic | \< 5 | Inpatient | 0.039 | 0.009 | 0.128 | 0.049 | 0.012 | 0.150 |
| Kyrgyz Republic | \>= 5 | Inpatient | 0.039 | 0.009 | 0.127 | 0.050 | 0.011 | 0.156 |
| Kyrgyz Republic | \< 5 | Outpatient | 0.039 | 0.010 | 0.125 | 0.049 | 0.012 | 0.154 |
| Kyrgyz Republic | \>= 5 | Outpatient | 0.039 | 0.010 | 0.129 | 0.050 | 0.012 | 0.160 |
| Lao PDR | \< 5 | Inpatient | 0.016 | 0.003 | 0.041 | 0.022 | 0.003 | 0.067 |
| Lao PDR | \>= 5 | Inpatient | 0.016 | 0.003 | 0.043 | 0.023 | 0.003 | 0.070 |
| Lao PDR | \< 5 | Outpatient | 0.016 | 0.003 | 0.039 | 0.022 | 0.004 | 0.065 |
| Lao PDR | \>= 5 | Outpatient | 0.016 | 0.003 | 0.040 | 0.022 | 0.004 | 0.067 |
| Latvia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Latvia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Latvia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Latvia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Lebanon | \< 5 | Inpatient | 0.045 | 0.012 | 0.116 | 0.061 | 0.016 | 0.160 |
| Lebanon | \>= 5 | Inpatient | 0.046 | 0.011 | 0.120 | 0.062 | 0.015 | 0.162 |
| Lebanon | \< 5 | Outpatient | 0.045 | 0.012 | 0.113 | 0.061 | 0.016 | 0.156 |
| Lebanon | \>= 5 | Outpatient | 0.046 | 0.012 | 0.115 | 0.062 | 0.016 | 0.157 |
| Lesotho | \< 5 | Inpatient | 0.036 | 0.017 | 0.066 | 0.047 | 0.018 | 0.098 |
| Lesotho | \>= 5 | Inpatient | 0.037 | 0.016 | 0.072 | 0.048 | 0.017 | 0.106 |
| Lesotho | \< 5 | Outpatient | 0.036 | 0.018 | 0.062 | 0.047 | 0.018 | 0.100 |
| Lesotho | \>= 5 | Outpatient | 0.037 | 0.017 | 0.067 | 0.048 | 0.017 | 0.106 |
| Liberia | \< 5 | Inpatient | 0.032 | 0.014 | 0.059 | 0.042 | 0.015 | 0.089 |
| Liberia | \>= 5 | Inpatient | 0.032 | 0.013 | 0.065 | 0.043 | 0.014 | 0.097 |
| Liberia | \< 5 | Outpatient | 0.032 | 0.015 | 0.057 | 0.042 | 0.016 | 0.088 |
| Liberia | \>= 5 | Outpatient | 0.032 | 0.014 | 0.063 | 0.043 | 0.015 | 0.097 |
| Libya | \< 5 | Inpatient | 0.042 | 0.015 | 0.093 | 0.055 | 0.017 | 0.133 |
| Libya | \>= 5 | Inpatient | 0.043 | 0.015 | 0.099 | 0.056 | 0.016 | 0.142 |
| Libya | \< 5 | Outpatient | 0.042 | 0.016 | 0.090 | 0.055 | 0.017 | 0.132 |
| Libya | \>= 5 | Outpatient | 0.042 | 0.016 | 0.094 | 0.056 | 0.016 | 0.140 |
| Lithuania | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Lithuania | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Lithuania | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Lithuania | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Luxembourg | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Luxembourg | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Luxembourg | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Luxembourg | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Madagascar | \< 5 | Inpatient | 0.032 | 0.014 | 0.059 | 0.042 | 0.015 | 0.089 |
| Madagascar | \>= 5 | Inpatient | 0.032 | 0.013 | 0.065 | 0.043 | 0.014 | 0.097 |
| Madagascar | \< 5 | Outpatient | 0.032 | 0.015 | 0.057 | 0.042 | 0.016 | 0.088 |
| Madagascar | \>= 5 | Outpatient | 0.032 | 0.014 | 0.063 | 0.043 | 0.015 | 0.097 |
| Malawi | \< 5 | Inpatient | 0.029 | 0.007 | 0.066 | 0.038 | 0.009 | 0.097 |
| Malawi | \>= 5 | Inpatient | 0.030 | 0.007 | 0.072 | 0.039 | 0.008 | 0.104 |
| Malawi | \< 5 | Outpatient | 0.029 | 0.007 | 0.066 | 0.038 | 0.009 | 0.099 |
| Malawi | \>= 5 | Outpatient | 0.030 | 0.007 | 0.071 | 0.039 | 0.008 | 0.105 |
| Malaysia | \< 5 | Inpatient | 0.012 | 0.002 | 0.034 | 0.017 | 0.003 | 0.056 |
| Malaysia | \>= 5 | Inpatient | 0.012 | 0.002 | 0.034 | 0.017 | 0.003 | 0.056 |
| Malaysia | \< 5 | Outpatient | 0.012 | 0.003 | 0.032 | 0.017 | 0.003 | 0.055 |
| Malaysia | \>= 5 | Outpatient | 0.012 | 0.003 | 0.033 | 0.017 | 0.003 | 0.054 |
| Maldives | \< 5 | Inpatient | 0.020 | 0.005 | 0.045 | 0.029 | 0.006 | 0.079 |
| Maldives | \>= 5 | Inpatient | 0.020 | 0.006 | 0.048 | 0.029 | 0.006 | 0.082 |
| Maldives | \< 5 | Outpatient | 0.020 | 0.006 | 0.044 | 0.029 | 0.006 | 0.077 |
| Maldives | \>= 5 | Outpatient | 0.020 | 0.006 | 0.045 | 0.029 | 0.006 | 0.079 |
| Mali | \< 5 | Inpatient | 0.034 | 0.012 | 0.079 | 0.046 | 0.013 | 0.118 |
| Mali | \>= 5 | Inpatient | 0.035 | 0.011 | 0.083 | 0.047 | 0.012 | 0.125 |
| Mali | \< 5 | Outpatient | 0.034 | 0.012 | 0.076 | 0.046 | 0.013 | 0.118 |
| Mali | \>= 5 | Outpatient | 0.035 | 0.012 | 0.081 | 0.047 | 0.012 | 0.126 |
| Malta | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Malta | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Malta | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Malta | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Marshall Islands | \< 5 | Inpatient | 0.011 | 0.002 | 0.033 | 0.016 | 0.003 | 0.051 |
| Marshall Islands | \>= 5 | Inpatient | 0.012 | 0.002 | 0.034 | 0.016 | 0.003 | 0.052 |
| Marshall Islands | \< 5 | Outpatient | 0.011 | 0.003 | 0.031 | 0.016 | 0.003 | 0.050 |
| Marshall Islands | \>= 5 | Outpatient | 0.012 | 0.003 | 0.032 | 0.016 | 0.003 | 0.051 |
| Mauritania | \< 5 | Inpatient | 0.036 | 0.017 | 0.066 | 0.047 | 0.018 | 0.098 |
| Mauritania | \>= 5 | Inpatient | 0.037 | 0.016 | 0.072 | 0.048 | 0.017 | 0.106 |
| Mauritania | \< 5 | Outpatient | 0.036 | 0.018 | 0.062 | 0.047 | 0.019 | 0.098 |
| Mauritania | \>= 5 | Outpatient | 0.037 | 0.017 | 0.067 | 0.048 | 0.018 | 0.105 |
| Mauritius | \< 5 | Inpatient | 0.036 | 0.013 | 0.082 | 0.047 | 0.015 | 0.113 |
| Mauritius | \>= 5 | Inpatient | 0.037 | 0.013 | 0.088 | 0.048 | 0.014 | 0.121 |
| Mauritius | \< 5 | Outpatient | 0.036 | 0.014 | 0.079 | 0.047 | 0.015 | 0.111 |
| Mauritius | \>= 5 | Outpatient | 0.037 | 0.013 | 0.085 | 0.048 | 0.014 | 0.120 |
| Mexico | \< 5 | Inpatient | 0.031 | 0.008 | 0.100 | 0.044 | 0.010 | 0.147 |
| Mexico | \>= 5 | Inpatient | 0.032 | 0.008 | 0.104 | 0.045 | 0.009 | 0.153 |
| Mexico | \< 5 | Outpatient | 0.031 | 0.008 | 0.098 | 0.044 | 0.010 | 0.144 |
| Mexico | \>= 5 | Outpatient | 0.032 | 0.008 | 0.103 | 0.044 | 0.010 | 0.149 |
| Micronesia, Fed. Sts. | \< 5 | Inpatient | 0.016 | 0.002 | 0.044 | 0.021 | 0.003 | 0.061 |
| Micronesia, Fed. Sts. | \>= 5 | Inpatient | 0.016 | 0.002 | 0.047 | 0.022 | 0.003 | 0.064 |
| Micronesia, Fed. Sts. | \< 5 | Outpatient | 0.016 | 0.002 | 0.043 | 0.021 | 0.004 | 0.059 |
| Micronesia, Fed. Sts. | \>= 5 | Outpatient | 0.016 | 0.002 | 0.045 | 0.022 | 0.004 | 0.061 |
| Moldova | \< 5 | Inpatient | 0.043 | 0.010 | 0.138 | 0.055 | 0.013 | 0.166 |
| Moldova | \>= 5 | Inpatient | 0.044 | 0.009 | 0.143 | 0.056 | 0.012 | 0.173 |
| Moldova | \< 5 | Outpatient | 0.043 | 0.010 | 0.138 | 0.055 | 0.013 | 0.168 |
| Moldova | \>= 5 | Outpatient | 0.044 | 0.010 | 0.141 | 0.056 | 0.013 | 0.173 |
| Monaco | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Monaco | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Monaco | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Monaco | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Mongolia | \< 5 | Inpatient | 0.016 | 0.003 | 0.041 | 0.022 | 0.003 | 0.067 |
| Mongolia | \>= 5 | Inpatient | 0.016 | 0.003 | 0.043 | 0.023 | 0.003 | 0.070 |
| Mongolia | \< 5 | Outpatient | 0.016 | 0.003 | 0.039 | 0.022 | 0.004 | 0.065 |
| Mongolia | \>= 5 | Outpatient | 0.016 | 0.003 | 0.040 | 0.022 | 0.004 | 0.067 |
| Montenegro | \< 5 | Inpatient | 0.043 | 0.010 | 0.138 | 0.057 | 0.013 | 0.179 |
| Montenegro | \>= 5 | Inpatient | 0.044 | 0.009 | 0.143 | 0.058 | 0.012 | 0.183 |
| Montenegro | \< 5 | Outpatient | 0.043 | 0.010 | 0.138 | 0.057 | 0.013 | 0.176 |
| Montenegro | \>= 5 | Outpatient | 0.044 | 0.010 | 0.141 | 0.058 | 0.013 | 0.180 |
| Morocco | \< 5 | Inpatient | 0.042 | 0.015 | 0.093 | 0.056 | 0.015 | 0.147 |
| Morocco | \>= 5 | Inpatient | 0.043 | 0.015 | 0.099 | 0.057 | 0.014 | 0.155 |
| Morocco | \< 5 | Outpatient | 0.042 | 0.016 | 0.090 | 0.056 | 0.015 | 0.144 |
| Morocco | \>= 5 | Outpatient | 0.042 | 0.016 | 0.094 | 0.057 | 0.014 | 0.153 |
| Mozambique | \< 5 | Inpatient | 0.037 | 0.013 | 0.085 | 0.050 | 0.015 | 0.127 |
| Mozambique | \>= 5 | Inpatient | 0.038 | 0.013 | 0.091 | 0.051 | 0.014 | 0.135 |
| Mozambique | \< 5 | Outpatient | 0.037 | 0.014 | 0.083 | 0.050 | 0.015 | 0.128 |
| Mozambique | \>= 5 | Outpatient | 0.038 | 0.013 | 0.089 | 0.051 | 0.014 | 0.136 |
| Myanmar | \< 5 | Inpatient | 0.022 | 0.009 | 0.045 | 0.030 | 0.010 | 0.071 |
| Myanmar | \>= 5 | Inpatient | 0.023 | 0.009 | 0.048 | 0.031 | 0.009 | 0.076 |
| Myanmar | \< 5 | Outpatient | 0.022 | 0.010 | 0.042 | 0.030 | 0.010 | 0.069 |
| Myanmar | \>= 5 | Outpatient | 0.023 | 0.009 | 0.045 | 0.031 | 0.010 | 0.074 |
| Namibia | \< 5 | Inpatient | 0.036 | 0.013 | 0.082 | 0.047 | 0.016 | 0.108 |
| Namibia | \>= 5 | Inpatient | 0.037 | 0.013 | 0.088 | 0.048 | 0.015 | 0.117 |
| Namibia | \< 5 | Outpatient | 0.036 | 0.014 | 0.079 | 0.047 | 0.016 | 0.107 |
| Namibia | \>= 5 | Outpatient | 0.037 | 0.013 | 0.085 | 0.048 | 0.015 | 0.115 |
| Nauru | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nauru | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nauru | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nauru | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nepal | \< 5 | Inpatient | 0.020 | 0.005 | 0.046 | 0.027 | 0.006 | 0.074 |
| Nepal | \>= 5 | Inpatient | 0.020 | 0.005 | 0.049 | 0.028 | 0.006 | 0.078 |
| Nepal | \< 5 | Outpatient | 0.020 | 0.006 | 0.045 | 0.028 | 0.007 | 0.071 |
| Nepal | \>= 5 | Outpatient | 0.020 | 0.005 | 0.047 | 0.028 | 0.006 | 0.075 |
| Netherlands | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Netherlands | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Netherlands | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Netherlands | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| New Zealand | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| New Zealand | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| New Zealand | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| New Zealand | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Nicaragua | \< 5 | Inpatient | 0.026 | 0.007 | 0.072 | 0.036 | 0.008 | 0.105 |
| Nicaragua | \>= 5 | Inpatient | 0.027 | 0.006 | 0.078 | 0.037 | 0.008 | 0.111 |
| Nicaragua | \< 5 | Outpatient | 0.026 | 0.007 | 0.071 | 0.036 | 0.009 | 0.103 |
| Nicaragua | \>= 5 | Outpatient | 0.027 | 0.007 | 0.075 | 0.037 | 0.008 | 0.109 |
| Niger | \< 5 | Inpatient | 0.032 | 0.014 | 0.059 | 0.042 | 0.014 | 0.093 |
| Niger | \>= 5 | Inpatient | 0.032 | 0.013 | 0.065 | 0.043 | 0.013 | 0.101 |
| Niger | \< 5 | Outpatient | 0.032 | 0.015 | 0.057 | 0.042 | 0.014 | 0.093 |
| Niger | \>= 5 | Outpatient | 0.032 | 0.014 | 0.063 | 0.043 | 0.014 | 0.102 |
| Nigeria | \< 5 | Inpatient | 0.056 | 0.020 | 0.148 | 0.078 | 0.021 | 0.226 |
| Nigeria | \>= 5 | Inpatient | 0.057 | 0.019 | 0.153 | 0.079 | 0.020 | 0.233 |
| Nigeria | \< 5 | Outpatient | 0.056 | 0.021 | 0.145 | 0.078 | 0.022 | 0.226 |
| Nigeria | \>= 5 | Outpatient | 0.057 | 0.020 | 0.150 | 0.079 | 0.021 | 0.233 |
| Niue | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Niue | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Niue | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Niue | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| North Macedonia | \< 5 | Inpatient | 0.043 | 0.010 | 0.138 | 0.055 | 0.013 | 0.167 |
| North Macedonia | \>= 5 | Inpatient | 0.044 | 0.009 | 0.143 | 0.056 | 0.012 | 0.173 |
| North Macedonia | \< 5 | Outpatient | 0.043 | 0.010 | 0.138 | 0.055 | 0.013 | 0.169 |
| North Macedonia | \>= 5 | Outpatient | 0.044 | 0.010 | 0.141 | 0.056 | 0.013 | 0.175 |
| Norway | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Norway | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Norway | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Norway | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Oman | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Oman | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Oman | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Oman | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Pakistan | \< 5 | Inpatient | 0.035 | 0.010 | 0.083 | 0.047 | 0.012 | 0.125 |
| Pakistan | \>= 5 | Inpatient | 0.036 | 0.010 | 0.089 | 0.048 | 0.011 | 0.134 |
| Pakistan | \< 5 | Outpatient | 0.035 | 0.011 | 0.080 | 0.047 | 0.012 | 0.124 |
| Pakistan | \>= 5 | Outpatient | 0.036 | 0.011 | 0.085 | 0.048 | 0.011 | 0.133 |
| Palau | \< 5 | Inpatient | 0.012 | 0.002 | 0.037 | 0.016 | 0.003 | 0.053 |
| Palau | \>= 5 | Inpatient | 0.012 | 0.002 | 0.038 | 0.016 | 0.003 | 0.054 |
| Palau | \< 5 | Outpatient | 0.012 | 0.002 | 0.036 | 0.016 | 0.003 | 0.051 |
| Palau | \>= 5 | Outpatient | 0.012 | 0.002 | 0.037 | 0.016 | 0.003 | 0.053 |
| Panama | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Panama | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Panama | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Panama | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Papua New Guinea | \< 5 | Inpatient | 0.016 | 0.003 | 0.041 | 0.022 | 0.003 | 0.067 |
| Papua New Guinea | \>= 5 | Inpatient | 0.016 | 0.003 | 0.043 | 0.023 | 0.003 | 0.070 |
| Papua New Guinea | \< 5 | Outpatient | 0.016 | 0.003 | 0.039 | 0.022 | 0.004 | 0.065 |
| Papua New Guinea | \>= 5 | Outpatient | 0.016 | 0.003 | 0.040 | 0.022 | 0.004 | 0.067 |
| Paraguay | \< 5 | Inpatient | 0.022 | 0.009 | 0.043 | 0.030 | 0.010 | 0.071 |
| Paraguay | \>= 5 | Inpatient | 0.022 | 0.009 | 0.047 | 0.031 | 0.009 | 0.076 |
| Paraguay | \< 5 | Outpatient | 0.022 | 0.010 | 0.041 | 0.030 | 0.010 | 0.068 |
| Paraguay | \>= 5 | Outpatient | 0.022 | 0.009 | 0.045 | 0.031 | 0.010 | 0.074 |
| Peru | \< 5 | Inpatient | 0.017 | 0.005 | 0.041 | 0.024 | 0.006 | 0.064 |
| Peru | \>= 5 | Inpatient | 0.018 | 0.004 | 0.045 | 0.025 | 0.006 | 0.068 |
| Peru | \< 5 | Outpatient | 0.017 | 0.005 | 0.040 | 0.024 | 0.006 | 0.061 |
| Peru | \>= 5 | Outpatient | 0.018 | 0.005 | 0.042 | 0.025 | 0.006 | 0.066 |
| Philippines | \< 5 | Inpatient | 0.016 | 0.003 | 0.041 | 0.022 | 0.003 | 0.067 |
| Philippines | \>= 5 | Inpatient | 0.016 | 0.003 | 0.043 | 0.023 | 0.003 | 0.070 |
| Philippines | \< 5 | Outpatient | 0.016 | 0.003 | 0.039 | 0.022 | 0.004 | 0.065 |
| Philippines | \>= 5 | Outpatient | 0.016 | 0.003 | 0.040 | 0.022 | 0.004 | 0.067 |
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
| Russian Federation | \< 5 | Inpatient | 0.043 | 0.010 | 0.138 | 0.057 | 0.013 | 0.179 |
| Russian Federation | \>= 5 | Inpatient | 0.044 | 0.009 | 0.143 | 0.058 | 0.012 | 0.183 |
| Russian Federation | \< 5 | Outpatient | 0.043 | 0.010 | 0.138 | 0.057 | 0.013 | 0.176 |
| Russian Federation | \>= 5 | Outpatient | 0.044 | 0.010 | 0.141 | 0.058 | 0.013 | 0.180 |
| Rwanda | \< 5 | Inpatient | 0.032 | 0.014 | 0.059 | 0.042 | 0.014 | 0.096 |
| Rwanda | \>= 5 | Inpatient | 0.032 | 0.013 | 0.065 | 0.043 | 0.013 | 0.103 |
| Rwanda | \< 5 | Outpatient | 0.032 | 0.015 | 0.057 | 0.043 | 0.014 | 0.097 |
| Rwanda | \>= 5 | Outpatient | 0.032 | 0.014 | 0.063 | 0.044 | 0.013 | 0.104 |
| Samoa | \< 5 | Inpatient | 0.016 | 0.003 | 0.041 | 0.022 | 0.003 | 0.067 |
| Samoa | \>= 5 | Inpatient | 0.016 | 0.003 | 0.043 | 0.023 | 0.003 | 0.070 |
| Samoa | \< 5 | Outpatient | 0.016 | 0.003 | 0.039 | 0.022 | 0.004 | 0.065 |
| Samoa | \>= 5 | Outpatient | 0.016 | 0.003 | 0.040 | 0.022 | 0.004 | 0.067 |
| San Marino | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| San Marino | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| San Marino | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| San Marino | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| São Tomé and Principe | \< 5 | Inpatient | 0.036 | 0.017 | 0.066 | 0.047 | 0.018 | 0.098 |
| São Tomé and Principe | \>= 5 | Inpatient | 0.037 | 0.016 | 0.072 | 0.048 | 0.017 | 0.106 |
| São Tomé and Principe | \< 5 | Outpatient | 0.036 | 0.018 | 0.062 | 0.047 | 0.018 | 0.100 |
| São Tomé and Principe | \>= 5 | Outpatient | 0.037 | 0.017 | 0.067 | 0.048 | 0.017 | 0.106 |
| Saudi Arabia | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Saudi Arabia | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Saudi Arabia | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Saudi Arabia | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Senegal | \< 5 | Inpatient | 0.032 | 0.007 | 0.074 | 0.043 | 0.009 | 0.113 |
| Senegal | \>= 5 | Inpatient | 0.033 | 0.007 | 0.080 | 0.044 | 0.008 | 0.120 |
| Senegal | \< 5 | Outpatient | 0.032 | 0.007 | 0.072 | 0.043 | 0.009 | 0.113 |
| Senegal | \>= 5 | Outpatient | 0.033 | 0.007 | 0.077 | 0.044 | 0.009 | 0.120 |
| Serbia | \< 5 | Inpatient | 0.043 | 0.010 | 0.138 | 0.057 | 0.013 | 0.179 |
| Serbia | \>= 5 | Inpatient | 0.044 | 0.009 | 0.143 | 0.058 | 0.012 | 0.183 |
| Serbia | \< 5 | Outpatient | 0.043 | 0.010 | 0.138 | 0.057 | 0.013 | 0.176 |
| Serbia | \>= 5 | Outpatient | 0.044 | 0.010 | 0.141 | 0.058 | 0.013 | 0.180 |
| Seychelles | \< 5 | Inpatient | 0.036 | 0.013 | 0.082 | 0.047 | 0.014 | 0.115 |
| Seychelles | \>= 5 | Inpatient | 0.037 | 0.013 | 0.088 | 0.048 | 0.013 | 0.123 |
| Seychelles | \< 5 | Outpatient | 0.036 | 0.014 | 0.079 | 0.047 | 0.014 | 0.114 |
| Seychelles | \>= 5 | Outpatient | 0.037 | 0.013 | 0.085 | 0.048 | 0.014 | 0.122 |
| Sierra Leone | \< 5 | Inpatient | 0.032 | 0.014 | 0.059 | 0.042 | 0.014 | 0.095 |
| Sierra Leone | \>= 5 | Inpatient | 0.032 | 0.013 | 0.065 | 0.043 | 0.013 | 0.103 |
| Sierra Leone | \< 5 | Outpatient | 0.032 | 0.015 | 0.057 | 0.042 | 0.014 | 0.097 |
| Sierra Leone | \>= 5 | Outpatient | 0.032 | 0.014 | 0.063 | 0.044 | 0.013 | 0.104 |
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
| Solomon Islands | \< 5 | Inpatient | 0.016 | 0.003 | 0.041 | 0.021 | 0.003 | 0.061 |
| Solomon Islands | \>= 5 | Inpatient | 0.016 | 0.003 | 0.043 | 0.022 | 0.003 | 0.064 |
| Solomon Islands | \< 5 | Outpatient | 0.016 | 0.003 | 0.039 | 0.021 | 0.004 | 0.059 |
| Solomon Islands | \>= 5 | Outpatient | 0.016 | 0.003 | 0.040 | 0.022 | 0.004 | 0.061 |
| Somalia | \< 5 | Inpatient | 0.037 | 0.013 | 0.090 | 0.050 | 0.014 | 0.129 |
| Somalia | \>= 5 | Inpatient | 0.038 | 0.012 | 0.095 | 0.051 | 0.013 | 0.137 |
| Somalia | \< 5 | Outpatient | 0.037 | 0.013 | 0.089 | 0.050 | 0.014 | 0.127 |
| Somalia | \>= 5 | Outpatient | 0.037 | 0.013 | 0.091 | 0.051 | 0.014 | 0.130 |
| South Africa | \< 5 | Inpatient | 0.034 | 0.009 | 0.085 | 0.044 | 0.012 | 0.114 |
| South Africa | \>= 5 | Inpatient | 0.035 | 0.008 | 0.092 | 0.046 | 0.011 | 0.121 |
| South Africa | \< 5 | Outpatient | 0.034 | 0.009 | 0.082 | 0.044 | 0.012 | 0.111 |
| South Africa | \>= 5 | Outpatient | 0.034 | 0.009 | 0.087 | 0.046 | 0.012 | 0.117 |
| South Sudan | \< 5 | Inpatient | 0.032 | 0.014 | 0.059 | 0.044 | 0.015 | 0.100 |
| South Sudan | \>= 5 | Inpatient | 0.032 | 0.013 | 0.065 | 0.045 | 0.014 | 0.107 |
| South Sudan | \< 5 | Outpatient | 0.032 | 0.015 | 0.057 | 0.044 | 0.016 | 0.096 |
| South Sudan | \>= 5 | Outpatient | 0.032 | 0.014 | 0.063 | 0.045 | 0.015 | 0.103 |
| Spain | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Spain | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Spain | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Spain | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Sri Lanka | \< 5 | Inpatient | 0.022 | 0.009 | 0.045 | 0.032 | 0.009 | 0.078 |
| Sri Lanka | \>= 5 | Inpatient | 0.023 | 0.009 | 0.048 | 0.032 | 0.009 | 0.083 |
| Sri Lanka | \< 5 | Outpatient | 0.022 | 0.010 | 0.042 | 0.032 | 0.010 | 0.075 |
| Sri Lanka | \>= 5 | Outpatient | 0.023 | 0.009 | 0.045 | 0.032 | 0.009 | 0.079 |
| St. Kitts and Nevis | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| St. Kitts and Nevis | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| St. Kitts and Nevis | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| St. Kitts and Nevis | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| St. Lucia | \< 5 | Inpatient | 0.023 | 0.008 | 0.052 | 0.034 | 0.008 | 0.094 |
| St. Lucia | \>= 5 | Inpatient | 0.024 | 0.008 | 0.055 | 0.035 | 0.008 | 0.098 |
| St. Lucia | \< 5 | Outpatient | 0.023 | 0.009 | 0.048 | 0.034 | 0.008 | 0.090 |
| St. Lucia | \>= 5 | Outpatient | 0.024 | 0.009 | 0.051 | 0.035 | 0.008 | 0.095 |
| St. Vincent and the Grenadines | \< 5 | Inpatient | 0.023 | 0.008 | 0.052 | 0.034 | 0.008 | 0.094 |
| St. Vincent and the Grenadines | \>= 5 | Inpatient | 0.024 | 0.008 | 0.055 | 0.035 | 0.008 | 0.098 |
| St. Vincent and the Grenadines | \< 5 | Outpatient | 0.023 | 0.009 | 0.048 | 0.034 | 0.008 | 0.090 |
| St. Vincent and the Grenadines | \>= 5 | Outpatient | 0.024 | 0.009 | 0.051 | 0.035 | 0.008 | 0.095 |
| Sudan | \< 5 | Inpatient | 0.037 | 0.013 | 0.090 | 0.049 | 0.013 | 0.133 |
| Sudan | \>= 5 | Inpatient | 0.038 | 0.012 | 0.095 | 0.050 | 0.012 | 0.139 |
| Sudan | \< 5 | Outpatient | 0.037 | 0.013 | 0.089 | 0.049 | 0.013 | 0.135 |
| Sudan | \>= 5 | Outpatient | 0.037 | 0.013 | 0.091 | 0.050 | 0.012 | 0.140 |
| Suriname | \< 5 | Inpatient | 0.023 | 0.008 | 0.052 | 0.034 | 0.008 | 0.094 |
| Suriname | \>= 5 | Inpatient | 0.024 | 0.008 | 0.055 | 0.035 | 0.008 | 0.098 |
| Suriname | \< 5 | Outpatient | 0.023 | 0.009 | 0.048 | 0.034 | 0.008 | 0.090 |
| Suriname | \>= 5 | Outpatient | 0.024 | 0.009 | 0.051 | 0.035 | 0.008 | 0.095 |
| Sweden | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Sweden | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Sweden | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Sweden | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Switzerland | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Switzerland | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Switzerland | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Switzerland | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Syrian Arab Republic | \< 5 | Inpatient | 0.037 | 0.013 | 0.090 | 0.050 | 0.014 | 0.129 |
| Syrian Arab Republic | \>= 5 | Inpatient | 0.038 | 0.012 | 0.095 | 0.051 | 0.013 | 0.137 |
| Syrian Arab Republic | \< 5 | Outpatient | 0.037 | 0.013 | 0.089 | 0.050 | 0.014 | 0.127 |
| Syrian Arab Republic | \>= 5 | Outpatient | 0.037 | 0.013 | 0.091 | 0.051 | 0.014 | 0.130 |
| Tajikistan | \< 5 | Inpatient | 0.039 | 0.009 | 0.128 | 0.050 | 0.010 | 0.163 |
| Tajikistan | \>= 5 | Inpatient | 0.039 | 0.009 | 0.127 | 0.051 | 0.010 | 0.167 |
| Tajikistan | \< 5 | Outpatient | 0.039 | 0.010 | 0.125 | 0.050 | 0.010 | 0.167 |
| Tajikistan | \>= 5 | Outpatient | 0.039 | 0.010 | 0.129 | 0.051 | 0.010 | 0.171 |
| Tanzania | \< 5 | Inpatient | 0.043 | 0.016 | 0.100 | 0.056 | 0.017 | 0.143 |
| Tanzania | \>= 5 | Inpatient | 0.044 | 0.015 | 0.106 | 0.058 | 0.017 | 0.150 |
| Tanzania | \< 5 | Outpatient | 0.043 | 0.016 | 0.099 | 0.057 | 0.017 | 0.143 |
| Tanzania | \>= 5 | Outpatient | 0.044 | 0.016 | 0.103 | 0.058 | 0.017 | 0.151 |
| Thailand | \< 5 | Inpatient | 0.015 | 0.003 | 0.040 | 0.021 | 0.004 | 0.065 |
| Thailand | \>= 5 | Inpatient | 0.015 | 0.003 | 0.042 | 0.022 | 0.004 | 0.067 |
| Thailand | \< 5 | Outpatient | 0.015 | 0.003 | 0.038 | 0.021 | 0.004 | 0.065 |
| Thailand | \>= 5 | Outpatient | 0.015 | 0.003 | 0.040 | 0.022 | 0.004 | 0.066 |
| Timor-Leste | \< 5 | Inpatient | 0.022 | 0.009 | 0.045 | 0.032 | 0.009 | 0.078 |
| Timor-Leste | \>= 5 | Inpatient | 0.023 | 0.009 | 0.048 | 0.032 | 0.009 | 0.083 |
| Timor-Leste | \< 5 | Outpatient | 0.022 | 0.010 | 0.042 | 0.032 | 0.010 | 0.075 |
| Timor-Leste | \>= 5 | Outpatient | 0.023 | 0.009 | 0.045 | 0.032 | 0.009 | 0.079 |
| Togo | \< 5 | Inpatient | 0.032 | 0.014 | 0.059 | 0.042 | 0.015 | 0.092 |
| Togo | \>= 5 | Inpatient | 0.032 | 0.013 | 0.065 | 0.043 | 0.014 | 0.099 |
| Togo | \< 5 | Outpatient | 0.032 | 0.015 | 0.057 | 0.042 | 0.015 | 0.092 |
| Togo | \>= 5 | Outpatient | 0.032 | 0.014 | 0.063 | 0.043 | 0.014 | 0.101 |
| Tonga | \< 5 | Inpatient | 0.012 | 0.002 | 0.034 | 0.017 | 0.003 | 0.056 |
| Tonga | \>= 5 | Inpatient | 0.012 | 0.002 | 0.034 | 0.017 | 0.003 | 0.056 |
| Tonga | \< 5 | Outpatient | 0.012 | 0.003 | 0.032 | 0.017 | 0.003 | 0.055 |
| Tonga | \>= 5 | Outpatient | 0.012 | 0.003 | 0.033 | 0.017 | 0.003 | 0.054 |
| Trinidad and Tobago | \< 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Trinidad and Tobago | \>= 5 | Inpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Trinidad and Tobago | \< 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Trinidad and Tobago | \>= 5 | Outpatient | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| Tunisia | \< 5 | Inpatient | 0.042 | 0.015 | 0.093 | 0.057 | 0.018 | 0.135 |
| Tunisia | \>= 5 | Inpatient | 0.043 | 0.015 | 0.099 | 0.058 | 0.017 | 0.142 |
| Tunisia | \< 5 | Outpatient | 0.042 | 0.016 | 0.090 | 0.057 | 0.019 | 0.130 |
| Tunisia | \>= 5 | Outpatient | 0.042 | 0.016 | 0.094 | 0.058 | 0.018 | 0.137 |
| Türkiye | \< 5 | Inpatient | 0.055 | 0.010 | 0.180 | 0.073 | 0.014 | 0.226 |
| Türkiye | \>= 5 | Inpatient | 0.056 | 0.010 | 0.181 | 0.074 | 0.013 | 0.229 |
| Türkiye | \< 5 | Outpatient | 0.055 | 0.011 | 0.179 | 0.073 | 0.014 | 0.224 |
| Türkiye | \>= 5 | Outpatient | 0.056 | 0.011 | 0.178 | 0.074 | 0.014 | 0.229 |
| Turkmenistan | \< 5 | Inpatient | 0.043 | 0.010 | 0.138 | 0.055 | 0.011 | 0.179 |
| Turkmenistan | \>= 5 | Inpatient | 0.044 | 0.009 | 0.143 | 0.056 | 0.011 | 0.182 |
| Turkmenistan | \< 5 | Outpatient | 0.043 | 0.010 | 0.138 | 0.056 | 0.011 | 0.179 |
| Turkmenistan | \>= 5 | Outpatient | 0.044 | 0.010 | 0.141 | 0.057 | 0.011 | 0.186 |
| Tuvalu | \< 5 | Inpatient | 0.012 | 0.002 | 0.034 | 0.017 | 0.003 | 0.056 |
| Tuvalu | \>= 5 | Inpatient | 0.012 | 0.002 | 0.034 | 0.017 | 0.003 | 0.056 |
| Tuvalu | \< 5 | Outpatient | 0.012 | 0.003 | 0.032 | 0.017 | 0.003 | 0.055 |
| Tuvalu | \>= 5 | Outpatient | 0.012 | 0.003 | 0.033 | 0.017 | 0.003 | 0.054 |
| Uganda | \< 5 | Inpatient | 0.032 | 0.014 | 0.059 | 0.042 | 0.014 | 0.095 |
| Uganda | \>= 5 | Inpatient | 0.032 | 0.013 | 0.065 | 0.043 | 0.013 | 0.102 |
| Uganda | \< 5 | Outpatient | 0.032 | 0.015 | 0.057 | 0.042 | 0.014 | 0.096 |
| Uganda | \>= 5 | Outpatient | 0.032 | 0.014 | 0.063 | 0.044 | 0.013 | 0.104 |
| Ukraine | \< 5 | Inpatient | 0.039 | 0.009 | 0.128 | 0.051 | 0.012 | 0.166 |
| Ukraine | \>= 5 | Inpatient | 0.039 | 0.009 | 0.127 | 0.052 | 0.011 | 0.168 |
| Ukraine | \< 5 | Outpatient | 0.039 | 0.010 | 0.125 | 0.051 | 0.012 | 0.162 |
| Ukraine | \>= 5 | Outpatient | 0.039 | 0.010 | 0.129 | 0.052 | 0.012 | 0.169 |
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
| Uzbekistan | \< 5 | Inpatient | 0.039 | 0.009 | 0.128 | 0.049 | 0.010 | 0.160 |
| Uzbekistan | \>= 5 | Inpatient | 0.039 | 0.009 | 0.127 | 0.050 | 0.010 | 0.164 |
| Uzbekistan | \< 5 | Outpatient | 0.039 | 0.010 | 0.125 | 0.050 | 0.011 | 0.165 |
| Uzbekistan | \>= 5 | Outpatient | 0.039 | 0.010 | 0.129 | 0.051 | 0.010 | 0.167 |
| Vanuatu | \< 5 | Inpatient | 0.016 | 0.003 | 0.041 | 0.022 | 0.003 | 0.067 |
| Vanuatu | \>= 5 | Inpatient | 0.016 | 0.003 | 0.043 | 0.023 | 0.003 | 0.070 |
| Vanuatu | \< 5 | Outpatient | 0.016 | 0.003 | 0.039 | 0.022 | 0.004 | 0.065 |
| Vanuatu | \>= 5 | Outpatient | 0.016 | 0.003 | 0.040 | 0.022 | 0.004 | 0.067 |
| Venezuela | \< 5 | Inpatient | 0.025 | 0.005 | 0.073 | 0.038 | 0.005 | 0.130 |
| Venezuela | \>= 5 | Inpatient | 0.026 | 0.005 | 0.077 | 0.039 | 0.005 | 0.137 |
| Venezuela | \< 5 | Outpatient | 0.025 | 0.006 | 0.071 | 0.038 | 0.006 | 0.126 |
| Venezuela | \>= 5 | Outpatient | 0.026 | 0.005 | 0.075 | 0.039 | 0.006 | 0.132 |
| Vietnam | \< 5 | Inpatient | 0.016 | 0.003 | 0.041 | 0.022 | 0.003 | 0.067 |
| Vietnam | \>= 5 | Inpatient | 0.016 | 0.003 | 0.043 | 0.023 | 0.003 | 0.070 |
| Vietnam | \< 5 | Outpatient | 0.016 | 0.003 | 0.039 | 0.022 | 0.004 | 0.065 |
| Vietnam | \>= 5 | Outpatient | 0.016 | 0.003 | 0.040 | 0.022 | 0.004 | 0.067 |
| Yemen, Rep. | \< 5 | Inpatient | 0.037 | 0.013 | 0.090 | 0.048 | 0.014 | 0.123 |
| Yemen, Rep. | \>= 5 | Inpatient | 0.038 | 0.012 | 0.095 | 0.049 | 0.014 | 0.130 |
| Yemen, Rep. | \< 5 | Outpatient | 0.037 | 0.013 | 0.089 | 0.049 | 0.015 | 0.122 |
| Yemen, Rep. | \>= 5 | Outpatient | 0.037 | 0.013 | 0.091 | 0.050 | 0.014 | 0.128 |
| Zambia | \< 5 | Inpatient | 0.033 | 0.009 | 0.075 | 0.044 | 0.010 | 0.112 |
| Zambia | \>= 5 | Inpatient | 0.034 | 0.009 | 0.079 | 0.045 | 0.010 | 0.117 |
| Zambia | \< 5 | Outpatient | 0.033 | 0.009 | 0.073 | 0.044 | 0.010 | 0.112 |
| Zambia | \>= 5 | Outpatient | 0.034 | 0.009 | 0.078 | 0.045 | 0.010 | 0.121 |
| Zimbabwe | \< 5 | Inpatient | 0.039 | 0.011 | 0.094 | 0.051 | 0.014 | 0.131 |
| Zimbabwe | \>= 5 | Inpatient | 0.040 | 0.011 | 0.099 | 0.053 | 0.013 | 0.140 |
| Zimbabwe | \< 5 | Outpatient | 0.039 | 0.012 | 0.091 | 0.052 | 0.014 | 0.132 |
| Zimbabwe | \>= 5 | Outpatient | 0.040 | 0.012 | 0.095 | 0.053 | 0.013 | 0.140 |

Estimated `r params$Pathogen` proportion by country, 2010 vs 2020

# Session info

``` r
sessioninfo::session_info()
```

    ## Warning in system2("quarto", "-V", stdout = TRUE, env = paste0("TMPDIR=", : running command '"quarto"
    ## TMPDIR=C:/Users/fbbu6966/AppData/Local/Temp/Rtmp0kgjaE/file38f021f825ad -V' had status 1

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
    ##  quarto   ERROR: Unknown command "TMPDIR=C:/Users/fbbu6966/AppData/Local/Temp/Rtmp0kgjaE/file38f021f825ad". Did you mean command "update"? @ C:\\Users\\fbbu6966\\AppData\\Local\\Programs\\RStudio\\RESOUR~1\\app\\bin\\quarto\\bin\\quarto.exe
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
