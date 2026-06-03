Diarrheal Entamoeba • Models fit
================
fbbu6966
2025-11-27

- [Settings](#settings)
- [Parameters](#parameters)
- [Data](#data)
- [Time trend](#time-trend)
- [BRMS](#brms)
  - [Reduced model fit](#reduced-model-fit)
  - [Model summary](#model-summary)
- [Session info](#session-info)

# Settings

``` r
## required packages ----
library(bd)
library(brms)
library(ggplot2)
library(metafor)
library(readxl)
library(rmarkdown)
library(rms)
library(tidyr)
library(kableExtra)
```

    ## 
    ## Attaching package: 'kableExtra'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     group_rows

``` r
library(stringr)

## global options ----
knitr::opts_chunk$set(fig.width = 10)
do.call(file.remove, list(list.files(params$PlotDir, full.names = TRUE)))
```

    ## logical(0)

# Parameters

| Parameters                       | Values             |
|:---------------------------------|:-------------------|
| Pathogen                         | Entamoeba          |
| Number of iteration              | 5000               |
| Warmup                           | 3000               |
| Delta value                      | 0.95               |
| Maximum tree-depth               | 15                 |
| Levels                           | Countries, Studies |
| Random effect on each data point | No                 |
| Stronger priors specified        | Normal(0,1)        |

Parameters of the model tested

# Data

``` r
## import data
```

``` r
source("01-data-v20251127.R")
```

    ## Warning in system2("quarto", "-V", stdout = TRUE, env = paste0("TMPDIR=", : running command '"quarto"
    ## TMPDIR=C:/Users/fbbu6966/AppData/Local/Temp/Rtmpgxvndi/file2db42d741384 -V' had status 1

``` r
es <- all_dta[[params$i]]
es$DTP_ID <- as.character(seq(1:length(es$SOURCE_ID)))
es$FLAG <-
  factor(es$FLAG, 
         levels=c(0,1,2,3,4,5,6, 7),
         labels=c("Keep data", "Data part of non WHO member states", "No WHO REG2 given",
                  "Year before 1990", "yi can't be calcualted", "TF choice to remove", 
                  "Excluded by preliminary checks", "Excluded in data cleaning"))

VacCov <- read_excel("../rota coverage.xlsx", sheet="Cov", range=cell_rows(5:206))
VacCov <- subset(VacCov, !is.na(`2023`))
VacCov <- VacCov %>% 
  mutate(COUNTRY2 = case_when(
    `Row Labels` == "Democratic Republic Of The Congo" ~ "congo (the democratic republic of the)",
    `Row Labels` == "Congo" ~ "congo (the)",
    `Row Labels` == "Côte D’ivoire" ~ "côte d'ivoire",
    `Row Labels` == "Dominican Republic" ~ "dominican republic (the)",
    `Row Labels` == "Gambia" ~ "gambia (the)",
    `Row Labels` == "Democratic People’s Republic Of Korea" ~ "korea (the democratic people's republic of)",
    `Row Labels` == "Republic Of Korea" ~ "korea (the republic of)",
    `Row Labels` == "Lao People’s Democratic Republic" ~ "lao people's dem. republic",
    `Row Labels` == "Micronesia (Federated States Of )" ~ "micronesia (fed. states of)",
    `Row Labels` == "Netherlands (Kingdom Of The)" ~ "netherlands",
    `Row Labels` == "São Tomé And Principe" ~ "sao tome and principe",
    `Row Labels` == "Türkiye" ~ "turkiye",
    `Row Labels` == "United Kingdom Of Great Britain And Northern Ireland" ~ "united kingdom",
    .default = tolower(`Row Labels`)
  ))
countries <- FERG2:::countries[,c("ISO3","COUNTRY")]
countries$COUNTRY2 <- tolower(countries$COUNTRY)
VacCov <- merge(VacCov, countries)
VacCov$`1990` <- 0
VacCov$`1991` <- 0
VacCov$`1992` <- 0
VacCov$`1993` <- 0
VacCov$`1994` <- 0
VacCov$`1995` <- 0
VacCov$`1996` <- 0
VacCov$`1997` <- 0
VacCov$`1998` <- 0
VacCov$`1999` <- 0
VacCov$`2000` <- 0
VacCov$`2001` <- 0
VacCov$`2002` <- 0
VacCov$`2003` <- 0
VacCov$`2004` <- 0
VacCov$`2005` <- 0
VacCov <- VacCov[,c(21, 23:38, 3:20)]
VacCov_long <-
  pivot_longer(VacCov, cols=c(2:35), names_to="YEAR", values_to="COVERAGE")
VacCov_long$COVERAGE <-
  if_else(is.na(VacCov_long$COVERAGE), 0, VacCov_long$COVERAGE)

es$YEAR_COV <- round(es$YEAR)
es <- merge(es, VacCov_long, all.x=TRUE, by.x=c("ISO3", "YEAR_COV"), by.y=c("ISO3", "YEAR"))
es$YEAR_COV <- NULL

# define reference strain
es$STRAIN <- factor(es$STRAIN)
es$STRAIN <- relevel(es$STRAIN, "histolytica")
print(table(es$STRAIN))
```

    ## 
    ##        histolytica histolytica_dispar 
    ##                361                 24

``` r
# reference DX
print(table(es$DIAGNOSIS))
```

    ## 
    ##     Culture   Molecular Immunoassay  Microscopy   Other/UNK 
    ##           0         230          72          83           0

``` r
print(table(es$REFERENCE))
```

    ## 
    ## Reference     Other 
    ##       230       155

``` r
saveRDS(VacCov_long, paste0(params$Dir,"/vaccine_coverage.RDS"))
saveRDS(es, paste0(params$Dir,"/es.RDS"))
```

# Time trend

``` r
## No time trend
fit_rma <-
  rma.mv(yi = yi ~ 1 + SYNDROMTYPE + AGE + REFERENCE + COVERAGE + STRAIN, 
         V = vi, 
         data = subset(es, as.integer(FLAG) == 1),
         random = ~ 1 | REG2 / SUB2 / COUNTRY / ID, 
         test = "t", 
         method = "REML",
         control=list(rel.tol=1e-8))

summary(fit_rma)
```

    ## 
    ## Multivariate Meta-Analysis Model (k = 351; method: REML)
    ## 
    ##    logLik   Deviance        AIC        BIC       AICc   
    ## -694.0198  1388.0395  1412.0395  1458.0923  1412.9850   
    ## 
    ## Variance Components:
    ## 
    ##             estim    sqrt  nlvls  fixed                factor 
    ## sigma^2.1  0.0000  0.0004      6     no                  REG2 
    ## sigma^2.2  0.0000  0.0001     12     no             REG2/SUB2 
    ## sigma^2.3  0.0000  0.0014     39     no     REG2/SUB2/COUNTRY 
    ## sigma^2.4  2.1324  1.4603     83     no  REG2/SUB2/COUNTRY/ID 
    ## 
    ## Test for Residual Heterogeneity:
    ## QE(df = 343) = 5603.5863, p-val < .0001
    ## 
    ## Test of Moderators (coefficients 2:8):
    ## F(df1 = 7, df2 = 343) = 19.8514, p-val < .0001
    ## 
    ## Model Results:
    ## 
    ##                           estimate      se      tval   df    pval    ci.lb    ci.ub      
    ## intrcpt                    -4.4517  0.2873  -15.4937  343  <.0001  -5.0169  -3.8866  *** 
    ## SYNDROMTYPEInpatient        0.6532  0.1520    4.2984  343  <.0001   0.3543   0.9521  *** 
    ## SYNDROMTYPEOutpatient      -0.2568  0.0725   -3.5437  343  0.0004  -0.3993  -0.1143  *** 
    ## AGEAge below 5             -0.5581  0.0650   -8.5838  343  <.0001  -0.6859  -0.4302  *** 
    ## AGEMixed ages               0.0420  0.0626    0.6706  343  0.5029  -0.0811   0.1651      
    ## REFERENCEOther              1.1564  0.3542    3.2651  343  0.0012   0.4598   1.8531   ** 
    ## COVERAGE                    0.0029  0.0050    0.5683  343  0.5702  -0.0070   0.0127      
    ## STRAINhistolytica_dispar    0.6753  0.5428    1.2441  343  0.2143  -0.3923   1.7430      
    ## 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
## Linear time trend
fit_rma_glb <-
  rma.mv(yi = yi ~ 1 + YEAR + SYNDROMTYPE + AGE + REFERENCE + COVERAGE + STRAIN, 
         V = vi, 
         data = subset(es, as.integer(FLAG) == 1),
         random = ~ 1 | REG2 / SUB2 / COUNTRY / ID, 
         test = "t", 
         method = "REML",
         control = list(rel.tol = 1e-8))

summary(fit_rma_glb)
```

    ## 
    ## Multivariate Meta-Analysis Model (k = 351; method: REML)
    ## 
    ##    logLik   Deviance        AIC        BIC       AICc   
    ## -691.9745  1383.9490  1409.9490  1459.8015  1411.0588   
    ## 
    ## Variance Components:
    ## 
    ##             estim    sqrt  nlvls  fixed                factor 
    ## sigma^2.1  0.0000  0.0002      6     no                  REG2 
    ## sigma^2.2  0.0000  0.0001     12     no             REG2/SUB2 
    ## sigma^2.3  0.0000  0.0012     39     no     REG2/SUB2/COUNTRY 
    ## sigma^2.4  2.1475  1.4654     83     no  REG2/SUB2/COUNTRY/ID 
    ## 
    ## Test for Residual Heterogeneity:
    ## QE(df = 342) = 5281.5054, p-val < .0001
    ## 
    ## Test of Moderators (coefficients 2:9):
    ## F(df1 = 8, df2 = 342) = 17.4182, p-val < .0001
    ## 
    ## Model Results:
    ## 
    ##                           estimate       se     tval   df    pval     ci.lb     ci.ub      
    ## intrcpt                    39.4957  63.1605   0.6253  342  0.5322  -84.7362  163.7276      
    ## YEAR                       -0.0219   0.0314  -0.6958  342  0.4870   -0.0836    0.0399      
    ## SYNDROMTYPEInpatient        0.6580   0.1521   4.3253  342  <.0001    0.3588    0.9573  *** 
    ## SYNDROMTYPEOutpatient      -0.2574   0.0725  -3.5513  342  0.0004   -0.3999   -0.1148  *** 
    ## AGEAge below 5             -0.5581   0.0650  -8.5847  342  <.0001   -0.6860   -0.4303  *** 
    ## AGEMixed ages               0.0413   0.0626   0.6595  342  0.5100   -0.0818    0.1644      
    ## REFERENCEOther              1.0447   0.3899   2.6792  342  0.0077    0.2777    1.8116   ** 
    ## COVERAGE                    0.0042   0.0054   0.7783  342  0.4370   -0.0064    0.0148      
    ## STRAINhistolytica_dispar    0.7195   0.5483   1.3122  342  0.1903   -0.3590    1.7980      
    ## 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(fit_rma_glb, fit_rma)
```

    ## Warning: REML comparisons not meaningful for models with different fixed effects
    ## (use 'refit=TRUE' to refit both models based on ML estimation).

    ## 
    ##         df       AIC       BIC      AICc    logLik    LRT   pval        QE 
    ## Full    13 1409.9490 1459.8015 1411.0588 -691.9745               5281.5054 
    ## Reduced 12 1412.0395 1458.0923 1412.9850 -694.0198 4.0905 0.0431 5603.5863

``` r
png(paste0(params$PlotDir, "/r_lineartimetrend.png"), width = 960)
ggplot(subset(es, as.integer(FLAG) == 1), aes(x = YEAR, y = yi)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw()
```

    ## `geom_smooth()` using formula = 'y ~ x'

``` r
dev.off()
```

png 2

``` r
setwd(params$Dir)
image <- paste0("02-fit-rota-strain-enta_files/figure-gfm/r_lineartimetrend.png")
cat("![](",image,")")
```

![](02-fit-rota-strain-enta_files/figure-gfm/r_lineartimetrend.png)

``` r
## Non-linear global time trend
```

``` r
png(paste0(params$PlotDir, "/r_nonlineartimetrend.png"), width = 960)
ggplot(subset(es, as.integer(FLAG) == 1), aes(x = YEAR, y = yi)) +
  geom_point() +
  geom_smooth(method = "loess") +
  theme_bw()
```

    ## `geom_smooth()` using formula = 'y ~ x'

``` r
dev.off()
```

png 2

``` r
setwd(params$Dir)
image <- paste0("02-fit-rota-strain-enta_files/figure-gfm/r_nonlineartimetrend.png")
cat("![](",image,")")
```

![](02-fit-rota-strain-enta_files/figure-gfm/r_nonlineartimetrend.png)

# BRMS

## Reduced model fit

``` r
fit_brms_reg_s <-
  brm(yi | se(sei) ~
        1 + YEAR + SYNDROMTYPE + AGE + REFERENCE + COVERAGE + STRAIN +
        (1 | REG2) +
        (1 | REG2:SUB2) +
        (1 | REG2:SUB2:COUNTRY) +
        (1 | REG2:SUB2:COUNTRY:ID) +
        (1 | REG2:SUB2:COUNTRY:ID:DTP_ID),
      chains = 5, iter = params$Iterations, warmup = params$Warmup,
      prior = prior(normal(0,1), class = sd),
      control = controls,
      cores = 5,
      data = subset(es, as.integer(FLAG) == 1),
      open_progress = FALSE,
      seed = 7)
```

    ## Compiling Stan program...

    ## Start sampling

    ## Warning: There were 9 divergent transitions after warmup. See
    ## https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
    ## to find out why this is a problem and how to eliminate them.

    ## Warning: Examine the pairs() plot to diagnose sampling problems

## Model summary

``` r
summary(fit_brms_reg_s)
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
est <- colnames(as_draws_df(fit_brms_reg_s))
selected_est <- est[grepl("^(sd_|b_|sigma)", est)]
length <- length(selected_est)
for (i in 1:length){
  png(paste0(params$PlotDir, "/r_summary-", i, ".png"), width = 960)
  plot(fit_brms_reg_s, variable = selected_est[i])
  dev.off()
}
```

``` r
setwd(params$Dir)
for (i in 1:length){
  image <- paste0("02-fit-rota-strain-enta_files/figure-gfm/r_summary-",i,".png")
  cat("![](",image,")")
}
```

![](02-fit-rota-strain-enta_files/figure-gfm/r_summary-1.png)![](02-fit-rota-strain-enta_files/figure-gfm/r_summary-2.png)![](02-fit-rota-strain-enta_files/figure-gfm/r_summary-3.png)![](02-fit-rota-strain-enta_files/figure-gfm/r_summary-4.png)![](02-fit-rota-strain-enta_files/figure-gfm/r_summary-5.png)![](02-fit-rota-strain-enta_files/figure-gfm/r_summary-6.png)![](02-fit-rota-strain-enta_files/figure-gfm/r_summary-7.png)![](02-fit-rota-strain-enta_files/figure-gfm/r_summary-8.png)![](02-fit-rota-strain-enta_files/figure-gfm/r_summary-9.png)![](02-fit-rota-strain-enta_files/figure-gfm/r_summary-10.png)![](02-fit-rota-strain-enta_files/figure-gfm/r_summary-11.png)![](02-fit-rota-strain-enta_files/figure-gfm/r_summary-12.png)![](02-fit-rota-strain-enta_files/figure-gfm/r_summary-13.png)![](02-fit-rota-strain-enta_files/figure-gfm/r_summary-14.png)![](02-fit-rota-strain-enta_files/figure-gfm/r_summary-15.png)

``` r
P <- plot(conditional_effects(fit_brms_reg_s), points = TRUE, ask=FALSE)
```

![](02-fit-rota-strain-enta_files/figure-gfm/r_conditional_effects-1.png)<!-- -->![](02-fit-rota-strain-enta_files/figure-gfm/r_conditional_effects-2.png)<!-- -->![](02-fit-rota-strain-enta_files/figure-gfm/r_conditional_effects-3.png)<!-- -->![](02-fit-rota-strain-enta_files/figure-gfm/r_conditional_effects-4.png)<!-- -->![](02-fit-rota-strain-enta_files/figure-gfm/r_conditional_effects-5.png)<!-- -->![](02-fit-rota-strain-enta_files/figure-gfm/r_conditional_effects-6.png)<!-- -->

``` r
length <- length(P)
for (i in 1:length){
  png(paste0(params$PlotDir, "/r_conditional_effects-", i, ".png"), width = 960)
  plot(P[[i]])
  dev.off()
}
```

``` r
# setwd(params$Dir)
for (i in 1:length){
  image <- paste0("02-fit-rota-strain-enta_files/figure-gfm/r_conditional_effects-",i,".png")
  cat("![](",image,")")
}
```

![](02-fit-rota-strain-enta_files/figure-gfm/r_conditional_effects-1.png)![](02-fit-rota-strain-enta_files/figure-gfm/r_conditional_effects-2.png)![](02-fit-rota-strain-enta_files/figure-gfm/r_conditional_effects-3.png)![](02-fit-rota-strain-enta_files/figure-gfm/r_conditional_effects-4.png)![](02-fit-rota-strain-enta_files/figure-gfm/r_conditional_effects-5.png)![](02-fit-rota-strain-enta_files/figure-gfm/r_conditional_effects-6.png)

``` r
saveRDS(fit_brms_reg_s, file = paste0(params$Dir, "/fit_brms_reg_s.rds"))

## show model code
```

``` r
stancode(fit_brms_reg_s)
```

    ## // generated with brms 2.22.0
    ## functions {
    ## }
    ## data {
    ##   int<lower=1> N;  // total number of observations
    ##   vector[N] Y;  // response variable
    ##   vector<lower=0>[N] se;  // known sampling error
    ##   int<lower=1> K;  // number of population-level effects
    ##   matrix[N, K] X;  // population-level design matrix
    ##   int<lower=1> Kc;  // number of population-level effects after centering
    ##   // data for group-level effects of ID 1
    ##   int<lower=1> N_1;  // number of grouping levels
    ##   int<lower=1> M_1;  // number of coefficients per level
    ##   array[N] int<lower=1> J_1;  // grouping indicator per observation
    ##   // group-level predictor values
    ##   vector[N] Z_1_1;
    ##   // data for group-level effects of ID 2
    ##   int<lower=1> N_2;  // number of grouping levels
    ##   int<lower=1> M_2;  // number of coefficients per level
    ##   array[N] int<lower=1> J_2;  // grouping indicator per observation
    ##   // group-level predictor values
    ##   vector[N] Z_2_1;
    ##   // data for group-level effects of ID 3
    ##   int<lower=1> N_3;  // number of grouping levels
    ##   int<lower=1> M_3;  // number of coefficients per level
    ##   array[N] int<lower=1> J_3;  // grouping indicator per observation
    ##   // group-level predictor values
    ##   vector[N] Z_3_1;
    ##   // data for group-level effects of ID 4
    ##   int<lower=1> N_4;  // number of grouping levels
    ##   int<lower=1> M_4;  // number of coefficients per level
    ##   array[N] int<lower=1> J_4;  // grouping indicator per observation
    ##   // group-level predictor values
    ##   vector[N] Z_4_1;
    ##   // data for group-level effects of ID 5
    ##   int<lower=1> N_5;  // number of grouping levels
    ##   int<lower=1> M_5;  // number of coefficients per level
    ##   array[N] int<lower=1> J_5;  // grouping indicator per observation
    ##   // group-level predictor values
    ##   vector[N] Z_5_1;
    ##   int prior_only;  // should the likelihood be ignored?
    ## }
    ## transformed data {
    ##   vector<lower=0>[N] se2 = square(se);
    ##   matrix[N, Kc] Xc;  // centered version of X without an intercept
    ##   vector[Kc] means_X;  // column means of X before centering
    ##   for (i in 2:K) {
    ##     means_X[i - 1] = mean(X[, i]);
    ##     Xc[, i - 1] = X[, i] - means_X[i - 1];
    ##   }
    ## }
    ## parameters {
    ##   vector[Kc] b;  // regression coefficients
    ##   real Intercept;  // temporary intercept for centered predictors
    ##   vector<lower=0>[M_1] sd_1;  // group-level standard deviations
    ##   array[M_1] vector[N_1] z_1;  // standardized group-level effects
    ##   vector<lower=0>[M_2] sd_2;  // group-level standard deviations
    ##   array[M_2] vector[N_2] z_2;  // standardized group-level effects
    ##   vector<lower=0>[M_3] sd_3;  // group-level standard deviations
    ##   array[M_3] vector[N_3] z_3;  // standardized group-level effects
    ##   vector<lower=0>[M_4] sd_4;  // group-level standard deviations
    ##   array[M_4] vector[N_4] z_4;  // standardized group-level effects
    ##   vector<lower=0>[M_5] sd_5;  // group-level standard deviations
    ##   array[M_5] vector[N_5] z_5;  // standardized group-level effects
    ## }
    ## transformed parameters {
    ##   real sigma = 0;  // dispersion parameter
    ##   vector[N_1] r_1_1;  // actual group-level effects
    ##   vector[N_2] r_2_1;  // actual group-level effects
    ##   vector[N_3] r_3_1;  // actual group-level effects
    ##   vector[N_4] r_4_1;  // actual group-level effects
    ##   vector[N_5] r_5_1;  // actual group-level effects
    ##   real lprior = 0;  // prior contributions to the log posterior
    ##   r_1_1 = (sd_1[1] * (z_1[1]));
    ##   r_2_1 = (sd_2[1] * (z_2[1]));
    ##   r_3_1 = (sd_3[1] * (z_3[1]));
    ##   r_4_1 = (sd_4[1] * (z_4[1]));
    ##   r_5_1 = (sd_5[1] * (z_5[1]));
    ##   lprior += student_t_lpdf(Intercept | 3, -4.6, 2.5);
    ##   lprior += normal_lpdf(sd_1 | 0, 1)
    ##     - 1 * normal_lccdf(0 | 0, 1);
    ##   lprior += normal_lpdf(sd_2 | 0, 1)
    ##     - 1 * normal_lccdf(0 | 0, 1);
    ##   lprior += normal_lpdf(sd_3 | 0, 1)
    ##     - 1 * normal_lccdf(0 | 0, 1);
    ##   lprior += normal_lpdf(sd_4 | 0, 1)
    ##     - 1 * normal_lccdf(0 | 0, 1);
    ##   lprior += normal_lpdf(sd_5 | 0, 1)
    ##     - 1 * normal_lccdf(0 | 0, 1);
    ## }
    ## model {
    ##   // likelihood including constants
    ##   if (!prior_only) {
    ##     // initialize linear predictor term
    ##     vector[N] mu = rep_vector(0.0, N);
    ##     mu += Intercept + Xc * b;
    ##     for (n in 1:N) {
    ##       // add more terms to the linear predictor
    ##       mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_2_1[J_2[n]] * Z_2_1[n] + r_3_1[J_3[n]] * Z_3_1[n] + r_4_1[J_4[n]] * Z_4_1[n] + r_5_1[J_5[n]] * Z_5_1[n];
    ##     }
    ##     target += normal_lpdf(Y | mu, se);
    ##   }
    ##   // priors including constants
    ##   target += lprior;
    ##   target += std_normal_lpdf(z_1[1]);
    ##   target += std_normal_lpdf(z_2[1]);
    ##   target += std_normal_lpdf(z_3[1]);
    ##   target += std_normal_lpdf(z_4[1]);
    ##   target += std_normal_lpdf(z_5[1]);
    ## }
    ## generated quantities {
    ##   // actual population-level intercept
    ##   real b_Intercept = Intercept - dot_product(means_X, b);
    ## }

# Session info

``` r
sessioninfo::session_info()
```

    ## Warning in system2("quarto", "-V", stdout = TRUE, env = paste0("TMPDIR=", : running command '"quarto"
    ## TMPDIR=C:/Users/fbbu6966/AppData/Local/Temp/Rtmpgxvndi/file2db4147449f3 -V' had status 1

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
    ##  quarto   ERROR: Unknown command "TMPDIR=C:/Users/fbbu6966/AppData/Local/Temp/Rtmpgxvndi/file2db4147449f3". Did you mean command "render"? @ C:\\Users\\fbbu6966\\AppData\\Local\\Programs\\RStudio\\RESOUR~1\\app\\bin\\quarto\\bin\\quarto.exe
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
