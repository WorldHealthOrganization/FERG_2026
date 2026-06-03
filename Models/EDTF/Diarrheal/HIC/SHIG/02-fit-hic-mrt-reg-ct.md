Diarrheal Mortality Shigella • Model fit
================
fbbu6966
2025-12-08

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

## global options ----
knitr::opts_chunk$set(fig.width = 10)
do.call(file.remove, list(list.files(params$PlotDir, full.names = TRUE)))
```

    ##  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [21] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE

# Parameters

| Parameters                | Values      |
|:--------------------------|:------------|
| Pathogen                  | Shigella    |
| Number of iteration       | 5000        |
| Warmup                    | 3000        |
| Delta value               | 0.9         |
| Maximum tree-depth        | 15          |
| Levels                    | REG_CT      |
| Stronger priors specified | Normal(0,1) |

Parameters of the model tested

# Data

``` r
## import data
```

``` r
source("01-data-hic-mrt-v20251208.R")
```

    ## Warning in RColorBrewer::brewer.pal(max_freq, "Greens"): minimal value for n is 3, returning requested palette with 3 different levels

    ## Warning in RColorBrewer::brewer.pal(max_freq, "Greens"): minimal value for n is 3, returning requested palette with 3 different levels

    ## Warning in RColorBrewer::brewer.pal(max_freq, "Greens"): minimal value for n is 3, returning requested palette with 3 different levels

    ## Warning in RColorBrewer::brewer.pal(max_freq, "Greens"): minimal value for n is 3, returning requested palette with 3 different levels

    ## Warning in system2("quarto", "-V", stdout = TRUE, env = paste0("TMPDIR=", : running command '"quarto"
    ## TMPDIR=C:/Users/fbbu6966/AppData/Local/Temp/RtmpwhXRBk/file15646e0914c2 -V' had status 1

``` r
es<-es[[params$i]]
es$DTP_ID <- as.character(seq(1:length(es$SOURCE_ID)))
es$FLAG <-
  factor(es$FLAG, 
         levels = c(0, 1, 2, 3, 4, 5, 6, 7),
         labels = c("Keep data", "Data part of non WHO member states", "No WHO REG2 given",
                    "Year before 1990", "yi can't be calcualted", "TF choice to remove", 
                    "Excluded by preliminary checks", "Excluded in data cleaning"))
saveRDS(es, paste0(params$Dir, sprintf("/es_hic_mrt_%s.rds", bd::today())))
```

# Time trend

``` r
## No time trend
fit_rma <-
  rma.mv(yi = yi , 
         V = vi, 
         data = subset(es, as.integer(FLAG) == 1),
         # random = ~ 1 | REG2 / ID, 
         random = ~ 1 | ID, 
         test = "t", 
         method = "REML")
summary(fit_rma)
```

    ## 
    ## Multivariate Meta-Analysis Model (k = 11; method: REML)
    ## 
    ##   logLik  Deviance       AIC       BIC      AICc   
    ## -16.8395   33.6790   37.6790   38.2841   39.3932   
    ## 
    ## Variance Components:
    ## 
    ##             estim    sqrt  nlvls  fixed  factor 
    ## sigma^2    0.8484  0.9211     11     no      ID 
    ## 
    ## Test for Heterogeneity:
    ## Q(df = 10) = 33.6891, p-val = 0.0002
    ## 
    ## Model Results:
    ## 
    ## estimate      se      tval  df    pval    ci.lb    ci.ub      
    ##  -5.0590  0.4160  -12.1600  10  <.0001  -5.9859  -4.1320  *** 
    ## 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
## Linear time trend
if(nrow(subset(es, as.integer(FLAG) == 1)) > 2){
  fit_rma_glb <-
    rma.mv(yi = yi ~ 1 + YEAR, 
           V = vi, 
           data = subset(es, as.integer(FLAG) == 1),
           # random = ~ 1 | REG2 / ID, 
           random = ~ 1 | ID, 
           test = "t", 
           method = "REML",
           control = list(rel.tol = 1e-8))
  summary(fit_rma_glb)
}
```

    ## 
    ## Multivariate Meta-Analysis Model (k = 11; method: REML)
    ## 
    ##   logLik  Deviance       AIC       BIC      AICc   
    ## -12.2889   24.5777   30.5777   31.1694   35.3777   
    ## 
    ## Variance Components:
    ## 
    ##             estim    sqrt  nlvls  fixed  factor 
    ## sigma^2    0.1219  0.3491     11     no      ID 
    ## 
    ## Test for Residual Heterogeneity:
    ## QE(df = 9) = 7.7384, p-val = 0.5607
    ## 
    ## Test of Moderators (coefficient 2):
    ## F(df1 = 1, df2 = 9) = 11.9165, p-val = 0.0073
    ## 
    ## Model Results:
    ## 
    ##          estimate       se     tval  df    pval    ci.lb     ci.ub     
    ## intrcpt  254.1894  74.9241   3.3926   9  0.0080  84.6993  423.6796  ** 
    ## YEAR      -0.1291   0.0374  -3.4520   9  0.0073  -0.2137   -0.0445  ** 
    ## 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
if (nrow(subset(es, as.integer(FLAG) == 1)) > 2){
  anova(fit_rma_glb, fit_rma)
}
```

    ## Warning: REML comparisons not meaningful for models with different fixed effects
    ## (use 'refit=TRUE' to refit both models based on ML estimation).

    ## 
    ##         df     AIC     BIC    AICc   logLik    LRT   pval      QE 
    ## Full     3 30.5777 31.1694 35.3777 -12.2889                7.7384 
    ## Reduced  2 37.6790 38.2841 39.3932 -16.8395 9.1013 0.0026 33.6891

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
image <- file.path(pdir, "r_lineartimetrend.png")
cat("![](",image,")")
```

![](02-fit-hic-mrt-reg-ct_files/figure-gfm/r_lineartimetrend.png)

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
image <- file.path(pdir, "r_nonlineartimetrend.png")
cat("![](",image,")")
```

![](02-fit-hic-mrt-reg-ct_files/figure-gfm/r_nonlineartimetrend.png)

# BRMS

## Reduced model fit

``` r
if (params$Level == "REG_YR"){
  fit_brms_reg_s <-
    brm(yi | se(sei) ~ 1 + YEAR + (1 | REG2) + (1 | REG2:ID),
        chains = 5, iter = params$Iterations, warmup = params$Warmup,
        prior = prior(normal(0,1), class = sd),
        control = controls,
        cores = 5,
        data = subset(es, as.integer(FLAG) == 1),
        open_progress = FALSE,
        seed = 7)
  saveRDS(fit_brms_reg_s, file = paste0(params$Dir, "/fit_brms_mrt_reg_yr.rds"))
  
} else if (params$Level == "GLB_YR"){
  fit_brms_reg_s <-
    brm(yi | se(sei) ~ 1 + YEAR + (1 | ID),
        chains = 5, iter = params$Iterations, warmup = params$Warmup,
        prior = prior(normal(0,1), class = sd),
        control = controls,
        cores = 5,
        data = subset(es, as.integer(FLAG) == 1),
        open_progress = FALSE,
        seed = 7)
  saveRDS(fit_brms_reg_s, file = paste0(params$Dir, "/fit_brms_mrt_glb_yr.rds"))
  
} else if (params$Level == "REG_CT") {
  fit_brms_reg_s <-
    brm(yi | se(sei) ~ 1 + (1 | REG2) + (1 | REG2:ID),
        chains = 5, iter = params$Iterations, warmup = params$Warmup,
        prior = prior(normal(0,1), class = sd),
        control = controls,
        cores = 5,
        data = subset(es, as.integer(FLAG) == 1),
        open_progress = FALSE,
        seed = 7)
  saveRDS(fit_brms_reg_s, file = paste0(params$Dir, "/fit_brms_mrt_reg_ct.rds"))
  
} else if (params$Level == "GLB_CT"){
  fit_brms_reg_s <-
    brm(yi | se(sei) ~ 1 + (1 | ID),
        chains = 5, iter = params$Iterations, warmup = params$Warmup,
        prior = prior(normal(0,1), class = sd),
        control = controls,
        cores = 5,
        data = subset(es, as.integer(FLAG) == 1),
        open_progress = FALSE,
        seed = 7)
  saveRDS(fit_brms_reg_s, file = paste0(params$Dir, "/fit_brms_mrt_glb_ct.rds"))
}
```

    ## Compiling Stan program...

    ## Start sampling

    ## Warning: There were 16 divergent transitions after warmup. See
    ## https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
    ## to find out why this is a problem and how to eliminate them.

    ## Warning: Examine the pairs() plot to diagnose sampling problems

## Model summary

``` r
summary(fit_brms_reg_s)
```

    ## Warning: There were 16 divergent transitions after warmup. Increasing adapt_delta above 0.9 may help. See
    ## http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup

    ##  Family: gaussian 
    ##   Links: mu = identity; sigma = identity 
    ## Formula: yi | se(sei) ~ 1 + (1 | REG2) + (1 | REG2:ID) 
    ##    Data: subset(es, as.integer(FLAG) == 1) (Number of observations: 11) 
    ##   Draws: 5 chains, each with iter = 5000; warmup = 3000; thin = 1;
    ##          total post-warmup draws = 10000
    ## 
    ## Multilevel Hyperparameters:
    ## ~REG2 (Number of levels: 3) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.77      0.49     0.05     1.90 1.00     3778     2989
    ## 
    ## ~REG2:ID (Number of levels: 11) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     0.73      0.35     0.21     1.53 1.00     4303     5008
    ## 
    ## Regression Coefficients:
    ##           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept    -5.04      0.62    -6.30    -3.83 1.00     4766     4535
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
for (i in 1:length) {
  image <- file.path(pdir, sprintf("r_summary-%s.png", i))
  cat("![](",image,")")
}
```

![](02-fit-hic-mrt-reg-ct_files/figure-gfm/r_summary-1.png)![](02-fit-hic-mrt-reg-ct_files/figure-gfm/r_summary-2.png)![](02-fit-hic-mrt-reg-ct_files/figure-gfm/r_summary-3.png)![](02-fit-hic-mrt-reg-ct_files/figure-gfm/r_summary-4.png)

``` r
if (grepl("YR", params$Level)) {
  P <- plot(conditional_effects(fit_brms_reg_s), points = TRUE, ask=FALSE)
  length <- length(P)
  for (i in 1:length) {
    png(paste0(params$PlotDir, "/r_conditional_effects-", i, ".png"), width = 960)
    plot(P[[i]])
    dev.off()
  }
}
```

``` r
# setwd(params$Dir)
if (grepl("YR", params$Level)) {
  for (i in 1:length){
    image <- file.path(pdir, sprintf("r_conditional_effects-%s.png", i))
    cat("![](",image,")")
  }
}

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
    ##   int prior_only;  // should the likelihood be ignored?
    ## }
    ## transformed data {
    ##   vector<lower=0>[N] se2 = square(se);
    ## }
    ## parameters {
    ##   real Intercept;  // temporary intercept for centered predictors
    ##   vector<lower=0>[M_1] sd_1;  // group-level standard deviations
    ##   array[M_1] vector[N_1] z_1;  // standardized group-level effects
    ##   vector<lower=0>[M_2] sd_2;  // group-level standard deviations
    ##   array[M_2] vector[N_2] z_2;  // standardized group-level effects
    ## }
    ## transformed parameters {
    ##   real sigma = 0;  // dispersion parameter
    ##   vector[N_1] r_1_1;  // actual group-level effects
    ##   vector[N_2] r_2_1;  // actual group-level effects
    ##   real lprior = 0;  // prior contributions to the log posterior
    ##   r_1_1 = (sd_1[1] * (z_1[1]));
    ##   r_2_1 = (sd_2[1] * (z_2[1]));
    ##   lprior += student_t_lpdf(Intercept | 3, -5.3, 2.5);
    ##   lprior += normal_lpdf(sd_1 | 0, 1)
    ##     - 1 * normal_lccdf(0 | 0, 1);
    ##   lprior += normal_lpdf(sd_2 | 0, 1)
    ##     - 1 * normal_lccdf(0 | 0, 1);
    ## }
    ## model {
    ##   // likelihood including constants
    ##   if (!prior_only) {
    ##     // initialize linear predictor term
    ##     vector[N] mu = rep_vector(0.0, N);
    ##     mu += Intercept;
    ##     for (n in 1:N) {
    ##       // add more terms to the linear predictor
    ##       mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_2_1[J_2[n]] * Z_2_1[n];
    ##     }
    ##     target += normal_lpdf(Y | mu, se);
    ##   }
    ##   // priors including constants
    ##   target += lprior;
    ##   target += std_normal_lpdf(z_1[1]);
    ##   target += std_normal_lpdf(z_2[1]);
    ## }
    ## generated quantities {
    ##   // actual population-level intercept
    ##   real b_Intercept = Intercept;
    ## }

# Session info

``` r
sessioninfo::session_info()
```

    ## Warning in system2("quarto", "-V", stdout = TRUE, env = paste0("TMPDIR=", : running command '"quarto"
    ## TMPDIR=C:/Users/fbbu6966/AppData/Local/Temp/RtmpwhXRBk/file15642576453e -V' had status 1

    ## ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────
    ##  setting  value
    ##  version  R version 4.5.0 (2025-04-11 ucrt)
    ##  os       Windows 10 x64 (build 19045)
    ##  system   x86_64, mingw32
    ##  ui       RStudio
    ##  language (EN)
    ##  collate  English_United States.utf8
    ##  ctype    English_United States.utf8
    ##  tz       Europe/Brussels
    ##  date     2025-12-08
    ##  rstudio  2025.05.0+496 Mariposa Orchid (desktop)
    ##  pandoc   3.4 @ C:/Users/fbbu6966/AppData/Local/Programs/RStudio/resources/app/bin/quarto/bin/tools/ (via rmarkdown)
    ##  quarto   ERROR: Unknown command "TMPDIR=C:/Users/fbbu6966/AppData/Local/Temp/RtmpwhXRBk/file15642576453e". Did you mean command "update"? @ C:\\Users\\fbbu6966\\AppData\\Local\\Programs\\RStudio\\RESOUR~1\\app\\bin\\quarto\\bin\\quarto.exe
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────
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
    ## ──────────────────────────────────────────────────────────────────────────────────────────────────────────
