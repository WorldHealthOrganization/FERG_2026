### FERG2 ESTIMATES - GENERIC SETTINGS
### 20/08/2025

## required packages
library(dalymod)
library(data.table)
library(dplyr)
library(FERG2)

## check package version
cat(sQuote("dalymod"), as.character(packageVersion("dalymod")), "\n")
cat(sQuote("FERG2"), as.character(packageVersion("FERG2")), "\n")

## helper functions
today <- function() format(Sys.Date(), "%Y%m%d")

## year range
all_yrs <- 2000:2021

## define age levels
age_levels <-
  c("[0,1)", "[1,5)", "[5,10)", "[10,15)", "[15,20)", "[20,25)",
    "[25,30)", "[30,35)", "[35,40)", "[40,45)", "[45,50)", "[50,55)",
    "[55,60)", "[60,65)", "[65,70)", "[70,75)", "[75,80)",
    "[80,85)", "[85,Inf)")

## WHO RLE
## https://cdn.who.int/media/docs/default-source/gho-documents/
## global-health-estimates/ghe2021_daly_methods.pdf?sfvrsn=690b16c3_2
rle <-
  expand.grid(
    AGE = age_levels,
    SEX = c("Male", "Female"))
rle$RLE <- 
  c(mean(c(92.65, 92.21)),
    89.74, 85.25, 80.26, 75.27, 70.28, 65.31, 60.34, 55.38, 50.43,
    45.51, 40.61, 35.74, 30.92, 26.21, 21.62, 17.19, 13.08,  7.28)

## population
## .. how does this work for age==0?
pop <- FERG2:::pop
pop$AGEGRP <-
  cut(pop$AGE, c(0, 1, seq(5, 85, 5), Inf),
      include.lowest = FALSE, right = FALSE)
# table(pop$AGEGRP)
# pop_agg <- aggregate(POP ~ AGEGRP + SEX + ISO3 + YEAR, data = pop, FUN = sum)
setDT(pop)
pop_agg <- pop[, .(POP = sum(POP)), .(AGEGRP, SEX, ISO3, YEAR)]
pop_agg <- as.data.frame(pop_agg)
names(pop_agg)[names(pop_agg) == "AGEGRP"] <- "AGE" 
# head(pop_agg)
# levels(pop_agg$AGE)

## local life expectancy
lle <- as.data.frame(FERG2:::life_exp)
lle$YEAR <- as.numeric(lle$YEAR)
lle$AGE <- as.numeric(lle$AGE) # contains '100+'
lle$LE <- as.numeric(lle$LE)
head(lle)
mid <- c(0.5, 3, seq(7.5, 87.5, 5))

# Split data frame into subsets by groups
groups <- split(lle, list(lle$ISO3, lle$YEAR), drop = TRUE)

# Apply approx to each subset and store the results
lle_lst <- lapply(groups, function(subset) {
  interpolated <- approx(x = subset$AGE, y = subset$LE, xout = mid)
  data.frame(
    ISO3 = unique(subset$ISO3),
    YEAR = unique(subset$YEAR),
    AGE = interpolated$x,
    LLE = interpolated$y
  )
})

# Combine all results into a single data frame
lle <- do.call(rbind, lle_lst)
rownames(lle) <- NULL
lle$AGE <- age_levels
head(lle)
tail(lle)

# calculate RLE by broad age groups
pop_agg_le <- merge(merge(pop_agg, rle), lle)
head(pop_agg_le)

get_weighted_rle <-
  function(age) {
    out <-
    pop_agg_le |>
      filter(AGE %in% dalymod:::split_age_string(age)) |>
      group_by(ISO3, YEAR) |>
      summarise(RLE = weighted.mean(RLE, POP), .groups = "drop")
    class(out) <- "data.frame"
    return(out)
  }
get_weighted_lle <-
  function(age) {
    out <-
      pop_agg_le |>
      filter(AGE %in% dalymod:::split_age_string(age)) |>
      group_by(ISO3, YEAR) |>
      summarise(LLE = weighted.mean(LLE, POP), .groups = "drop")
    class(out) <- "data.frame"
    return(out)
  }

yll_ch <- get_weighted_rle("<5")
yll_ad <- get_weighted_rle("5+")

## age and sex groupings
age_agg <-
  list("<5" = c("[0,1)", "[1,5)"),
       "5+" = c("[5,10)", "[10,15)", "[15,20)", "[20,25)", "[25,30)", 
                "[30,35)", "[35,40)", "[40,45)", "[45,50)", "[50,55)",
                "[55,60)", "[60,65)", "[65,70)", "[70,75)", "[75,80)",
                "[80,85)", "[85,Inf)"))
age_agg_all <-
  list("ALL AGES" = names(age_agg))
age_all <-
  list("ALL AGES" = age_levels)
sex_agg <-
  list("Both sexes" = c("Male", "Female"))

## country groupings
glb <- list(GLOBAL = FERG2:::countries$ISO3)

reg_names <- unique(FERG2:::countries$REG2)
reg <- list()
for (i in seq_along(reg_names))
  reg[[reg_names[i]]] <- subset(FERG2:::countries, REG2 == reg_names[i])$ISO3

sub_names <- unique(FERG2:::countries$SUB2)
sub <- list()
for (i in seq_along(sub_names))
  sub[[sub_names[i]]] <- subset(FERG2:::countries, SUB2 == sub_names[i])$ISO3


## get Source Attribution estimates
get_sa <-
  function(hazard, n_samples = NULL, source = "Food", folder = "Source-attribution-estimates-final") {
    # define subregions
    sub <- sub_space <- sort(unique(FERG2:::countries$SUB2))
    sub_space <- gsub("AFR", "AFR ", sub_space)
    sub_space <- gsub("AMR", "AMR ", sub_space)
    sub_space <- gsub("EMR", "EMR ", sub_space)
    sub_space <- gsub("EUR", "EUR ", sub_space)
    sub_space <- gsub("SEAR", "SEAR ", sub_space)
    sub_space <- gsub("WPR", "WPR ", sub_space)
    
    # compile output
    out <- vector("list", length(sub))
    names(out) <- sub_space
    
    # if 100% foodborne
    if (hazard == 1) {
      for (s in sub_space)
        out[[s]] <- 1
      names(out) <- gsub(" ", "", names(out))
      return(out)
    }
    
    # .. otherwise
    for (s in sub_space) {
      # check if file exists
      f <- sprintf(
        "SA/%1$s/%2$s/%2$s_%3$s_SAMPLES.rds",
        folder, hazard, s)
      
      # use estimate if file exists; otherwise force to zero
      if (file.exists(f)) {
        out[[s]] <- 0.01 * readRDS(f)[, source]
        
      } else {
        warning("File not found for ", sQuote(s), ". Estimate forced to zero.")
        out[[s]] <- 0
      }
    }
    names(out) <- gsub(" ", "", names(out))
    
    # replace NA by 1
    out <- lapply(out, function(x) {x[is.na(x)] <- 1;return(x)})
    
    # sample values if needed
    if (!is.null(n_samples)) {
      id <- sample(seq(max(lengths(out))), n_samples, replace = TRUE)
      out <- lapply(out, function(x) x[id])
      out <- lapply(out, function(x) {x[is.na(x)] <- 0;return(x)})
    }
    
    # return output
    return(out)
  }

## render Rmd and include date in output filename
render_today <-
  function(input, ...) {
    rmarkdown::render(
      input,
      output_file = paste0(
        gsub("\\.Rmd", "", input, ignore.case = TRUE),
        "-",
        today(),
        ".html"), 
      envir = globalenv(),
      ...)
  }

## summarise dalycalc_agg_age
# .. country/sub/reg/glb x age/all
dalyout <-
  function(daly_calc_agg) {
    # aggregates
    daly_calc_agg_age <-
      dalycalc_aggregate_agesex(daly_calc_agg, age_agg, sex_agg)
    daly_calc_agg_all <-
      dalycalc_aggregate_agesex(daly_calc_agg, age_all, sex_agg)

    out_yr_nr <-
      out_yr_rt <-
      rbind(
        dalycalc_summary(daly_calc_agg_age),
        dalycalc_summary(dalycalc_aggregate_country(daly_calc_agg_age, sub)),
        dalycalc_summary(dalycalc_aggregate_country(daly_calc_agg_age, reg)),
        dalycalc_summary(dalycalc_aggregate_country(daly_calc_agg_age, glb)),
        dalycalc_summary(daly_calc_agg_all),
        dalycalc_summary(dalycalc_aggregate_country(daly_calc_agg_all, sub)),
        dalycalc_summary(dalycalc_aggregate_country(daly_calc_agg_all, reg)),
        dalycalc_summary(dalycalc_aggregate_country(daly_calc_agg_all, glb)))

    out_yr_nr$METRIC[grepl("_NR", out_yr_nr$MEASURE)] <- "NR"
    out_yr_nr$METRIC[grepl("_INC", out_yr_nr$MEASURE)] <- "CS"
    out_yr_nr$MEASURE <- gsub("_.*", "", out_yr_nr$MEASURE)
    # str(out_yr_nr)
    
    out_yr_rt <- subset(out_yr_rt, grepl("_NR", out_yr_rt$MEASURE))
    out_yr_rt$VAL_MEAN <- 1e5 * out_yr_rt$VAL_MEAN / out_yr_rt$POP
    out_yr_rt$VAL_MEDIAN <- 1e5 * out_yr_rt$VAL_MEDIAN / out_yr_rt$POP
    out_yr_rt$VAL_LWR <- 1e5 * out_yr_rt$VAL_LWR / out_yr_rt$POP
    out_yr_rt$VAL_UPR <- 1e5 * out_yr_rt$VAL_UPR / out_yr_rt$POP
    out_yr_rt$MEASURE <- gsub("_.*", "", out_yr_rt$MEASURE)
    out_yr_rt$METRIC <- "RT"
    # str(out_yr_rt)
    
    out_yr <- rbind(out_yr_nr, out_yr_rt)
    out_yr <- out_yr[, c(1, 10, 2:9)]
    # str(out_yr)
    
    return(out_yr)
  }

# .. country/sub/reg/glb x age_agg/all
dalyout_agg <-
  function(daly_calc_agg) {
    # aggregates
    daly_calc_agg_all <-
      dalycalc_aggregate_agesex(daly_calc_agg, age_agg_all)
    
    out_yr_nr <-
      out_yr_rt <-
      rbind(
        dalycalc_summary(daly_calc_agg),
        dalycalc_summary(dalycalc_aggregate_country(daly_calc_agg, sub)),
        dalycalc_summary(dalycalc_aggregate_country(daly_calc_agg, reg)),
        dalycalc_summary(dalycalc_aggregate_country(daly_calc_agg, glb)),
        dalycalc_summary(daly_calc_agg_all),
        dalycalc_summary(dalycalc_aggregate_country(daly_calc_agg_all, sub)),
        dalycalc_summary(dalycalc_aggregate_country(daly_calc_agg_all, reg)),
        dalycalc_summary(dalycalc_aggregate_country(daly_calc_agg_all, glb)))
    
    out_yr_nr$METRIC[grepl("_NR", out_yr_nr$MEASURE)] <- "NR"
    out_yr_nr$METRIC[grepl("_INC", out_yr_nr$MEASURE)] <- "CS"
    out_yr_nr$MEASURE <- gsub("_.*", "", out_yr_nr$MEASURE)
    # str(out_yr_nr)
    
    out_yr_rt <- subset(out_yr_rt, grepl("_NR", out_yr_rt$MEASURE))
    out_yr_rt$VAL_MEAN <- 1e5 * out_yr_rt$VAL_MEAN / out_yr_rt$POP
    out_yr_rt$VAL_MEDIAN <- 1e5 * out_yr_rt$VAL_MEDIAN / out_yr_rt$POP
    out_yr_rt$VAL_LWR <- 1e5 * out_yr_rt$VAL_LWR / out_yr_rt$POP
    out_yr_rt$VAL_UPR <- 1e5 * out_yr_rt$VAL_UPR / out_yr_rt$POP
    out_yr_rt$MEASURE <- gsub("_.*", "", out_yr_rt$MEASURE)
    out_yr_rt$METRIC <- "RT"
    # str(out_yr_rt)
    
    out_yr <- rbind(out_yr_nr, out_yr_rt)
    out_yr <- out_yr[, c(1, 10, 2:9)]
    # str(out_yr)
    
    return(out_yr)
  }

# summarize samples
summarize <-
  function(x, year = "2021") {
    out <- t(sapply(subset(x, YEAR == year)$SAMPLES, mean_ci))
    rownames(out) <- subset(x, YEAR == year)$COUNTRY
    knitr::kable(out)
  }

# get age distribution for HIC/diarrhea
split_prop_age <-
  function(prop, age, country, year) {
    setDT(pop_agg)
    setkey(pop_agg, ISO3, YEAR, AGE)
    p <- pop_agg[.(country, as.numeric(year), dalymod:::split_age_string(age))]
    p$POP_PROP <- p$POP / sum(p$POP)
    p$PROP <- prop * p$POP_PROP
    return(data.frame(p))
  }
split_prop_both <-
  function(age, country, year) {
    do.call(
      "rbind",
      mapply(split_prop_age,
             age, names(age),
             MoreArgs = list(country = country, year = year),
             SIMPLIFY = FALSE))
  }
get_prob_ch <-
  function(x, country, year) {
    p_agesex <- split_prop_both(x, country, as.numeric(year))
    sum(subset(p_agesex, AGE %in% c("[0,1)", "[1,5)"))$PROP)
  }
get_prob_ch_country <-
  function(x, country) {
    prob_ch <-
      data.frame(
        PROB = c(get_prob_ch(x[[1]], country, names(x[1])),
                 get_prob_ch(x[[2]], country, names(x[2])),
                 get_prob_ch(x[[3]], country, names(x[3]))),
        YEAR = as.numeric(names(x)))
    prob_ch_full <-
      data.frame(YEAR = all_yrs)
    prob_ch_full$PROB <-
      approx(prob_ch$YEAR, prob_ch$PROB, prob_ch_full$YEAR)$y
    prob_ch_full$COUNTRY <- country
    return(prob_ch_full)
  }
get_prob_all_country <-
  function(x, country) {
    prob_all <-
      rbind(split_prop_both(x[[1]], country, names(x[1])),
            split_prop_both(x[[2]], country, names(x[2])),
            split_prop_both(x[[3]], country, names(x[3])))
    
    agegrp <- names(x[[1]])
    agegrp_map <- lapply(agegrp, dalymod:::split_age_string)
    agegrp_map <- data.frame(
      AGEGRP = rep(agegrp, lengths(agegrp_map)),
      AGE = unlist(agegrp_map, use.names = FALSE))

    prob_all <- merge(prob_all, agegrp_map)
    prob_all <- aggregate(PROP~AGEGRP+ISO3+YEAR, prob_all, sum)
    
    prob_all_full <-
      group_by(prob_all, AGEGRP) |>
      reframe(P = as.data.frame(approx(YEAR, PROP, all_yrs))) |>
      as.data.frame()
    prob_all_full$YEAR <- prob_all_full$P$x
    prob_all_full$PROB <- prob_all_full$P$y
    prob_all_full$P <- NULL
    prob_all_full$COUNTRY <- country
    
    return(prob_all_full)
  }

## Settings regarding saving sim files as directly saving them gives issues
# Settings
dir_copy <- function(hazard){
  dirsettings <- list()
  dirsettings$id_sim <- paste0(hazard, "_SIM")
  dirsettings$id_sim_food <- paste0(hazard, "_SIM_FOOD")
  dirsettings$local_temp_sim <- file.path(tempdir(),  dirsettings$id_sim)
  dirsettings$local_temp_sim_food <- file.path(tempdir(),  dirsettings$id_sim_food)
  dir.create( dirsettings$local_temp_sim, recursive = TRUE, showWarnings = FALSE)
  dir.create( dirsettings$local_temp_sim_food, recursive = TRUE, showWarnings = FALSE)
  dirsettings$executedir <- getwd()
  return(dirsettings)
}


# sim files are saved on a folder on the computer before moving it due to connection issues
safe_copy <- function(src_dir, dest_dir, tries = 3) {
  files <- list.files(src_dir, full.names = TRUE)
  for (i in seq_len(tries)) {
    ok <- all(file.copy(files, dest_dir, overwrite = TRUE, recursive = TRUE))
    if (ok) return(TRUE)
    Sys.sleep(2)
  }
  FALSE
}
