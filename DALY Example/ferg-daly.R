## ---- ferg-daly
### FERG2 ESTIMATES - GENERIC DALY SCRIPT
### 10/05/2025: add 'n_samples' to get_sa()
### 17/05/2025: introduce 'per case' metric

## setup
source("../ferg-settings.R")
dirsettings <- dir_copy(hazard)

## process DALY model
daly_mod <- dalymod(file, n_samples)
# str(daly_mod)
saveRDS(daly_mod, sprintf("%s-dalymod-%s.rds", hazard, bd::today()))

## get source attribution
# if(!exists("folder_sa")){
#   folder_sa <- "20250517"
# }
# sa <- get_sa(hazard_sa, folder_sa, n_samples)
sa <- get_sa(hazard_sa, n_samples)
knitr::kable(t(sapply(sa, bd::mean_ci)))

## calculate DALYs per year
out <- NULL

for (year in all_yrs) {
  # population for year
  pop_year <- subset(pop_agg, YEAR == year)
  
  # calculate DALYs
  daly_calc <- dalycalc(daly_mod, year, pop_year, rle)
  
  # aggregate nodes
  daly_calc_agg <- dalycalc_aggregate_nodes(daly_calc)
  
  # apply source attribution
  daly_calc_agg_food <- daly_calc_agg
  for (i in seq_along(daly_calc_agg_food)) {
    reg2 <- subset(FERG2:::countries, ISO3 == names(daly_calc_agg)[i])$SUB2
    daly_calc_agg_food[i] <-
      dalycalc_mult(daly_calc_agg_food[i], sa[[reg2]])
  }
  
  # aggregate age-sex
  daly_calc_agg_age <-
    dalycalc_aggregate_agesex(daly_calc_agg, age_agg, sex_agg)
  daly_calc_agg_age_food <-
    dalycalc_aggregate_agesex(daly_calc_agg_food, age_agg, sex_agg)
  
  # save simulations
  saveRDS(
    daly_calc_agg_age,
    sprintf("%s/%s-dalycalc-%s.rds", dirsettings$local_temp_sim, hazard, year), compress = FALSE)
  saveRDS(
    daly_calc_agg_age_food,
    sprintf("%s/%s-dalycalc-food-%s.rds", dirsettings$local_temp_sim_food, hazard, year), compress = FALSE)
  
  # summaries
  # .. country/sub/reg/glb x age/all
  daly_calc_agg_all <-
    dalycalc_aggregate_agesex(daly_calc_agg, age_all, sex_agg)
  daly_calc_agg_all_food <-
    dalycalc_aggregate_agesex(daly_calc_agg_food, age_all, sex_agg)
  
  # . absolute numbers
  out_yr_nr_all <-
    out_yr_rt_all <-
    rbind(
      dalycalc_summary(daly_calc_agg_age),
      dalycalc_summary(dalycalc_aggregate_country(daly_calc_agg_age, sub)),
      dalycalc_summary(dalycalc_aggregate_country(daly_calc_agg_age, reg)),
      dalycalc_summary(dalycalc_aggregate_country(daly_calc_agg_age, glb)),
      dalycalc_summary(daly_calc_agg_all),
      dalycalc_summary(dalycalc_aggregate_country(daly_calc_agg_all, sub)),
      dalycalc_summary(dalycalc_aggregate_country(daly_calc_agg_all, reg)),
      dalycalc_summary(dalycalc_aggregate_country(daly_calc_agg_all, glb)))
  out_yr_nr_all$SOURCE <- out_yr_rt_all$SOURCE <- "ALL"
  
  out_yr_nr_food <-
    out_yr_rt_food <-
    rbind(
      dalycalc_summary(daly_calc_agg_age_food),
      dalycalc_summary(dalycalc_aggregate_country(daly_calc_agg_age_food, sub)),
      dalycalc_summary(dalycalc_aggregate_country(daly_calc_agg_age_food, reg)),
      dalycalc_summary(dalycalc_aggregate_country(daly_calc_agg_age_food, glb)),
      dalycalc_summary(daly_calc_agg_all_food),
      dalycalc_summary(dalycalc_aggregate_country(daly_calc_agg_all_food, sub)),
      dalycalc_summary(dalycalc_aggregate_country(daly_calc_agg_all_food, reg)),
      dalycalc_summary(dalycalc_aggregate_country(daly_calc_agg_all_food, glb)))
  out_yr_nr_food$SOURCE <- out_yr_rt_food$SOURCE <- "FOOD"
  
  out_yr_nr <- rbind(out_yr_nr_all, out_yr_nr_food)
  out_yr_rt <- rbind(out_yr_rt_all, out_yr_rt_food)
  
  out_yr_nr$METRIC[grepl("_NR", out_yr_nr$MEASURE)] <- "NR"
  out_yr_nr$METRIC[grepl("_INC", out_yr_nr$MEASURE)] <- "CS"
  out_yr_nr$MEASURE <- gsub("_.*", "", out_yr_nr$MEASURE)
  # str(out_yr_nr)
  
  # . rates per 100k
  out_yr_rt <- subset(out_yr_rt, grepl("_NR", out_yr_rt$MEASURE))
  out_yr_rt$VAL_MEAN <- 1e5 * out_yr_rt$VAL_MEAN / out_yr_rt$POP
  out_yr_rt$VAL_MEDIAN <- 1e5 * out_yr_rt$VAL_MEDIAN / out_yr_rt$POP
  out_yr_rt$VAL_LWR <- 1e5 * out_yr_rt$VAL_LWR / out_yr_rt$POP
  out_yr_rt$VAL_UPR <- 1e5 * out_yr_rt$VAL_UPR / out_yr_rt$POP
  out_yr_rt$MEASURE <- gsub("_.*", "", out_yr_rt$MEASURE)
  out_yr_rt$METRIC <- "RT"
  # str(out_yr_rt)
  
  # . compile all
  out_yr <- rbind(out_yr_nr, out_yr_rt)
  out_yr <- out_yr[, c(1, 11, 10, 2:9)]
  out_yr <- cbind(YEAR = year, out_yr)
  # str(out_yr)
  
  out <- rbind(out, out_yr)
}

# Move sim folders from temporary folder to output folder
ok_sim <- safe_copy(dirsettings$local_temp_sim, paste0(dirsettings$executedir, "/SIM"))
if (!ok_sim) warning("Could not copy results from ", dirsettings$local_temp_sim, " to ", paste0(dirsettings$executedir, "/SIM"))
ok_sim_food <- safe_copy(dirsettings$local_temp_sim_food, paste0(dirsettings$executedir, "/SIM_FOOD"))
if (!ok_sim_food) warning("Could not copy results from ", dirsettings$local_temp_sim_food, " to ", paste0(dirsettings$executedir, "/SIM_FOOD"))

# Clean temporary folder
unlink(dirsettings$local_temp_sim, recursive = TRUE, force = TRUE)
unlink(dirsettings$local_temp_sim_food, recursive = TRUE, force = TRUE)

out$MEASURE <- factor(out$MEASURE)
out$METRIC <- factor(out$METRIC)
out$COUNTRY <- factor(out$COUNTRY)
# saveRDS(out, sprintf("%s-dalyout-%s.rds", hazard, bd::today()))
# writexl::write_xlsx(out, sprintf("%s-dalyout-%s.xlsx", hazard, bd::today()))
# str(out)

DT::datatable(
  rownames = FALSE,
  subset(out, AGE == "ALL AGES" & COUNTRY == "GLOBAL" & METRIC == "NR")) |>
  DT::formatRound(columns = 8:12, digits = c(0, rep(3, 4)))
