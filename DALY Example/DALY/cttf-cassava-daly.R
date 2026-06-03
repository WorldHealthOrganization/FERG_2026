#' ---
#' title: "cttf-cassava-daly"
#' output:
#'  github_document:
#'    toc: true
#'    toc_depth: 2
#'    html_preview: true
#' date: "`r Sys.Date()`"
#' ---

#' # Settings
hazard <- "cttf-cassava"
hazard_sa <- 1
file <- "ferg2-daly-cassava-20250415.xlsx"
n_samples <- 1e3
set.seed(264)
source("../ferg-settings.R")

#' # DALY
source("../ferg-daly.R") 
