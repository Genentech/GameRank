#' ---
#' title: "GameRank Algorithm (Vignette)"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{GameRank}
#'   %\VignetteEngine{knitr::rmarkdown}
#'   %\VignetteEncoding{UTF-8}
#' ---
# /* Hit Ctrl + Shift + K in Rstudio to generate html_document for quick view. Note: Vignette YAML header doesn't work then. */
# /* Run knitr::spin("vignettes/01_GameRank.R", format="Rmd") to turn this into a Vignette */

#+ setup, include=FALSE
rm( list=ls() )
library( dplyr )
library( tidyr )
library( purrr )

library( survival )
library( survminer )

devtools::load_all("~/GameRank/")

# setting outputs 
dir_varsel <- "out/vig5/"
file_varchk <- system.file( "data", "dlbcl_surv_varchk.rds", package = "GameRank" )
file_rnd <- file.path( dir_varsel, "dlbcl_surv_rndsel.rds")
file_gmr <- file.path( dir_varsel, "dlbcl_surv_gmrsel.rds")

# load data
load( "data/tcga_dlbcl_cna_cnv.Rdata" )
dat <- dat %>% filter( grepl( "01A-11D", dat$SampleID_cna ) ) # Filter only one half of samples 
vck <- readRDS( file = file_varchk )

# Brush up Overall Survival endpoint by imputing missing days by 1+max days survival
max_days <- 1 + max( dat$days_to_death, na.rm=TRUE )
dat <- dat %>% mutate( days_to_death = ifelse( is.na(days_to_death), max_days, days_to_death ) )

# Prepare variable selection
re <- vck %>% filter( is_response ) %>% pull( variable )
va <- vck %>% filter( "Perfect"==check_missing ) %>% arrange( desc(mutual_information), desc(entropy), rot.p ) %>% head( 200 ) %>% pull( variable ) %>% unique
dat <- dat %>% dplyr::select( -setdiff( lst_vars, va ) ) # At this point, it makes sense to remove all unused variables for efficiency

# Quick check on OS before starting
ft <- survfit( formula( sprintf( "%s ~  1", re ) ), dat )
ggsurvplot( ft, dat )
u <- 1 * 365L # Set landmark time point to be 2-year OS

# All samples become validation samples, they'll be used during bootstrap validation
ds <- matrix( 2, nrow=nrow(dat), ncol=1 )

# Run random variable selection for Cox PH models
fn_bt_eval_cox <- ff_fn_eval_bootstrap( fn_train_cox, fn_eval_cox, 20L )
rnd <- random_selection( dat, re, va, fn_train_dummy, fn_bt_eval_cox, 2L, ds, FALSE, 40L, u = u )
rnd

gmr <- game_rank( dat, re, va, fn_train_dummy, fn_bt_eval_cox, 2L, dsi = c(2L,2L), FALSE, 4L, 5L, 5L, u = u )
gmr


# Make sure the model is identifyable, ie. doesn't include NA coefficients
mod <- fn_train_cox( dat, re, rnd$variable_selections[[1]] ) 
mod

# Save results and render reports
saveRDS( rnd, file = file_rnd )
saveRDS( gmr, file = file_gmr )
render_random_summary( rnd, dir_varsel )
render_game_rank_summary( gmr, dir_varsel )

rnd$agg_results %>% filter( opt )
gmr$game_rank_selection

# Little hack to satisfy calibration output report, we just provide 3 x the same dataset as training, validation and test for plotting
dat_rpt <- bind_rows( mutate( dat, ds = 1 ), mutate( dat, ds = 2 ) ) %>% bind_rows( mutate( dat, ds = 3 ) )
ds_fin <- matrix( dat_rpt %>% pull( ds ), ncol = 1 )

render_model_calibration_cox( ds_fin, dat_rpt, re, rnd$variable_selections[[1]], k = 1, u = u, dir_varsel )
render_model_calibration_cox( ds_fin, dat_rpt, re, gmr$game_rank_selection, k = 1, u = u, dir_varsel )

