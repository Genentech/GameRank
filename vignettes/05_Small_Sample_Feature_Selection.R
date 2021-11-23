#' ---
#' title: "Small Sample Feature Selection (Vignette)"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{GameRank}
#'   %\VignetteEngine{knitr::rmarkdown}
#'   %\VignetteEncoding{UTF-8}
#' ---
# /* Hit Ctrl + Shift + K in Rstudio to generate html_document for quick view. Note: Vignette YAML header doesn't work then. */
# /* Run knitr::spin("vignettes/01_GameRank.R", format="Rmd") to turn this into a Vignette */

#' 
#' In this vignette, we'll take a look into the case where sample size does not permit a full fledged training:validation:test split evaluation.
#' Recommendations for those situations are resorting to resampling methods for estimating generalization performance. Amongst the most recommended ones is
#' the bootstrap [REF]. Here we'll show how the optimism-adjusted bootstrap estimate can be used in combination with our wrapper methods.
#' 
#' The main element is the avoidance of a training split through by-passing model generation. The model validation step then will make use of the complete dataset
#' and apply the bootstrap method to evaluate the objective function. Using by applying the full data model together with a model trained on the bootstrap data, generalization 
#' error and its downward bias can be estimated and used to obtain a bias-adjusted performance.
#' 
#' For this showcase, we'll use the DLBCL data from TCGA that only comprises 48 tumor samples. As training function, we'll use fn_train_dummy that just returns a string with
#' class 'DummyTrain'. The validation function however will be wrapped into a function that calls the right model generation and model validation functions within a bootstrap procedure on the complete dataset. This comes in extremely costly in terms of computation time.
#' 

#+ eval=TRUE
rm( list=ls() )
library( dplyr )
library( tidyr )
library( purrr )

library( survival )
library( survminer )

devtools::load_all("~/GameRank/")

# setting outputs 
dir_varsel <- "~/GameRank/out/vig5/"
file_varchk <- file.path( "~/GameRank/data", "dlbcl_surv_varchk.rds" )
file_rnd <- file.path( dir_varsel, "dlbcl_surv_rndsel.rds")
file_rnd <- file.path( "~/GameRank/data", "dlbcl_surv_rndsel.rds")
file_gmr <- file.path( dir_varsel, "dlbcl_surv_gmrsel.rds")
file_gmr <- file.path( "~/GameRank/data", "dlbcl_surv_gmrsel.rds")

# load data
load( "~/GameRank/data/tcga_dlbcl_cna_cnv.Rdata" )
dat <- dat %>% filter( grepl( "01A-11D", dat$SampleID_cna ) ) # Filter only one half of samples 

#' Brush up Overall Survival endpoint by imputing missing days by 1+max days survival
max_days <- 1 + max( dat$days_to_death, na.rm=TRUE )
dat <- dat %>% mutate( days_to_death = ifelse( is.na(days_to_death), max_days, days_to_death ) )
# re <- "Surv( days_to_death, vital_status )"
# va <- lst_vars
# vck <- check_variables( dat, re, va )
# save( vck, file = file_varchk )

#+ eval=TRUE
vck <- readRDS( file = file_varchk )

#' Prepare variable selection
re <- vck %>% filter( is_response ) %>% pull( variable )
va <- vck %>% filter( "Perfect"==check_missing ) %>% arrange( desc(mutual_information), desc(entropy), rot.p ) %>% head( 50 ) %>% pull( variable ) %>% unique
dat <- dat %>% dplyr::select( -setdiff( lst_vars, va ) ) # At this point, it makes sense to remove all unused variables for efficiency

#' Quick check on OS before starting
#+ fig.width=7
ft <- survfit( formula( sprintf( "%s ~  1", re ) ), dat )
ggsurvplot( ft, dat )
u <- 1 * 365L # Set landmark time point to be 2-year OS

#' Set the split-sample matrix, where all samples become validation samples, as they'll be used
#' during bootstrap validation
ds <- matrix( 2, nrow=nrow(dat), ncol=1 )

#' Run random variable selection for Cox PH models
#+ eval=FALSE
fn_bt_eval_cox <- ff_fn_eval_bootstrap( fn_train_cox, fn_eval_cox, 20L )
start_time <- Sys.time()
rnd <- random_selection( dat, re, va, fn_train_dummy, fn_bt_eval_cox, 2L, ds, FALSE, 40L, u = u )
end_time <- Sys.time()
dt <- difftime( end_time, start_time )
cat( sprintf( "Random search with bootstrap ran from %s to %s (%1.2f %s) \n", start_time, end_time, dt, units(dt) ) )
saveRDS( rnd, file = file_rnd )
rnd

#' Run GameRank selection for Cox PH models. Please note that instead of the split-sample matrix _dsi_, the split-sample indices,
#' which are used to randomly split observations into training:validation splits per match, are set to 2s [_c(2L,2L)_] indicating
#' validation samples only. _Note: one can shift the fraction of training:validation per match by adding 1s or 2s to the dsi array._
#+ eval=FALSE
start_time <- Sys.time()
gmr <- game_rank( dat, re, va, fn_train_dummy, fn_bt_eval_cox, 2L, dsi = c(2L,2L), FALSE, 6L, 5L, 3L, u = u )
end_time <- Sys.time()
dt <- difftime( end_time, start_time )
cat( sprintf( "GameRank with bootstrap ran from %s to %s (%1.2f %s) \n", start_time, end_time, dt, units(dt) ) )
saveRDS( gmr, file = file_gmr )
gmr
# Save results and render reports
render_random_summary( rnd, dir_varsel )
render_game_rank_summary( gmr, dir_varsel )

#' As the previous steps run a while, let's fast-forward and load pre-computed results.
#+ eval=TRUE
rnd <- readRDS( file = file_rnd )
gmr <- readRDS( file = file_gmr )

#' Make sure the model is identifiable, ie. doesn't include NA coefficients
modr <- fn_train_cox( dat, re, rnd$variable_selections[[1]] ) 
modr
modg <- fn_train_cox( dat, re, gmr$game_rank_selection ) 
modg

#' Little hack to satisfy calibration output report, we just provide 3 x the same dataset as training, validation and test for plotting
#+ eval=FALSE
dat_rpt <- bind_rows( mutate( dat, ds = 1 ), mutate( dat, ds = 2 ) ) %>% bind_rows( mutate( dat, ds = 3 ) )
ds_fin <- matrix( dat_rpt %>% pull( ds ), ncol = 1 )
render_model_calibration_cox( ds_fin, dat_rpt, re, rnd$variable_selections[[1]], k = 1, u = u, dir_varsel, output_file = "rnd_summary" )
render_model_calibration_cox( ds_fin, dat_rpt, re, gmr$game_rank_selection, k = 1, u = u, dir_varsel, output_file = "gmr_summary" )

#' For small sample sizes we can only visualize the overall calibration for both selections:
#+ eval=TRUE, fig.width=7
gplot_predictions_cox( dat, re, gmr$game_rank_selection, modg, u )

# eval=TRUE, fig.width=7
# gplot_predictions_cox( dat, re, rnd$variable_selections[[1]], modr, u )
