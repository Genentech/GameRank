#' ---
#' title: "Feature Selection (Vignette)"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{GameRank}
#'   %\VignetteEngine{knitr::rmarkdown}
#'   %\VignetteEncoding{UTF-8}
#' ---
# /* Hit Ctrl + Shift + K in Rstudio to generate html_document for quick view. Note: Vignette YAML header doesn't work then. */
# /* Run knitr::spin("vignettes/01_GameRank.R", format="Rmd") to turn this into a Vignette */

#' 
#' After having constructed better features and selected a reasonable subset, the task is now to start building and evaluating predictive models.
#' Most wrapper algorithms for feature selection make use of two functions: one to generate a model from data, and one to evaluate its generalization performance.
#' Usually this requires a split of the data into a training and a validation part (split set evaluation) of sufficient sizes. We'll consider a small sample
#' scenario in a separate vignette.
#' 
#' Since wrappers will provide selections that are biased to optimal performance on the validation set, a third split is necessary: the hold-out test set. This set
#' is used just once to measure the performance of the finally selected model. It is advisable to ensure via a sample size calculation for the width of confidence intervals
#' that this set is sufficiently large to obtain a reasonable precision.
#' To mitigate the wrapper being biased towards a single validation set, we provide the option to define _multiple parallel training:validation splits_. At each step the wrapper
#' will evaluate multiple models created on varying training:validation splits and average those validation performances for feature evaluation. 
#' These multiple parallel splits can be defined by a matrix with columns providing a 1 or 2 per sample row to denote training (=1) or validation (=2) usage of that observation.
#' 
#' As a first attempt, we propose to explore the feature space by the most robust wrapper algorithm: the random search wrapper. Given a target size $m$ of the feature combination
#' it generates a random search matrix of disjoint feature combinations for evaluation. Then it evaluates the full or just a subsample, if a maximum number of evaluations 
#' would be exceeded, of these combinations and returns a database of selection performances and the best one. Using the random search is easy to get an first idea of how well 
#' models perform in predicting the outcome furthermore selecting the best combinations provides a refined list of candidate features together with those that were not used yet.
#' 
#' In this example we'll search for a predictive Cox Proportional Hazards model predicting the 2-year Overall Survival probability. We'll estimate calibration, that is the closeness 
#' of predicted versus observed probability, following the method from [REF]. The wrapper will be run to minimize the absolute difference between predicted and 
#' observed 1-year survival probability.
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
dir_varsel <- "~/GameRank/out/vig3/"
file_varchk <- file.path( "~/GameRank/data", "ucec_surv_varchk.rds" )
file_rnd <- file.path( dir_varsel, "ucec_surv_rndsel.rds")
file_rnd <- file.path( "~/GameRank/data", "ucec_surv_rndsel.rds")

# load data
load( "~/GameRank/data/tcga_ucec_cna_cnv.Rdata" )
dat <- dat %>% filter( Sample %in% c("01","02") )  # Filter only tumor samples 
vck <- readRDS( file = file_varchk )

# Brush up Overall Survival endpoint by imputing missing days by 1+max days survival
max_days <- 1 + max( dat$days_to_death, na.rm=TRUE )
dat <- dat %>% mutate( days_to_death = ifelse( is.na(days_to_death), max_days, days_to_death ) )

# Prepare variable selection
re <- vck %>% filter( is_response ) %>% pull( variable )
va <- vck %>% 
  filter( "Perfect"==check_missing ) %>% 
  arrange( desc(mutual_information), desc(entropy), rot.p ) %>% 
  head( 200 ) %>%
  pull( variable ) %>%
  unique

#' After some setup code, let's have a quick look at the survival curve.
#+ fig.width=7
# Quick check on OS before starting
ft <- survfit( formula( sprintf( "%s ~  1", re ) ), dat )
ggsurvplot( ft, dat )
u <- 2 * 365L # Set landmark time point to be 2-year OS

#' Next we'll call the first feature selection wrapper algorithm available: the random search.
#' This algorithm is the most robust and error tolerant one. It generates a sampling matrix
#' of n_eval rows each containing a disjoint random feature combination for evaluation.
#' 
# Randomly split data into training (60%), validation (20%) and test (20%)
rr <- rep_len( c(1,1,1,1,2,2,3,3), length.out = nrow(dat) ) # 1=training, 2=validation, 3=test
rr <- rr[ order( runif( length( rr ) )) ]
table( rr )
df_trval <- dat[which( rr %in% c(1,2) ),]
df_tst   <- dat[which( rr %in% c(3) ),]
ds <- matrix( rr[ which( rr %in% c(1,2) ) ], ncol=1 )
ds_fin <- matrix( rr, ncol=1 )

#' The following code runs for a while, thus we skip it by loading its results.
#+ eval=FALSE
# Run random variable selection for Cox PH models
start_time <- Sys.time()
rnd <- random_selection( df_trval, re, va, fn_train_cox, fn_eval_cox, 5L, ds, FALSE, 100L, u = u )
end_time <- Sys.time()
dt <- difftime( end_time, start_time )
cat( sprintf( "Random search feature selection ran from %s to %s (%1.2f %s) \n", start_time, end_time, dt, units(dt) ) )
saveRDS( rnd, file = file_rnd )
# Save results and render reports
render_random_summary( rnd, dir_varsel )
render_model_calibration_cox( ds_fin, dat, re, rnd$variable_selections[[1]], k = 1, u = u, dir_varsel )

#' Load results from random_selection run:
#+ eval=TRUE
rnd <- readRDS( file = file_rnd )

#' All wrappers return a list with their high-level input parameters, time used for computation and
#' then datasets for each individual evaluation (that is per parallel split), aggregated evaluations
#' and finally the best obtained feature selection.
#' 
#' The best selection found consists of the following features
sel <- unlist( rnd$variable_selections )
sel

#' Let's get the model:
set.seed(23)
mod <- fn_train_cox( dat[which(ds_fin[,1] %in% c(1,2)),], re, sel )
summary( mod )

#' 
#' Let's plot it's calibration on training, validation and test:
#+ fig.width=7
gplot_predictions_cox( dat[which(1==ds_fin[,1]),], re, sel, mod, u )

#+ fig.width=7
gplot_predictions_cox( dat[which(2==ds_fin[,1]),], re, sel, mod, u )

#+ fig.width=7
gplot_predictions_cox( dat[which(3==ds_fin[,1]),], re, sel, mod, u )

#' Note: Please note that non-linear models like the logistic regression or the
#' Cox PH models are non-convex optimization problems and are *not* guaranteed to
#' always converge to the same result.
#' Model robustness should also be a criterion for subsequent choice but it
#' was not part of the random variable selection wrapper.
#' _Rerunning these graphs a couple of times will give varying results except for
#' selections that lead to a unique global optimum._
#' 
#' The best results (possibly aggregated across multiple parallel evaluations) are:
rnd$agg_results %>% arrange( mean_validation ) %>% head

#' The bare dataset of evaluations, including a column for the split defined in rnd$ds used, is:
rnd$results %>% head 

#' There is a list of other wrapper algorithms that also use the same interface as random_selection:
#' 
#' - forward,
#' - backward,
#' - bidirectional, and 
#' - lrsearch
#' 
#' All of them make use of one or more training:validation splits defined by a ds-matrix consisting of columns with 1s and 2s'.
