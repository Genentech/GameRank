#' ---
#' title: "Model Evaluation Algorithm (Vignette)"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{GameRank}
#'   %\VignetteEngine{knitr::rmarkdown}
#'   %\VignetteEncoding{UTF-8}
#' ---
# /* Hit Ctrl + Shift + K in Rstudio to generate html_document for quick view. Note: Vignette YAML header doesn't work then. */
# /* Run knitr::spin("vignettes/01_GameRank.R", format="Rmd") to turn this into a Vignette */


#' 
#' Once the best feature combination and the best trained model is determine, this model needs to be checked. There are two steps needed: first, 
#' plots for it's calibration have to be done for training, validation and test splits, and second, the data used to generate it has to checked for 
#' influential observations, ie those that more than others determine the model's parameters. These can be identified by the jackknife approach, that
#'  is by comparing models and their predictions that were trained by leaving an single observation out to the model obtained by the complete dataset. 
#'  This approach is found in the _dffits_ and _dfbetas_ functions for 'lm' and 'glm' classes in R. Here we'll show a generalized approach.
#'  


#' First let's load the some program settings
#+ eval=TRUE
rm( list=ls() )
library( dplyr )
library( tidyr )
library( purrr )

library( survival )
library( survminer )

devtools::load_all("~/GameRank/")

# setting outpu1ts 
dir_varsel <- "~/GameRank/out/vig4/"
file_varchk <- file.path( "~/GameRank/data", "ucec_surv_varchk.rds" )
file_gmr <- file.path( dir_varsel, "ucec_surv_gmr.rds")
file_gmr <- file.path( "~/GameRank/data", "ucec_surv_gmr.rds")
file_mck <- file.path( dir_varsel, "ucec_surv_modchk.rds")
file_mck <- file.path( "~/GameRank/data", "ucec_surv_modchk.rds")

# load data
load( "~/GameRank/data/tcga_ucec_cna_cnv.Rdata" )
dat <- dat %>% filter( Sample %in% c("01","02") )  # Filter only tumor samples 
vck <- readRDS( file = file_varchk )

# Brush up Overall Survival endpoint by imputing missing days by 1+max days survival
max_days <- 1 + max( dat$days_to_death, na.rm=TRUE )
dat <- dat %>% mutate( days_to_death = ifelse( is.na(days_to_death), max_days, days_to_death ) )

# Prepare variable selection
re <- vck %>% 
  filter( is_response ) %>% 
  pull( variable )
va <- vck %>% 
  filter( "Perfect"==check_missing ) %>% 
  arrange( desc(mutual_information), desc(entropy), rot.p ) %>% 
  head( 100 ) %>% 
  pull( variable ) %>% 
  unique

#' Let's have a quick check on OS before starting
#+ fig.width=7
ft <- survfit( formula( sprintf( "%s ~  1", re ) ), dat )
ggsurvplot( ft, dat )
u <- 1 * 365L # Set landmark time point to be 2-year OS

#' Randomly split data into training (60%), validation (20%) and test (20%)
rr <- rep_len( c(1,1,1,1,2,2,3,3), length.out = nrow(dat) )
rr <- rr[ order( runif( length( rr ) )) ]
table( rr )
df_trval <- dat[which( rr %in% c(1,2) ),]
df_tst   <- dat[which( rr %in% c(3) ),]
ds <- matrix( rr[ which( rr %in% c(1,2) ) ], ncol=1 )
ds_fin <- matrix( rr, ncol=1 )

#' Run the GameRank wrapper algorithm to obtain a variable selection. As this part runs
#' a while, we will skip forward and load some pre-computed results.
#+ eval=FALSE
start_time <- Sys.time()
gmr <- game_rank( df_trval, re, va, fn_train_cox, fn_eval_cox, 5L, dsi = c(1L,2L), FALSE, 7L, 7L, 7L, u = u )
end_time <- Sys.time()
dt <- difftime( end_time, start_time )
cat( sprintf( "GameRank feature selection ran from %s to %s (%1.2f %s) \n", start_time, end_time, dt, units(dt) ) )
saveRDS( gmr, file = file_gmr )
render_game_rank_summary( gmr, dir_varsel )
render_model_calibration_cox( ds_fin, dat, re, gmr$game_rank_selection, 1, u, dir_varsel )

#' Loading the result and creating a model, we see that it has a NA as coefficient.
#+ eval=TRUE
gmr <- readRDS( file = file_gmr ) 
mod <- fn_train_cox( df_trval, re, gmr$game_rank_selection )
mod$coefficients[which(is.na(mod$coefficients))] <- 0 # If coefficients are NA, we can fix them by setting to 0
mod

#' Plotting model calibration for training+validation
#+ fig.width=7
gplot_predictions_cox( df_trval, re, gmr$game_rank_selection, mod, u  )

#' Plotting model calibration for hold-out test.
#+ fig.width=7
gplot_predictions_cox( df_tst, re, gmr$game_rank_selection, mod, u  )

#' # Check for influential data points
#'
#' for lm and glm models there are basic R functions, like dffits, cooks.model and influence.measures to evaluate 
#' for influential observations. Down we propose an attempt for Cox PHP models.
#' 
#+ eval=FALSE
start_time <- Sys.time()
tips <- influential_observations( df_trval, re, gmr$game_rank_selection, fn_train_cox, fn_eval_cox, fn_infl_coefficients, fn_predict_cox, u )
end_time <- Sys.time()
dt <- difftime( end_time, start_time )
cat( sprintf( "Determination of influential points ran from %s to %s (%1.2f %s) \n", start_time, end_time, dt, units(dt) ) )
render_influence_summary( tips, dir_varsel )
saveRDS( tips, file = file_mck )

#' Since the influential_observations functions runs a while, we'll load some precomputed example.
#+ eval=TRUE
tips <- readRDS( file = file_mck )
tips %>% filter( is_influential ) %>% dplyr::select( row, is_influential, is_influential_co, tidyr::everything() )

#' There are 4 influential observations, on influencing dffit and two of the dfbetas.

#' Let's have a look and plot the distributions of the dffits and dfbeta values:
#+ fig.width=7, fig.height=14
val_cols <- c("deval", "dffit", grep( "_dfbeta$", colnames(tips), value=TRUE ) )
tips %>%
  filter( !is.na(row) ) %>%
  dplyr::select( all_of( c("row", val_cols ) ) ) %>%
  tidyr::pivot_longer( cols = val_cols, names_to = "measure" ) %>%
  filter( !is.na(value) ) %>%
  ggplot( aes( x=value, y=..density.., color=measure, group=measure ) ) +
  facet_wrap( c("measure"), ncol=2 ) +
  geom_density( bw = "ucv" ) +
  geom_rug( aes( x=value, y=NA ), sides = "t" )

#' All distributions look quite centered with a view outlier exceptions.

