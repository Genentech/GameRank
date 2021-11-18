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
dir_varsel <- "out/vig4/"
file_varchk <- system.file( "data", "ucec_surv_varchk.rds", package = "GameRank" )
file_gmr <- file.path( dir_varsel, "ucec_surv_gmr.rds")
file_mck <- file.path( dir_varsel, "ucec_surv_modchk.rds")

# load data
load( "data/tcga_ucec_cna_cnv.Rdata" )
dat <- dat %>% filter( Sample %in% c("01","02") )  # Filter only tumor samples 
vck <- readRDS( file = file_varchk )

# Brush up Overall Survival endpoint by imputing missing days by 1+max days survival
max_days <- 1 + max( dat$days_to_death, na.rm=TRUE )
dat <- dat %>% mutate( days_to_death = ifelse( is.na(days_to_death), max_days, days_to_death ) )

# Prepare variable selection
re <- vck %>% filter( is_response ) %>% pull( variable )
va <- vck %>% filter( "Perfect"==check_missing ) %>% arrange( desc(mutual_information), desc(entropy), rot.p ) %>% head( 100 ) %>% pull( variable ) %>% unique

# Quick check on OS before starting
ft <- survfit( formula( sprintf( "%s ~  1", re ) ), dat )
ggsurvplot( ft, dat )
u <- 1 * 365L # Set landmark time point to be 2-year OS

# Randomly split data into training (60%), validation (20%) and test (20%)
rr <- rep_len( c(1,1,1,1,2,2,3,3), length.out = nrow(dat) )
rr <- rr[ order( runif( length( rr ) )) ]
table( rr )
df_trval <- dat[which( rr %in% c(1,2) ),]
df_tst   <- dat[which( rr %in% c(3) ),]
ds <- matrix( rr[ which( rr %in% c(1,2) ) ], ncol=1 )
ds_fin <- matrix( rr, ncol=1 )

gmr <- game_rank( df_trval, re, va, fn_train_cox, fn_eval_cox, 5L, dsi = c(1L,2L), FALSE, 7L, 7L, 7L, u = u )
gmr

mod <- fn_train_cox( df_trval, re, gmr$game_rank_selection )
mod

gplot_predictions_cox( df_trval, re, gmr$game_rank_selection, mod, u  )
gplot_predictions_cox( df_tst, re, gmr$game_rank_selection, mod, u  )

# Check for influential data points
#
# for lm and glm models there are basic R functions, like dffits, cooks.model and influence.measures to evaluate 
# for influential observations. Down we propose an attempt for Cox PHP models.
# 
tips <- influential_observations( df_trval, re, gmr$game_rank_selection, fn_train_cox, fn_eval_cox, fn_infl_coefficients, fn_predict_cox, u )
tips %>% filter( is_influential )

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

saveRDS( gmr, file = file_gmr )
saveRDS( tips, file = file_mck )
render_game_rank_summary( gmr, dir_varsel )
render_model_calibration_cox( ds_fin, dat, re, gmr$game_rank_selection, 1, u, dir_varsel )
render_influence_summary( tips, dir_varsel )
