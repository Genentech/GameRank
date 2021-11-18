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
dir_varsel <- "out/vig3/"
file_varchk <- system.file( "data", "ucec_surv_varchk.rds", package = "GameRank" )
file_rnd <- file.path( dir_varsel, "ucec_surv_rndsel.rds")

# load data
load( "data/tcga_ucec_cna_cnv.Rdata" )
dat <- dat %>% filter( Sample %in% c("01","02") )  # Filter only tumor samples 
vck <- readRDS( file = file_varchk )

# Brush up Overall Survival endpoint by imputing missing days by 1+max days survival
max_days <- 1 + max( dat$days_to_death, na.rm=TRUE )
dat <- dat %>% mutate( days_to_death = ifelse( is.na(days_to_death), max_days, days_to_death ) )

# Prepare variable selection
re <- vck %>% filter( is_response ) %>% pull( variable )
va <- vck %>% filter( "Perfect"==check_missing ) %>% arrange( desc(mutual_information), desc(entropy), rot.p ) %>% head( 200 ) %>% pull( variable ) %>% unique

# Quick check on OS before starting
ft <- survfit( formula( sprintf( "%s ~  1", re ) ), dat )
ggsurvplot( ft, dat )
u <- 2 * 365L # Set landmark time point to be 2-year OS

# Randomly split data into training (60%), validation (20%) and test (20%)
rr <- rep_len( c(1,1,1,1,2,2,3,3), length.out = nrow(dat) )
rr <- rr[ order( runif( length( rr ) )) ]
table( rr )
df_trval <- dat[which( rr %in% c(1,2) ),]
df_tst   <- dat[which( rr %in% c(3) ),]
ds <- matrix( rr[ which( rr %in% c(1,2) ) ], ncol=1 )
ds_fin <- matrix( rr, ncol=1 )

# Run random variable selection for Cox PH models
rnd <- random_selection( df_trval, re, va, fn_train_cox, fn_eval_cox, 5L, ds, FALSE, 100L, u = u )
rnd

# Save results and render reports
saveRDS( rnd, file = file_rnd )
render_random_summary( rnd, dir_varsel )
render_model_calibration_cox( ds_fin, dat, re, rnd$variable_selections[[1]], k = 1, u = u, dir_varsel )




