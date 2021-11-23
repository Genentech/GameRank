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

#' 
#' As part of this package, we release a novel wrapper algorithm for variable ranking. This vignette is completely devoted to this algorithm and how it can be used.
#' In the introduction of the first vignette, we mentioned that wrappers do make use of a validation split to estimate and measure the objective: the prediction error.
#' However through repeated iteration, wrappers will become biased towards _this_ validation set, and therefore we proposed to instead use a set of _parallel_ training:validation splits of the _same_ data to average out the error for choosing the split.
#' 
#' The algorithm shown now, named GameRank, doesn't require a split into training and validation. Instead it generates its own training:validation splits while evaluating 
#' pairs of random feature combinations against each other. During each round it just records how often which combination performed better compared to another, similar to 
#' football or basket ball teams playing matches. After a sufficient number of matches have been recorded, GameRank uses a maximum-likelihood ranking model to estimate the
#' team contribution of each individual feature. The basic assumption is that the sum of the individual contributions determines the team strengths that translates into the 
#' winning probability of one team over the other. Since this model is a maximum-likelihood model, it comes in with all the large sample advantages of this method, that it is
#' asymptotic consistency. Since training:validation splits are not kept fixed but rather sampled afresh each time, GameRank avoids becoming biased against a fixed combination.
#' Together with the possibility to parallelize generation of feature set combinations data across computing nodes, this may render GameRank a large-scale wrapper algorithm.
#' 
#' Here is how it is used.
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
dir_varsel <- "~/GameRank/out/vig6/"
file_varchk <- file.path( "~/GameRank/data", "ucec_surv_varchk.rds" )
file_gmr <- file.path( dir_varsel, "ucec_surv_gmrsel.rds")
file_gmr <- file.path( "~/GameRank/data", "ucec_surv_gmrsel.rds")

# load data
load( "~/GameRank/data/tcga_ucec_cna_cnv.Rdata" )
dat <- dat %>% filter( Sample %in% c("01","02") ) # Filter only tumor samples 
vck <- readRDS( file = file_varchk )

# Brush up Overall Survival endpoint by imputing missing days by 1+max days survival
max_days <- 1 + max( dat$days_to_death, na.rm=TRUE )
dat <- dat %>% mutate( days_to_death = ifelse( is.na(days_to_death), max_days, days_to_death ) )

#' Prepare variable selection
re <- vck %>% filter( is_response ) %>% pull( variable )
va <- vck %>% filter( "Perfect"==check_missing ) %>% arrange( desc(mutual_information), desc(entropy), rot.p ) %>% head( 200 ) %>% pull( variable ) %>% unique

#' Quick check on OS before starting
#+ fig.width=7
ft <- survfit( formula( sprintf( "%s ~  1", re ) ), dat )
ggsurvplot( ft, dat )
u <- 2 * 365L # Set landmark time point to be 2-year OS

#' Randomly split data into training (60%), validation (20%) and test (20%)
rr <- rep_len( c(1,1,1,1,2,2,3,3), length.out = nrow(dat) )
rr <- rr[ order( runif( length( rr ) )) ]
table( rr )
df_trval <- dat[which( rr %in% c(1,2) ),]
df_tst   <- dat[which( rr %in% c(3) ),]
ds <- matrix( rr[ which( rr %in% c(1,2) ) ], ncol=1 )
ds_fin <- matrix( rr, ncol=1 )

#' Run GameRank variable selection for Cox PH models on the combined training+validation splits.
#+ eval=FALSE
gmr <- game_rank( df_trval, re, va, fn_train_cox, fn_eval_cox, 5L, c(1L,2L), FALSE, 7L, 10L, 5L, u = u )
# Save results and render reports
render_game_rank_summary( gmr, dir_varsel )
render_model_calibration_cox( ds_fin, dat, re, gmr$variable_selections[[1]], k = 1, u = u, dir_varsel )
saveRDS( gmr, file = file_gmr )

#+ eval=TRUE
gmr <- readRDS( file = file_gmr )

#' Key algorithm parameters are the team_size, that is how many variables are chosen per match and compared against another
#' selection of equal size. Also important is the number of round and the number of matches every variable needs to 
#' participate in. These parameters influence the match matrix.
gmr$team_size
gmr$rounds
gmr$min_matches_per_var

#' GameRank also returns the time it spends in each stage: generating the sample matrix, evaluating matches,
#' computing the maximum-likelihood estimator to rank variables, and the overall time.
dt <- difftime( gmr$end, gmr$start, units="auto" )
cat(sprintf("GameRank ran from %s to %s for %1.4f [%s]. \n", gmr$start, gmr$end, dt, units(dt) ))
dt <- difftime( gmr$match_matrix_time, gmr$start, units="auto" )
cat(sprintf("Time to sample match matrix was %1.4f [%s]. \n", dt, units(dt) ))
dt <- difftime( gmr$match_played_time, gmr$match_matrix_time, units="auto" )
cat(sprintf("Time to evaluate matches was %1.4f [%s]. \n", dt, units(dt) ))
dt <- difftime( gmr$fit_time, gmr$match_played_time, units="auto" )
cat(sprintf("Time to fit maximum-likelihood ranking model was %1.4f [%s]. \n", dt, units(dt) ))

#' The final selection determined by the parameter m is returned
gmr$game_rank_selection

#' A table with the individual variable strengths, their standard error and whether they would be selected (if they have
#' a positive contribution to the team):
gmr$variable_ranking

#' The match matrix can be obtained also: The n.pos and n.neg columns denote how often the (+1)-team and how ofthen the (-1)-team won. The
#' +1s and -1s indicate which variable was in the respective team:
gmr$match_results[1:10,1:10]
