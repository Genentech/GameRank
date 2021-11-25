#' ---
#' title: "GameRank Vignette"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{GameRank}
#'   %\VignetteEngine{knitr::rmarkdown}
#'   %\VignetteEncoding{UTF-8}
#' ---
# /* Hit Ctrl + Shift + K in Rstudio to generate html_document for quick view. Note: Vignette YAML header doesn't work then. */
# /* Run knitr::spin("vignettes/01_GameRank.R", format="Rmd") to turn this into a Vignette */

#+ setup, eval=TRUE, include=FALSE
rm( list=ls() )
library( dplyr )
library( tidyr )
library( purrr )
library( ggplot2 )

devtools::load_all("~/GameRank/")

#' # Introduction
#' This vignette will show how to perform wrapper-based feature selection using the GameRank package.
#' In an usual feature selection scenario the following individual steps will be performed
#' 
#' 1. Feature screening - Evaluating how much information each feature contains about the outcome
#' 2. Feature construction - If the data doesn't comprise good features, try to construct better ones
#' 3. Feature selection - Apply variable selection methods to determine best combination, here: we'll use wrappers
#' 4. Model evaluation - Check performance of final model on hold-out data
#' 5. Modle exploitation - Using the model [not discussed here]
#' 
#' Within this vignette, we'll use the following toy dataset
summary( toy_data )
vars <- grep( "the_|rnd", colnames(toy_data), value=TRUE )
resp <- "resp"

#' # 1. Feature screening
#' 
#' A first step in predictive modelling is variable screening. This is to check for missing data,
#' outliers and - most importantly - variables that bear information about the response variable.
#'
#' GameRank provides a one-stop function for that: check_variables function.
#' 
vck <- check_variables( toy_data, resp, vars )
vck %>% summary
vck %>% filter( !is_response ) %>% arrange( desc(entropy) )

#' A look into the variable entropy and mutual information with respect to the response is a good idea to
#' identify variables that are constant or bear low information 
#+ fig.width=7
vck %>% 
  ggplot(aes(x=entropy, y=mutual_information) ) +
  geom_point()
vck %>% filter( 0.006 <= mutual_information )


#' It is also important to understand if features are multi-modal. GameRank provides a function for this task: check_multimodality.
#' 
#' GameRank determines multi-modality as follows: Let k be the number of mixture components ranging from 1 to k. The algorithm
#' fits first m_fits GMM models with k components using flexmix. Only models that converge are retained. For each k the minimum AIC
#' is determined together with the number of converging models. All k with less than min_fits_converged models are removed. The
#' k for which the minimum AIC is attained is then chosen. In case of ties for this AIC the minimum number of components are chosen.
#' The first model obtaining these k and AIC is then used to determine cut-points if it has more than one component.
#' 
#' Cut-points are determined by the standard root finding algorithm, determining the points where the prior x mixture component distribution
#' are equal.
#' 
#' It is obvious that the chance for detecting multi-modal distributions depends on the available data and thus distributions reported
#' are not multi-modal or multi-modal distributions may go undected. Thus a visual review of some variable distributions is certainly
#' advisable.
#' 
toy_data %>%
  ggplot( aes( x=the_multi, y=..density.. ) ) +
  geom_histogram( bins = 100, alpha=0.5 ) +
  geom_density( bw = "ucv" ) 

mumo <- check_multimodality( dat = toy_data, resp = resp, vars = vars[1:9],n_comp = 3, m_fits = 25,  min_fits_converged = 20 )

mumo_vars <- mumo$transforms %>% keep( function(x) !is.null( pluck( x, "transformed_var" ) ) ) %>% map( "transformed_var" ) %>% as.character
mumo$data[,mumo_vars]

mumo$transforms$the_multi$aic_aggregate
mumo$transforms$the_multi$best_model
flexmix::parameters(mumo$transforms$the_multi$best_model)
flexmix::prior(mumo$transforms$the_multi$best_model)
pcuts <- mumo$transforms$the_multi$cut_points
mumo$transforms$the_multi$transformed_var

#' Lets create categorical variable and add it to the set.
toy_data <- toy_data %>% bind_cols( mumo$data[,mumo$transforms$the_multi$transformed_var] )
vars <- c(vars, mumo$transforms$the_multi$transformed_var )

#' # 2. Feature construction
#' Another useful task is to see if standard transforms improve Normality of the features. The one-stop facility in GameRank
#' tries square root, cube root, log and z-score transformations. Those that increase the Shapiro-Wilk W-statistics are considered
#' useful and added to the dataset.
#' 
smp <- simple_transforms( toy_data, vars = vars )
tfs <- smp$transformations %>% Reduce( bind_rows, ., NULL )
tfs %>% group_by( variable ) %>% filter( max(W)==W )
tfs %>% pull( transform ) %>% table 

#' Lets add some transformed variables.
svars <- tfs %>% group_by( variable ) %>% filter( max(W)==W ) %>% filter( "identity"!=transform )
svars
toy_data <- toy_data %>% bind_cols( smp$data[,svars$transformed_var] )
vars <- c(vars, svars$transformed_var )

#' Another feature construction approach that can be tried a second round is searching for Power-Transformations via the Box-Cox transformation.
#' However we will skip this here. Please take a look at the example code for box_cox_normal and box_cox_binomial.
#' 

#' 
#' # 3. Feature selection
#' Now, let's run two feature selection algorithms, the bidirectional search that applies forward and backward selection and the GameRank algorithm.
#' First, we'll split the dataset into thirds: one for training the model, one for validating it and one final hold-out dataset.
#' 
rr <- rep_len( c(1L,2L,3L), length.out = nrow(toy_data) ) 
rr <- rr[ order( runif( length(rr) )  )]
df_test <- toy_data[which(3==rr),]
df_sel  <- toy_data[which(rr %in% c(1L,2L)),]
ds <- prepare_splits( ds = 1L, dat = df_sel, resp = resp, vars = vars, fn_train = fn_train_binomial, fn_eval = fn_eval_binomial )

bds <- bidirectional( dat = df_sel, resp = resp, vars = vars, fn_train = fn_train_binomial, fn_eval = fn_eval_binomial_auroc, m = 6L, ds = ds, maximize = TRUE )
bds$variable_selections
bds$agg_results %>% arrange( desc(mean_validation) )
bds$agg_results %>% arrange( desc(mean_validation) ) %>% filter( opt )

#' Now let's try GameRank. A note, since GameRank doesn't use a validation set, the dsi parameter receives and index vector of 1s and 2s that is then
#' repeated to the length of the dataset and thus denotes the relative proportions of training to validation split per round. In small sample scenarios
#' it contains just 2s and the fn_eval function performs a bootstrap.
#' 
gmr <- game_rank( dat = df_sel, resp = resp, vars = vars, fn_train = fn_train_binomial, fn_eval = fn_eval_binomial_auroc, m = 6L, dsi = c(1L,2L), maximize = TRUE, 
                  team_size = 3L, rounds = 10L, min_matches_per_var = 5L )
gmr$variable_ranking %>% as.data.frame
gmr$game_rank_selection

#' 
#' # 4. Model evaluation
#' Having obtained a proposed variable selection, a model needs to be run and evaluated for calibration, ie if it's predictions correlate with
#' the observed outcome. This is easy for regression problems, for probability predictions or survival predictions it involves estimating the
#' observed distribution.
#' 
#+ fig.width=7
bds_fsel <- bds %>% purrr::pluck( "variable_selections" ) %>% purrr::pluck( 1L)
mod_bds <- fn_train_binomial( dat = df_sel, resp = resp, selection = bds_fsel  )
mod_bds
gplot_predictions_binomial( dat = df_sel, resp = resp, selection = bds_fsel, mod = mod_bds )

#+ fig.width=7
gmr_fsel <- gmr %>% purrr::pluck( "game_rank_selection" ) 
mod_gmr <- fn_train_binomial( dat = df_sel, resp = resp, selection = gmr_fsel  )
mod_gmr
gplot_predictions_binomial( dat = df_sel, resp = resp, selection = bds_fsel, mod = mod_gmr )

