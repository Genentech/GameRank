#
# Bidirectional Search wrapper algorithm
#

#' @import tibble dplyr formula.tools

#' @title Bidirectional search algorithm
#' 
#' @description Performs forward and backward selection steps per iteration to converge to
#' the same, consistent variable selection set.
#' 
#' @details The Bidirectional Search algorithm runs as follows:
#' \code{ \cr
#' 1) Start SFS with Yf = {} \cr
#' 2) Start SBS with Yb = X \cr
#' 3) Select the best feature \cr
#'         x+ = arg max_{x not in Yfk; x in Ybk} J( Yfk + x ) \cr
#'         Yfk+1 = Yfk + x+ \cr
#' 4) Remove the worst feature \cr
#'         x- = arg max_{x in Ybk; x not in Yfk+1} J( Ybk - x ) \cr
#'         Ybk+1 = Ybk - x-; \cr
#'         k = k + 1 \cr
#' 5) Go to 3 \cr
#' }
#' 
#' @param fo Only for call with formula as first argument. Extracts lhs ~ rhs into resp and vars, and calls backward( dat, resp, vars, ... )
#' @param dat data.frame or tibble comprising data for model generation and validation.
#' @param resp Character string defining the response (lhs) of the model formula.
#' @param vars Character vector defining the list of variables for selection. Those are concatenated by '+' 
#' as the right hand side (rhs) of the modelling formula.
#' @param fn_train Function with signature function( dat, resp, selection, ... ) that returns a model or NULL in any other case on the given data dat.
#' @param fn_eval Function with signature function( dat, resp, selection, ... ) that returns a real number or NA in any other case, e.g. when model is NULL.
#' @param m Size of final partition size. 
#' Note this parameter is used for stopping only. Best selection will be determined by whole set of evaluated selections, i.e., can be larger than m.
#' @param ds Definition of (parallel) training:validation splits
#'  - a matrix with d columns containing 1s and 2s, where 1 denotes sample is used for training the model and 2 denotes sample used for validation.
#'    The average of all d training:validation results is used for selection.
#'  - an integer number determing the number of random training:validation splits that should be generated. The sampling will ensure a sufficient number
#'    of complete cases in the training split.
#' @param maximize A logic value determining if fn_eval is maximized (set to TRUE) or minimized (set to FALSE).
#' @param ... An other arguments passed to fn_train or fn_eval during calls, e.g. maybe 'u = 365' for Survival evaluations specifying the landmark day.
#'
#' @return List with elements
#' \describe{
#'  \item{response}{As from input parameters}
#'  \item{variables}{As from input parameters}
#'  \item{m}{As from input parameters}
#'  \item{splits}{As from input parameters}
#'  \item{maximize}{As from input parameters}
#'  \item{start}{Start time of core algorithm loop}
#'  \item{end}{End time of core algorithm loop}
#'  \item{variable_selections}{Best selections overall (regardless of m)}
#'  \item{results}{Dataset with one record per train:validation evaluation}
#'  \item{agg_results}{Dataset with averaged performance over splits}
#' }
#' @name bidirectional
NULL

#' @rdname bidirectional
#' @examples 
#' vars <- grep( "the_|rnd", colnames(toy_data), value=TRUE )
#' resp <- "resp"
#' res <- bidirectional( toy_data, resp, vars, fn_train_binomial, fn_eval_binomial_auroc, 4L, 1L, TRUE )
#' res$variable_selections
#' res$agg_results %>% filter( opt ) %>% arrange( desc(mean_validation) )
#' @export
bidirectional <- function( dat, resp, vars, 
                           fn_train = fn_train_binomial,
                           fn_eval  = fn_eval_binomial,
                           m = NULL,
                           ds = 5L, 
                           maximize = TRUE,
                           ... )
{
  # Check parameters
  stopifnot( is.data.frame(dat) | is_tibble(dat) )
  stopifnot( is.character(resp) )
  stopifnot( is.character(vars)  & 1 < length(vars) )
  stopifnot( is.logical(maximize) )
  
  if( is.null(m) ) { stop( "Please provide number m of features to select.\n" ) }
  
  # Obtain evaluation splits
  ds <- prepare_splits( ds, dat, resp, vars, fn_train, fn_eval, ... )
  
  start_time <- Sys.time()
  # Bidirectional Search algorithm:
  # 1) Start SFS with Yf = {}
  # 2) Start SBS with Yb = X
  # 3) Select the best feature
  #         x+ = arg max_{x not in Yfk; x in Ybk} J( Yfk + x )
  #         Yfk+1 = Yfk + x+
  # 4) Remove the worst feature
  #         x- = arg max_{x in Ybk; x not in Yfk+1} J( Ybk - x )
  #         Ybk+1 = Ybk - x-;
  #    k = k + 1
  # 5) Go to 3
  #
  df_evl <- NULL
  agg_evl <- NULL
  Yf <- c( "1" )
  Yb <- vars # Contains features to select
  k <- 1
  while( !setequal( Yb, intersect(Yb,Yf) ) & (length(intersect( Yf, Yb )) <= m) ) {
    # Forward step
    if( 0 < length(setdiff(Yb,Yf))) {
      best_vars <- eval_add_vars( ds, dat, resp, lst_vars, fn_train, fn_eval, maximize, Yf, setdiff( Yb, Yf ), ...  )
      if(!is.null(best_vars)) {
        df_evl <- bind_rows( df_evl, purrr::pluck( best_vars,'df_evl') %>% mutate( k = k ) )
        agg_evl <- bind_rows( agg_evl, purrr::pluck( best_vars, 'agg_evl' ) %>% mutate( k = k ) )
        bs <- purrr::pluck( best_vars, 'best_selections' )
        if( 0 < length( bs ) ) {
          Y <- purrr::pluck( bs,1L )
          Yf <- union( Y, Yf )  
        }
      }
    }
    # Backward step
    if( 0 < length(setdiff(Yb,Yf))) {
      best_vars <- eval_remove_vars( ds, dat, resp, lst_vars, fn_train, fn_eval, maximize, Yb, setdiff( Yb, Yf ), ...  )
      if(!is.null(best_vars)) {
        df_evl <- bind_rows( df_evl, purrr::pluck( best_vars,'df_evl') %>% mutate( k = k ) )
        agg_evl <- bind_rows( agg_evl, purrr::pluck( best_vars, 'agg_evl' ) %>% mutate( k = k ) )
        bs <- purrr::pluck( best_vars, 'best_selections' )
        if( 0 < length(bs) ) {
          Y <- purrr::pluck( bs,1L )
          Yb <- intersect( Yb, Y ) # Y is Yb removed by one variable
        }
      }
    }
    k <- k + 1
  } # while
  end_time <- Sys.time()
  
  # Determine best selection(s)
  agg <- agg_evals( df_evl, NULL, maximize )
  best_selections <- best_selection( agg )
  
  ret <- list( 
    algorithm = "bidirectional",
    # Input data
    # data = dat,
    response = resp,
    variables = vars,
    m = m,
    
    # Input parameters
    splits = ds, 
    maximize = maximize,
    
    # Time
    start = start_time, end = end_time,
    
    # Results
    variable_selections = best_selections,
    results = df_evl,
    agg_results = agg_evl
  )
  return( ret )
}

#' @rdname bidirectional
#' @export
bidirectional.formula <- function( fo, dat, 
                                   fn_train = fn_train_binomial,
                                   fn_eval  = fn_eval_binomial,
                                   m = NULL,
                                   ds = 5L, 
                                   maximize = TRUE,
                                   ... )
{
  # Check inputs
  stopifnot( "formula"==class( fo ) ) 
  stopifnot( is.data.frame(dat) | is_tibble(dat) )
  stopifnot( is.function(fn_train) )
  stopifnot( is.function(fn_eval) )
  
  resp <- as.character(lhs.vars( fo ) )
  vars <- as.character(rhs.vars( fo ) )
  
  bidirectional( dat = dat,
                 resp = resp, 
                 vars = vars,
                 m = m,
                 fn_train = fn_train,
                 fn_eval  = fn_eval,
                 ds = ds, 
                 maximize = maximize,
                 ... )
}