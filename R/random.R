#
#  Random sampling wrapper algorithm
#

#' @import tibble dplyr 
#' @importFrom formula.tools lhs.vars rhs.vars

#' @title Random search algorithm
#' 
#' @description First constructs a random sampling plan of disjoint 
#' feature selections that are covering all provided features for a predefined 
#' number of evaluations. Then run training:validation evaluation for each 
#' selection, and returns the best combination.
#' 
#' @details The random search algorithm runs as follows:
#' \preformatted{ 
#' 1. Build a sampling matrix of disjoint and unique feature combinations 
#' 2. For each combination evaluate the model per training:validation splits 
#' 3. Return the best selection
#' }
#' 
#' @param fo Only for call with formula as first argument. Extracts lhs ~ rhs 
#' into resp and vars, and calls backward( dat, resp, vars, ... )
#' @param dat data.frame or tibble comprising data for model generation 
#' and validation.
#' @param resp Character string defining the response (lhs) of the 
#' model formula.
#' @param vars Character vector defining the list of variables for selection. 
#' Those are concatenated by '+' 
#' as the right hand side (rhs) of the modeling formula.
#' @param fn_train Function with signature function( dat, resp, selection, ... )
#' that returns a model or NULL in any other case on the given data dat.
#' @param fn_eval Function with signature function( dat, resp, selection, ... ) 
#' that returns a real number or NA in any other case, e.g. when model is NULL.
#' @param m Size of final partition size. 
#' Note this parameter is used for stopping only. Best selection will be 
#' determined by whole set of evaluated selections, i.e., can be larger than m.
#' @param ds Definition of (parallel) training:validation splits
#'  - a matrix with d columns containing 1s and 2s, where 1 denotes sample is 
#'    used for training the model and 2 denotes sample used for validation.
#'    The average of all d training:validation results is used for selection.
#'  - an integer number determining the number of random training:validation 
#'    splits that should be generated. The sampling will ensure a 
#'    sufficient number of complete cases in the training split.
#' @param maximize A logic value determining if fn_eval is maximized (set to 
#' TRUE) or minimized (set to FALSE).
#' @param nevals Number of training:validation evaluations.
#' @param ... An other arguments passed to fn_train or fn_eval during calls, 
#' e.g. maybe 'u = 365' for Survival evaluations specifying the landmark day.
#'
#' @return List with elements
#' \describe{
#'  \item{response}{As from input parameters}
#'  \item{variables}{As from input parameters}
#'  \item{m}{As from input parameters}
#'  \item{splits}{As from input parameters}
#'  \item{maximize}{As from input parameters}
#'  \item{nevals}{As from input parameters}
#'  \item{start}{Start time of core algorithm loop}
#'  \item{end}{End time of core algorithm loop}
#'  \item{variable_selections}{Best selections overall (regardless of m)}
#'  \item{results}{Dataset with one record per train:validation evaluation}
#'  \item{agg_results}{Dataset with averaged performance over splits}
#' }
#' 
#' @examples 
#' vars <- grep( "the_|rnd", colnames(toy_data), value=TRUE )
#' resp <- "resp"
#' res <- random_selection( toy_data, resp, vars, 
#'                          fn_train_binomial, fn_eval_binomial_auroc, 
#'                          4L, 1L, TRUE )
#' res$variable_selections
#' res$agg_results %>% filter( opt ) %>% arrange( desc(mean_validation) )
#' 
#' @name random
NULL



build_sample_matrix <- function( vars, m, nevals ) {
  
  sm <- matrix( NA_character_, nrow=0, ncol=m )
  while( nrow(sm) < nevals ) {
    idx <- seq_len( length(vars) )
    idx <- idx[ order( stats::runif( length(idx) )) ]
    while( m < length(idx) ) {
      r <- sort( vars[ idx[ seq_len(m) ] ] )
      idx <- idx[-seq_len(m)]
      sm <- rbind( sm, r )
    }
    sm <- unique( sm )
  }
  sm <- sm[seq_len(nevals),]
  rownames(sm) <- NULL
  colnames(sm) <- NULL
  
  return( sm ) 
}




#' @rdname random
#' @export
random_selection <- function( dat,
                              resp,
                              vars,
                              fn_train = fn_train_binomial,
                              fn_eval  = fn_eval_binomial,
                              m = NULL,
                              ds = 5L, 
                              maximize = TRUE,
                              nevals = 100L,
                              ... )
{
  # Check inputs ----
  stopifnot( is.data.frame(dat) | is_tibble(dat) )
  stopifnot( is.character(vars) & 1 < length(vars) )
  stopifnot( is.character(vars) )
  
  stopifnot( is.function(fn_train) ) 
  stopifnot( is.function(fn_eval) ) 
  
  if( is.null(m) ) { stop( "Please provide number m of features to select.\n" ) }
  
  stopifnot( is.integer(nevals) )
  stopifnot( is.logical(maximize) )
  
  message( sprintf("random_selection: Starting random search algorithm for %d observations with %d features (m=%d, maximize=%d).", nrow(dat), length(vars), m, maximize ))
  # Obtain evaluation splits
  ds <- prepare_splits( ds, dat, resp, vars, fn_train, fn_eval, ... )
  
  start_time <- Sys.time()
  # Prepare sampling plan with the following properties:
  #  - Each variable is considered at least min_sample_per_var times
  #  - Variables are evaluated in groups of m sets
  # Let a be # variables, b be group size, c be min # each var evaluated then
  samps <- build_sample_matrix( vars = vars, m = m, nevals = nevals )
  
  # Evaluate random combinations
  df_evl <- NULL
  for( r in seq_len( nrow(samps) ) ) {
    message( sprintf( "Evaluating selection %d.", r )  )
    sel <- as.character( samps[r, ] )
    evl <- mutate( eval_splits( ds, dat, resp, sel, fn_train, fn_eval, ... ), row = r )
    df_evl <- bind_rows( df_evl, evl )
  }
  end_time <- Sys.time()
  
  # Determine best selection(s)
  df_agg <- agg_evals( df_evl, "row", maximize )
  best_selections <- best_selection( df_agg )
  
  ret <- list( 
    algorithm = "random",
    # Input data
    # data = dat,
    response = resp,
    variables = vars,
    m = m,
    
    # Input parameters
    splits = ds, 
    nevals = nevals,
    maximize = maximize,
    
    # Time
    start = start_time, end = end_time,
    
    # Results
    variable_selections = best_selections,
    results = df_evl,
    agg_results = df_agg
  )
  return( ret )
}
  
#' @rdname random
#' @export
random_selection.formula <- function( fo,
                              dat,
                              fn_train = fn_train_binomial,
                              fn_eval  = fn_eval_binomial,
                              m = 5L,
                              ds = 5L, 
                              maximize = TRUE,
                              nevals = 100L,
                              ... )
{
  # Check inputs
  stopifnot( "formula"==class( fo ) ) 
  stopifnot( is.data.frame(dat) | is_tibble(dat) )
  stopifnot( is.function(fn_train) )
  stopifnot( is.function(fn_eval) )
  
  resp <- as.character(formula.tools::lhs.vars( fo ) )
  vars <- as.character(formula.tools::rhs.vars( fo ) )
  
  random_selection( dat = dat,
            resp = resp, 
            vars = vars,
            fn_train = fn_train,
            fn_eval  = fn_eval,
            m = m,
            ds = ds, 
            maximize = maximize,
            nevals = nevals,
            ... )
}