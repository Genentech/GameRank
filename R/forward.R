#
# Forward Selection wrapper algorithm
#

#' @import tibble dplyr formula.tools

#' @title Sequential Forward Selection algorithm
#' 
#' @description Starts with the empty set of available features and adds the best feature
#' during each iteration until determined partition size is reached. In case of ties, all partitions
#' are explored.
#' 
#' @details The Sequential Forward Selection algorithm runs as follows:
#' \code{ \cr
#'  1) Start with the empty set Y0 = {} \cr
#'  2) Select the next best feature x+ = arg max J( Yk + x ) for x not in Yk \cr
#'  3) Update Yk+1 = Yk + x+; k = k + 1 \cr
#'  4) Go to 2 \cr
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
#' @name forward
NULL

#' @rdname forward
#' @export
forward <- function(  dat, resp, vars, 
                      fn_train = fn_train_binomial,
                      fn_eval  = fn_eval_binomial,
                      m = NULL,
                      ds = 5L, 
                      maximize = TRUE,
                      ... ) {
  # Check parameters
  stopifnot( is.data.frame(dat) | is_tibble(dat) )
  stopifnot( is.character(resp) )
  stopifnot( is.character(vars)  & 1 < length(vars) )
  stopifnot( is.logical(maximize) )
  
  if( is.null(m) ) { stop( "Please provide number m of features to select.\n" ) }
  
  # Obtain evaluation splits, if necessary
  ds <- prepare_splits( ds, dat, resp, vars, fn_train, fn_eval, ... )

  start_time <- Sys.time()
  # Sequential Forward Selection algorithm:
  # 1) Start with the empty set Y0 = {}
  # 2) Select the next best feature x+ = arg max J( Yk + x ) for x not in Yk
  # 3) Update Yk+1 = Yk + x+; k = k + 1
  # 4) Go to 2
  # -> Implemented as depth first search, in case of ties or no variables
  #    In case of ties, all combinations are generated, in case of no best variable
  #    the first one is added
  
  df_evl <- NULL
  agg_evl <- NULL
  Y <- c("1")  # Start with offset, which is anyway included in each model
  queue <- list()
  queue[[length(queue)+1]] <- Y
  k <- 1
  while( 0 < length(queue) ) {
    # Pull first element from search queue
    Y <- queue[[1]]
    queue[[1]] <- NULL
    
    if( length(Y) <= m ) {
      best_vars <- eval_add_vars( ds, dat, resp, lst_vars, fn_train, fn_eval, maximize, Y, setdiff( vars, Y ), ...  )
      if( !is.null(best_vars) ) {
        df_evl <- bind_rows( df_evl, best_vars[['df_evl']] %>% mutate( k = k ) )
        agg_evl <- bind_rows( agg_evl, best_vars[['agg_evl']] %>% mutate( k = k ) )
        
        bs <- best_vars[['best_selections']]
        if( 0==length(bs) ) bs <- list( as.character( union( Y, utils::head( setdiff(vars,Y),1 ) ) ) )
        queue <- append( queue, bs )
        cat( sprintf("No of partitions %d in search queue \n", length( queue ) ) )
      }
    }
    cat( sprintf( "Finished iteration %d \n", k ) )
    k <- k + 1
  } # while
  end_time <- Sys.time()
  # browser()
  
  # Determine best selection(s)
  agg <- agg_evals( df_evl, NULL, maximize )
  best_selections <- best_selection( agg )
  
  ret <- list( 
    algorithm = "forward",
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
} # forward (END)

#' @rdname forward
#' @export
forward.formula <- function( fo, dat, 
                             fn_train = fn_train_binomial,
                             fn_eval  = fn_eval_binomial,
                             m = NULL,
                             ds = 5L, 
                             maximize = TRUE,
                             ...  ) 
{
  # Check inputs
  stopifnot( "formula"==class( fo ) ) 
  stopifnot( is.data.frame(dat) | is_tibble(dat) )
  stopifnot( is.function(fn_train) )
  stopifnot( is.function(fn_eval) )
  
  resp <- as.character(lhs.vars( fo ) )
  vars <- as.character(rhs.vars( fo ) )
  
  forward(  dat = dat,
            resp = resp, 
            vars = vars,
            fn_train = fn_train,
            fn_eval  = fn_eval,
            m = NULL,
            ds = ds, 
            maximize = maximize,
            ... )
}
