#
# Plus-L, Minus-R Search
#

#' @import tibble dplyr 
#' @importFrom formula.tools lhs.vars rhs.vars

#' @title Plus-L, Minus-R search algorithm
#' 
#' @description Performs L forward and R backward selection steps per iteration
#' until a desired partition size is reached. Depending on L<R or L>R it starts 
#' with either the full set of variables or the empty set of variables.
#' 
#' @details The Plus-L, Minus-R algorithm runs as follows:
#' \preformatted{ 
#' 1. if L > R then  
#'      Y0 = {} 
#'    else  
#'      Y0 = X; go to step 3 
#'    fi 
#' 2. repeat L times 
#'      x+ = arg max_{x not in Yk} J( Yk + x ) 
#'      Yk+1 = Yk + x+; k = k + 1 
#' 3. repeat R times 
#'      x- = arg max_{x in Yk} J( Yk - x ) 
#'      Yk+1 = Yk - x-; k = k + 1 
#' 4. go to step 2 
#' }
#' 
#' @param fo Only for call with formula as first argument. Extracts lhs ~ rhs 
#' into resp and vars, and calls backward( dat, resp, vars, ... )
#' @param dat data.frame or tibble comprising data for model generation and 
#' validation.
#' @param resp Character string defining the response (lhs) of the 
#' model formula.
#' @param vars Character vector defining the list of variables for selection. 
#' Those are concatenated by '+' as the right hand side (rhs) of the 
#' modeling formula.
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
#'    splits that should be generated. The sampling will ensure a sufficient 
#'    number of complete cases in the training split.
#' @param maximize A logic value determining if fn_eval is maximized (set to 
#' TRUE) or minimized (set to FALSE).
#' @param L Number of forward steps per iteration.
#' @param R Number of backward steps per iteration.
#' @param kmax Limit number of iterations
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
#'  \item{L}{As from input parameters}
#'  \item{R}{As from input parameters}
#'  \item{kmax}{As from input parameters}
#'  \item{start}{Start time of core algorithm loop}
#'  \item{end}{End time of core algorithm loop}
#'  \item{variable_selections}{Best selections overall (regardless of m)}
#'  \item{results}{Dataset with one record per train:validation evaluation}
#'  \item{agg_results}{Dataset with averaged performance over splits}
#' }
#' @name lrsearch
NULL

#' @rdname lrsearch
#' @examples 
#' vars <- grep( "the_|rnd", colnames(toy_data), value=TRUE )
#' resp <- "resp"
#' res <- lrsearch( toy_data, resp, vars, 
#'                  fn_train_binomial, fn_eval_binomial_auroc, 
#'                  4L, 1L, TRUE, 3L, 5L )
#' res$variable_selections
#' res$agg_results %>% filter( opt ) %>% arrange( desc(mean_validation) )
#' @export
lrsearch <- function(  dat, resp, vars, 
                       fn_train = fn_train_binomial,
                       fn_eval  = fn_eval_binomial,
                       m = NULL,
                       ds = 5L, 
                       maximize = TRUE,
                       L, R, 
                       kmax = 1000,
                       ... )
{
  # Check parameters
  stopifnot( is.data.frame(dat) | is_tibble(dat) )
  stopifnot( is.character(resp) )
  stopifnot( is.character(vars)  & 1 < length(vars) )
  stopifnot( is.logical(maximize) )
  stopifnot( is.integer(R) )
  stopifnot( is.integer(L) )
  stopifnot( R!=L )
  stopifnot( 0 < kmax )
  
  if( is.null(m) ) { stop( "Please provide number m of features to select.\n" ) }
  
  # Obtain evaluation splits, if necessary
  ds <- prepare_splits( ds, dat, resp, vars, fn_train, fn_eval, ... )
  
  message( sprintf("lrsearch: Starting Plus-L,Minus-R search algorithm for %d observations with %d features (m=%d, maximize=%d).", nrow(dat), length(vars), m, maximize ))
  start_time <- Sys.time()
  # Plus-L, Minus-R algorithm:
  # 1. if L > R then 
  #      Y0 = {}
  #    else 
  #      Y0 = X; go to step 3
  #    fi
  # 2. repeat L times
  #      x+ = arg max_{x not in Yk} J( Yk + x )
  #      Yk+1 = Yk + x+; k = k + 1
  # 3. repeat R times
  #      x- = arg max_{x in Yk} J( Yk - x )
  #      Yk+1 = Yk - x-; k = k + 1
  # 4. go to step 2
  # 
  df_evl <- NULL
  agg_evl <- NULL
  Y <- NULL
  d <- 0 # direction counter, is always counted towards 0. To avoid 'go to's from pseudocode.
  if( L > R ) {
    Y <- c("1")
    d <- -L
  } else {
    Y <- vars
    d <- +R
  }
  
  k <- 0
  while( (k <= kmax) & (m!=length(Y)) & (d != 0) ) {
    # L forward steps
    while( (k <= kmax) & (m!=length(Y)) & (d < 0) ) {
      best_vars <- eval_add_vars( ds, dat, resp, vars, fn_train, fn_eval, maximize, Y, setdiff( vars, Y ), ...  )
      if(!is.null(best_vars))  {
        df_evl <- bind_rows( df_evl, best_vars[['df_evl']] %>% mutate( k = k ) )
        agg_evl <- bind_rows( agg_evl, best_vars[['agg_evl']] %>% mutate( k = k ) )
        bs <- purrr::pluck( best_vars, 'best_selections' )
        Y <- purrr::pluck( bs, 1 )
        d <- d + 1
      }
    }
    d <- +R
    k <- k + 1
    
    # R backward steps
    while( (k <= kmax) & (m!=length(Y)) & (0 < d) ) {
      best_vars <- eval_remove_vars( ds, dat, resp, vars, fn_train, fn_eval, maximize, Y, Y, ...  )
      if(!is.null(best_vars)) {
        df_evl <- bind_rows( df_evl, best_vars[['df_evl']] %>% mutate( k = k ) )
        agg_evl <- bind_rows( agg_evl, best_vars[['agg_evl']] %>% mutate( k = k ) )
        bs <- purrr::pluck( best_vars, 'best_selections' )
        Y <- purrr::pluck( bs, 1 )
        d <- d - 1
      }
    }  
    d <- -L
    k <- k + 1
    
    # iterate to next +R/-L round
  }
  end_time <- Sys.time()
  
  # Determine best selection(s)
  agg <- agg_evals( df_evl, NULL, maximize )
  best_selections <- best_selection( agg )
  
  ret <- list( 
    algorithm = "lrsearch",
    # Input data
    # data = dat,
    response = resp,
    variables = vars,
    m = m,
    
    # Input parameters
    splits = ds, 
    L = L, R = R,
    kmax = kmax,
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

#' @rdname lrsearch
#' @export
lrsearch.formula <- function( fo, dat, 
                              fn_train = fn_train_binomial,
                              fn_eval  = fn_eval_binomial,
                              m = NULL,
                              ds = 5L, 
                              maximize = TRUE,
                              L, R,
                              kmax = 1000,
                              ...  ) 
{
  # Check inputs
  stopifnot( "formula"==class( fo ) ) 
  stopifnot( is.data.frame(dat) | is_tibble(dat) )
  stopifnot( is.function(fn_train) )
  stopifnot( is.function(fn_eval) )
  
  resp <- as.character(formula.tools::lhs.vars( fo ) )
  vars <- as.character(formula.tools::rhs.vars( fo ) )
  
  lrsearch(  dat = dat,
            resp = resp, 
            vars = vars,
            fn_train = fn_train,
            fn_eval  = fn_eval,
            m = m,
            ds = ds, 
            maximize = maximize,
            L = L, R = R, 
            kmax = kmax,
            ... )
}
