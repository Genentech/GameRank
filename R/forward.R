#
# Forward Selection wrapper algorithm
#

#
# Return the best variable to add; in case of ties can be multiple, in case of ties; if none is found
# return the first of 'vars' in ret[['best_vars']]; also return evaluations and aggregate evaluations.
#

# next_forward <- function( dat, resp, vars, fn_train, fn_eval, ds, maximize, current_selection, ...  ) {
#   # browser()
#   remaining_vars <- setdiff( vars, current_selection )
#   cat( "Evaluating:")
#   l_evl <- lapply( remaining_vars, function( nvar, ... ) {
#     cat( " ", nvar, " " )
#     ret <- eval_splits( ds, dat, resp, union( current_selection, nvar ), fn_train, fn_eval, ... )
#     ret$added <- nvar
#     return( ret )
#   }, ... )
#   df_evl <- Reduce( bind_rows, l_evl, NULL )
#   cat("\n")
#   
#   agg_evl <- df_evl %>% 
#     group_by( ch_selection, added, selection ) %>%
#     summarise( mean_train = mean( eval_train, na.rm=TRUE ),
#                mean_validation = mean( eval_validation, na.rm=TRUE ),
#                mean_bias = mean( bias, na.rm=TRUE ) ) %>%
#     ungroup
#   
#   ret <- list( df_evl = df_evl, agg_evl = agg_evl )
#   if( all( is.na(agg_evl$mean_validation) ) ) {
#     ret[['best_vars']] <- c( remaining_vars[1] )
#     cat( sprintf( "No best found, adding first (%s) \n", ret[['best_vars']] ) )
#   } else if( maximize ) {
#     ret[['best_vars']] <- agg_evl %>% 
#       filter( !is.na(mean_validation) ) %>% 
#       filter( max( mean_validation, na.rm=TRUE ) <= mean_validation ) %>%
#       pull( added )
#     cat( sprintf( "Best found, adding variable (%s) \n", ret[['best_vars']] ) )
#   } else if( !maximize ) {
#     ret[['best_vars']] <- agg_evl %>% 
#       filter( !is.na(mean_validation) ) %>% 
#       filter( mean_validation <= min( mean_validation, na.rm=TRUE ) ) %>%
#       pull( added )
#     cat( sprintf( "Best found, adding variable (%s) \n", ret[['best_vars']] ) )
#   }
#   return( ret )
# }

#
# Conduct one forward wrapper algorithm
#
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
      df_evl <- bind_rows( df_evl, best_vars[['df_evl']] %>% mutate( k = k ) )
      agg_evl <- bind_rows( agg_evl, best_vars[['agg_evl']] %>% mutate( k = k ) )
      
      queue <- append( queue, best_vars[['best_selections']] )
      cat( sprintf("No of partitions %d in search queue \n", length( queue ) ) )
    }
    cat( sprintf( "Finished iteration %d \n", k ) )
    k <- k + 1
  } # while

  # Determine best selection(s)
  agg <- agg_evals( df_evl, NULL, maximize )
  best_selections <- best_selection( agg )
  
  ret <- list( # Input data
    # data = dat,
    response = resp,
    variables = vars,
    m = m,
    
    # Input parameters
    splits = ds, 
    maximize = maximize,
    
    # Results
    variable_selections = best_selections,
    results = df_evl,
    agg_results = agg_evl
  )
  return( ret )
} # forward (END)

forward.formula <- function( fo, dat, 
                             fn_train = fn_train_binomial,
                             fn_eval  = fn_eval_binomial,
                             m = NULL,
                             ds = 5L, 
                             maximize = TRUE,
                             ...  ) 
{
  # Check inputs
  stopifnot( is.formula( fo ) ) 
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
