#
# Backward Selection wrapper algorithm
#

#
# Return the best variable to add; in case of ties can be multiple, in case of ties; if none is found
# return the first of 'vars' in ret[['best_vars']]; also return evaluations and aggregate evaluations.
#
next_backward <- function( dat, resp, vars, fn_train, fn_eval, ds, maximize, current_selection, ...  ) {
  # browser()
  remaining_vars <- current_selection
  cat( "Evaluating:")
  l_evl <- lapply( remaining_vars, function( nvar, ... ) {
    cat( " ", nvar, " " )
    ret <- eval_splits( ds, dat, resp, setdiff( current_selection, nvar ), fn_train, fn_eval, ... )
    ret$removed <- nvar
    return( ret )
  }, ... )
  df_evl <- Reduce( bind_rows, l_evl, NULL )
  cat("\n")
  
  agg_evl <- df_evl %>% 
    group_by( ch_selection, removed, selection ) %>%
    summarise( mean_train = mean( eval_train, na.rm=TRUE ),
               mean_validation = mean( eval_validation, na.rm=TRUE ),
               mean_bias = mean( bias, na.rm=TRUE ) ) %>%
    ungroup
  
  ret <- list( df_evl = df_evl, agg_evl = agg_evl )
  if( all( is.na(agg_evl$mean_validation) ) ) {
    ret[['best_vars']] <- c( remaining_vars[length(remaining_vars)] )
    cat( sprintf( "No best found, removing last (%s) \n", ret[['best_vars']] ) )
  } else if( maximize ) {
    ret[['best_vars']] <- agg_evl %>% 
      filter( !is.na(mean_validation) ) %>% 
      filter( max( mean_validation, na.rm=TRUE ) <= mean_validation ) %>%
      pull( removed )
    cat( sprintf( "Best found, removal variable (%s) \n", ret[['best_vars']] ) )
  } else if( !maximize ) {
    ret[['best_vars']] <- agg_evl %>% 
      filter( !is.na(mean_validation) ) %>% 
      filter( mean_validation <= min( mean_validation, na.rm=TRUE ) ) %>%
      pull( removed )
    cat( sprintf( "Best found, removal variable (%s) \n", ret[['best_vars']] ) )
  }
  return( ret )
}

#
# Conduct one backward wrapper algorithm
#
backward <- function( dat, resp, vars, 
                      fn_train = fn_train_binomial,
                      fn_eval  = fn_eval_binomial,
                      m = m,
                      ds = 5L, 
                      maximize = TRUE, 
                      ... ) {
  # Check parameters
  stopifnot( is.data.frame(dat) | is_tibble(dat) )
  stopifnot( is.character(resp) )
  stopifnot( is.character(vars)  & 1 < length(vars) )
  stopifnot( is.logical(maximize) )
  
  if( is.null(m) ) { stop( "Please provide number m of features to select.\n" ) }
  
  # Obtain evaluation splits
  ds <- prepare_splits( ds, dat, resp, vars, fn_train, fn_eval, ... )
  
  # Sequential Backward Selection algorithm:
  # 1) Start with the full set Y0 = X
  # 2) Remove the worst feature x- = arg max_{x in Yk} J( Yk - x)
  # 3) Update Yk+1 = Yk - x-; k = k + 1
  # 4) Go to 2)
  
  df_evl <- NULL
  agg_evl <- NULL
  Y <- vars  # Start with all variables included
  queue <- list()
  queue[[length(queue)+1]] <- Y
  k <- 1
  while( 0 < length(queue) ) {
    # Pull first element from search queue
    Y <- queue[[1]]
    queue[[1]] <- NULL
    
    if( m <= length(Y) ) {
      best_vars <- next_backward( dat, resp, vars, fn_train, fn_eval, ds, maximize, Y, ... )
      df_evl <- bind_rows( df_evl, best_vars[['df_evl']] %>% mutate( step = k ) )
      agg_evl <- bind_rows( agg_evl, best_vars[['agg_evl']] %>% mutate( step = k ) )
      
      cat( "Removing vars: ")
      cat( best_vars[['best_vars']] )
      for( bvar in best_vars[['best_vars']] ) {
        queue[[length(queue)+1]] <- setdiff( Y, bvar )
      }
      cat("\n")
      cat( sprintf("No of partitions %d in search queue \n", length( queue ) ) )
    }
    cat( sprintf( "Finished iteration %d \n", k ) )
    k <- k + 1
  } # while
  
  # Determine best selection(s)
  best_results <- NULL
  if( maximize ) {
    best_results <- agg_evl %>% 
      filter( !is.na(mean_validation) ) %>% 
      filter( max( mean_validation, na.rm=TRUE ) <= mean_validation )
  } else if( !maximize ) {
    best_results <- agg_evl %>% 
      filter( !is.na(mean_validation) ) %>% 
      filter( mean_validation <= min( mean_validation, na.rm=TRUE ) )
  }
  best_selections <- best_results %>% pull( selection )
  
  ret <- list( # Input data
    # data = dat,
    response = resp,
    variables = vars,
    m = m,
    
    # Input parameters
    splits = ds, 
    min_partition = min_partition,
    maximize = maximize,
    
    # Results
    variable_selections = best_selections,
    selection_results = best_results,
    results = df_evl,
    agg_results = agg_evl
  )
  return( ret )
} # backward (END)

backward.formula <- function( fo, dat, 
                              fn_train = fn_train_binomial,
                              fn_eval  = fn_eval_binomial,
                              m = NULL,
                              ds = 5L, 
                              maximize = TRUE,
                              min_partition = NULL,
                              ...  ) 
{
  # Check inputs
  stopifnot( is.formula( fo ) ) 
  stopifnot( is.data.frame(dat) | is_tibble(dat) )
  stopifnot( is.function(fn_train) )
  stopifnot( is.function(fn_eval) )
  
  resp <- as.character(lhs.vars( fo ) )
  vars <- as.character(rhs.vars( fo ) )
  
  backward( dat = dat,
            resp = resp, 
            vars = vars,
            fn_train = fn_train,
            fn_eval  = fn_eval,
            m = m,
            ds = ds, 
            maximize = maximize,
            ... )
}