#
# Plus-L, Minus-R Search
#


next_lrforward <- function( dat, resp, vars, fn_train, fn_eval, ds, maximize, current_selection, ...  ) {
  # browser()
  remaining_vars <- setdiff( vars, current_selection )
  cat( "Evaluating:")
  l_evl <- lapply( remaining_vars, function( nvar, ... ) {
    cat( " ", nvar, " " )
    ret <- eval_splits( ds, dat, resp, union( current_selection, nvar ), fn_train, fn_eval, ... )
    ret$added <- nvar
    return( ret )
  }, ... )
  df_evl <- Reduce( bind_rows, l_evl, NULL )
  cat("\n")
  
  agg_evl <- df_evl %>% 
    group_by( ch_selection, added, selection ) %>%
    summarise( mean_train = mean( eval_train, na.rm=TRUE ),
               mean_validation = mean( eval_validation, na.rm=TRUE ),
               mean_bias = mean( bias, na.rm=TRUE ) ) %>%
    ungroup
  
  ret <- list( df_evl = df_evl, agg_evl = agg_evl )
  if( all( is.na(agg_evl$mean_validation) ) ) {
    ret[['best_vars']] <- c( remaining_vars[1] )
    cat( sprintf( "No best found, adding first (%s) \n", ret[['best_vars']] ) )
  } else if( maximize ) {
    ret[['best_vars']] <- agg_evl %>% 
      filter( !is.na(mean_validation) ) %>% 
      filter( max( mean_validation, na.rm=TRUE ) <= mean_validation ) %>%
      pull( added )
    cat( sprintf( "Best found, adding variable (%s) \n", ret[['best_vars']] ) )
  } else if( !maximize ) {
    ret[['best_vars']] <- agg_evl %>% 
      filter( !is.na(mean_validation) ) %>% 
      filter( mean_validation <= min( mean_validation, na.rm=TRUE ) ) %>%
      pull( added )
    cat( sprintf( "Best found, adding variable (%s) \n", ret[['best_vars']] ) )
  }
  return( ret )
}

next_lrbackward <- function( dat, resp, vars, fn_train, fn_eval, ds, maximize, current_selection, ...  ) {
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
      best_vars <- next_lrforward( dat, resp, vars, fn_train, fn_eval, ds, maximize, Y, ... )
      df_evl  <- bind_rows(  df_evl, best_vars[['df_evl']] %>% mutate( step = k, d = d ) )
      agg_evl <- bind_rows( agg_evl, best_vars[['agg_evl']] %>% mutate( step = k, d = d ) )
      cat( "Adding first of vars: ")
      cat( best_vars[['best_vars']] )
      bvar <- best_vars[['best_vars']]
      Y <- union( Y, bvar[1] )
      cat("\n")
      d <- d + 1
    }
    d <- +R
    k <- k + 1
    
    # R backward steps
    while( (k <= kmax) & (m!=length(Y)) & (0 < d) ) {
      best_vars <- next_lrbackward( dat, resp, vars, fn_train, fn_eval, ds, maximize, Y, ... )
      df_evl  <- bind_rows(  df_evl, best_vars[['df_evl']] %>% mutate( step = k, d = d ) )
      agg_evl <- bind_rows( agg_evl, best_vars[['agg_evl']] %>% mutate( step = k, d = d ) )
      cat( "Removing first of vars: ")
      cat( best_vars[['best_vars']] )
      bvar <- best_vars[['best_vars']]
      Y <- setdiff( Y, bvar[1] )
      cat("\n")
      d <- d - 1
    }  
    d <- -L
    k <- k + 1
    
    # iterate to next +R/-L round
  }
  
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
    L = L, R = R,
    kmax = kmax,
    
    # Results
    variable_selections = best_selections,
    selection_results = best_results,
    results = df_evl,
    agg_results = agg_evl
  )
  return( ret )
}



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
  stopifnot( is.formula( fo ) ) 
  stopifnot( is.data.frame(dat) | is_tibble(dat) )
  stopifnot( is.function(fn_train) )
  stopifnot( is.function(fn_eval) )
  
  resp <- as.character(lhs.vars( fo ) )
  vars <- as.character(rhs.vars( fo ) )
  
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
