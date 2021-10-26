#
# Bidirectional Search wrapper algorithm
#


next_bds_forward <- function( dat, resp, vars, fn_train, fn_eval, ds, maximize, Yf, Yb, ...  ) {
  # browser()
  remaining_vars <- setdiff( Yb, Yf )
  cat( "Evaluating:")
  l_evl <- lapply( remaining_vars, function( nvar, ... ) {
    cat( " ", nvar, " " )
    ret <- eval_splits( ds, dat, resp, union( Yf, nvar ), fn_train, fn_eval, ... )
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

next_bds_backward <- function( dat, resp, vars, fn_train, fn_eval, ds, maximize, Yf, Yb, ...  ) {
  # browser()
  remaining_vars <- setdiff( Yb, Yf )
  cat( "Evaluating:")
  l_evl <- lapply( remaining_vars, function( nvar, ... ) {
    cat( " ", nvar, " " )
    ret <- eval_splits( ds, dat, resp, setdiff( Yb, nvar ), fn_train, fn_eval, ... )
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
  while( m <= length( setdiff( Yb, Yf ) ) ) {
    
    # Forward step
    best_vars <- next_bds_forward( dat, resp, vars, fn_train, fn_eval, ds, maximize, Yf, Yb, ... )
    df_evl <- bind_rows( df_evl, best_vars[['df_evl']] %>% mutate( step = k, type = "forward" ) )
    agg_evl <- bind_rows( agg_evl, best_vars[['agg_evl']] %>% mutate( step = k, type = "forward" ) )
    
    cat( "Adding vars: ")
    cat( best_vars[['best_vars']] )
    Yf <- union( Yf, best_vars[['best_vars']] )
    cat("\n")
    
    # Backward step
    best_vars <- next_bds_backward( dat, resp, vars, fn_train, fn_eval, ds, maximize, Yf, Yb, ... )
    df_evl <- bind_rows( df_evl, best_vars[['df_evl']] %>% mutate( step = k, type = "backward" ) )
    agg_evl <- bind_rows( agg_evl, best_vars[['agg_evl']] %>% mutate( step = k, type = "backward" ) )
    
    cat( "Removing vars: ")
    cat( best_vars[['best_vars']] )
    Yb <- setdiff( Yb, best_vars[['best_vars']] )
    cat("\n")
    
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
    maximize = maximize,
    
    # Results
    variable_selections = best_selections,
    selection_results = best_results,
    results = df_evl,
    agg_results = agg_evl
  )
  return( ret )
}

bidirectional.formula <- function( fo, dat, 
                                   fn_train = fn_train_binomial,
                                   fn_eval  = fn_eval_binomial,
                                   m = NULL,
                                   ds = 5L, 
                                   maximize = TRUE,
                                   ... )
{
  # Check inputs
  stopifnot( is.formula( fo ) ) 
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