#
# Sequential Floating Forward Selection (SFFS) and
# Sequential Floating Backward Selection (SFBS)
#


# Sequential Floating Forward Search ----
# sffs_forward <- function( ds, dat, resp, vars, fn_train, fn_eval, Y, maximize, ... ) {
#   
#   remaining_vars <- setdiff( vars, Y )
#   df_evl <- map_dfr( .x=remaining_vars, .f=function( ff ) mutate( eval_splits( ds, dat, resp, vars, fn_train, fn_eval, union( Y,ff ), ... ), added = ff )  )
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
# sffs_backward <- function( ds, dat, resp, vars, fn_train, fn_eval, Y, J, maximize, ... ) {
#   
#   df_evl <- map_dfr( .x=Y, .f=function( ff ) mutate( eval_splits( ds, dat, resp, vars, fn_train, fn_eval, setdiff( Y,ff ), ... ), removed = ff )  )
#   
#   agg_evl <- df_evl %>% 
#     group_by( ch_selection, removed, selection ) %>% 
#     summarise( mean_train = mean( eval_train, na.rm=TRUE ),
#                mean_validation = mean( eval_validation, na.rm=TRUE ),
#                mean_bias = mean( bias, na.rm=TRUE ) ) %>%
#     ungroup
#   
#   ret <- list( df_evl = df_evl, agg_evl = agg_evl )
#   if( maximize ) {
#     ret[['best_vars']] <- agg_evl %>% 
#       filter( !is.na(mean_validation) ) %>% 
#       filter( (J < mean_validation) &  max( mean_validation, na.rm=TRUE ) <= mean_validation ) %>%
#       pull( removed )
#     cat( sprintf( "Best found, removing variable (%s) \n", ret[['best_vars']] ) )
#   } else if( !maximize ) {
#     ret[['best_vars']] <- agg_evl %>% 
#       filter( !is.na(mean_validation) ) %>% 
#       filter( (J > mean_validation) & mean_validation <= min( mean_validation, na.rm=TRUE ) ) %>%
#       pull( removed )
#     cat( sprintf( "Best found, removing variable (%s) \n", ret[['best_vars']] ) )
#   }
#   return( ret )
# }
# 

sffs <- function( dat, resp, vars, 
                  fn_train = fn_train_binomial,
                  fn_eval  = fn_eval_binomial,
                  m = NULL,
                  ds = 5L, 
                  maximize = TRUE,
                  kmax = 100,
                  ... ) 
{
  # Check parameters
  stopifnot( is.data.frame(dat) | is_tibble(dat) )
  stopifnot( is.character(resp) )
  stopifnot( is.character(vars)  & 1 < length(vars) )
  stopifnot( is.logical(maximize) )
  
  if( is.null(m) ) { stop( "Please provide number m of features to select.\n" ) }
  
  # Obtain evaluation splits, if necessary
  ds <- prepare_splits( ds, dat, resp, vars, fn_train, fn_eval, ... )
  
  # Ref.: "Floating Search Methods for Feature Selection with Nonmonotonic Criterion Functions"; Pudil, Ferri, Novovicova, Kittler, 1994, IEEE
  # Sequential Floating Forward Selection algorithm:
  # 1. Y = {}
  # 2. Select the best feature
  #      x+ = arg max_{x not in Yk} J( Yk + x )
  #      Yk = Yk + x+; k = k + 1
  # 3. Select the worst feature* 
  #    [Note some book-keeping is necessary to avoid infinite loops]
  #      x- = arg max_{x in Yk} J( Yk - x )
  # 4. if J( Yk - x- ) > J( Yk ) then
  #      Yk+1 = Yk - x-; k = k + 1
  #      go to step 3
  #    else
  #      go to step 2
  # 
  df_evl <- NULL
  agg_evl <- NULL
  Y <- c("1")
  k <- 1
  while( (k <= kmax) & length(Y) <= m ) {
    best_vars <- eval_add_vars( ds, dat, resp, lst_vars, fn_train, fn_eval, maximize, Y, setdiff( vars, Y ), ...  )
    df_evl <- bind_rows( df_evl, best_vars[['df_evl']] %>% mutate( k = k ) )
    agg_evl <- bind_rows( agg_evl, best_vars[['agg_evl']] %>% mutate( k = k ) )
    Y <- best_vars[['agg_evl']] %>% filter( opt ) %>% pull( selection ) %>% pluck(1)
    J <- best_vars[['agg_evl']] %>% filter( opt ) %>% pull( mean_validation ) %>% pluck(1)
    k <- k + 1
    
    browser()
    removing <- TRUE
    while( removing ) {
      removing <- FALSE
      best_vars <- eval_remove_vars( ds, dat, resp, lst_vars, fn_train, fn_eval, maximize, Y, Y, ...  )
      df_evl <- bind_rows( df_evl, best_vars[['df_evl']] %>% mutate( k = k ) )
      agg_evl <- bind_rows( agg_evl, best_vars[['agg_evl']] %>% mutate( k = k ) )
      
      if( maximize ) {
        rr <- best_vars[['agg_evl']] %>% filter( J < mean_validation ) %>% filter( mean_validation == max(mean_validation, na.rm=TRUE ) )
        if( 0 < nrow(rr) ) {
          Y <- rr$selection[[1]]
          removing <- TRUE
          k <- k + 1
        }
      } else if( !maximize ) {
        rr <- best_vars[['agg_evl']] %>% filter( mean_validation < J ) %>% filter( mean_validation == max(mean_validation, na.rm=TRUE ) )
        if( 0 < nrow(rr) ) {
          Y <- rr$selection[[1]]
          removing <- TRUE
          k <- k + 1
        }
      }
      
    } # while
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
    kmax = kmax,
    
    # Results
    variable_selections = best_selections,
    results = df_evl,
    agg_results = agg_evl
  )
  return( ret )
}

sffs.formula <- function( fo, dat, 
                          fn_train = fn_train_binomial,
                          fn_eval  = fn_eval_binomial,
                          m = NULL,
                          ds = 5L, 
                          maximize = TRUE,
                          kmax = 100,
                          ... ) 
{
  # Check inputs
  stopifnot( is.formula( fo ) ) 
  stopifnot( is.data.frame(dat) | is_tibble(dat) )
  stopifnot( is.function(fn_train) )
  stopifnot( is.function(fn_eval) )
  
  resp <- as.character(lhs.vars( fo ) )
  vars <- as.character(rhs.vars( fo ) )
  
  sffs( dat = dat,
        resp = resp, 
        vars = vars,
        fn_train = fn_train,
        fn_eval  = fn_eval,
        m = NULL,
        ds = ds, 
        maximize = maximize,
        kmax = kmax,
        ... )
}

# Sequential Backward Forward Search ----
sfbs <- function( dat, resp, vars, 
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
  
  # Obtain evaluation splits, if necessary
  ds <- prepare_splits( ds, dat, resp, vars, fn_train, fn_eval, ... )
  
  sfbs_backward <- function( ds, dat, resp, vars, fn_train, fn_eval, Y, maximize, ... ) {
    
    df_evl <- map_dfr( .x=Y, .f=function( ff ) mutate( eval_splits( ds, dat, resp, vars, fn_train, fn_eval, setdiff( Y,ff ), ... ), removed = ff )  )
    
    agg_evl <- df_evl %>% 
      group_by( ch_selection, removed, selection ) %>% 
      summarise( mean_train = mean( eval_train, na.rm=TRUE ),
                 mean_validation = mean( eval_validation, na.rm=TRUE ),
                 mean_bias = mean( bias, na.rm=TRUE ) ) %>%
      ungroup
    
    ret <- list( df_evl = df_evl, agg_evl = agg_evl )
    if( all( is.na(agg_evl$mean_validation) ) ) {
      ret[['best_vars']] <- c( remaining_vars[1] )
      cat( sprintf( "No best found, removing first (%s) \n", ret[['best_vars']] ) )
    } else if( maximize ) {
      ret[['best_vars']] <- agg_evl %>% 
        filter( !is.na(mean_validation) ) %>% 
        filter( max( mean_validation, na.rm=TRUE ) <= mean_validation ) %>%
        pull( removed )
      cat( sprintf( "Best found, removing variable (%s) \n", ret[['best_vars']] ) )
    } else if( !maximize ) {
      ret[['best_vars']] <- agg_evl %>% 
        filter( !is.na(mean_validation) ) %>% 
        filter( mean_validation <= min( mean_validation, na.rm=TRUE ) ) %>%
        pull( removed )
      cat( sprintf( "Best found, removing variable (%s) \n", ret[['best_vars']] ) )
    }
    return( ret )
  }
  
  sfbs_forward <- function( ds, dat, resp, vars, fn_train, fn_eval, Y, J, maximize, ... ) {
    
    df_evl <- map_dfr( .x=setdiff(vars,Y), .f=function( ff ) mutate( eval_splits( ds, dat, resp, vars, fn_train, fn_eval, union( Y,ff ), ... ), added = ff )  )
    
    agg_evl <- df_evl %>% 
      group_by( ch_selection, added, selection ) %>% 
      summarise( mean_train = mean( eval_train, na.rm=TRUE ),
                 mean_validation = mean( eval_validation, na.rm=TRUE ),
                 mean_bias = mean( bias, na.rm=TRUE ) ) %>%
      ungroup
    
    ret <- list( df_evl = df_evl, agg_evl = agg_evl )
    if( maximize ) {
      ret[['best_vars']] <- agg_evl %>% 
        filter( !is.na(mean_validation) ) %>% 
        filter( (J < mean_validation) &  max( mean_validation, na.rm=TRUE ) <= mean_validation ) %>%
        pull( added )
      cat( sprintf( "Best found, adding variable (%s) \n", ret[['best_vars']] ) )
    } else if( !maximize ) {
      ret[['best_vars']] <- agg_evl %>% 
        filter( !is.na(mean_validation) ) %>% 
        filter( (J > mean_validation) & mean_validation <= min( mean_validation, na.rm=TRUE ) ) %>%
        pull( added )
      cat( sprintf( "Best found, adding variable (%s) \n", ret[['best_vars']] ) )
    }
    return( ret )
  }
  
  
  # Ref.: "Floating Search Methods for Feature Selection with Nonmonotonic Criterion Functions"; Pudil, Ferri, Novovicova, Kittler, 1994, IEEE
  # Sequential Floating Forward Selection algorithm:
  # 1. Y = {}
  # 2. Select the best feature
  #      x- = arg max_{x in Yk} J( Yk - x )
  #      Yk = Yk - x-; k = k + 1
  # 3. Select the worst feature* 
  #    [Note some book-keeping is necessary to avoid infinite loops]
  #      x+ = arg max_{x not in Yk} J( Yk + x )
  # 4. if J( Yk + x+ ) > J( Yk ) then
  #      Yk+1 = Yk + x+; k = k + 1
  #      go to step 3
  #    else
  #      go to step 2
  # 
  Y <- vars
  k <- 1
  while( length(Y) <= m & !any( duplicated( agg_evl$ch_selection ) ) ) {
    
    best_vars = sfbs_backward( dat, resp, vars, fn_train, fn_eval, Y, maximize, ... )
    df_evl <- bind_rows( df_evl, best_vars[['df_evl']] %>% mutate( step = k ) )
    agg_evl <- bind_rows( agg_evl, best_vars[['agg_evl']] %>% mutate( step = k ) )
    nvars <- best_vars[['best_vars']]
    
    cat( "Adding first of vars: ")
    cat( nvars )
    cat( "\n" )
    Y <- union( Y, nvars[1] )
    J <- agg_evl %>% filter( added == nvars[1] ) %>% pull( mean_validation )
    cat( sprintf( "Best backward evaluation is %1.4f \n", J ) )
    k <- k + 1
    
    adding <- TRUE
    while( adding ) {
      adding <- FALSE
      best_vars <- sfbs_forward( ds, dat, resp, vars, fn_train, fn_eval, Y, J, maximize, ... )
      df_evl <- bind_rows( df_evl, best_vars[['df_evl']] %>% mutate( step = k ) )
      agg_evl <- bind_rows( agg_evl, best_vars[['agg_evl']] %>% mutate( step = k ) )
      nvars <- best_vars[['best_vars']]
      if( 0 < length(nvars) ) {
        cat( sprintf( "Adding %s \n", nvars[1] ) )
        Y <- setdiff( Y, nvars[1] )
        adding <- TRUE
        k <- k + 1
      } # if
    } # while
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

sfbs.formula <- function( fo, dat, 
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
  
  sfbs( dat = dat,
        resp = resp, 
        vars = vars,
        fn_train = fn_train,
        fn_eval  = fn_eval,
        m = NULL,
        ds = ds, 
        maximize = maximize,
        ... )
}

