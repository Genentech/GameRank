#
#  Random sampling wrapper algorithm
#
random_selection <- function( dat,
                              resp,
                              vars,
                              fn_train = fn_train_binomial,
                              fn_eval  = fn_eval_binomial,
                              ds = 5L, 
                              ksize = 5L,
                              nevals = 100L,
                              maximize = TRUE,
                              ... )
{
  # Check inputs ----
  stopifnot( is.data.frame(dat) | is_tibble(dat) )
  stopifnot( is.character(vars) & 1 < length(vars) )
  stopifnot( is.character(vars) )
  
  stopifnot( is.function(fn_train) ) 
  stopifnot( is.function(fn_eval) ) 
  
  stopifnot( is.integer(ksize) )
  stopifnot( is.integer(nevals) )
  stopifnot( is.logical(maximize) )
  
  # Obtain evaluation splits
  ds <- prepare_splits( ds, dat, resp, vars, fn_train, fn_eval, ... )
  
  # Prepare sampling plan with the following properties:
  #  - Each variable is considered at least min_sample_per_var times
  #  - Variables are evaluated in groups of ksize sets
  # Let a be # variables, b be group size, c be min # each var evaluated then
  samps <- comboSample( v = vars, m = ksize, repetition = FALSE, n = nevals )
  
  # Evaluate random combinations
  results <- NULL
  for( r in 1:nrow(samps) ) {
    cat( sprintf( "Evaluating selection %d \n", r )  )
    sel <- as.character( samps[r, ] )
    evl <- eval_splits( ds, dat, resp, sel, fn_train, fn_eval, ... )
    evl$row <- r
    results <- bind_rows( results, evl )
  }
  
  agg <- results %>%
    group_by( ch_selection, row ) %>%
    summarise( mean_train = mean( eval_train ),
               mean_validation = mean( eval_validation ) ) %>%
    ungroup %>%
    arrange( desc(mean_validation), row )
  
  if( maximize ) {
    best <- agg %>% filter( mean_validation >= max( mean_validation) )  
  } else {
    best <- agg %>% filter( mean_validation <= min( mean_validation) )  
  }
  
  best_results <- results %>% 
    filter( ch_selection %in% best$ch_selection ) 
  best_selections <- best_results %>%
    pull( selection ) %>%
    unique
  
  ret <- list( # Input data
    # data = dat,
    response = resp,
    variables = vars,
    
    # Input parameters
    ksize = ksize,
    nevals = nevals,
    maximize = maximize,
    
    # Results
    variable_selections = best_selections,
    selection_results = best_results,
    results = results
  )
  return( ret )
}
  

random_selection.formula <- function( fo,
                              dat,
                              fn_train = fn_train_binomial,
                              fn_eval  = fn_eval_binomial,
                              ds = 5L, 
                              ksize = 5L,
                              nevals = 100L,
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
  
  random_selection( dat = dat,
            resp = resp, 
            vars = vars,
            fn_train = fn_train,
            fn_eval  = fn_eval,
            ds = ds, 
            ksize = ksize,
            nevals = nevals,
            maximize = maximize,
            ... )
}