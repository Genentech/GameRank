
library( formula.tools )
# Formula interface
forward_selection_formula <- function( dat, # data.frame or tibble
                                       fo,  # formula of response ~ all variables for selection
                                       # Training function
                                       fn_train = function( fo, dat, ... ) { 
                                         glm( fo, "binomial", dat ) 
                                       },
                                       # Evaluation function
                                       fn_eval  = function( fo, dat, mod, ... ) { 
                                         prfp <- performance( 
                                           prediction( predictions = predict( mod, newdata = dat, type="response"), 
                                                       labels = model.frame(fo,dat)[,1] ),
                                           "auc" )
                                         auc <- as.numeric( prfp@y.values )
                                         return( auc )
                                       },
                                       # can be a matrix with 1s (development) and 2s (validation) splits of nrow(dat) rows, an integer denoting the number of 
                                       # of predefined independent splits
                                       ds = 5L,
                                       n_free = 5L,
                                       max_iter = 1E6,
                                       ... )
{
  resp <- as.character( lhs.vars( fo ) )
  nvars <- as.character( rhs.vars( fo ) )
  forward_selection( dat,
                     resp, nvars,
                     fn_train, fn_eval, ds, n_free, max_iter,
                     ... )
}

forward_selection <- function( dat, # data.frame or tibble
                               # fo,  # formula of response ~ all variables for selection
                               resp,
                               nvars,
                               # Training function
                               fn_train = function( fo, dat, ... ) { 
                                 glm( fo, "binomial", dat ) 
                               },
                               # Evaluation function
                               fn_eval  = function( fo, dat, mod, ... ) { 
                                 prfp <- performance( 
                                   prediction( predictions = predict( mod, newdata = dat, type="response"), 
                                               labels = model.frame(fo,dat)[,1] ),
                                   "auc" )
                                 auc <- as.numeric( prfp@y.values )
                                 return( auc )
                               },
                               # can be a matrix with 1s (development) and 2s (validation) splits of nrow(dat) rows, an integer denoting the number of 
                               # of predefined independent splits
                               ds = 5L,
                               n_free = 5L,
                               max_iter = 1E6,
                               ...
) 
{
  stopifnot( is.data.frame(dat) | is_tibble(dat) )
  MIN_DATA_N <- 30
  stopifnot( nrow(dat) > MIN_DATA_N ) # Do not run this algorithm if too few data!
  stopifnot( is.character(resp) )
  stopifnot( is.character(nvars) )
  
  # stopifnot( is_formula(fo) )
  # resp <- as.character( lhs.vars( fo ) )
  # nvars <- as.character( rhs.vars( fo ) )
  fo <- formula( sprintf( "%s ~ %s", resp, paste(nvars, collapse = " + ") ) )
  
  # Backward selection algorithm ----
  # 1. - Randomly split into training and development (maybe prepare k independent splits)
  splits <- NULL
  if( is.matrix(ds) && nrow(ds)==nrow(dat) && ncol(ds) > 0 && setequal(unique(as.numeric(ds)),c(1,2)) ) {
    splits <- ds
  } else if( is.integer(ds) ) {
    # By default, we use a third for validation and two thirds for model building in this case
    splits <- matrix( rep_len( c(1,1,2), length.out = nrow(dat) * ds ), ncol=ds, nrow=nrow(dat), byrow = FALSE )
    for( i in 1:ncol(splits) ) { 
      # Note: It is key to obtain initial models for each split; that is why we invest 2/3rds for model building
      success <- FALSE
      j <- 0
      while( !success ) {
        cat( sprintf( "Searching initial model (split %d, round %d) \n", i, j ) )
        splits[, i ] <- splits[ order( runif( length( splits[,i] ) ) ) , i ]   
        mo <- tryCatch({
          mo <- fn_train( fo, dat[which(1==splits[,i]),], ... )
          mo
        }, error = function(e) NA )
        success = (!is.na(mo))
        j <- j + 1
      }
    }
  } else {
    cat( "Parameter ds:")
    cat( "is.integer : ", is.integer(ds), "\n")
    cat( "is.matrix : ", is.matrix(ds), "\n")
    cat( "nrow equal # data :", (nrow(ds)==nrow(dat)), "\n")
    cat( "ncols > 0  :", (ncol(ds) > 0), "\n")
    cat( "only 1s/2s :", (setequal(unique(as.numeric(ds)),c(1,2)) ), "\n")
  }
  stopifnot(!is.null(splits) || !is.matrix(splits))
  
  #    - Initialize set of variables
  rest <- nvars
  selection <- "1"
  #    - Initialize iteration counter
  t <- 0
  added <- TRUE
  perf <- NULL
  csel <- paste( selection, collapse = " + " )
  fo <- formula( sprintf( "%s ~ %s", resp, csel ) )
  for( j in 1:ncol(splits) ) {
    try({
      mo <- fn_train( fo, dat[which(1==splits[,j]),], ... )
      dp <-  fn_eval( fo, dat[which(1==splits[,j]),], mo, ... )
      vp <-  fn_eval( fo, dat[which(2==splits[,j]),], mo, ... )
      re <- tibble( train_performance = dp, validation_performance = vp, split = j, round = t, added_var = "", selection = as.character(csel) ) %>%
        mutate( bias = sqrt( (train_performance - validation_performance)^2 ) )
      perf <- dplyr::bind_rows(perf,re)
    })
  }
  t <- 1
  
  # 2. WHILE variables were added DO
  addition_log <- NULL
  while( (t < max_iter ) & 0 < length(rest) & added ) {
    added <- FALSE
    cat( sprintf( "Iteration %d: # of variables in current selection %s \n", t, length(selection) ) )
    
    # Get base performance selection for baseline 
    bperf <- NULL
    csel <- paste( selection, collapse = " + " )
    fo <- formula( sprintf( "%s ~ %s", resp, csel ) )
    for( j in 1:ncol(splits) ) {
      try({
        mo <- fn_train( fo, dat[which(1==splits[,j]),], ... )
        dp <-  fn_eval( fo, dat[which(1==splits[,j]),], mo, ... )
        vp <-  fn_eval( fo, dat[which(2==splits[,j]),], mo, ... )
        re <- tibble( train_performance = dp, validation_performance = vp, split = j, round = t, added_var = "", selection = as.character(csel) ) %>%
          mutate( bias = sqrt( (train_performance - validation_performance)^2 ) )
        bperf <- dplyr::bind_rows(bperf,re)
      })
    }
    bperf <- bperf %>% 
      mutate( last_bias = bias,
              last_train_performance = train_performance,
              last_validation_performance = validation_performance )
    
    #      - evaluate removal of each variable in all splits through independent models
    rperf <- NULL
    for( var in rest ) {
      cat( sprintf( "Evaluating var %s \n", var ) )
      nsel <- union( selection, var )
      csel <- paste( nsel, collapse = " + " )
      fo <- formula( sprintf( "%s ~ %s", resp, csel ) )
      for( j in 1:ncol(splits) ) {
        try({
          mo <- fn_train( fo, dat[which(1==splits[,j]),], ... )
          dp <-  fn_eval( fo, dat[which(1==splits[,j]),], mo, ... )
          vp <-  fn_eval( fo, dat[which(2==splits[,j]),], mo, ... )
          re <- tibble( train_performance = dp, validation_performance = vp, split = j, round = t, added_var = var, selection = as.character(csel) ) %>%
            mutate( bias = sqrt( (train_performance - validation_performance)^2 ) )
          rperf <- dplyr::bind_rows(rperf,re)
        })
      } # for
    } # for
    rperf <- rperf %>% 
      dplyr::left_join( bperf %>% dplyr::select( split, last_bias, last_train_performance, last_validation_performance ), "split" ) %>%
      mutate( diff_bias = bias - last_bias, 
              diff_train_performance = train_performance - last_train_performance,
              diff_validation_performance = validation_performance - last_validation_performance )
    
    perf <- perf %>% bind_rows( rperf )
    
    #      - check that bias is reduced to previous selection
    aperf <- rperf %>%
      group_by( added_var, selection, round ) %>%
      summarise( m_diff_bias = mean( diff_bias )  ) %>%
      ungroup %>% 
      # filter( m_diff_bias < 0 ) %>%
      filter( m_diff_bias == min( m_diff_bias ) )
    
    #      - remove variable from set that best improves bias between training and validation
    if( 0 < nrow(aperf) ) {
      cat( sprintf( "Adding variable %s with mean bias reduction %1.4f \n", 
                    aperf %>% pull( added_var ),
                    aperf %>% pull( m_diff_bias ) ) )
      nv <- aperf %>% pull( added_var )
      selection <- union( selection, nv )
      rest <- setdiff( rest, nv )
      addition_log <- addition_log %>% bind_rows( aperf %>% mutate( t = t ) )
      added <- TRUE
      if( 0 < aperf %>% pull( m_diff_bias ) ) {
        if( 0 < n_free ) {
          n_free <- n_free - 1
          cat( sprintf( "Reducing free, sub-optimal trails left to %d \n", n_free ) )
        } else {
          cat( "No more free, sub-optimal trails left, stopping selection \n")
          added <- FALSE # No more free, sub-optimal trials, so stop iteration
        }
      }
    } else {
      added <- FALSE
    }
    
    #      - if anything removed, continue with next selection turn
    t <- t + 1
  } # while
  #    OD 
  
  # Building return 
  cat( "Compiling results \n" )
  ret <- list(
    data = dat,
    formula = fo,
    response = resp,
    all_vars = nvars,
    
    ds = ds,
    max_iter = max_iter,
    
    splits = splits,
    variable_selection = selection,
    addition_log = addition_log,
    
    logs = perf
  )
  return( ret )
}