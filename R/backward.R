#
# Backward Selection wrapper algorithm
#
source( "R/utils.R" )

#
# Conduct one backward wrapper algorithm
#
backward <- function( dat, resp, vars, 
                      fn_train = fn_train_binomial,
                      fn_eval  = fn_eval_binomial,
                      ds = 5L, 
                      max_evals = 1E8,
                      maximize = TRUE, 
                      if_na_take_last_k = 1,
                      ... ) {
  # Check parameters
  stopifnot( is.data.frame(dat) | is_tibble(dat) )
  stopifnot( is.character(resp) )
  stopifnot( is.character(vars)  & 1 < length(vars) )
  stopifnot( is.numeric(max_evals) )
  stopifnot( is.logical(maximize) )
  
  # Obtain evaluation splits
  ds <- prepare_splits( ds, dat, resp, vars, fn_train, fn_eval, ... )
  
  ranks <- unique( tibble( variable = vars, rank = 1:length(vars) ) )
  
  col_sep <- ","
  results <- NULL
  queue <- list( list( selection = vars ) )
  t <- max_evals
  while( 0 <= t & 1 <= length(queue) ) {
    partition <- queue[[1]]
    queue[[1]] <- NULL
    
    cur_vars <- unique( partition[["selection"]] )
    cat( sprintf( "Current selection " ) )
    cat( cur_vars, "\n" )
    if( 2<=length(cur_vars) ) {
      # a) Evaluate removal per variable this round
      rnd <- NULL
      for( var in cur_vars ) {
        sel <- setdiff( cur_vars, var )
        evl <- eval_splits( ds, dat, resp, sel, fn_train, fn_eval, ... )
        evl$removed <- var
        evl$n_evaluation <- t
        t <- t - 1 
        rnd <- dplyr::bind_rows(rnd,evl)
        results <- dplyr::bind_rows(results,evl)
        cat( sprintf( "Evaluationg removal of %s (%1.4f) \n", var, mean( evl$eval_train, na.rm=TRUE ) ) )
      }
      # b) Determine partitions for next rounds
      # browser()
      agg <- rnd %>%
        group_by( removed ) %>%
        summarise( mean_train = mean( eval_train ),
                   mean_validation = mean( eval_validation, na.rm=TRUE ) ) %>%
        ungroup %>%
        arrange( desc(mean_validation) ) %>% 
        dplyr::left_join( ranks, c("removed"="variable") )
      
      if( all( is.na( agg$mean_validation ) ) ) {
        nxt <- agg %>% filter( max(rank) - if_na_take_last_k < rank )
      } else if( maximize ) {
        nxt <- agg %>% filter( mean_validation >= max( mean_validation, na.rm=TRUE) )  
      } else {
        nxt <- agg %>% filter( mean_validation <= min( mean_validation, na.rm=TRUE) )  
      }
      
      cat( sprintf( "Adding %d new paritions to queue (size %d) \n", nrow(nxt), length(queue) ))
      cat( sprintf( "Adding paritions for removed variable %s \n", nxt$removed ))
      for( var in nxt$removed ) {
        nxt_part <- setdiff( cur_vars, var )
        queue[[length(queue)+1]] <- list( selection = unique( nxt_part ) ) # Adding partition
      }
    }
    cat( sprintf( "Completed round of evaluations, length of queue %d \n", length(queue) ) )
  } # while
  
  agg <- results %>%
    group_by( ch_selection, size ) %>%
    summarise( mean_train = mean( eval_train, na.rm=TRUE ),
               mean_validation = mean( eval_validation, na.rm=TRUE ) ) %>%
    ungroup %>%
    arrange( desc(mean_validation), size )
  
  if( maximize ) {
    best <- agg %>% filter( mean_validation >= max( mean_validation, na.rm=TRUE) )  
  } else {
    best <- agg %>% filter( mean_validation <= min( mean_validation, na.rm=TRUE) )  
  }
  
  best_results <- results %>% 
    filter( ch_selection %in% best$ch_selection ) 
  best_selections <- best_results %>%
    pull( selection ) %>%
    unique
  cat( "Best selection(s) found is/are : \n" )
  print( best_results )
  
  ret <- list( # Input data
    # data = dat,
    response = resp,
    variables = vars,
    
    # Input parameters
    splits = ds, 
    max_evals = max_evals,
    maximize = maximize,
    
    # Results
    variable_selections = best_selections,
    selection_results = best_results,
    results = results
  )
  return( ret )
} # backward (END)

backward.formula <- function( fo, dat, 
                              fn_train = fn_train_binomial,
                              fn_eval  = fn_eval_binomial,
                              ds = 5L, 
                              max_evals = 1E8,
                              maximize = TRUE,
                              if_na_take_last_k = 1,
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
            ds = ds, 
            max_evals = max_evals,
            maximize = maximize,
            if_na_take_last_k = if_na_take_last_k,
            ... )
}