#
# Forward Selection wrapper algorithm
#
source( "R/utils.R" )

#
# Conduct one forward wrapper algorithm
#
forward <- function(  dat, resp, vars, 
                      fn_train = fn_train_binomial,
                      fn_eval  = fn_eval_binomial,
                      ds = 5L, 
                      max_evals = 1E8,
                      maximize = TRUE,
                      ... ) {
  # Check parameters
  stopifnot( is.data.frame(dat) | is_tibble(dat) )
  stopifnot( is.character(resp) )
  stopifnot( is.character(vars)  & 1 < length(vars) )
  stopifnot( is.logical(maximize) )
  
  # Obtain evaluation splits
  ds <- prepare_splits( ds, dat, resp, vars, fn_train, fn_eval, ... )

  col_sep <- ","
  results <- NULL
  queue <- list( list( selection = c() ) )
  t <- max_evals
  while( 0 <= t & 1 <= length(queue) ) {
    partition <- queue[[1]]
    queue[[1]] <- NULL
    
    part_sel <- partition[["selection"]]
    cur_vars <- setdiff( vars, part_sel )
    if( 1<=length(cur_vars) ) {
      # a) Evaluate removal per variable this round
      rnd <- NULL
      for( var in cur_vars ) {
        sel <- union( var, part_sel )
        evl <- eval_splits( ds, dat, resp, sel, fn_train, fn_eval, ... )
        evl$added <- var
        evl$n_evaluation <- t
        t <- t - 1 
        rnd <- dplyr::bind_rows(rnd,evl)
        results <- dplyr::bind_rows(results,evl)
      }
      # b) Determine partitions for next rounds
      agg <- rnd %>%
        group_by( added ) %>%
        summarise( mean_train = mean( eval_train ),
                   mean_validation = mean( eval_validation ) ) %>%
        ungroup %>%
        arrange( desc(mean_validation) )
      
      if( maximize ) {
        nxt <- agg %>% filter( mean_validation >= max( mean_validation) )  
      } else {
        nxt <- agg %>% filter( mean_validation <= min( mean_validation) )  
      }
      
      cat( sprintf( "Addint %d new paritions to queue (size %d) \n", nrow(nxt), length(queue) ))
      cat( sprintf( "Adding paritions for added variables %s \n", nxt$added ))
      for( var in nxt$added ) {
        nxt_part <- union( var, part_sel )
        queue[[length(queue)+1]] <- list( selection = nxt_part ) # Adding partition
      }
    }
    cat( sprintf( "Completed round of evaluations, length of queue %d, remaining evaluations %d \n", length(queue), t ) )
  } # while
  
  agg <- results %>%
    group_by( ch_selection, size ) %>%
    summarise( mean_train = mean( eval_train ),
               mean_validation = mean( eval_validation ) ) %>%
    ungroup %>%
    arrange( desc(mean_validation), size )
  
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
    splits = ds, 
    max_evals = max_evals,
    maximize = maximize,
    
    # Results
    variable_selections = best_selections,
    selection_results = best_results,
    results = results
  )
  return( ret )
} # forward (END)

forward.formula <- function( fo, dat, 
                              fn_train = fn_train_binomial,
                              fn_eval  = fn_eval_binomial,
                              ds = 5L, 
                              max_evals = 1E8,
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
            ds = ds, 
            max_evals = max_evals,
            maximize = maximize,
            ... )
}