#
# Utility evaluation
#
#' @import labelled dplyr


#' @title Split evaluation function
#' 
#' @description Function receives a matrix with one or more training:validation splits, defined by 1 as training and
#' 2 as validation. This matrix has nrow(dat) rows and evaluates fn_train and fn_eval for each parallel defined split
#' for a given selection. All wraper algorithms except GroupRank make use of this function.
#' 
#' @param ds Matrix with n rows and k columns containing only 1s and 2s. n must equal nrow(dat) and 1<=k.
#' @param dat data.frame or tibble comprising data for model generation and validation.
#' @param resp Character string defining the response (lhs) of the model formula.
#' @param selection Character vector defining the list of variables for selection. Those are concatenated by '+' 
#' as the right hand side (rhs) of the modelling formula.
#' @param fn_train Function with signature function( dat, resp, selection, ... ) that returns a model or NULL in any other case on the given data dat.
#' @param fn_eval Function with signature function( dat, resp, selection, ... ) that returns a real number or NA in any other case, e.g. when model is NULL.
#' @param ... Further arguments passed to fn_train and fn_eval.
#' @param var_sep Separating character used to build unique concatenated string for selection, which is saved in column ch_selection of the output.
#' 
#' @return tibble with elements
#' \describe{
#'  \item{selection}{list entry with input selection}
#'  \item{ch_selection}{character string with input selection separated by var_sep}
#'  \item{m}{number of features in selection}
#'  \item{response}{As from input parameters}
#'  \item{split}{index into ds columns determining split evaluated}
#'  \item{eval_train}{result of model evaluation on training}
#'  \item{eval_validation}{result of model evaluation on validation (this is used for optimization)}
#'  \item{bias}{square root of squared difference between training and validation values}
#' }
#' 
eval_splits <- function( ds, dat, resp, selection, fn_train, fn_eval, ..., var_sep = "," ) 
{
  ret <- NULL
  ret <- purrr::map_dfr( 1:ncol(ds), function(k) {
      mod  <- fn_train( dat[which(1==ds[,k]),], resp, selection, ... )   # Obtain model from data in 1-fold
      evl1 <- fn_eval(  dat[which(1==ds[,k]),], resp, selection, mod, ... ) # Evaluate model on data in 1-fold
      evl2 <- fn_eval(  dat[which(2==ds[,k]),], resp, selection, mod, ... ) # Evaluate model on data in 2-fold
      selection <- base::sort( selection )
      txt_selection <- base::paste( selection, collapse = var_sep )
      cnt <- base::length(selection)
      rr <- tibble( selection = list( selection ),
                    ch_selection = txt_selection,
                    m = cnt,
                    response = resp, 
                    split = k,
                    eval_train = evl1,
                    eval_validation = evl2,
                    bias = sqrt( (evl1 - evl2)^2 ) )
  })
  # for( k in 1:ncol(ds) ) {
  #   mod  <- fn_train( dat[which(1==ds[,k]),], resp, selection, ... )   # Obtain model from data in 1-fold
  #   evl1 <- fn_eval(  dat[which(1==ds[,k]),], resp, selection, mod, ... ) # Evaluate model on data in 1-fold
  #   evl2 <- fn_eval(  dat[which(2==ds[,k]),], resp, selection, mod, ... ) # Evaluate model on data in 2-fold
  #   selection <- sort( selection )
  #   txt_selection <- paste( selection, collapse = var_sep )
  #   cnt <- length(selection)
  #   rr <- tibble( selection = list( selection ),
  #                 ch_selection = txt_selection,
  #                 m = cnt,
  #                 response = resp, 
  #                 split = k,
  #                 eval_train = evl1,
  #                 eval_validation = evl2,
  #                 bias = sqrt( (evl1 - evl2)^2 ) )
  #   ret <- bind_rows( ret, rr )
  # } # for
  return( ret )
} # eval_splits (END)

#' @title Aggregation function for eval_splits results
#' 
#' @description Aggregates results from eval_splits call by considering additional grouping variables and
#' whether a maximization or minimization is required. Flags optimal selection across all results.
#' 
#' @param df_eval Results from one or more eval_splits calls.
#' @param var Additional grouping variable for results aggregation, e.g. may be variable 'added' or 'removed' during forward or backward selection.
#' @param maximize Logical value indicating if optimal result is largest (TRUE) or smallest result (FALSE)
#' 
#' @return tibble with aggregated (averaged) training and validation errors
#' \describe{
#' \item{ch_selection}{Column ch_selection from eval_splits call}
#' \item{var}{Columns passed as var parameter}
#' \item{selection}{List column of selections, that is from eval_splits call}
#' \item{m}{Size of aggregated paritions}
#' \item{mean_train}{Averaged eval_train results, ignoring NAs}
#' \item{mean_validation}{Averaged eval_validation results, ignoring NAs}
#' \item{mean_bias}{Averaged bias results, ignoring NAs}
#' }
#'
agg_evals <- function( df_evl, var, maximize ) {
  agg <- df_evl %>%
    group_by( across( c("ch_selection",  as.character( var ), "selection") ) ) %>%
    summarise( m = unique(m),
               mean_train = mean( eval_train, na.rm=TRUE ),
               mean_validation = mean( eval_validation, na.rm=TRUE ),
               mean_bias = mean( bias, na.rm=TRUE ) ) %>%
    ungroup
  if( maximize ) {
    agg <- agg %>% dplyr::mutate( opt = ( mean_validation >= max(mean_validation, na.rm=TRUE) ) )
  } else if( !maximize ) {
    agg <- agg %>% dplyr::mutate( opt = ( mean_validation <= min(mean_validation, na.rm=TRUE) ) )
  } else {
    agg <- agg %>% dplyr::mutate( opt = NA )
  }

  return( agg )
}

#' @title Helper function to extract best selections
#' 
#' @param agg Result from agg_evals call
#' 
#' @return selection column for rows flagged as opt (=TRUE)
best_selection <- function( agg ) {
  sel <- agg %>% filter( !is.na(opt) & (TRUE==opt) ) %>% pull( selection )
  return( sel )
}

#' @title Generic helper function that evaluates a set of variables for extending a selection
#' 
#' @param ds Definition of (parallel) training:validation splits by a matrix with d columns containing 1s and 2s, where 1 denotes sample is used for training the model and 2 denotes sample used for validation.
#' @param dat data.frame or tibble comprising data for model generation and validation.
#' @param resp Character string defining the response (lhs) of the model formula.
#' @param vars Character vector defining the list of variables for selection. Those are concatenated by '+' 
#' as the right hand side (rhs) of the modelling formula.
#' @param fn_train Function with signature function( dat, resp, selection, ... ) that returns a model or NULL in any other case on the given data dat.
#' @param fn_eval Function with signature function( dat, resp, selection, ... ) that returns a real number or NA in any other case, e.g. when model is NULL.
#' @param maximize A logic value determining if fn_eval is maximized (set to TRUE) or minimized (set to FALSE).
#' @param selection Character vector of currently selected variables.
#' @param add_vars Character vector of variables that should be used to extend and evaluate selection.
#' @param ... An other arguments passed to fn_train or fn_eval during calls, e.g. maybe 'u = 365' for Survival evaluations specifying the landmark day.
#' @param var_sep Character used to generated standardized string for selection by concatenation
#' 
#' @return List with elements
#' \describe{
#' \item{df_evl}{tibble with all results from eval_split calls being concatenated}
#' \item{agg_evl}{tibble with aggregated evaluations by (ch_)selection and additional column 'added'. Best selections are flagged by 'opt'} 
#' \item{best_selections}{List of selections that were identified best on validation split.}
#' }
#' 
eval_add_vars <- function( ds, dat, resp, vars, fn_train, fn_eval, maximize, selection, add_vars = NULL, ..., var_sep = ","  ) {
  if( is.null(add_vars) ) add_vars <- base::setdiff( vars, selection )
  if( 0 == base::length(add_vars) ) return( NULL )
  
  df_evl <- purrr::map_dfr( .x=add_vars, .f=function(vv) mutate( eval_splits(ds,dat,resp,union(selection,vv),fn_train,fn_eval, ... ), added = vv ) )
  agg_evl <- agg_evals( df_evl, "added", maximize )
  best_selections <- best_selection( agg_evl )
  
  return( list( df_evl = df_evl, agg_evl = agg_evl, best_selections = best_selections ) )
}

#' @title Generic helper function that evaluates a set of variables for reducing a selection
#' 
#' @param ds Definition of (parallel) training:validation splits by a matrix with d columns containing 1s and 2s, where 1 denotes sample is used for training the model and 2 denotes sample used for validation.
#' @param dat data.frame or tibble comprising data for model generation and validation.
#' @param resp Character string defining the response (lhs) of the model formula.
#' @param vars Character vector defining the list of variables for selection. Those are concatenated by '+' 
#' as the right hand side (rhs) of the modelling formula.
#' @param fn_train Function with signature function( dat, resp, selection, ... ) that returns a model or NULL in any other case on the given data dat.
#' @param fn_eval Function with signature function( dat, resp, selection, ... ) that returns a real number or NA in any other case, e.g. when model is NULL.
#' @param maximize A logic value determining if fn_eval is maximized (set to TRUE) or minimized (set to FALSE).
#' @param selection Character vector of currently selected variables.
#' @param remove_vars Character vector of variables that should be used to reduce and evaluate selection.
#' @param ... An other arguments passed to fn_train or fn_eval during calls, e.g. maybe 'u = 365' for Survival evaluations specifying the landmark day.
#' @param var_sep Character used to generated standardized string for selection by concatenation
#' 
#' @return List with elements
#' \describe{
#' \item{df_evl}{tibble with all results from eval_split calls being concatenated}
#' \item{agg_evl}{tibble with aggregated evaluations by (ch_)selection and additional column 'removed'. Best selections are flagged by 'opt'} 
#' \item{best_selections}{List of selections that were identified best on validation split.}
#' }
#' 
eval_remove_vars <- function( ds, dat, resp, vars, fn_train, fn_eval, maximize, selection, remove_vars = NULL, ..., var_sep = ","  ) {
  if( is.null(remove_vars) ) remove_vars <- selection
  if( 0 == length(remove_vars) ) return( NULL )
  
  df_evl <- purrr::map_dfr( .x=remove_vars, .f=function(vv) mutate( eval_splits(ds,dat,resp,setdiff(selection,vv),fn_train,fn_eval, ... ), removed = vv ) )
  agg_evl <- agg_evals( df_evl, "removed", maximize )
  best_selections <- best_selection( agg_evl )
  
  return( list( df_evl = df_evl, agg_evl = agg_evl, best_selections = best_selections ) )
}


