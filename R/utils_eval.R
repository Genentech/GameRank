#
# Utility evaluation
#

#' @title Split evaluation function
#' 
#' @description 
#' 
#' @export
eval_splits <- function( ds, dat, resp, selection, fn_train, fn_eval, ..., var_sep = "," ) 
{
  ret <- NULL
  for( k in 1:ncol(ds) ) {
    mod  <- fn_train( dat[which(1==ds[,k]),], resp, selection, ... )   # Obtain model from data in 1-fold
    evl1 <- fn_eval(  dat[which(1==ds[,k]),], resp, selection, mod, ... ) # Evaluate model on data in 1-fold
    evl2 <- fn_eval(  dat[which(2==ds[,k]),], resp, selection, mod, ... ) # Evaluate model on data in 2-fold
    selection <- sort( selection )
    txt_selection <- paste( selection, collapse = var_sep )
    cnt <- length(selection)
    rr <- tibble( selection = list( selection ),
                  ch_selection = txt_selection,
                  m = cnt,
                  response = resp, 
                  split = k,
                  eval_train = evl1,
                  eval_validation = evl2,
                  bias = sqrt( (evl1 - evl2)^2 ) )
    ret <- bind_rows( ret, rr )
  } # for
  return( ret )
} # eval_splits (END)


agg_evals <- function( df_evl, var, maximize ) {
  agg <- df_evl %>%
    group_by( across( c("ch_selection",  as.character( var ), "selection") ) ) %>%
    summarise( m = unique(m),
               mean_train = mean( eval_train, na.rm=TRUE ),
               mean_validation = mean( eval_validation, na.rm=TRUE ),
               mean_bias = mean( bias, na.rm=TRUE ) ) %>%
    ungroup
  if( maximize ) {
    agg <- agg %>% mutate( opt = ( mean_validation >= max(mean_validation, na.rm=TRUE) ) )
  } else if( !maximize ) {
    agg <- agg %>% mutate( opt = ( mean_validation <= min(mean_validation, na.rm=TRUE) ) )
  } else {
    agg <- agg %>% mutate( opt = NA )
  }

  return( agg )
}

best_selection <- function( agg ) {
  sel <- agg %>% filter( !is.na(opt) & (TRUE==opt) ) %>% pull( selection )
  return( sel )
}

eval_add_vars <- function( ds, dat, resp, vars, fn_train, fn_eval, maximize, selection, add_vars = NULL, ..., var_sep = ","  ) {
  if( is.null(add_vars) ) add_vars <- setdiff( vars, selection )
  
  df_evl <- purrr::map_dfr( .x=add_vars, .f=function(vv) mutate( eval_splits(ds,dat,resp,union(selection,vv),fn_train,fn_eval, ... ), added = vv ) )
  agg_evl <- agg_evals( df_evl, "added", maximize )
  best_selections <- best_selection( agg_evl )
  
  return( list( df_evl = df_evl, agg_evl = agg_evl, best_selections = best_selections ) )
}

eval_remove_vars <- function( ds, dat, resp, vars, fn_train, fn_eval, maximize, selection, remove_vars = NULL, ..., var_sep = ","  ) {
  if( is.null(remove_vars) ) remove_vars <- selection
  
  df_evl <- purrr::map_dfr( .x=remove_vars, .f=function(vv) mutate( eval_splits(ds,dat,resp,setdiff(selection,vv),fn_train,fn_eval, ... ), removed = vv ) )
  agg_evl <- agg_evals( df_evl, "removed", maximize )
  best_selections <- best_selection( agg_evl )
  
  return( list( df_evl = df_evl, agg_evl = agg_evl, best_selections = best_selections ) )
}


