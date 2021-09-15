#
# File with util functions
#

# Check splits helper function
check_split <- function( ds, dat, resp, vars, fn_train, fn_eval, ... ) {
  ret <- tryCatch({
    mod <- fn_train( dat[which(1==ds),], resp, vars, ... )   # Obtain model from data in 1-fold
    evl <- fn_eval(  dat[which(2==ds),], resp, vars, mod, ... ) # Evaluate model on data in 2-fold
    TRUE
  }, error = function(e) FALSE )
  return( ret )
}

# Prepare splits helper function
prepare_splits <- function( ds, dat, resp, vars, fn_train, fn_eval, ... ) {
  rr <- c(1,1,2) # Split into 2/3 development (=1) and 1/3 validation (=2)
  
  # 1) If ds is integer, check and convert to 1-column matrix
  if( is.integer(ds) & length(ds)==nrow(dat) ) {
    ds <- matrix( ds, ncol=1, nrow=nrow(dat) )
  }
  # 2) if ds is matrix with enough rows and columns check and return it
  if( is.matrix(ds) && nrow(ds)==nrow(dat) && 1<=ncol(ds) && setequal( c(1,2), as.integer(unique(unlist(ds))) ) ) {
    check <- TRUE
    for( c in 1:ncols(ds) ) {
      check <- check & check_split( ds[,c], dat, resp, vars, fn_train, fn_eval, ... )  
    }
    if( check )  return( ds )
    stop( "Not all splits are fit for evaluation with all variables! \n" )
  }
  # 3) 
  if( is.integer(ds) ) {
    cat( sprintf( "Generating %d splits \n", ds ) )
    mm <- matrix( NA, ncol = ds, nrow = nrow(dat) ) 
    aa <- rep_len( rr, length.out = nrow(dat) )
    for( k in 1:ds ) {
      aa <- aa[ order( runif( length(aa) ) ) ]
      run <- check_split( aa, dat, resp, vars, fn_train, fn_eval, ... )
      while( !run ) {
        aa <- aa[ order( runif( length(aa) ) ) ]
        run <- check_split( aa, dat, resp, vars, fn_train, fn_eval, ... )
      }
      mm[,k] <- aa
    }
    return( mm )
  }
  stop( "Could not generate valid splits \n" )
}

eval_splits <- function( ds, dat, resp, selection, fn_train, fn_eval, ... ) {
  var_sep <- ","
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
                  size = cnt,
                  response = resp, 
                  split = k,
                  eval_train = evl1,
                  eval_validation = evl2,
                  bias = sqrt( (evl1 - evl2)^2 ) )
    ret <- bind_rows( ret, rr )
  } # for
  return( ret )
} # eval_splits (END)
