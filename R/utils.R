#
# File with util functions
#

# Construct m random splits that work to generate fn_train models if possible
build_splits <- function( m, dat, resp, vars, fn_train, fn_eval, ... ) {
  rr <- c( 1,1,2 ) # Use 2/3 train and 1/3 validation
  
  # Check number of complete cases
  min_cc <- 25
  idx_cc <- which( complete.cases( dat[,vars] ) )
  if( 0==length(idx_cc) | length(idx_cc) < min_cc ) { 
    warning( sprintf( "Not enough complete cases! (%d < %d)", length(idx_cc), min_cc ) )
    return( NULL ) 
  }

  # Check if complete cases do cover complete domain for each categorical variable
  for( var in vars ) {
    if( is.character(dat[,var]) | is.factor(dat[,var])) {
      if( ! setequal( unique(dat[idx_cc,var]), setdiff( unique( dat[,var] ), NA ) ) ) {
        cat( sprintf( "Variable %s does not cover full domain by complete cases.", var ) )
        return( NULL )
      }
    }
  }

  idx_ncc <- setdiff( 1:nrow(dat), idx_cc )
  idx <- c(idx_cc,idx_ncc)
  rr_cc  <- rep_len( rr, length(idx_cc) )
  rr_ncc <- rep_len( rr, length(idx_ncc) )

  ret <- matrix( NA, ncol=0, nrow=nrow(dat) )
  k <- 100
  while( 0<=k & ncol(ret) < m ) {
    # Independently sample from complete cases and non-complete cases
    rr_cc  <- rr_cc[order(runif(length(rr_cc )))]
    rr_ncc <- rr_cc[order(runif(length(rr_ncc)))]
    rr <- c(rr_cc, rr_ncc)
    
    dd <- matrix( NA, ncol=1, nrow=nrow(dat) )
    dd[idx[which(rr==1)],1] <- 1
    dd[idx[which(rr==2)],1] <- 2
    
    ck <- check_split( dd[,1], dat, resp, vars, fn_train, fn_eval, ... )
    if( ck ) {
      # If split is good, add it
      ret <- cbind( ret, dd )
    }
    k <- k - 1
  }
  return( ret )
}


# Check splits helper function
# 
check_split <- function( ds, dat, resp, vars, fn_train, fn_eval, ... ) {
  mo <- fn_train( dat[which(1==ds),], resp, vars, ... )   # Obtain model from data in 1-fold
  return( !is.null(mo) )
  # Remark: Convergence is *not* checked!
}

# Prepare splits helper function
prepare_splits <- function( ds, dat, resp, vars, fn_train, fn_eval, ... ) {
  rr <- c(1,1,2) # Split into 2/3 development (=1) and 1/3 validation (=2)
  
  # 1) If ds is integer, check and convert to 1-column matrix
  if( is.integer(ds) & length(ds)==nrow(dat) ) {
    ds <- matrix( ds, ncol=1, nrow=nrow(dat) )
    ret <- check_split( ds[,1], dat, resp, vars, fn_train, fn_eval, ... )
    if( ret ) return( ds ) else return( NULL )
  }
  # 2) if ds is matrix with enough rows and columns check and return it
  if( is.matrix(ds) && nrow(ds)==nrow(dat) && 1<=ncol(ds) && setequal( c(1,2), as.integer(unique(unlist(ds))) ) ) {
    ret <- TRUE
    for( k in 1:ncol(ds) ) {
      ret <- ret & check_split( ds[,k], dat, resp, vars, fn_train, fn_eval, ... )  
    }
    if( ret ) return( ds ) else return( NULL )
  }
  # 3) 
  if( is.integer(ds) ) {
    cat( sprintf( "Generating %d splits \n", ds ) )
    ds <- build_splits( ds, dat, resp, vars, fn_train, fn_eval, ... )
    if( !is.null(ds) )  return( ds )
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
