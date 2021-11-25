#
# File with util functions
#

# Construct m random splits that work to generate fn_train models if possible

#' @title Helper function to generate splits that are likely to produce models
#' 
#' @description The function should generate m random splits for dat, resp, vars such
#' that fn_train is likely to produce models. E.g. if there are a lot of missing values
#' in the dataset and splitting reduces these in an unfavorable way, model parameters for lm
#' or glm might not be estimable. That happens likely if the number of complete.cases is too low.
#' Therefore this function first separates the complete and incomplete cases, and then
#' randomly samples and distributes them into training and validation splits. This ensures that
#' both splits include 2/3rds and 1/3rd complete cases each.
#' 
#' @param m Number of parallel random splits to generate
#' @param dat data.frame or tibble comprising data for model generation and validation.
#' @param resp Character string defining the response (lhs) of the model formula.
#' @param vars Character vector defining the list of variables for selection. Those are concatenated by '+' 
#' as the right hand side (rhs) of the modelling formula.
#' @param fn_train Function with signature function( dat, resp, selection, ... ) that returns a model or NULL in any other case on the given data dat.
#' @param fn_eval Function with signature function( dat, resp, selection, ... ) that returns a real number or NA in any other case, e.g. when model is NULL.
#' @param ... An other arguments passed to fn_train or fn_eval during calls, e.g. maybe 'u = 365' for Survival evaluations specifying the landmark day.
#' 
#' @return Matrix with nrow(dat) rows and m columns comprising 1s and 2s to indicate if sample is part of training or validation split.
#' 
#' @export
build_splits <- function( m, dat, resp, vars, fn_train, fn_eval, ... ) {
  rr <- c( 1,1,2 ) # Use 2/3 train and 1/3 validation
  
  fo <- stats::formula( sprintf( "%s ~ %s ", resp, paste( vars, collapse = " + ")) )
  mf <- stats::model.frame( formula = fo, data = dat, na.action = stats::na.pass ) 
  
  # Check number of complete cases
  min_cc <- 25
  idx_cc <- which( stats::complete.cases( mf ) )
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
    rr_cc  <- rr_cc[order(stats::runif(length(rr_cc )))]
    rr_ncc <- rr_cc[order(stats::runif(length(rr_ncc)))]
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
  if( is.matrix(ds) && nrow(ds)==nrow(dat) && 1<=ncol(ds) && all(unique(unlist(ds) %in% c(1,2)) ) ) {
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


