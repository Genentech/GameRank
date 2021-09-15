
#
# Write custom summary parts for TLS reporting
#

summary.logical <- function( ... ) {
  
  x.all <- as.logical( as.list( ... ) ) 
  
  x <- x.all[which(!is.na(x.all))]
  if( 0 == length(x) ) {
    ret <- c( 
      N = as.integer( length( x.all ) ),
      n = as.integer( length( which(!is.na(x.all)) ) )
    )
    return( ret )
  }
  
  N = as.integer( length( x.all ) )
  n.true <- length( which( x.all ) )
  n.false <- length( which( !x.all ) )
  ret <- c( 
    N = N,
    n = as.integer( length( which(!is.na(x.all)) ) ),
    "TRUE" = sprintf( "%d (%1.1f%%)", n.true, (100*n.true/N) ),
    "FALSE" = sprintf( "%d (%1.1f%%)", n.false, (100*n.false/N) )
  )
  return( ret ) 
}

summary.numeric <- function( ... ) {
  
  x.all <- as.numeric( as.list( ... ) ) 
  
  x <- x.all[which(!is.na(x.all))]
  if( 0 == length(x) ) {
    ret <- c( 
      N = as.integer( length( x.all ) ),
      n = as.integer( length( which(!is.na(x.all)) ) ),
      
      mean = NA,
      sd = NA,
      
      min = NA,
      q1 = NA,
      median = NA,
      q3 = NA,
      max = NA
    )
    return( ret )
  }
  
  ret <- c(
    N = as.integer( length( x.all ) ),
    n = as.integer( length( which(!is.na(x.all)) ) ),
    
    mean = as.numeric( mean( x, na.rm=TRUE ) ),
    sd = as.numeric( sd( x, na.rm=TRUE ) ),
    
    min = as.numeric( min( x, na.rm=TRUE ) ),
    q1 = as.numeric( quantile( x, 0.25, na.rm=TRUE ) ),
    median = as.numeric( median( x, na.rm=TRUE ) ),
    q3 = as.numeric( quantile( x, 0.75, na.rm=TRUE ) ),
    max = as.numeric( max( x, na.rm=TRUE ) )
  )  
  
  return( ret )
}



summary.integer <- function( ... ) {
  
  x.all <- as.numeric( as.list( ... ) ) 
  
  x <- x.all[which(!is.na(x.all))]
  if( 0 == length(x) ) {
    ret <- c( 
      N = as.integer( length( x.all ) ),
      n = as.integer( length( which(!is.na(x.all)) ) ),
      
      mean = NA,
      sd = NA,
      
      min = NA,
      q1 = NA,
      median = NA,
      q3 = NA,
      max = NA
    )
    return( ret )
  }
  
  ret <- c(
    N = as.integer( length( x.all ) ),
    n = as.integer( length( which(!is.na(x.all)) ) ),
    
    mean = as.numeric( mean( x, na.rm=TRUE ) ),
    sd = as.numeric( sd( x, na.rm=TRUE ) ),
    
    min = as.numeric( min( x, na.rm=TRUE ) ),
    q1 = as.numeric( quantile( x, 0.25, na.rm=TRUE ) ),
    median = as.numeric( median( x, na.rm=TRUE ) ),
    q3 = as.numeric( quantile( x, 0.75, na.rm=TRUE ) ),
    max = as.numeric( max( x, na.rm=TRUE ) )
  )  
  
  return( ret )
}
