# Density estimators
#
# Cross-validated bin width for histograms, Larry Wassermann, Chapter 6.2ff, p.127
# Carsten Henneges

Jhat <- function( x, m ) {
  n <- length( x )
  h <- 1 / m
  hh <- hist( x=x, breaks = m, plot=FALSE )
  pp <- hh$counts / n
  J <- 2.0 / ( h * (m - 1.0) ) - (m + 1.0)/( h * (m - 1.0) ) * sum( pp^2 )
  return( J )
}

Jloss_ucv <- function( x ) {
  n <- length( x )
  mmax <- sqrt( n )
  hh <- 1:mmax
  jj <- sapply( hh, FUN=function( mm ) Jhat( x, mm ) )
  oo <- which.min( jj ) 
  return( list( bins = hh, loss = jj, bins.opt = hh[ oo ], loss.opt = jj[ oo ] ) )
} 

nbins_ucv <-  function( x ) {
  n <- length( x )
  mmax <- sqrt( n )
  hh <- 1:mmax
  jj <- sapply( hh, FUN=function( mm ) Jhat( x, mm ) )
  oo <- which.min( jj ) 
  return( hh[ oo ] )
} 


bins_ucv <-  function( x ) {
  n <- length( x )
  mmax <- sqrt( n )
  hh <- 1:mmax
  jj <- sapply( hh, FUN=function( mm ) Jhat( x, mm ) )
  oo <- which.min( jj ) 
  bb <- hist( x=x, breaks = hh[ oo ], plot=FALSE )
  return( bb$breaks )
} 

