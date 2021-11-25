# Density estimators
#
# Cross-validated bin width for histograms, Larry Wassermann, Chapter 6.2ff, p.127
# Carsten Henneges

Jhat <- function( x, m ) {
  n <- length( x )
  h <- 1 / m
  hh <- graphics::hist( x=x, breaks = m, plot=FALSE )
  pp <- hh$counts / n
  J <- 2.0 / ( h * (m - 1.0) ) - (m + 1.0)/( h * (m - 1.0) ) * sum( pp^2 )
  return( J )
}

Jloss_ucv <- function( x ) {
  n <- length( x )
  mmax <- sqrt( n )
  hh <- seq_len( mmax )
  jj <- vapply( hh, FUN=function( mm ) Jhat( x, mm ), 1.0 )
  oo <- which.min( jj ) 
  return( list( bins = hh, loss = jj, bins.opt = hh[ oo ], loss.opt = jj[ oo ] ) )
} 

nbins_ucv <-  function( x ) {
  n <- length( x )
  mmax <- sqrt( n )
  hh <- seq_len( mmax )
  jj <- vapply( hh, FUN=function( mm ) Jhat( x, mm ), 1.0 )
  oo <- which.min( jj ) 
  return( hh[ oo ] )
} 

#' @title Helper function to determine the optimal number of breaks for a histogram
#' 
#' @description This function determines the optimal set of breaks for a 
#' histogram for numerical sample
#' by leave-one-out cross-validation minimizing the empirical risk.
#' 
#' @param x A numeric vector of variables for which the set of breaks should 
#' be computed.
#' 
#' @return A numeric vector of cut-points points that will define an optimal 
#' histogram density estimator for the sample.
#' 
#' @examples 
#' library( ggplot2 )
#' xv <- rnorm( 100 )
#' bk <- bins_ucv( xv )
#' # Plot histogram with properly selected smoothed bin widths
#' ggplot( aes( x=xv, y=..density.. ), data = data.frame( xv = xv ) ) + 
#'    geom_histogram( breaks = bk )
#' @export
bins_ucv <-  function( x ) {
  n <- length( x )
  mmax <- sqrt( n )
  hh <- seq_len( mmax )
  jj <- vapply( hh, FUN=function( mm ) Jhat( x, mm ), 1.0 )
  oo <- which.min( jj ) 
  bb <- graphics::hist( x=x, breaks = hh[ oo ], plot=FALSE )
  return( bb$breaks )
} 

