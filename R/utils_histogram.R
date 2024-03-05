# Density estimators
#
# Cross-validated bin width for histograms, Larry Wassermann, Chapter 6.2ff, p.127
# Carsten Henneges

Jhat <- function( x, m ) {
  
  # Code from https://www.stat.cmu.edu/~larry/all-of-statistics/=Rprograms/density.examples.r  
  ### histogram cross validation function
  # n <- length(x)
  # a <- min(x)
  # b <- max(x)
  # k <- 100
  # nbins <- seq(1,n,length=k)  ###number of bins
  # nbins <- round(nbins)
  # h <- (b-a)/nbins            ###width of bins
  # risk <- rep(0,k)
  # for(i in 1:k){
  #   ###get counts N_j
  #   br <- seq(a,b,length=nbins[i]+1)
  #   N <- hist(x,breaks=br,plot=F)$counts
  #   risk[i] <- sum(N^2)/(n^2*h[i])  - (2/(h[i]*n*(n-1)))*sum(N*(N-1))
  # }
  # hbest <- h[risk==min(risk)]
  # hbest <- hbest[1]  ###in case of tie take first (smallest) one
  # mbest <- (b-a)/hbest   ###optimal number of bins
  # list(risk=risk,nbins=nbins,h=h,mbest=mbest)
  
  n <- length( x )
  mx <- max( x )
  mn <- min( x )
  h <- (mx - mn) / m
  hh <- graphics::hist( x=x, 
                        # if breaks is set to the number of bins, it will be 
                        # processed by pretty() which doesn't give the right 
                        # buckets, thus we'll have to use seq() here
                        breaks = seq( min(x,na.rm=TRUE),
                                      max(x,na.rm=TRUE), 
                                      length.out = 1+m ), 
                        include.lowest=TRUE,
                        plot=FALSE )
  pp <- hh$counts / n
  # J <- 2.0 / ( h * (m - 1.0) ) - (m + 1.0)/( h * (m - 1.0) ) * sum( pp^2 )
  J <- sum(pp^2)/(n^2*h)  - (2/(h*n*(n-1)))*sum(pp*(pp-1)) # CHE/2024-02-21: Finally fixed.
  
  return( J )
}

Jloss_ucv <- function( x ) {
  n <- length( x )
  mmax <- sqrt( n )
  hh <- seq_len( mmax )
  jj <- vapply( hh, FUN=function( mm ) Jhat( x, mm ), 1.0 )
  oo <- which.min( jj ) 
  return( list( bins = hh, loss = jj, bins.opt = hh[ oo ], 
                loss.opt = jj[ oo ] ) )
} 

nbins_ucv <-  function( x ) {
  n <- length( x )
  mmax <- sqrt( n )
  hh <- seq_len( mmax )
  jj <- vapply( hh, FUN=function( mm ) Jhat( x, mm ), 1.0 )
  oo <- which.min( jj ) 
  return( hh[ oo ] )
} 

#' @title Helper function to determine the optimal number of breaks 
#' for a histogram
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
  bb <- graphics::hist( x=x, 
                        breaks = seq( min(x,na.rm=TRUE),
                                      max(x,na.rm=TRUE), 
                                      length.out = 1+hh[ oo ] ), 
                        include.lowest=TRUE,
                        plot=FALSE )
  return( bb$breaks )
} 

