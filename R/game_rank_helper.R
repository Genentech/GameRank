
#
# Helper file for game_rank.R comprising the maximum likelihood function
# and it's gradient.
#

#' @title  Evaluating Group Rank Maximum Likelihood estimator 
#' @description 
#' Define group rank negative log-likelihood function from
#' Huang et al., 2008, p.10, Eq.31
#' @export
ll_gr_helper <- function( vs, res_matches ) {
    Tp <- ( ( res_matches[,-c(1,2)] > 0 ) %*% vs )
    Tm <- ( ( res_matches[,-c(1,2)] < 0 ) %*% vs )
    LL <- exp( Tp + Tm - ( res_matches$n.pos - res_matches$n.neg ) ) / 
      ( exp( Tp - ( res_matches$n.pos - res_matches$n.neg ) ) + exp( Tm ) )^2
    ret <- -sum( log( LL ) )
    return( ret )
  } # function (END)


#' @title  Define group rank gradient for negative log-likelihood function 
#' @examples 
#' require(GameRank)
#' require(numDeriv)
#' # Some example match data frame
#' res_matches <- as.data.frame( matrix( c(
#' 1,4, +1, +1, -1, -1,  0,  0,  0,
#' 3,2, +1, +1,  0, -1, -1,  0,  0,
#' 0,5, +1, -1, -1, +1,  0,  0,  0,
#' 5,0, +1, +1,  0,  0,  0, -1, -1,
#' 2,3, -1, +1,  0, +1,  0,  0, -1
#' ),
#' nrow = 5, ncol = 9, byrow = TRUE,
#' dimnames = list( NULL, c("n.pos","n.neg", sprintf( "f%d", 1:7 ) ) ) ) )
#' res_matches
#' 
#' vs <- c( 0.1, 0.2, 0.3, -0.2, -0.1, 0.2, -0.2 )
#' vs <- runif( 7 )
#' names(vs) <- colnames(res_matches)[-c(1,2)]
#' 
#' ll_gr_helper( vs, res_matches )
#' 
#' numericDeriv( quote( ll_gr_helper( vs, res_matches ) ), theta = "vs", rho = as.environment( list( vs=vs, res_matches = res_matches, ll_gr_helper = ll_gr_helper  ) ) )
#' numDeriv::grad( func = function(vs) ll_gr_helper( vs, res_matches ), x = vs )
#' ll_gr_grad_helper( vs, res_matches )
#' 
#' @export
ll_gr_grad_helper <- function( vs, res_matches ) {
    Tp <- ( ( res_matches[,-c(1,2)] > 0 ) %*% vs )
    Tm <- ( ( res_matches[,-c(1,2)] < 0 ) %*% vs )
    
    pp <- exp( Tp + res_matches[,2] )/ 
      ( exp( Tp + res_matches[,2] ) + exp( Tm + res_matches[,1] ) )
    # CHE/2024-09-09: Spotted error in derivative. This line must use Tm and 
    # row 33 must use pn. Also extracted log-likelihood and its gradient into
    # helper file.
    pn <- exp( Tm + res_matches[,1] )/ 
      ( exp( Tp + res_matches[,2] ) + exp( Tm + res_matches[,1] ) )
    gr <- vapply( names(vs), FUN=function(co) {
      ms <- sum( abs( res_matches[,co] ) )
      pps <- sum( pp[which(res_matches[,co] > 0)] )
      pns <- sum( pn[which(res_matches[,co] < 0)] )
      return( -ms + 2 * (pps + pns) )
    }, 1.0)
    names( gr ) <- names( vs )
    return( gr )
  } # function (END)

