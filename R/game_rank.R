#
# File including the GameRank algorithm
# Author: Carsten Henneges, hennegc1@gene.com
#
# GameRank is a maximum-likelihood based wrapper algorithm for feature selection.
# It is based on the algorithm developed by Huang et al., JMLR, 2008 to rank
# individuals by group comparisons.
#


#' @title GameRank algorithm
#' 
#' @description GameRank is an two phase algorithm that first builds a comparison of feature set pairs and evaluates them
#' on repeated 50:50 random training:validation splits. In the second phase it uses those evaluations to estimate a maximum
#' likelihood ranking model for feature combinations.
#' 
#' @details The Game Rank algorithm runs as follows:
#' \code{ \cr
#' 1. Build a sampling matrix of disjoint and unique feature combinations
#' 2. For each combination evaluate the model per training:validation splits.
#' 3. Return the best selection.
#' }
#' 
#' @param fo Only for call with formula as first argument. Extracts lhs ~ rhs into resp and vars, and calls backward( dat, resp, vars, ... )
#' @param dat data.frame or tibble comprising data for model generation and validation.
#' @param resp Character string defining the response (lhs) of the model formula.
#' @param vars Character vector defining the list of variables for selection. Those are concatenated by '+' 
#' as the right hand side (rhs) of the modelling formula.
#' @param fn_train Function with signature function( dat, resp, selection, ... ) that returns a model or NULL in any other case on the given data dat.
#' @param fn_eval Function with signature function( dat, resp, selection, ... ) that returns a real number or NA in any other case, e.g. when model is NULL.
#' @param m Size of final partition size. 
#' @param team_size Selection sizes that are evaluated during match phase against each other
#' @param rounds Number of rounds each pair of selections is evaluated on randomly selected train/eval sets
#' @param min_matches_per_var Minimum number of matches each variable needs to be part of before the match phase can end
#' @param opt_method Should be either 'BFGS' or 'CG', the latter in case Hessian-based optimization is infeasible for the match data
#' @param max_iter Maximum number of optimization iterations for group rank model fit, passed to optim( ... ).
#' @param ... An other arguments passed to fn_train or fn_eval during calls, e.g. maybe 'u = 365' for Survival evaluations specifying the landmark day.
#'
#' @return List with elements
#' \describe{
#'  \item{response}{As from input parameters}
#'  \item{variables}{As from input parameters}
#'  \item{m}{As from input parameters}
#'  \item{maximize}{As from input parameters}
#'  \item{team_size}{As from input parameters}
#'  \item{rounds}{As from input parameters}
#'  \item{min_matches_per_var}{As from input parameters}
#'  \item{opt_method}{As from input parameters}
#'  \item{max_iter}{As from input parameters}
#'  \item{start}{Start time of core algorithm loop}
#'  \item{end}{End time of core algorithm loop}
#'  \item{match_matrix_time}{Time after start when the match matrix has been sampled.}
#'  \item{match_played_time}{Time after start when the match comparisons have been completed.}
#'  \item{fit_time}{Time after start when the maximum-likelihood model has been fit, excluding time to compute the score vector and Hessian matrix.}
#'  \item{match_results}{The match matrix with scoring results that was used to determine the feature selection comparisons.}
#'  \item{variable_ranking}{Tibble with variable ranking, including standard errors}
#'  \item{game_rank_selection}{Best m variables as final selection per maximum likelihood estimate.}

#'  \item{optimization_result}{Result of optimization call.}
#'  \item{solution}{Optimization solution vector.}
#'  \item{score_vector}{Gradient at optimal solution (score vector)}
#'  \item{inv_hessian}{Inverse Hessian matrix at optimal solution. Can be used for Delta method.}
#' }
#' @name game_rank
NULL

#
# Algorithm to generate a match matrix that enumerates random feature selection pairs
# with predefined size and ensures that each feature is part of a minimum number of evaluations.
#
# The match matrix bears +1 for the (+) team and -1 for the (-) team, and 0 for all others.
# An index vector into the list of variables is created and then chunked in sizes of team_size
# to alternatingly define (+) and (-) teams being added to the match matrix until it is empty and 
# then is refilled. The process continues until each feature is evaluated min_matches_per_var times.
#
build_match_matrix <- function( sel, team_size, min_matches_per_var ) {
  
  sel[1:length(sel)] <- 0 # Set all elements to 0
  match_matrix <- matrix( NA, nrow=0, ncol = length(sel) )
  colnames(match_matrix) <- names(sel)
  while( is.null(match_matrix) | !all( min_matches_per_var < colSums( abs( match_matrix ) )  ) ) {
    idx <- 1:length(sel)
    idx <- idx[ order( runif( length(idx) ) ) ]
    while( 2*team_size < length(idx)  ) {
      idxp <- idx[1:team_size]
      idx <- idx[-c(1:team_size)]
      idxn <- idx[1:team_size]
      idx <- idx[-c(1:team_size)]
      Tpm <- sel
      Tpm[idxp] <- +1L
      Tpm[idxn] <- -1L
      match_matrix <- rbind( match_matrix, Tpm )
    }
  }

  return( match_matrix )  
}

# GameRank ----
#' @rdname game_rank
#' @export
game_rank <- function( dat,
                       resp,
                       vars,
                       fn_train = fn_train_binomial,
                       fn_eval  = fn_eval_binomial,
                       m = NULL,
                       maximize = TRUE,
                       team_size = 5L,
                       rounds = 50L,
                       min_matches_per_var = 10L,
                       opt_method = "BFGS",
                       max_iter = 1E8L, 
                       ...
                       )
{
  # Check inputs ----
  stopifnot( is.data.frame(dat) | is_tibble(dat) )
  stopifnot( is.character(vars) & 1 < length(vars) )
  stopifnot( is.character(vars) )
  
  stopifnot( is.function(fn_train) ) 
  stopifnot( is.function(fn_eval) ) 
  
  if( is.null(m) ) { stop( "Please provide number m of features to select.\n" ) }
  stopifnot( is.integer(team_size) )
  stopifnot( is.integer(rounds) ) 
  stopifnot( is.integer(min_matches_per_var) )
  stopifnot( opt_method %in% c("BFGS","CG") )
  stopifnot( is.integer(max_iter) )
  stopifnot( is.logical(maximize) )
  
  start_time <- Sys.time()
  # Initialize variables ----
  sel <- rep_len( 0.0, length.out = length(vars) )
  names(sel) <- vars
  # Initialize team vector, set first to +1, second to -1 and rest to 0
  # Tpm <- sel
  # Tpm[1:team_size] <- +1L
  # Tpm[(team_size+1):(2*team_size)] <- -1L
  MM <- build_match_matrix( sel = sel, team_size = team_size, min_matches_per_var = min_matches_per_var )
  match_matrix_time <- Sys.time()
  
  # Initialize row selection vector: 1 = development, 2 = evaluation in a round
  ds <- rep_len( c(1,2), length.out = nrow(dat) )
  
  # Define team evaluation function (internal) ----
  eval_team <- function( dat, resp, selection, ds, ... ) {
    ret <- tryCatch( {
      mod <- fn_train( dat[which(1==ds),], resp, selection, ... )   # Obtain model from data in 1-fold
      evl <- fn_eval(  dat[which(2==ds),], resp, selection, mod, ... ) # Evaluate model on data in 2-fold
      evl
    }, error = function(e) NA )
    return( ret )
  }
  
  # Run GameRank matches ----
  nmm <- nrow(MM)
  cat( sprintf( "Comparing variable selections (# matches %d)--- \n", nmm ) )
  t <- 1
  res_matches <- NULL
  while( t <= nmm ) {
    # Tpm[ 1:length(Tpm) ] <- Tpm[ order( runif( length( Tpm ) ) ) ] # Shuffle teams
    Tpm <- MM[t,]
    
    fop <- names( which( Tpm > 0 ) )
    fon <- names( which( Tpm < 0 ) )
    
    np <- 0L
    nn <- 0L
    for( r in 1:rounds ) {
      ds[ 1:length(ds) ] <- ds[ order( runif( length( ds ) ) ) ] # Shuffle rows to development and validation
      
      scop <- eval_team( dat, resp, fop, ds, ... )
      scon <- eval_team( dat, resp, fon, ds, ... )
      
      # Compare results
      if( !is.na(scop) & !is.na(scon) ) {
        # Note: larger scores are better!
        if( maximize ) {
          # Maximize
          if(  scop < scon ) {
            nn <- nn + 1 # (-) team wins
          } else if( scon < scop ) {
            np <- np + 1 # (+) team wins
          }
        } else {
          # Minimize
          if(  scop > scon ) {
            nn <- nn + 1 # (-) team wins
          } else if( scon > scop ) {
            np <- np + 1 # (+) team wins
          }
        }
        
      } else if( !is.na(scop) & is.na(scon) ) {
        np <- np + 1 # (+) team wins, as it only produces results
      } else if( !is.na(scon) & is.na(scop) ) {
        nn <- nn + 1 # (-) team wins, as it only produces results
      } # if
    } # for (END) 
    
    res_match <- cbind( data.frame( n.pos = np, n.neg = nn), t( Tpm ) )
    res_matches <- bind_rows( res_matches, res_match )
    
    cat( sprintf( "Iteration %4d of %4d -- (+) : (-) scored %4d : %4d \n", t, nmm, np, nn ) )
    t <- t + 1 # Iteration counter
  } # while (END)
  match_played_time <- Sys.time()
  
  # Evaluating Group Rank Maximum Likelihood estimator 
  # Define group rank negative log-likelihood function ----
  # Huang et al., 2008, p.10, Eq.31
  ll_gr <- local({
    function( vs ) {
      Tp <- ( ( res_matches[,-c(1:2)] > 0 ) %*% vs )
      Tm <- ( ( res_matches[,-c(1:2)] < 0 ) %*% vs )
      LL <- exp( Tp + Tm - ( res_matches$n.pos - res_matches$n.neg ) ) / ( exp( Tp - ( res_matches$n.pos - res_matches$n.neg ) ) + exp( Tm ) )^2
      ret <- -sum( log( LL ) )
      return( ret )
    } # function (END)
  }, envir = new.env())
  # Define group rank gradient for negative log-likelihood function ----
  ll_gr_grad <- local({
    function( vs ) {
      Tp <- ( ( res_matches[,-c(1:2)] > 0 ) %*% vs )
      Tm <- ( ( res_matches[,-c(1:2)] < 0 ) %*% vs )
      
      pp <- exp( Tp + res_matches[,2] )/ ( exp( Tp + res_matches[,2] ) + exp( Tm + res_matches[,1] ) )
      pn <- exp( Tp + res_matches[,1] )/ ( exp( Tp + res_matches[,2] ) + exp( Tm + res_matches[,1] ) )
      gr <- sapply( names(vs), FUN=function(co) {
        ms <- sum( abs( res_matches[,co] ) )
        pps <- sum( pp[which(res_matches[,co] > 0)] )
        pns <- sum( pp[which(res_matches[,co] < 0)] )
        return( -ms + 2 * (pps + pns) )
      })
      names( gr ) <- names( vs )
      return( gr )
    } # function (END)
  }, envir = new.env() )
  
  # Fitting Group Rank model ----
  cat( "Optimizing maximum likelihood \n" )
  oo <- optim( par = sel, fn = ll_gr, gr = ll_gr_grad, method = opt_method, control = list( fnscale = +1L, maxit = max_iter ) )
  oo
  fit_time <- Sys.time()
  
  cat( "Calculating score vector \n" )
  gg <- ll_gr_grad( oo$par )  
  gg
  
  cat( "Calculating Hessian matrix \n" )
  hh <- jacobian( func = ll_gr_grad, x = oo$par )
  hh
  
  vv <- chol2inv( hh )
  rownames(vv) <- colnames(vv) <- colnames(hh)
  vv
  end_time <- Sys.time()
  
  # Compiling results  ----
  cat( "Compiling results \n" )
  vsel_result <- tibble( variable = names( oo$par ),
                         vs = as.numeric( oo$par),
                         vs.var = diag( vv ) ) %>%
    mutate( selected = (vs > 0) ) %>%
    arrange( desc( vs ), vs.var )
  
  var_selection <- vsel_result$variable[1:m]
  
  # Return list object with parameters and results
  ret <- list(
    # Input data
    # data = dat,
    response = resp,
    variables = vars,
    m = m,
    
    # Input parameters
    team_size = team_size,
    rounds = rounds,
    min_matches_per_var = min_matches_per_var,
    opt_method = opt_method,
    max_iter = max_iter,
    maximize = maximize,
    
    
    # Time
    start = start_time, end = end_time,
    match_matrix_time = match_matrix_time, match_played_time = match_played_time, fit_time = fit_time,
    
    # Results
    match_results = res_matches,
    variable_ranking = vsel_result,
    game_rank_selection = var_selection,
    
    optimization_result = oo,
    solution = oo$par,
    score_vector = gg,
    inv_hessian = hh
  )
  class( ret ) <- c( "GameRank", class(ret) )
  return( ret )
} # game_rank (END)

# GameRank formula interface ----
#' @rdname game_rank
#' @export
game_rank.formula <- function( fo, dat, 
                               fn_train = fn_train_binomial,
                               fn_eval  = fn_eval_binomial,
                               m = NULL,
                               maximize = TRUE,
                               team_size = 5L,
                               rounds = 50L,
                               min_matches_per_var = 10L,
                               opt_method = "BFGS",
                               max_iter = 1E8L, 
                               ... ) 
{
  # Check inputs
  stopifnot( is.formula( fo ) ) 
  stopifnot( is.data.frame(dat) | is_tibble(dat) )
  stopifnot( is.function(fn_train) )
  stopifnot( is.function(fn_eval) )
  
  resp <- as.character(lhs.vars( fo ) )
  vars <- as.character(rhs.vars( fo ) )
  
  game_rank( dat = dat,
             resp = resp,
             vars = vars,
             fn_train = fn_train,
             fn_eval  = fn_eval,
             m = m,
             maximize = maximize,
             team_size = team_size,
             rounds = rounds,
             min_matches_per_var = min_matches_per_var,
             opt_method = opt_method,
             max_iter = max_iter,
             ... )
} # game_rank.formula (END)


#
# Given a Match matrix, we want to estimate (differential) team score and their
# standard errors using the Delta method. We also construct (1-alpha)% confidence
# intervals using Normal approximation.
#
# TODO Debug and check!
#
estimate_T_matches <- function( Tm, vs, HH, alpha = 0.05 ) {
  stopifnot( is.matrix(Tm) & 0==length( setdiff( unique(as.numeric(Tm)), c(-1L,0L,+1L) ) ) )
  stopifnot( is.numeric(vs) )
  stopifnot( (length(vs)==nrow(HH)) & (length(vs)==ncol(HH)) )
  
  # TODO: Check this implementation for correctness!
  Tscore <- function( m ) { as.numeric( m %*% vs ) } 
  se_Tscore <- function( m ) {
    gg <- grad( func = Tscore, x = m )
    # hh <- hessian( func = Tscore, x = m )
    as.numeric( gg %*% HH %*% gg )
  } 
  
  zz <- qnorm( 1 - alpha / 2.0 )
  rr <- matrix( NA, ncol = 4, nrow = nrow(Tm) )
  colnames(rr) <- c("dT","dT.se","dT.LCL","dT.UCL")
  for( i in 1:nrow(rr) ) {
    rr[i,"dT"] <- Tscore( Tm[i,] )
    rr[i,"dT.se"] <- se_Tscore( Tm[i,] )
  }
  rr[,"dT.LCL"] <- rr[,"dT"] - zz * rr[,"dT.se"] 
  rr[,"dT.UCL"] <- rr[,"dT"] + zz * rr[,"dT.se"]
  return( rr )  
}

estimate_T_matches.GameRank <- function( gmr_result, Tm ) {
  estimate_T_matches( Tm = Tm, vs = gmr_result$solution ) 
}

