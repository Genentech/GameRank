#
# File including the GameRank algorithm
# Author: Carsten Henneges, hennegc1@gene.com
#
# GameRank is a maximum-likelihood based wrapper algorithm for feature selection.
# It is based on the algorithm developed by Huang et al., JMLR, 2008 to rank
# individuals by group comparisons.
#

#' @import tibble dplyr 
#' @importFrom formula.tools lhs.vars rhs.vars
#' @import numDeriv
#' @importFrom rlang .data

#' @title GameRank algorithm
#' 
#' @description GameRank is an two phase algorithm that first builds a 
#' comparison of feature set pairs and evaluates them
#' on repeated 50:50 random training:validation splits. In the second phase it 
#' uses those evaluations to estimate a maximum likelihood ranking model for 
#' feature combinations.
#' 
#' @details The Game Rank algorithm runs as follows:
#' \preformatted{ 
#' 1. Build a sampling matrix of disjoint and unique feature combinations 
#'    for two teams
#' 2. For each combination evaluate the model per training:validation splits
#'    and recorde a score for the winning team
#' 3. Estimate the feature contribution to team performance using the group
#'    ranking model from Huang et al. (2008, JMLR)
#' 3. Return the feature ranking and best selection
#' }
#' 
#' @param fo Only for call with formula as first argument. Extracts lhs ~ rhs 
#' into resp and vars, and calls backward( dat, resp, vars, ... )
#' @param dat data.frame or tibble comprising data for model generation and 
#' validation.
#' @param resp Character string defining the response (lhs) of the 
#' model formula.
#' @param vars Character vector defining the list of variables for selection. 
#' Those are concatenated by '+' as the right hand side (rhs) of the 
#' modeling formula.
#' @param fn_train Function with signature function( dat, resp, selection, ... )
#' that returns a model or NULL in any other case on the given data dat.
#' @param fn_eval Function with signature function( dat, resp, selection, ... ) 
#' that returns a real number or NA in any other case, e.g. when model is NULL.
#' @param m Size of final partition size. 
#' @param dsi Integer array of 1s and 2s to determine the proportion of training 
#' vs validation splits during match rounds. Default is c(1,2) to 
#' use 50:50 splits.  Can be set to c(2,2) if only validation should be used,
#' e.g. in the case of bootstrapping the validation performance.
#' @param maximize A logic value determining if fn_eval is maximized (set to 
#' TRUE) or minimized (set to FALSE).
#' @param team_size Selection sizes that are evaluated during match phase 
#' against each other
#' @param rounds Number of rounds each pair of selections is evaluated on 
#' randomly selected train/eval sets
#' @param min_matches_per_var Minimum number of matches each variable needs to 
#' be part of before the match phase can end
#' @param opt_method Should be either 'BFGS' or 'CG', the latter in case 
#' Hessian-based optimization is infeasible for the match data
#' @param max_iter Maximum number of optimization iterations for group rank 
#' model fit, passed to optim( ... ).
#' @param ... An other arguments passed to fn_train or fn_eval during calls, 
#' e.g. maybe 'u = 365' for Survival evaluations specifying the landmark day.
#'
#' @return List with elements
#' \describe{
#'  \item{response}{As from input parameters}
#'  \item{variables}{As from input parameters}
#'  \item{m}{As from input parameters}
#'  \item{dsi}{As from input parameters}
#'  \item{maximize}{As from input parameters}
#'  \item{team_size}{As from input parameters}
#'  \item{rounds}{As from input parameters}
#'  \item{min_matches_per_var}{As from input parameters}
#'  \item{opt_method}{As from input parameters}
#'  \item{max_iter}{As from input parameters}
#'  \item{start}{Start time of core algorithm loop}
#'  \item{end}{End time of core algorithm loop}
#'  \item{match_matrix_time}{Time after start when the match matrix has 
#'  been sampled.}
#'  \item{match_played_time}{Time after start when the match comparisons have 
#'  been completed.}
#'  \item{fit_time}{Time after start when the maximum-likelihood model has been 
#'  fit, excluding time to compute the score vector and Hessian matrix.}
#'  \item{match_results}{The match matrix with scoring results that was used to
#'  determine the feature selection comparisons.}
#'  \item{variable_ranking}{Tibble with variable ranking, including 
#'  standard errors}
#'  \item{game_rank_selection}{Best m variables as final selection per maximum 
#'  likelihood estimate.}

#'  \item{optimization_result}{Result of optimization call.}
#'  \item{solution}{Optimization solution vector.}
#'  \item{score_vector}{Gradient at optimal solution (score vector)}
#'  \item{inv_hessian}{Inverse Hessian matrix at optimal solution. Can be used 
#'  for Delta method.}
#' }
#' 
#' @name game_rank
NULL

#
# Algorithm to generate a match matrix that enumerates random feature selection
# pairs with predefined size and ensures that each feature is part of a minimum
# number of evaluations.
#
# The match matrix bears +1 for the (+) team and -1 for the (-) team, and 0 for 
# all others. An index vector into the list of variables is created and then
# chunked in sizes of team_size to alternatingly define (+) and (-) teams being 
# added to the match matrix until it is empty and then is refilled. The process 
# continues until each feature is evaluated min_matches_per_var times.
#
build_match_matrix <- function( sel, team_size, min_matches_per_var ) {
  
  sel[seq_along(sel)] <- 0 # Set all elements to 0
  match_matrix <- matrix( NA, nrow=0, ncol = length(sel) )
  colnames(match_matrix) <- names(sel)
  while( is.null(match_matrix) | 
         !all( min_matches_per_var < colSums( abs( match_matrix ) )  ) ) {
    idx <- seq_along(sel)
    idx <- idx[ order( stats::runif( length(idx) ) ) ]
    while( 2*team_size < length(idx)  ) {
      idxp <- idx[seq_len(team_size)]
      idx <- idx[-seq_len(team_size)]
      idxn <- idx[seq_len(team_size)]
      idx <- idx[-seq_len(team_size)]
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
#' @examples 
#' vars <- grep( "the_|rnd", colnames(toy_data), value=TRUE )
#' resp <- "resp"
#' res <- game_rank( toy_data, resp, vars, 
#'                   fn_train_binomial, fn_eval_binomial_auroc, 
#'                   4L, c(1L,2L), TRUE, 3L, 7L, 3L )
#' res$game_rank_selection
#' res$variable_ranking
#' @export
game_rank <- function( dat,
                       resp,
                       vars,
                       fn_train = fn_train_binomial,
                       fn_eval  = fn_eval_binomial,
                       m = NULL,
                       dsi = c(1L,2L),
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
  
  if(is.null(m)) { stop( "Please provide number m of features to select.\n" ) }
  stopifnot( is.vector(dsi) & is.integer(dsi) & all( dsi %in% c(1,2) ) )
  stopifnot( is.integer(team_size) )
  stopifnot( is.integer(rounds) ) 
  stopifnot( is.integer(min_matches_per_var) )
  stopifnot( opt_method %in% c("BFGS","CG") )
  stopifnot( is.integer(max_iter) )
  stopifnot( is.logical(maximize) )
  
  message( sprintf("game_rank: Starting GameRank selection algorithm for %d observations with %d features (m=%d, maximize=%d).", nrow(dat), length(vars), m, maximize ))
  message( sprintf("game_rank: Additional GameRank parameters team_size=%d, rounds=%d, min matches per var = %d.", team_size, rounds, min_matches_per_var ))
  start_time <- Sys.time()
  # Initialize variables ----
  sel <- rep_len( 0.0, length.out = length(vars) )
  names(sel) <- vars
  # Initialize team vector, set first to +1, second to -1 and rest to 0
  # Tpm <- sel
  # Tpm[1:team_size] <- +1L
  # Tpm[(team_size+1):(2*team_size)] <- -1L
  MM <- build_match_matrix( sel = sel, team_size = team_size, 
                            min_matches_per_var = min_matches_per_var )
  match_matrix_time <- Sys.time()
  
  # Initialize row selection vector: 1 = development, 2 = evaluation in a round
  ds <- rep_len( dsi, length.out = nrow(dat) )
  
  # Define team evaluation function (internal) ----
  eval_team <- function( dat, resp, selection, ds, ... ) {
    ret <- tryCatch( {
      # Note: Careful! R argument matching gives exact matches priority, then followed by substring-matching, ie.
      #       selection=..., and se=1E-4 in the ... will confuse, since finally positional matching is applied.
      #       CHE/2023-07-13.
      # Obtain model from data in 1-fold
      mod <- fn_train( dat=dat[which(1==ds),], resp=resp, selection=selection, ... )   
      # Evaluate model on data in 2-fold
      evl <- fn_eval(  dat=dat[which(2==ds),], resp=resp, selection=selection, mod=mod, ... ) 
      evl
    }, error = function(e) NA )
    return( ret )
  }
  
  # Run GameRank matches ----
  nmm <- nrow(MM)
  message( sprintf( "Comparing variable selections (# matches %d)--- ", nmm ) )
  t <- 1
  res_matches <- NULL
  while( t <= nmm ) {
    # Shuffled teams
    Tpm <- MM[t,]
    
    fop <- names( which( Tpm > 0 ) )
    fon <- names( which( Tpm < 0 ) )
    
    np <- 0L
    nn <- 0L
    for( r in seq_len(rounds) ) {
      # Shuffle rows to development and validation
      ds[ seq_along(ds) ] <- ds[ order( stats::runif( length( ds ) ) ) ] 
      
      scop <- eval_team( dat=dat, resp=resp, selection=fop, ds=ds, ... )
      scon <- eval_team( dat=dat, resp=resp, selection=fon, ds=ds, ... )
      
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
    
    message( sprintf( "Iteration %4d of %4d -- (+) : (-) scored %4d : %4d ", 
                      t, nmm, np, nn ) )
    t <- t + 1 # Iteration counter
  } # while (END)
  match_played_time <- Sys.time()
  
  # Evaluating Group Rank Maximum Likelihood estimator 
  # Define group rank negative log-likelihood function ----
  # Huang et al., 2008, p.10, Eq.31
  ll_gr <- local({
    function( vs ) {
      Tp <- ( ( res_matches[,-c(1,2)] > 0 ) %*% vs )
      Tm <- ( ( res_matches[,-c(1,2)] < 0 ) %*% vs )
      LL <- exp( Tp + Tm - ( res_matches$n.pos - res_matches$n.neg ) ) / 
        ( exp( Tp - ( res_matches$n.pos - res_matches$n.neg ) ) + exp( Tm ) )^2
      ret <- -sum( log( LL ) )
      return( ret )
    } # function (END)
  }, envir = new.env())
  # Define group rank gradient for negative log-likelihood function ----
  ll_gr_grad <- local({
    function( vs ) {
      Tp <- ( ( res_matches[,-c(1,2)] > 0 ) %*% vs )
      Tm <- ( ( res_matches[,-c(1,2)] < 0 ) %*% vs )
      
      pp <- exp( Tp + res_matches[,2] )/ 
        ( exp( Tp + res_matches[,2] ) + exp( Tm + res_matches[,1] ) )
      pn <- exp( Tp + res_matches[,1] )/ 
        ( exp( Tp + res_matches[,2] ) + exp( Tm + res_matches[,1] ) )
      gr <- vapply( names(vs), FUN=function(co) {
        ms <- sum( abs( res_matches[,co] ) )
        pps <- sum( pp[which(res_matches[,co] > 0)] )
        pns <- sum( pp[which(res_matches[,co] < 0)] )
        return( -ms + 2 * (pps + pns) )
      }, 1.0)
      names( gr ) <- names( vs )
      return( gr )
    } # function (END)
  }, envir = new.env() )
  
  # Fitting Group Rank model ----
  message( "Optimizing maximum likelihood " )
  oo <- stats::optim( par = sel, fn = ll_gr, gr = ll_gr_grad, 
                      method = opt_method, 
                      control = list( fnscale = +1L, maxit = max_iter ) )
  oo
  fit_time <- Sys.time()
  
  message( "Calculating score vector " )
  gg <- ll_gr_grad( oo$par )  
  gg
  
  message( "Calculating Hessian matrix " )
  hh <- NULL
  hh <- tryCatch({numDeriv::jacobian( func = ll_gr_grad, x = oo$par )}, error = function(ee) NULL )
  hh
  
  vv <- NULL
  vv <- tryCatch({
    vv <- solve( hh )
    rownames(vv) <- colnames(vv) <- colnames(hh)
    vv
  }, error = function(ee) NULL)
  vv
  end_time <- Sys.time()
  
  # Compiling results  ----
  message( "Compiling results " )
  vsel_result <- tibble( variable = names( oo$par ),
                         vs = as.numeric( oo$par) ) %>%
    mutate( selected = (.data$vs > 0) ) %>%
    arrange( desc( .data$vs ) )
  
  vsel_result <- tryCatch({
    tmp <- vsel_result %>% mutate( vs.var = diag( vv ) )
    tmp <- tmp %>% 
      arrange( desc( .data$vs ), .data$vs.var )
    tmp
  }, error = function(ee) vsel_result )
  
  var_selection <- vsel_result$variable[seq_len(m)]
  
  # Return list object with parameters and results
  ret <- list(
    algorithm = "game_rank",
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
    match_matrix_time = match_matrix_time, 
    match_played_time = match_played_time, fit_time = fit_time,
    
    # Results
    match_results = res_matches,
    variable_ranking = vsel_result,
    game_rank_selection = var_selection,
    
    optimization_result = oo,
    solution = oo$par,
    score_vector = gg,
    inv_hessian = hh,
    hessian = vv
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
                               dsi = c(1L,2L),
                               maximize = TRUE,
                               team_size = 5L,
                               rounds = 50L,
                               min_matches_per_var = 10L,
                               opt_method = "BFGS",
                               max_iter = 1E8L, 
                               ... ) 
{
  # Check inputs
  stopifnot( "formula"==class( fo ) ) 
  stopifnot( is.data.frame(dat) | is_tibble(dat) )
  stopifnot( is.function(fn_train) )
  stopifnot( is.function(fn_eval) )
  
  resp <- as.character(formula.tools::lhs.vars( fo ) )
  vars <- as.character(formula.tools::rhs.vars( fo ) )
  
  game_rank( dat = dat,
             resp = resp,
             vars = vars,
             fn_train = fn_train,
             fn_eval  = fn_eval,
             m = m,
             dsi = dsi,
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
estimate_T_matches <- function( Tm, vs, HH, alpha = 0.05 ) {
  stopifnot( is.matrix(Tm) & 
               0==length( setdiff( unique(as.numeric(Tm)), c(-1L,0L,+1L) ) ) )
  stopifnot( is.numeric(vs) )
  stopifnot( (length(vs)==nrow(HH)) & (length(vs)==ncol(HH)) )
  
  # TODO: Check this implementation for correctness!
  Tscore <- function( m ) { as.numeric( m %*% vs ) } 
  se_Tscore <- function( m ) {
    gg <-   numDeriv::grad( func = Tscore, x = m )
    # hh <- hessian( func = Tscore, x = m )
    as.numeric( gg %*% HH %*% gg )
  } 
  
  zz <- stats::qnorm( 1 - alpha / 2.0 )
  rr <- matrix( NA, ncol = 4, nrow = nrow(Tm) )
  colnames(rr) <- c("dT","dT.se","dT.LCL","dT.UCL")
  for( i in seq_len( nrow(rr) ) ) {
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

