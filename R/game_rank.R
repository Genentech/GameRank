#
# File including the GameRank algorithm
# Author: Carsten Henneges, hennegc1@gene.com
#
# GameRank is a maximum-likelihood based wrapper algorithm for feature selection.
# It is based on the algorithm developed by Huang et al., JMLR, 2008 to rank
# individuals by group comparisons.
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
#' @title GameRank - A maximum-likelihood based wrapper for variable selection
#' 
#' @param dat Is either a data.frame or tibble with the full dataset used for selection
#' @param resp Character string defining response variable (e.g. either a variable name or an formula expression, like 'Surv(time,censored)', used by fn_train and fn_eval)
#' @param vars Character strings with variables or variable expressions. Subselections will be fed into fn_train and fn_eval functions.
#' @param team_size Selection sizes that are evaluated during match phase against each other
#' @param rounds Number of rounds each pair of selections is evaluated on randomly selected train/eval sets
#' @param min_matches_per_var Minimum number of matches each variable needs to be part of before the match phase can end
#' @param opt_method Should be either 'BFGS' or 'CG', the latter in case Hessian-based optimization is infeasible for the match data
#' @param max_iter Maximum number of rounds of the match phase and for group rank model fit.
#' 
#' @returns List with selection results [TBD]
game_rank <- function( dat,
                       resp,
                       vars,
                       fn_train = fn_train_binomial,
                       fn_eval  = fn_eval_binomial,
                       m = NULL,
                       team_size = 5L,
                       rounds = 50L,
                       min_matches_per_var = 10L,
                       opt_method = "BFGS",
                       max_iter = as.integer(1E8), 
                       maximize = TRUE,
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
  
  # Initialize variables ----
  total_start_time = Sys.time()
  sel <- rep_len( 0.0, length.out = length(vars) )
  names(sel) <- vars
  # Initialize team vector, set first to +1, second to -1 and rest to 0
  # Tpm <- sel
  # Tpm[1:team_size] <- +1L
  # Tpm[(team_size+1):(2*team_size)] <- -1L
  MM <- build_match_matrix( sel = sel, team_size = team_size, min_matches_per_var = min_matches_per_var )
  
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
  cat( sprintf( "Comparing variable selections (# matches %d)--- \n", nrow(MM) ) )
  t <- 1
  res_matches <- NULL
  time_start <- Sys.time()
  # while( n_matches < max_iter && (is.null(res_matches) || !all( min_matches_per_var < colSums( abs( res_matches[,-c(1,2)] ) ) ) ) ) {
  while( t <= nrow(MM) ) {
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
    
    cat( sprintf( "Iteration %d -- (+) : (-) scored %d : %d \n", t, np, nn ) )
    t <- t + 1 # Iteration counter
  } # while (END)
  time_stop <- Sys.time()
  match_time <- difftime( time1 = time_stop, time2 = time_start, units = "secs" )
  cat( sprintf( "Match phase took %1.4s \n", as.double(match_time) ))
  
  # Evaluating Group Rank Maximumn Likelihood estimator 
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
  time_start <- Sys.time()
  oo <- optim( par = sel, fn = ll_gr, gr = ll_gr_grad, method = opt_method, control = list( fnscale = +1L, maxit = max_iter ) )
  oo
  
  cat( "Calculating score vector \n" )
  gg <- ll_gr_grad( oo$par )  
  gg
  
  cat( "Calculating Hessian matrix \n" )
  hh <- jacobian( func = ll_gr_grad, x = oo$par )
  hh
  
  vv <- chol2inv( hh )
  rownames(vv) <- colnames(vv) <- colnames(hh)
  vv
  time_stop <- Sys.time()
  grprnk_time <- difftime( time1 = time_stop, time2 = time_start, units = "secs" )
  cat( sprintf( "Fitting GroupRank model took %1.4s \n", as.double(grprnk_time) ))
  
  total_stop_time = Sys.time()
  total_time <- difftime( time1 = total_stop_time, time2 = total_start_time, units = "secs" )
  cat( sprintf( "Overall selection process took %1.4s \n", as.double(total_time) ))
  
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
    
    # Timings 
    match_time = match_time,
    grprnk_time = grprnk_time,
    total_time = total_time,
    
    # Results
    match_matrix = MM,
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
game_rank.formula <- function( fo, dat, 
                               fn_train = fn_train_binomial,
                               fn_eval  = fn_eval_binomial,
                               m = NULL,
                               team_size = 5L,
                               rounds = 50L,
                               min_matches_per_var = 10L,
                               opt_method = "BFGS",
                               max_iter = as.integer(1E8), 
                               maximize = TRUE,
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
             team_size = team_size,
             rounds = rounds,
             min_matches_per_var = min_matches_per_var,
             opt_method = opt_method,
             max_iter = max_iter,
             maximize = maximize,
             ... )
} # game_rank.formula (END)


#
# Given a Match matrix, we want to estimate (differential) team score and their
# standard errors using the Delta method. We also construct (1-alpha)% confidence
# intervals using Normal approximation.
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

