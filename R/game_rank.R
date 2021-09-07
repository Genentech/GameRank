
library( formula.tools )

# Formula interface
game_rank_formula <- function( dat, # data.frame or tibble
                       fo,  # formula of response ~ all variables for selection
                       # Training function
                       fn_train = function( fo, dat, ... ) { 
                         glm( fo, "binomial", dat ) 
                       },
                       # Evaluation function
                       fn_eval  = function( fo, dat, mod, ... ) { 
                         prfp <- performance( 
                           prediction( predictions = predict( mod, newdata = dat, type="response"), 
                                       labels = model.frame(fo,dat)[,1] ),
                           "auc" )
                         auc <- as.numeric( prfp@y.values )
                         return( auc )
                       },
                       tsize  = 5,  # Team Size
                       rounds = 25, # Rounds to evaluate per match
                       min_matches_per_var = 10, # minimum number of matches per variable
                       max_iter = 1E6, # Maximum number of iterations for matches and model fit
                       opt_method = "BFGS", # Alternative is "CG"
                       ... )
{
  resp <- as.character( lhs.vars( fo ) )
  nvars <- as.character( rhs.vars( fo ) )
  game_rank( dat,
             resp, nvars,
             fn_train, fn_eval, tsize, rounds, min_matches_per_var, max_iter, opt_method,
             ... )
}
  
  #' @title GameRank Algorithm - Variable Selection Wrapper using Random Group Comparisons and a Maximum-Likelihood Estimation
  #' @param dat data.frame or tibble comprising the whole encoded dataset
  #' @param fo Formula with LHS being response and RHS variables all variables from which to select. Tip: use as.formula( df )
  #' @param resp Character string with response variable (LHS for formula)
  #' @param nvars Character list with variables for selection (RHS terms for formula)
  #' @param fn_train Function with signature function( fo, dat, ... ) returning a model generated on a development fold
  #' @param fn_eval Function to score a model on an unknown evaluation set with signature function( fo, dat, mod, ... )
  #' @param tsize Team size determining how many variables are randomly selected and compared against each other. Default: 5
  #' @param rounds Number of rounds that each variable comparison is run on random shuffles. Default: 25
  #' @param min_matches_per_var How often does each variable have to be evaluated before the Maximum-Likelihood Estimates can be inferred. Default: 10
#' @param opt_method Parameter for optim(...) to fit the log-likelihood function for GroupRank. Default: "BFGS", switch to "CG" if Hessian becomes too large.
#' @param ... Additional parameters passed to fn_train and fn_eval for support.
#' 
#' @return list with
#' @export
#'
#' @examples
game_rank <- function( dat, # data.frame or tibble
                       # fo,  # formula of response ~ all variables for selection
                       resp, 
                       nvars,
                       # Training function
                       fn_train = function( fo, dat, ... ) { 
                         glm( fo, "binomial", dat ) 
                       },
                       # Evaluation function
                       fn_eval  = function( fo, dat, mod, ... ) { 
                         prfp <- performance( 
                           prediction( predictions = predict( mod, newdata = dat, type="response"), 
                                       labels = model.frame(fo,dat)[,1] ),
                           "auc" )
                         auc <- as.numeric( prfp@y.values )
                         return( auc )
                       },
                       tsize  = 5,  # Team Size
                       rounds = 25, # Rounds to evaluate per match
                       min_matches_per_var = 10, # minimum number of matches per variable
                       max_iter = 1E6, # Maximum number of iterations for matches and model fit
                       opt_method = "BFGS", # Alternative is "CG"
                       ...
) 
{
  
  # Check arguments ----
  stopifnot( is.data.frame(dat) | is_tibble(dat) )
  stopifnot( is.character(resp) )
  stopifnot( is.character(nvars) )
  
  # stopifnot( is_formula(fo) )
  fo <- formula( sprintf( "%s ~ %s", resp, paste(nvars, collapse = " + ") ) )
  
  # Initialize variable selection vector with 0.0
  var <- rep_len( 0.0, length(nvars) )
  names(var) <- nvars
  
  # Initialize team vector, set first team to +1, second to -1, rest to 0
  Tpm <- var
  Tpm[1:tsize] <- +1.0
  Tpm[(tsize+1):(2*tsize)] <- -1.0
  
  # Initialize row selection vector: 1 = development, 2 = evaluation
  ds <- rep_len( c(1,2), nrow(dat) )
  
  # Define team evaluation function ----
  # dat - data.frame 
  # fo  - model formula
  # ds  - vector of length nrow(dat) with 1=development, 2=evaluation data
  # returns numeric score, e.g. AUC etc with "x wins against y on the given data with the given variables <=> score(x) > score(y); NA is used if failure"
  eval_team <- function( dat, fo, ds, ... ) {
    sco <- tryCatch( {
      # mod <- glm( fo, "binomial", dat[which(1==ds),] )
      mod <- fn_train( fo, dat[which(1==ds),], ... )
      # prfp <- performance( prediction( predictions = predict( mod, newdata = dat[which(2==ds),], type="response"), 
      #                                  labels = dat[which(2==ds),resp] ), 
      #                      "auc" )
      # auc <- as.numeric( prfp@y.values )
      # auc
      evl <-  fn_eval( fo, dat[which(2==ds),], mod, ... )
      evl
    }, error = function(e) NA )
    return(sco)
  }
  
  # Run GameRank matches ----
  cat( "Comparing variable selections \n" )
  n_matches <- 1 # match counter
  res_matches <- NULL
  run <- TRUE
  while( run ) {
    Tpm[ 1:length( Tpm ) ] <- Tpm[ order( runif( length( Tpm ) ) ) ] # Shuffle teams
    # sum( Tpm )
    # which( ( Tpm < 0 ) )
    # which( ( Tpm > 0 ) )
    
    fop <- formula( sprintf( "%s ~ %s", resp, paste( names( which( Tpm > 0 ) ), collapse = " + " )) )
    fon <- formula( sprintf( "%s ~ %s", resp, paste( names( which( Tpm < 0 ) ), collapse = " + " )) )
    # print( fop )
    # print( fon )
    
    np <- 0
    nn <- 0
    for( r in 1:rounds ) {
      # Randomly split dataset into development and evaluation rows 
      ds <- ds[ order( runif( length(ds) ) ) ]
      
      # Evaluate positive Team
      scop <- eval_team( dat, fop, ds, ... )
      
      # Evaluate negative Team
      scon <- eval_team( dat, fon, ds, ... )
      
      # Compare results
      if( !is.na(scop) & !is.na(scon) ) {
        if( scop < scon ) {
          nn <- nn + 1 # (-) team wins
        } else if( scop > scon ) {
          np <- np + 1 # (+) team wins
        }
      } else if( !is.na(scop) &  is.na(scon) ) {
        np <- np + 1 # (+) team wins, as it only produces results
      } else if(  is.na(scop) & !is.na(scon) ) {
        nn <- nn + 1 # (-) team wins, as it only produces results
      }
    } # for
    
    # cat( "(+) : (-) = ", np, " : ", nn, "\n" )
    res_match <- cbind( data.frame( n.pos = np, n.neg = nn ), t( Tpm ) )
    res_matches <- bind_rows( res_matches, res_match ) 
    
    cat( "Iteration ", n_matches, " -- ",  "(+) : (-) scored ", np, " : ", nn, "\n" )
    n_matches <- n_matches + 1 # iteration counter
    run <- ((n_matches < max_iter) && (!all( colSums( abs( res_matches[,-c(1,2)])) > min_matches_per_var ) )) # Evaluate loop condition
  } # while
  
  # Evaluating Group Rank Maximum Likelihood estimator 
  # Define group rank negative log-likelihood function ----
  # Huang et al., 2008, p.10, Eq.31
  ll_gr <- local( {
    function( vs ) {
      Tp <- ( ( res_matches[,-c(1:2)] > 0 ) %*% vs )
      Tm <- ( ( res_matches[,-c(1:2)] < 0 ) %*% vs )
      LL <- exp( Tp + Tm - ( res_matches$n.pos - res_matches$n.neg ) ) / ( exp( Tp - (res_matches$n.pos - res_matches$n.neg) ) + exp( Tm ) )^2
      ret <- -sum( log( LL ) ) 
      return( ret )
    }
  }, envir = new.env() )
  
  # Define group rank gradient for negative log-likelihood function ----
  # Huang et al., 2008, p.10, Eq.33, without regularization term: ... + mu(exp(vs) - exp(-vs))
  ll_gr_grad <- local({
    function( vs ) {
      Tp <- ( ( res_matches[,-c(1:2)] > 0 ) %*% vs )
      Tm <- ( ( res_matches[,-c(1:2)] < 0 ) %*% vs )
      
      pp <- exp( Tp + res_matches[,2] ) / ( exp( Tp + res_matches[,2] ) + exp( Tm + res_matches[,1] ) )
      pn <- exp( Tm + res_matches[,1] ) / ( exp( Tp + res_matches[,2] ) + exp( Tm + res_matches[,1] ) )
      
      gr <- sapply( names(vs), FUN=function(co) {
        ms <- sum( abs(res_matches[,co]) )
        pps <- sum( pp[which(res_matches[,co] > 0)] )
        pns <- sum( pn[which(res_matches[,co] < 0)] )
        return( -ms + 2 * (pps + pns) )
      })
      names( gr ) <- names( vs )
      
      # browser()
      # gg = grad( ll_gr, vs )
      # ep = sum( abs( gr - gg ) )
      # if( ep > 1E-3 ) {
      #   cat( "Large difference between gradient ")
      #   cat( "Analytical Gradient ", gr, "\n" )
      #   cat( "Numerical Gradient  ", attr( gg, "gradient" ), "\n" )
      # } 
      
      return( gr )
    }
  }, envir = new.env() )
  
  cat( "Optimizing maximum likelihood \n" )
  oo <- optim( par = var, fn = ll_gr, gr = ll_gr_grad, method = opt_method, control = list( maxit = max_iter )  )
  oo
  
  cat( "Calculating score vector \n" )
  # gg <- grad( func = ll_gr, x=oo$par )
  # names(gg) <- names(oo$par)
  gg <- ll_gr_grad( oo$par )  
  gg
  
  cat( "Calculating hessian matrix \n" )
  # hh <- hessian( func=ll_gr, x=oo$par  )
  # rownames(hh) <- colnames(hh) <- names(oo$par)
  # hh
  hh <- jacobian( func = ll_gr_grad, x = oo$par )
  hh
  
  vv <- chol2inv( hh )
  rownames(vv) <- colnames(vv) <- colnames(hh)
  vv
  
  cat( "Compiling results \n" )
  vsel_result <- tibble( variable = names( oo$par ),
                         vs = as.numeric( oo$par ),
                         vs.var = diag( vv ) ) %>%
    mutate( selected = (vs > 0) ) %>%
    arrange( desc( vs ) )
  vsel_result %>% head
  
  # Return list object with parameters and results
  ret <- list(
    data = dat,
    formula = fo,
    response = resp,
    all_vars = nvars,
    
    team_size = tsize,
    rounds = rounds,
    min_matches_per_var = min_matches_per_var, 
    opt_method = opt_method,
    opt_result = oo,
    
    variable_ranking = vsel_result,
    
    n_matches = n_matches,
    match_results = res_matches,
    
    effects = oo$par,
    score_vector = gg,
    inv_hessian = hh
  )
  
  return( ret )
} # game_rank (END)

