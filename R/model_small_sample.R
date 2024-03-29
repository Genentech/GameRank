#
#  Utils - Small Sample Selection
#
# 

# ' @importFrom rlang .

#' @title Small Sample Variable Selection helper functions
#' 
#' @description Wrapper-based variable selection makes use of two functions, one
#' for model training and one for model evaluation. This roots in the assumption
#' of using the split set method for measuring generalization performance on 
#' true hold-out data. However, when you have a very small sample set this might 
#' seem problematic.
#' 
#' An way out then could be to estimate generalization performance using 
#' **cross-validation** or the **bootstrap**. In both approaches, the training 
#' function must be called on repeated subsamples (folds) such that the error
#' can be estimated by combining predictions from the excluded samples across 
#' all folds.
#' 
#' Here we'll provide two function factory methods, ie functions that take a 
#' training or a validation function and wrap it into one of these algorithms 
#' such that the call to the wrapper functions keeps the same interface. 
#' The idea will be to replace fn_train by a dummy method that just returns a 
#' dummy model. The fn_eval parameter will then consist of a function that
#' incorporates the right fn_train function and applies it within the estimation 
#' procedure.
#' 
#' \describe{
#' \item{fn_train_dummy}{Dummy train function that just returns and 
#' class "DummyModel" object.}
#' \item{ff_fn_eval_cross_validation}{Function generator (ff_) function that 
#' takes an fn_train and fn_eval function together with n_folds parameter, and 
#' returns a function that calls fn_train and fn_eval to return a n_fold
#' cross-validated estimate.}
#' }
#' @param dat data.frame or tibble rows from the full dataseet provided to the
#' wrapper that should be used for generating or evaluating models.
#' @param resp Response variable being the lhs of the model formula
#' @param selection Current selection for model generation or evaluation
#' @param ... Any other arguments passed to both types of functions, 
#' 
#' @return fn_train_... functions return a fitted model object or NULL if that
#' fails. fn_eval_... functions return a numeric real value measuring the 
#' validation performance on the given data, or NA if that fails.
#' 
#' @examples 
#' vars <- grep( "the_|rnd", colnames(toy_data), value=TRUE )
#' resp <- "resp"
#' ipos <- which(  toy_data$resp )
#' ineg <- which( !toy_data$resp )
#' smtoy <- rbind( 
#'   # draw 50 random observations from each response category
#'   toy_data[ ipos[1:25],], 
#'   toy_data[ ineg[1:25],]
#' )
#' smtoy
#' 
#' # ff_fn_eval_bootstrap and ff_fn_eval_cross_validation are function factories
#' # that return a created function which makes use of the supplied fn_train and
#' # fn_eval function during bootstrap or cross-validation.
#' fn_bt_eval_binomial_auroc <- ff_fn_eval_bootstrap( fn_train_binomial, 
#'                                                    fn_eval_binomial_auroc, 
#'                                                    25L )
#' 
#' # For small sample sizes, the GameRank index vector dsi needs to be set 
#' # to c(2L,2L) such that all observations are used for validation. Similarly, 
#' # the training:validation split matrix then also completely consists of 2s.
#' res <- game_rank( smtoy, resp, vars, 
#'                   fn_train_dummy, fn_bt_eval_binomial_auroc, 
#'                   4L, c(2L,2L), TRUE, 3L, 5L, 3L )
#' res$game_rank_selection
#' res$variable_ranking
#' 
#' 
#' @name model_small_sample
NULL

#' @rdname model_small_sample
#' @export
fn_train_dummy <- function( dat, resp, selection, ... ) {
  ret <- "Small sample Dummy model"
  class(ret) <- "DummyModel"
  return( ret )
}


#' @rdname model_small_sample
#' @param arg_fn_train Function to generate model on non-hold-out 
#' folds (fn_train)
#' @param arg_fn_eval Function for validation on hold-out fold (fn_eval)
#' @param arg_n_folds Number of cross-validation folds to use
#' @export
ff_fn_eval_cross_validation <- function( arg_fn_train = NULL, 
                                         arg_fn_eval = NULL, 
                                         arg_n_folds = 10L ) {
  stopifnot( !is.null(arg_fn_train) & !is.null(arg_fn_eval) )
  stopifnot( is.function(arg_fn_train) )
  stopifnot( is.function(arg_fn_eval) )
  stopifnot( is.integer(arg_n_folds) & (1 < arg_n_folds) )
  
  fu <- local({
    # Definition of constructed fn_eval 
    function( dat, resp, selection, mod, ... ) {
      if( is.null(mod) ) return( NA )
      if( "DummyModel"!=class(mod) ) return( NA )
      ret <- NA
      
      #
      # Reference:
      # Larry Wasssermann, "All of non-parametric statistics", Wiley, 2008
      # Chapter 3.1, "The Jackknife", p.27, equation 3.1
      #
      
      # 1) Estimate Tn
      mod <- fn_train( dat, resp, selection, ... )
      Tn <- fn_eval( dat, resp, selection, mod, ... )
      
      # 2) Estimate bias theta of Tn via jackknife / cross-validation
      # Perform n_fold Cross-Validation here. Note: to avoid over-fitting to the
      # folder, each time this function is called a new assignment to folds is
      # sampled. This may seem less robust but prevents variable selection to
      # favor a specific fixed fold configuration.
      n <- nrow(dat)
      ds <- rep_len( x = seq_len( n_folds ), length.out = n )
      ds <- ds[ order( stats::runif( length(ds) ) ) ]
      theta <- rep_len( NA, length.out = length(ds) )
      for( k in seq_len( n_folds ) ) {
        idx <- which( k == ds )
        cvmod <- fn_train( dat[-idx,], resp, selection, ... )
        cvevl <- fn_eval(  dat[ idx,], resp, selection, cvmod, ...  )
        theta[idx] <- cvevl
      }
      # 3) Adjust estimate Tn by bias theta
      bjack <- (n-1) * ( mean( theta, na.rm=TRUE ) - Tn )
      Tjack <- Tn - bjack
      
      ret <- Tjack 
      return( ret )
    } # END(function)
  }, envir = list( # Wrap function factory arguments into local environment
    fn_train = arg_fn_train, 
    fn_eval = arg_fn_eval, 
    n_folds = arg_n_folds 
    # END( envir argument )
  ) )
  return( fu ) # Return constructed function
}


#' @rdname model_small_sample
#' @param arg_fn_train Function to generate model on non-hold-out 
#' folds (fn_train)
#' @param arg_fn_eval Function for validation on hold-out fold (fn_eval)
#' @param arg_n_boots Number of bootstrapped folds to use
#' @export
ff_fn_eval_bootstrap <- function( arg_fn_train = NULL, 
                                  arg_fn_eval = NULL, 
                                  arg_n_boots = 25L ) {
  stopifnot( !is.null(arg_fn_train) & !is.null(arg_fn_eval) )
  stopifnot( is.function(arg_fn_train) )
  stopifnot( is.function(arg_fn_eval) )
  stopifnot( is.integer(arg_n_boots) & (1 < arg_n_boots) )
  
  fu <- local({
    # Definition of constructed fn_eval 
    function( dat, resp, selection, mod, ... ) {
      if( is.null(mod) ) return( NA )
      if( "DummyModel"!=class(mod) ) return( NA )
      if( 0==nrow(dat) ) return( NA )
      ret <- NA
      
      #
      # From Efron and Tibshirani, "An introduction to the bootstrap", 
      # Chapman & Hall/CRC, 1993, Ch. 17 "Cross-validation and other estimates 
      # of prediction error", p. 237ff,  Ch. 17.6-7 Bootstrap estimates of
      # prediction error / The .632 bootstrap estimator, p. 247f
      #
      didx <- seq_len( nrow(dat) )
      the_statistic <- function( dd,ii ) {
        # browser()
        mm  <- fn_train( dd[ii,], resp, selection, ... )
        vv  <- fn_eval( dd, resp, selection, mm, ... )
        vvb <- fn_eval( dd[ii,], resp, selection, mm, ... )
        jj  <- setdiff( didx, ii )
        vvo <- NA
        err0 <- NA
        if( 0 < length(jj) ) {
          vvo <- fn_eval( dd[jj,], resp, selection, mm, ... )  
          err0 <- vvo / length(jj)
        }
        
        st  <- c( bte_pred_error = vv,  # Bootstrap estimate of prediction error
                  apparent_error = vvb, # Apparent error
                  optimism = vv - vvb,  # Optimism
                  oobs_error = vvo,     # Out-of-bag sample error 
                  epsilon0   = err0     # Averaged OOB sample error (epsilon0)
        )
        return( st )
      } # the_statistics (END)
      
      bt <- boot::boot( data = dat,
                        statistic = the_statistic, 
                        R = n_boots,
                        sim = "ordinary",
                        stype = "i" )
      # Go with the optimism corrected bootstrap prediction error for now
      ret <- mean( bt$t[,which(names(bt$t0)=="bte_pred_error")], na.rm=TRUE ) + 
        mean( bt$t[,which(names(bt$t0)=="optimism")], na.rm=TRUE ) 
      return( ret )
    } # END(function)
  }, envir = list( # Wrap function factory arguments into local environment
    fn_train = arg_fn_train, 
    fn_eval = arg_fn_eval, 
    n_boots = arg_n_boots 
    # END( envir argument )
  ) )
  return( fu ) # Return constructed function
}
