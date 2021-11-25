#
# Support file providing some modeling and evaluation functions
#

#' @import flexmix rms

#' @title Functions for Model fitting and evaluation
#' 
#' @description The interface for model generating functions is
#' \strong{function( dat, resp, selection, ... )}
#' where dat is a data.frame or tipple sliced out of the overall
#' dataset according to the split matrix (1 = training) and (2 = validation).
#' Furthermore, resp denotes a character string that should be used as
#' response variable (on the lhs of a formula) and selection denotes
#' a vector of character strings that represents the feature combination to
#' be evaluated.
#' Model function should return NULL if the model cannot be fit or determined or
#' other reasons exist for why it cannot be generated.
#' 
#' The interface for evaluation functions is \strong{function( dat, resp, selection, mod, ... )}
#' with the same parameters as model generating functions but extended by the model mod. 
#' Evaluation functions should return a numeric value and NA in case of the model is NULL
#' or any other failure preventing evaluation.
#' 
#' Model functions ending with _imp modify the selection variables to perform a maximum-likelihood
#' imputation for missing values. Precisely, let \eqn{x} be the variable it is replaced by
#' \eqn{\phi(x) + \psi(x)} where \eqn{\phi(x) = x if x is not missing, and 0 otherwise} and
#' \eqn{\psi(x) = 1 if x is missing, and 0 otherwise}. This equates to a 0-1-Dummy coding for
#' missing values and the effect for \eqn{\psi} can be interpreted as the maximum-likelihood
#' imputation value to be used if \eqn{x} is missing.
#' 
#' 
#' @param dat data.frame or tibble rows from the full dataseet provided to the wrapper that should
#' be used for generating or evaluating models.
#' @param resp Response variable being the lhs of the model formula
#' @param selection Current selection for model generation or evaluation
#' @param mod For evaluation functions the model to be evaluated on dat
#' @param ... Any other arguments passed to both types of functions, e.g. 'u = 365' to define the
#' landmark day for survival probability evaluation.
#'  
#' @name model_functions
NULL

# Regression model ----
#' @rdname model_functions
#' @export
fn_train_normal <- function( dat, resp, selection, ... ) {
  mod <- NULL
  mod <- tryCatch({
    fo <- stats::formula( sprintf( "%s ~ %s", resp, paste( selection, collapse = " + " ) ) )
    mod <- stats::glm( formula = fo, family = stats::gaussian, data = dat )
  }, error = function( e ) NULL )
  return( mod )
}

#' @rdname model_functions
#' @export
fn_eval_normal <- function( dat, resp, selection, mod, ... ) {
  if( is.null(mod) ) return( NA )
  ret <- NA
  ret <- tryCatch({
    tibble( y = dat %>% pull( resp ), 
            yhat = as.numeric( predict(mod, newdata=dat, type="response" ) ) ) %>%
      mutate( sq = (yhat-y)^2 ) %>%
      pull( sq ) %>%
      mean
  }, error = function(e) NA )
  return( ret )
}

# Binomial model ----
#' @rdname model_functions
#' @export
fn_train_binomial <- function( dat, resp, selection, ... ) {
  mod <- NULL
  mod <- tryCatch({
    fo <- stats::formula( sprintf( "%s ~ %s", resp, paste( selection, collapse = " + " ) ) )
    mod <- stats::glm( formula = fo, family = stats::binomial, data = dat )
  }, error = function( e ) NULL )
  return( mod )
}

#' @rdname model_functions
#' @export
fn_eval_binomial <- function( dat, resp, selection, mod, ... ) {
  if( is.null(mod) ) return( NA )
  ret <- NA
  ret <- tryCatch({
    tibble( y = dat %>% pull( resp ), 
            yhat = as.numeric( predict(mod, newdata=dat, type="response" ) ) ) %>%
      mutate( di = abs( y-yhat ) ) %>% 
      pull( di ) %>%
      mean 
  }, error = function( e ) NA )
  return( ret )  
}

#' @rdname model_functions
#' @export
fn_eval_binomial_auroc <- function( dat, resp, selection, mod, ... ) {
  if( is.null(mod) ) return( NA )
  ret <- NA
  ret <- tryCatch({
    fo <- stats::as.formula( mod )
    prfp <- ROCR::performance( prediction = ROCR::prediction( predictions = predict( mod, newdata = dat, type = "response" ),
                                                              labels = as.character( stats::model.frame( formula = fo, data = dat )[,1] ) ),
                               measure = "auc"  )
    auc <- as.numeric( prfp@y.values )
    auc
  }, error = function(e) NA )
  return( ret )
}

# Survival model ----
#' @rdname model_functions
#' @export
fn_train_cox <- function( dat, resp, selection, ... ) {
  mod <- NULL
  mod <- tryCatch({
    fo <- stats::formula( sprintf( "%s ~ %s", resp, paste( selection, collapse = " + " ) ) )
    mod <- coxph( formula = fo, data = dat, x=TRUE, y=TRUE )
  }, error = function(e) NULL )
  return( mod )
}

#
# Following Austin et al., "Graphical calibration curves and the integration calibration index (ICI) for survival models, Statistics in Medicine, 2020
#
#' @rdname model_functions
#' @export
fn_eval_cox <- function( dat, resp, selection, mod, u = NULL, ... ) {
  if( is.null(mod) ) return( NA )
  stopifnot( is.numeric(u) & (0<=u) )
  ret <- NA
  ret <- tryCatch({
    mod$coefficients[which(is.na(mod$coefficients))] <- 0
    prd <- 1 - pec::predictSurvProb( mod, newdata = dat, times = u, ... )
    cll.prd <- log( -log(1-prd) ) 
    cal.cox <- coxph( stats::formula( sprintf( "%s ~ rms::rcs(cll.prd,3)", resp )), data = bind_cols( dat, tibble( cll.prd = cll.prd )), x=TRUE )
    # cal.cox <- coxph( formula( sprintf( "%s ~ rms::lsp(cll.prd,3)", resp )), data = bind_cols( dat, tibble( cll.prd = cll.prd )), x=TRUE )
    grd.cox <- seq( stats::quantile( prd, probs = 0.01, na.rm=TRUE ), 
                    stats::quantile( prd, probs = 0.99, na.rm=TRUE ), 
                    length = 100 )
    grd.prd.cll <- log(-log(1-grd.cox) )
    df.grd.cox <- data.frame( grd.cox, grd.prd.cll ) %>% stats::setNames(c("prd","cll.prd"))
    df.grd.cox$prd.cal <- 1 - pec::predictSurvProb( cal.cox, df.grd.cox, times = u )
    ret <- mean( abs( df.grd.cox$prd - df.grd.cox$prd.cal ) )
  }, error = function( e ) NA )
  return( ret )  
}

# Linear Discriminant Analysis ----
#' @rdname model_functions
#' @export
fn_train_lda <- function( dat, resp, selection, lda_fit_type = "mle", ... ) {
  mod <- NULL
  mod <- tryCatch({
    fo <- stats::formula( sprintf( "%s ~ %s", resp, paste( selection, collapse = " + " ) ) )
    mod <- MASS::lda( formula = fo, data = dat, method = lda_fit_type, ... )
  }, error = function(e) NULL )
  return( mod )
}

#' @rdname model_functions
#' @export
fn_eval_lda <- function( dat, resp, selection, mod, lda_pred_type = "plug-in", ... ) {
  if( is.null(mod) ) return( NA )
  ret <- NA
  ret <- tryCatch({
    tibble( y = dat %>% pull( resp ) %>% as.character, 
            yhat = as.character( predict(mod, newdata=dat, method = lda_pred_type, ... )$class ) ) %>%
      mutate( di = abs( y == yhat ) ) %>% 
      pull( di ) %>%
      mean 
  }, error = function( e ) NA )
  return( ret )
}

# Quadratic Discriminant Analysis ----
#' @rdname model_functions
#' @export
fn_train_qda <- function( dat, resp, selection, qda_fit_type = "mle", ... ) {
  mod <- NULL
  mod <- tryCatch({
    fo <- stats::formula( sprintf( "%s ~ %s", resp, paste( selection, collapse = " + " ) ) )
    mod <- MASS::qda( formula = fo, data = dat, method = qda_fit_type, ... )
  }, error = function(e) NULL )
  return( mod )
}

#' @rdname model_functions
#' @export
fn_eval_qda <- function( dat, resp, selection, mod, qda_pred_type = "plug-in", ... ) {
  if( is.null(mod) ) return( NA )
  ret <- NA
  ret <- tryCatch({
    tibble( y = dat %>% pull( resp ) %>% as.character, 
            yhat = as.character( predict(mod, newdata=dat, method = qda_pred_type, ... )$class ) ) %>%
      mutate( di = abs( y == yhat ) ) %>% 
      pull( di ) %>%
      mean 
  }, error = function( e ) NA )
  return( ret )
}