#
# Support file providing some modelling and evaluation functions
#

# Regression model ----
#' @rdname model_functions
#' @export
fn_train_normal_ml_imp <- function( dat, resp, selection, ... ) {
  fo <- stats::formula( sprintf( "%s ~ %s", 
                                 resp, paste( selection, collapse = " + " ) ) )
  fo <- rewrite_formula( fo, dat )
  mo <- stats::glm( formula = fo, family = stats::gaussian, data = dat )
  return( mo )
}

# Binomial model ----
#' @rdname model_functions
#' @export
fn_train_binomial_ml_imp <- function( dat, resp, selection, ... ) {
  fo <- stats::formula( sprintf( "%s ~ %s", 
                                 resp, paste( selection, collapse = " + " ) ) )
  fo <- rewrite_formula( fo, dat )
  mod <- stats::glm( formula = fo, family = stats::binomial, data = dat )
  return( mod )
}

# Survival model ----
#' @rdname model_functions
#' @export
fn_train_cox_ml_imp <- function( dat, resp, selection, ... ) {
  fo <- stats::formula( sprintf( "%s ~ %s", 
                                 resp, paste( selection, collapse = " + " ) ) )
  fo <- rewrite_formula( fo, dat )
  mod <- coxph( formula = fo, data = dat, x=TRUE, y=TRUE )
  return( mod )
}
