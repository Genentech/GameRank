#
# Support file providing some modelling and evaluation functions
#

# Regression model ----
fn_train_normal_ml_imp <- function( dat, resp, selection, ... ) {
  fo <- formula( sprintf( "%s ~ %s", resp, paste( selection, collapse = " + " ) ) )
  fo <- rewrite_formula( fo, dat )
  mo <- glm( formula = fo, family = gaussian, data = dat )
  return( mo )
}

# Binomial model ----
fn_train_binomial_ml_imp <- function( dat, resp, selection, ... ) {
  fo <- formula( sprintf( "%s ~ %s", resp, paste( selection, collapse = " + " ) ) )
  fo <- rewrite_formula( fo, dat )
  mod <- glm( formula = fo, family = binomial, data = dat )
  return( mod )
}

# fn_eval_binomial <- function( dat, resp, selection, mod, ... ) {
#   fo <- as.formula( mod )
#   prfp <- ROCR::performance( prediction = ROCR::prediction( predictions = predict( mod, newdata = dat, type = "response" ),
#                                                             labels = as.character( model.frame( formula = fo, data = dat )[,1] ) ), 
#                              measure = "auc"  )
#   auc <- as.numeric( prfp@y.values )
#   return( auc )
# }


# Survival model ----
fn_train_cox_ml_imp <- function( dat, resp, selection, ... ) {
  fo <- formula( sprintf( "%s ~ %s", resp, paste( selection, collapse = " + " ) ) )
  fo <- rewrite_formula( fo, dat )
  mod <- coxph( formula = fo, data = dat, x=TRUE, y=TRUE )
  return( mod )
}
