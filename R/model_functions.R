#
# Support file providing some modelling and evaluation functions
#

# Regression model ----
fn_train_normal <- function( dat, resp, selection, ... ) {
  fo <- formula( sprintf( "%s ~ %s", resp, paste( selection, collapse = " + " ) ) )
  mo <- glm( formula = fo, family = gaussian, data = dat )
  return( mo )
}

fn_eval_normal <- function( dat, resp, selection, mod, ... ) {
  ret <- tibble( y = dat %>% pull( resp ), 
                 yhat = as.numeric( predict(mod, newdata=dat, type="response" ) ) ) %>%
    mutate( sq = (yhat-y)^2 ) %>%
    pull( sq ) %>%
    mean

  return( ret )
}

# Binomial model ----
fn_train_binomial <- function( dat, resp, selection, ... ) {
  fo <- formula( sprintf( "%s ~ %s", resp, paste( selection, collapse = " + " ) ) )
  mod <- glm( formula = fo, family = binomial, data = dat )
  return( mod )
}

fn_eval_binomial <- function( dat, resp, selection, mod, ... ) {
  ret <- tibble( y = dat %>% pull( resp ), 
                 yhat = as.numeric( predict(mod, newdata=dat, type="response" ) ) ) %>%
    mutate( di = abs( y-yhat ) ) %>% 
    pull( di ) %>%
    mean 

    return( ret )  
}
'<.mad' <- function( a,b ) (as.numeric(b) < as.numeric(a)) # check b better than a [a<b]

# fn_eval_binomial <- function( dat, resp, selection, mod, ... ) {
#   fo <- as.formula( mod )
#   prfp <- ROCR::performance( prediction = ROCR::prediction( predictions = predict( mod, newdata = dat, type = "response" ),
#                                                             labels = as.character( model.frame( formula = fo, data = dat )[,1] ) ), 
#                              measure = "auc"  )
#   auc <- as.numeric( prfp@y.values )
#   return( auc )
# }


# Survival model ----
fn_train_cox <- function( dat, resp, selection, ... ) {
  fo <- formula( sprintf( "%s ~ %s", resp, paste( selection, collapse = " + " ) ) )
  mod <- coxph( formula = fo, data = dat, x=TRUE, y=TRUE )
  return( mod )
}

#
# Following Austin et al., "Graphical calibration curves and the integration calibration index (ICI) for survival models, Statistics in Medicine, 2020
#
fn_eval_cox <- function( dat, resp, selection, mod, u = NULL, ... ) {
  stopifnot( !is.null(u) & is.numeric(u) & (0<=u) )
  
  prd <- 1 - predictSurvProb( mod, newdata = dat, times = u, ... )
  cll.prd <- log( -log(1-prd) ) 
  cal.cox <- coxph( formula( sprintf( "%s ~ rcs(cll.prd,3)", resp )), data = bind_cols( dat, tibble( cll.prd = cll.prd )), x=TRUE )
  grd.cox <- seq( quantile( prd, probs = 0.01 ), 
                  quantile( prd, probs = 0.99 ), 
                  length = 100 )
  grd.prd.cll <- log(-log(1-grd.cox))
  df.grd.cox <- data.frame( grd.cox, grd.prd.cll ) %>% setNames(c("prd","cll.prd"))
  df.grd.cox$prd.cal <- 1 - predictSurvProb( cal.cox, df.grd.cox, times = u )
  ret <- mean( abs( df.grd.cox$prd - df.grd.cox$prd.cal ) )

  return( ret )  
}
