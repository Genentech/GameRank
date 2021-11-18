


#' @title Basic get_params function for lm and glm models calling coefficients
fn_infl_coefficients <- function( mod, ... ) coefficients(mod)

#' @title Linear predictor function for Cox PH models
fn_predict_cox <- function( mod, dat, ... ) predict( mod, newdata = dat, type = "lp" )

#' @title Determine influential points for the current selection
#' 
#' @param dat data.frame or tibble comprising data for model generation and validation.
#' @param resp Character string defining the response (lhs) of the model formula.
#' @param selection Character vector defining the list of variables for selection. Those are concatenated by '+' 
#' as the right hand side (rhs) of the modeling formula.
#' @param fn_train Function with signature function( dat, resp, selection, ... ) that returns a model or NULL in any other case on the given data dat.
#' @param fn_eval Function with signature function( dat, resp, selection, ... ) that returns a real number or NA in any other case, e.g. when model is NULL.
#' @param fn_get_params Function with signature function( model ) that returns an 1-row tibble of model parameters
#' @param fn_predict Function with signature function( model, data ) that returns an 1-column tibble of model predictions
#' 
#' @return Returns a tibble with columns for observation (row), evaluation without observation i (ei), difference to full dataset, difference 
#' between prediction from full model vs model without observation i (y0, yi and deval), as well as differences between parameters (_dfbeta).
influential_observations <- function( dat, resp, selection, 
                                      fn_train, fn_eval,
                                      fn_get_params, fn_predict, ... ) 
{
  
  # Reference model
  m0 <- fn_train( dat, resp, selection, ... )
  c0 <- fn_get_params( m0, ... )
  e0 <- fn_eval( dat, resp, selection, m0, ... )
  res <- tibble( row = NA_integer_, ei = e0, deval = NA_real_, yi = NA_real_, dffit = NA_real_ ) %>%
    bind_cols( as_tibble( matrix( c0, nrow=1 ) ) %>% setNames( sprintf( "%s_dfbeta", names(c0)) ) )
  
  ret <- map_dfr( 1:nrow(dat),
                  function(i) {
                    # browser()
                    mm <- fn_train( dat[-i,], resp, selection )
                    cc <- fn_get_params( mm )
                    ee <- fn_eval( dat, resp, selection, mm, ... )
                    
                    y0 <- fn_predict( m0, dat[i,], ... ) # predict( m0, newdata = dat[i,], type = "lp" )
                    yy <- fn_predict( mm, dat[i,], ... ) # predict( mm, newdata = dat[i,], type = "lp" )
                    
                    res <- tibble( row = i, ei = ee, deval = e0 - ei, yi = yy, dffit = y0 - yi ) %>%
                      bind_cols( as_tibble( matrix( c0 - cc, nrow=1 ) ) %>% setNames( sprintf( "%s_dfbeta", names(c0)) ) )
                  })
  
  ret$is_influential <- FALSE
  ret$is_influential_co <- ""
  for( co in c("deval","dffit", sprintf("%s_dfbeta",names(c0))) ) { 
    ret$is_influential[ which.max( abs(ret[[co]] ) ) ] <- TRUE
    ret$is_influential_co[ which.max( abs(ret[[co]] ) ) ] <- sprintf( "%s %s", ret$is_influential_co[ which.max( abs(ret[[co]] ) ) ], co )
  }
  
  ret <- bind_rows(res, ret)
  return( ret )
}
