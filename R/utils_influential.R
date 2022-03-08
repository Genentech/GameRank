
#' @import purrr
#' @importFrom rlang .data

#' @rdname utils_influential
#' @param mod Model to obtain coefficients from via stats::coefficients
#' @export
fn_infl_coefficients <- function( mod, ... ) stats::coefficients(mod)

#' @rdname utils_influential
#' @param mod Model to obtain predictions from via predict(...)
#' @export
fn_predict_cox <- function( mod, dat, ... ) predict( mod, newdata = dat, 
                                                     type = "lp" )

#' @rdname utils_influential
#' @param mod Model to obtain predictions from via predict(...)
#' @export
fn_predict_glm <- function( mod, dat, ... ) predict( mod, newdata = dat, 
                                                     type = "response" )

#' @title Determine influential points for the current selection
#' 
#' @param dat data.frame or tibble comprising data for model generation and 
#' validation.
#' @param resp Character string defining the response (lhs) of the model formula.
#' @param selection Character vector defining the list of variables for 
#' selection. Those are concatenated by '+' 
#' as the right hand side (rhs) of the modeling formula.
#' @param fn_train Function with signature function( dat, resp, selection, ... ) 
#' that returns a model or NULL in any other case on the given data dat.
#' @param fn_eval Function with signature function( dat, resp, selection, ... ) 
#' that returns a real number or NA in any other case, e.g. when model is NULL.
#' @param fn_get_params Function with signature function( model ) that returns 
#' an 1-row tibble of model parameters
#' @param fn_predict Function with signature function( model, data ) that 
#' returns an 1-column tibble of model predictions
#' @param c_out Constant for the Robust Outlier Test determining the scaling of 
#' IQR to determine non-outlier range from Q1 - c_out x IQR to Q3 + c_out x IQR. 
#' @param ... Further arguments passed to fn_train or fn_eval functions
#' 
#' @return Returns a tibble with columns for observation (row), evaluation 
#' without observation i (ei), difference to full dataset, difference 
#' between prediction from full model vs model without observation i (y0, yi and 
#' deval), as well as differences between parameters (_dfbeta).
#' 
#' @examples 
#' vars <- grep( "the_|rnd", colnames(toy_data), value=TRUE )
#' resp <- "resp"
#' gmr <- game_rank( dat = toy_data, resp = resp, vars = vars, 
#'                   fn_train = fn_train_binomial, fn_eval = fn_eval_binomial_auroc, 
#'                   m = 6L, dsi = c(1L,2L), maximize = TRUE, 
#'                   team_size = 3L, rounds = 10L, min_matches_per_var = 5L )
#' gmr$variable_ranking %>% as.data.frame
#' gmr_fsel <- gmr$game_rank_selection
#' 
#' ifo <- influential_observations( toy_data, resp, gmr_fsel, 
#'                                  fn_train_binomial, fn_eval_binomial_auroc, 
#'                                  fn_infl_coefficients, fn_predict_glm )
#' ifo 
#' 
#' @name utils_influential
#' @export
influential_observations <- function( dat, resp, selection, 
                                      fn_train, fn_eval,
                                      fn_get_params, fn_predict, 
                                      c_out = 1.5, ... ) 
{
  stopifnot( is.data.frame(dat) | is_tibble(dat) )
  stopifnot( is.character(selection) & 1 < length(selection) )
  stopifnot( is.character(resp) )
  
  stopifnot( is.function(fn_train) ) 
  stopifnot( is.function(fn_eval) ) 
  stopifnot( is.function(fn_get_params) ) 
  stopifnot( is.function(fn_predict) ) 
  stopifnot( is.numeric(c_out) & (0 < c_out) ) 
  
  message( sprintf("influential_observations: Determining influence statistics [dfbetas,dffits approach] for %d observations with %d variables.", nrow(dat), length(selection)) )
  # Reference model
  m0 <- fn_train( dat, resp, selection, ... )
  c0 <- fn_get_params( m0, ... )
  e0 <- fn_eval( dat, resp, selection, m0, ... )
  res <- tibble( row = NA_integer_, is_influential = NA, 
                 is_influential_co = NA_character_, 
                 ei = e0, deval = NA_real_, yi = NA_real_, dffit = NA_real_ ) %>%
    bind_cols( as_tibble( matrix( c0, nrow=1 ) ) %>% 
                 stats::setNames( sprintf( "%s_dfbeta", names(c0)) ) )
  
  ret <- purrr::map_dfr( seq_len( nrow(dat) ),
                  function(i) {
                    # browser()
                    mm <- fn_train( dat[-i,], resp, selection )
                    cc <- fn_get_params( mm )
                    ee <- fn_eval( dat, resp, selection, mm, ... )
                    
                    y0 <- fn_predict( m0, dat[i,], ... ) 
                    yy <- fn_predict( mm, dat[i,], ... )
                    
                    res <- tibble( row = i, ei = ee, deval = e0 - ee, 
                                   yi = yy, dffit = y0 - yy ) %>%
                      bind_cols( as_tibble( matrix( c0 - cc, nrow=1 ) ) %>% 
                                   stats::setNames( sprintf( "%s_dfbeta", 
                                                             names(c0)) ) )
                  })
  
  ret$is_influential <- FALSE
  ret$is_influential_co <- ""
  for( co in c("deval","dffit", sprintf("%s_dfbeta",names(c0))) ) { 
    cut_min <- stats::quantile( ret[[co]], probs = 0.25, na.rm = TRUE ) - 
      c_out * stats::IQR( ret[[co]], na.rm = TRUE )
    cut_max <- stats::quantile( ret[[co]], probs = 0.75, na.rm = TRUE ) + 
      c_out * stats::IQR( ret[[co]], na.rm = TRUE )
    idx <- which( ret[[co]] < cut_min | cut_max < ret[[co]] )
    ret$is_influential[ idx ] <- TRUE
    ret$is_influential_co[ idx ] <- sprintf( "%s %s", 
                                             ret$is_influential_co[ idx ], co )
  }
  
  ret <- bind_rows(res, ret) %>% 
    mutate( is_influential_co = trimws( .data$is_influential_co ))
  return( ret )
}
