
#' @import ggplot2 knitr rmarkdown broom

#' @title Render a standard template reports for a given objects.
#' 
#' @param obj Algorithm results output from forward, backward and other 
#' selection methods
#' @param output_dir Output directory for the report to be saved to.
#' @param dat data.frame or tibble rows from the full dataseet provided to the
#' wrapper that should be used for generating or evaluating models.
#' @param resp Response variable being the lhs of the model formula
#' @param selection Current selection for model generation or evaluation
#' @param mod For evaluation functions the model to be evaluated on dat
#' 
#' @return render_... functions don't have return values. They generated the 
#' required output documents.
#' 
#' @examples 
#' 
#' vars <- grep( "the_|rnd", colnames(toy_data), value=TRUE )
#' resp <- "resp"
#' gmr <- game_rank( dat = toy_data, resp = resp, vars = vars, 
#'                   fn_train = fn_train_binomial, fn_eval = fn_eval_binomial_auroc,
#'                   m = 6L, dsi = c(1L,2L), maximize = TRUE, 
#'                   team_size = 3L, rounds = 10L, min_matches_per_var = 5L )
#' gmr$variable_ranking %>% as.data.frame
#' gmr_fsel <- gmr$game_rank_selection
#'   
#' \donttest{  render_game_rank_summary( gmr, getwd() ) }
#' 
#' @name utils_reporting
NULL


#' @title Render a standard template report for a given object. (Internal function)
#' 
#' @param template_name File name of the standard template
#' @param output_dir Output directory to which to render
#' @param params The Rmd parameter list. Needs to match the expected 
#' parameter definitions of the Rmd template.
#' @param output_file Output filename (with or without suffix). Default NULL. 
#' If not NULL, it is appended to output_dir
#' @param envir Environment passed to rmarkdown::render containing R objects
#' for the Rmd file.
#' @param ... Further parameters passed to rmarkdown::render(...) function.
#' 
#' @rdname utils_reporting
render_std_template <- function( template_name, output_dir = NULL, params, output_file = NULL, envir = new.env(), ... ) {
  stopifnot(!is.null(output_dir))
  stopifnot( is.list(params))
  
  template_file <- system.file( "templates", template_name, package="GameRank" )
  rmarkdown::render( input = template_file, 
                     output_dir = output_dir,
                     output_file = output_file,
                     params = params, 
                     envir = envir,
                     ...  )
}

#' @rdname utils_reporting
#' @export
render_backward_summary <- function( obj, output_dir = NULL, output_file = NULL ) {
  render_std_template(  template_name = "rmd_backward_summary.Rmd", 
                        output_dir = output_dir, 
                        params = list( bwd = obj ), 
                        output_file = output_file, envir = new.env() )
}

#' @rdname utils_reporting
#' @export
render_forward_summary <- function( obj, output_dir = NULL, output_file = NULL ) {
  render_std_template(  template_name = "rmd_backward_summary.Rmd", 
                        output_dir = output_dir, 
                        params = list( fwd = obj ), 
                        output_file = output_file, envir = new.env() )
}

#' @rdname utils_reporting
#' @export
render_bidirectional_summary <- function( obj, output_dir = NULL, output_file = NULL ) {
  render_std_template(  template_name = "rmd_bidirectional_summary.Rmd", 
                        output_dir = output_dir, 
                        params = list( bds = obj ), 
                        output_file = output_file, envir = new.env() )
}

#' @rdname utils_reporting
#' @export
render_lrsearch_summary <- function( obj, output_dir = NULL, output_file = NULL ) {
  render_std_template(  template_name = "rmd_lrsearch_summary.Rmd", 
                        output_dir = output_dir, 
                        params = list( lrs = obj ), 
                        output_file = output_file, envir = new.env() )
}

#' @rdname utils_reporting
#' @export
render_random_summary <- function( obj, output_dir = NULL, output_file = NULL ) {
  render_std_template(  template_name = "rmd_random_summary.Rmd", 
                        output_dir = output_dir, 
                        params = list( rnd = obj ), 
                        output_file = output_file, envir = new.env() )
}

#' @rdname utils_reporting
#' @export
render_game_rank_summary <- function( obj, output_dir = NULL, output_file = NULL ) {
  render_std_template(  template_name = "rmd_game_rank_summary.Rmd", 
                        output_dir = output_dir, 
                        params = list( gmr = obj ), 
                        output_file = output_file, envir = new.env() )
}

#' @rdname utils_reporting
#' @export
render_variable_checks_summary <- function( obj, output_dir = NULL, output_file = NULL ) {
  render_std_template(  template_name = "rmd_variable_checks.Rmd", 
                        output_dir = output_dir, 
                        params = list( vck = obj ), 
                        output_file = output_file, envir = new.env() )
}

#' @rdname utils_reporting
#' @param ds Definition of (parallel) training:validation splits
#'  - a matrix with d columns containing 1s and 2s, where 1 denotes sample is 
#'    used for training the model and 2 denotes sample used for validation.
#'    The average of all d training:validation results is used for selection.
#'  - an integer number determining the number of random training:validation 
#'    splits that should be generated. The sampling will ensure a sufficient 
#'    number of complete cases in the training split.
#' @param k Column from ds to use for split into training:validation:test via
#' 1s, 2s and 3s.
#' @export
render_model_calibration_normal <- function( ds, dat, resp, selection, k, 
                                             output_dir = NULL, output_file = NULL ) {
  render_std_template(  template_name = "rmd_calibration_normal.Rmd", 
                        output_dir = output_dir, 
                        params = list( ds = ds, dat = dat, resp = resp, selection = selection, k = k ), 
                        output_file = output_file, envir = new.env() )
}


#' @rdname utils_reporting
#' @param k Column from ds to use for split into training:validation:test via
#' 1s, 2s and 3s.
#' @export
render_model_calibration_binomial <- function( ds, dat, resp, selection, k, 
                                               output_dir = NULL, output_file = NULL ) {
  render_std_template(  template_name = "rmd_calibration_binomial.Rmd", 
                        output_dir = output_dir, 
                        params = list( ds = ds, dat = dat, resp = resp, selection = selection, k = k ), 
                        output_file = output_file, envir = new.env() )
}


#' @rdname utils_reporting
#' @param k Column from ds to use for split into training:validation:test via
#' 1s, 2s and 3s.
#' @param u Landmark time point at which the survival probability is evaluated
#' @export
render_model_calibration_cox <- function( ds, dat, resp, selection, k, u,
                                          output_dir = NULL, output_file = NULL ) {
  render_std_template(  template_name = "rmd_calibration_cox.Rmd", 
                        output_dir = output_dir, 
                        params = list( ds = ds, dat = dat, resp = resp, selection = selection, k = k, u = u ), 
                        output_file = output_file, envir = new.env() )
}

#' @rdname utils_reporting
#' @export
render_influence_summary <- function( obj, output_dir = NULL, output_file = NULL ) {
  render_std_template(  template_name = "rmd_influence_summary.Rmd", 
                        output_dir = output_dir, 
                        params = list( tinf = obj ), 
                        output_file = output_file, envir = new.env() )
}

