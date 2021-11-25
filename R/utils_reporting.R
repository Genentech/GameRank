
#' @import ggplot2 gridExtra gridExtra knitr rmarkdown broom

#' @title Render a standard template reports for a given objects.
#' 
#' @param obj Algorithm results output from forward, backward and other selection methods
#' @param output_dir Output directory for the report to be saved to.
#' 
#' @return render_... functions don't have return values. They generated the required output documents.
#' 
#' @name utils_reporting
NULL


#' @title Render a standard template report for a given object. (Internal function)
#' 
#' @param template_name File name of the standard template
#' @param output_dir Output directory to which to render
#' @param params The Rmd parameter list. Needs to match the expected parameter definitions of the Rmd template.
#' @param output_file Output filename (with or without suffix). Default NULL. If not NULL, it is appended to output_dir
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
render_backward_summary <- function( bwd, output_dir = NULL, output_file = NULL ) {
  render_std_template(  template_name = "rmd_backward_summary.Rmd", 
                        output_dir = output_dir, 
                        params = list( bwd = bwd ), 
                        output_file = output_file, envir = new.env() )
}

#' @rdname utils_reporting
#' @export
render_forward_summary <- function( fwd, output_dir = NULL, output_file = NULL ) {
  render_std_template(  template_name = "rmd_backward_summary.Rmd", 
                        output_dir = output_dir, 
                        params = list( fwd = fwd ), 
                        output_file = output_file, envir = new.env() )
}

#' @rdname utils_reporting
#' @export
render_bidirectional_summary <- function( bds, output_dir = NULL, output_file = NULL ) {
  render_std_template(  template_name = "rmd_bidirectional_summary.Rmd", 
                        output_dir = output_dir, 
                        params = list( bds = bds ), 
                        output_file = output_file, envir = new.env() )
}

#' @rdname utils_reporting
#' @export
render_lrsearch_summary <- function( lrs, output_dir = NULL, output_file = NULL ) {
  render_std_template(  template_name = "rmd_lrsearch_summary.Rmd", 
                        output_dir = output_dir, 
                        params = list( lrs = lrs ), 
                        output_file = output_file, envir = new.env() )
}

#' @rdname utils_reporting
#' @export
render_random_summary <- function( rnd, output_dir = NULL, output_file = NULL ) {
  render_std_template(  template_name = "rmd_random_summary.Rmd", 
                        output_dir = output_dir, 
                        params = list( rnd = rnd ), 
                        output_file = output_file, envir = new.env() )
}

#' @rdname utils_reporting
#' @export
render_game_rank_summary <- function( gmr, output_dir = NULL, output_file = NULL ) {
  render_std_template(  template_name = "rmd_game_rank_summary.Rmd", 
                        output_dir = output_dir, 
                        params = list( gmr = gmr ), 
                        output_file = output_file, envir = new.env() )
}

#' @rdname utils_reporting
#' @export
render_variable_checks_summary <- function( vck, output_dir = NULL, output_file = NULL ) {
  render_std_template(  template_name = "rmd_variable_checks.Rmd", 
                        output_dir = output_dir, 
                        params = list( vck = vck ), 
                        output_file = output_file, envir = new.env() )
}

#' @rdname utils_reporting
#' @export
render_model_calibration_normal <- function( ds, dat, resp, selection, k, 
                                             output_dir = NULL, output_file = NULL ) {
  render_std_template(  template_name = "rmd_calibration_normal.Rmd", 
                        output_dir = output_dir, 
                        params = list( ds = ds, dat = dat, resp = resp, selection = selection, k = k ), 
                        output_file = output_file, envir = new.env() )
}


#' @rdname utils_reporting
#' @export
render_model_calibration_binomial <- function( ds, dat, resp, selection, k, 
                                               output_dir = NULL, output_file = NULL ) {
  render_std_template(  template_name = "rmd_calibration_binomial.Rmd", 
                        output_dir = output_dir, 
                        params = list( ds = ds, dat = dat, resp = resp, selection = selection, k = k ), 
                        output_file = output_file, envir = new.env() )
}


#' @rdname utils_reporting
#' @export
render_model_calibration_cox <- function( ds, dat, resp, selection, k, u,
                                          output_dir = NULL, output_file = NULL ) {
  render_std_template(  template_name = "rmd_calibration_cox.Rmd", 
                        output_dir = output_dir, 
                        params = list( ds = ds, dat = dat, resp = resp, selection = selection, k = k, u = u ), 
                        output_file = output_file, envir = new.env() )
}

render_influence_summary <- function( tinf, output_dir = NULL, output_file = NULL ) {
  render_std_template(  template_name = "rmd_influence_summary.Rmd", 
                        output_dir = output_dir, 
                        params = list( tinf = tinf ), 
                        output_file = output_file, envir = new.env() )
}

