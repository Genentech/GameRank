
#' @import ggplot2 ggExtra
#' @importFrom rlang .data

# Plot Normal Calibration
tbl_predictions_normal <- function( dat, resp, selection, mod, ... ) {
  stopifnot( !is.null(mod) )
  ret <- tryCatch({
    mf <- dat
    mf$prd <- predict( mod, newdata = dat, type="response", ... ) 
    mf$obs <- dat[,resp]
    res <- mf[,c("obs","prd")]
    res
  }, error = function( e ) NA )
  return( ret )
}

gplot_predictions_normal <- function( dat, resp, selection, mod, ... ) {
  res <- tbl_predictions_normal(  dat, resp, selection, mod, ... )
  ret <- ggplot( data = res, aes( x=.data$prd, y=.data$obs ) ) +
    geom_point() +
    geom_smooth( method = "loess", se = TRUE, color = "blue" ) +
    # theme_classic() +
    xlim( c(0.0,1.0 ) ) +
    ylim( c(0.0,1.0 ) ) +
    xlab( sprintf( "Predicted values" ) ) +
    ylab( "Observed values" ) 
  ggExtra::ggMarginal( ret, type="densigram", margins = "both" )
}

