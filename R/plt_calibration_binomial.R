# Plot Binomial Calibration
tbl_predictions_binomial <- function( dat, resp, selection, mod, ... ) {
  stopifnot( !is.null(mod) )
  ret <- tryCatch({
    mf <- dat
    mf$prd <- predict( mod, newdata = dat, type="response", ... ) 
    mf$obs <- as.logical( dat[,resp] )
    res <- mf[,c("obs","prd")]
    res
  }, error = function( e ) NA )
  return( ret )
}

gplot_predictions_binomial <- function( dat, resp, selection, mod, ... ) {
  res <- tbl_predictions_binomial(  dat, resp, selection, mod, ... )
  ret <- res %>%
    ggplot( aes( x=prd, y=obs ) ) +
    geom_point() +
    geom_smooth( method = "loess", se = TRUE, color = "blue" ) +
    # theme_classic() +
    xlim( c(0.0,1.0 ) ) +
    ylim( c(0.0,1.0 ) ) +
    xlab( sprintf( "Predicted Probability" ) ) +
    ylab( "Observed Probability" ) 
  ggExtra::ggMarginal( ret, type="densigram", margins = "x" )
}

