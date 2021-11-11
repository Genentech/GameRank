# Plot Survival Calibration
tbl_predictions_cox <- function( dat, resp, selection, mod, u, ... ) {
  stopifnot( !is.null(mod) )
  ret <- tryCatch({
    mf <- dat
    mf$prd <- as.numeric( 1 - pec::predictSurvProb( mod, newdata = dat, times = u, ... ) )
    mf$cll.prd <- log( -log(1-mf$prd) ) 
    cal.cox <- coxph( formula( sprintf( "%s ~ rms::rcs(cll.prd,3)", resp )), data = mf, x=TRUE )
    # cal.cox <- coxph( formula( sprintf( "%s ~ rms::lsp(cll.prd,3)", resp )), data = bind_cols( dat, tibble( cll.prd = cll.prd )), x=TRUE )
    # grd.cox <- seq( quantile( mf$prd, probs = 0.01, na.rm=TRUE ), 
    #                 quantile( mf$prd, probs = 0.99, na.rm=TRUE ), 
    #                 length = 100 )
    # grd.prd.cll <- log(-log(1-grd.cox) )
    # df.grd.cox <- data.frame( grd.cox, grd.prd.cll ) %>% setNames(c("prd","cll.prd"))
    # df.grd.cox$prd.cal <- 1 - predictSurvProb( cal.cox, df.grd.cox, times = u )
    mf$obs <- as.numeric( 1 - pec::predictSurvProb( cal.cox, newdata = mf, times = u ) )
    res <- model.frame( formula( sprintf("%s ~ %s + %s + %s", resp, "prd", "cll.prd", "obs") ), mf )
    res
  }, error = function( e ) NA )
  return( ret )
}

gplot_predictions_cox <- function( dat, resp, selection, mod, u, ... ) {
  res <- tbl_predictions_cox(  dat, resp, selection, mod, u, ... )
  lim <- max(res$prd,res$obs, na.rm=TRUE)
  co <- cor.test( x = res$prd, y = res$obs, method = "pearson" )
  epsi <- 0.0001
  msg <- sprintf( "Pearson Correlation %1.4f (%1.4f, %1.4f; %s)", co$estimate, co$conf.int[1], co$conf.int[2], ifelse( epsi < co$p.value, sprintf("p=%1.4f",co$p.value), "p<.0001" ) )
  ret <- res %>%
    ggplot( aes( x=prd, y=obs ) ) +
    geom_abline( slope = 1, intercept = 0, color = "gray" ) +
    geom_point() +
    geom_smooth( method = "loess", se = TRUE, color = "blue" ) +
    # theme_classic() +
    xlim( c(0,lim ) ) +
    ylim( c(0,lim ) ) +
    xlab( sprintf( "Predicted Survival Probability \n %s", msg ) ) +
    ylab( "Observed Survival Probability" ) 
  ggExtra::ggMarginal( ret, type="densigram" )
}


gplot_km_cox <- function( dat, resp, selection, mod, u, cutpoint = NULL, ... ) {
  plt <- dat
  plt$prd <- as.numeric( 1 - pec::predictSurvProb( mod, newdata = dat, times = u, ... ) )
  if( is.null(cutpoint) ) cutpoint <- median(plt$prd, na.rm=TRUE)
  plt$cut <- factor( cut( plt$prd, breaks = c(-Inf,cutpoint,+Inf), labels = c("Low","High"), include.lowest = TRUE  ) )
  fo <- formula( sprintf( "%s ~ cut", resp ) )
  fit <- survminer::surv_fit( fo, plt )
  ggsurvplot( fit, plt, conf.int = TRUE, risk.table = TRUE )
}
