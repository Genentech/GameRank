
#' @import ggplot2 ggExtra
#' @importFrom rlang .data

gplot_multimodal_variable <- function( dat, resp, vars, mumo, variable, ... ) {
  stopifnot( is.data.frame(dat) | is_tibble(dat) )
  stopifnot( is.character(resp) )
  stopifnot( is.character(vars)  & 1 < length(vars) )
  
  msel <- purrr::pluck( mumo$transforms, variable )
  tb_aic_raw <- purrr::pluck( msel, "aic_tab" ) 
  tb_aic_agg <- purrr::pluck( msel, "aic_aggregate" ) 
  flx_mod <- purrr::pluck( msel, "best_model" )
  cps <- purrr::pluck( msel, "cut_points" )
  stopifnot( !is.null(msel) )
  stopifnot( !is.null(tb_aic_raw) & !is.null(tb_aic_agg) )
  stopifnot( !is.null(flx_mod) & ("flexmix"==class(flx_mod)) )
  
  gp1 <- ggplot() +
    geom_line( data = tb_aic_agg, 
               mapping = aes( x=k, y=.data$min_aic ), color = "blue" ) + 
    geom_point( data = tb_aic_raw, 
                mapping = aes( x=k, y=.data$aic ) )
  
  gp2 <- ggplot() +
    geom_histogram( data = dat, 
                    mapping = aes_string( x=variable, y="..density.." ), 
                    breaks = bins_ucv( dat[[variable]] ), alpha=0.5 ) +
    geom_density( data = dat, 
                  mapping = aes_string( x=variable, y="..density.." ), 
                  bw = "ucv" ) 
  
  nn <- 500
  xs <- seq( min(dat[[variable]], na.rm=TRUE), 
             max(dat[[variable]], na.rm=TRUE), 
             length.out = nn )
  
  for( k in seq_len( ncol(flexmix::parameters(flx_mod)) ) ) {
    mn <- flexmix::parameters(flx_mod)[1,k]
    sig <- flexmix::parameters(flx_mod)[2,k]
    prio <- flexmix::prior(flx_mod)[k]
    plt <- tibble( xs = xs, prob = dnorm( x = xs, mean = mn, sd = sig ) * prio )
    gp2 <- gp2 + 
      geom_line( data = plt, mapping = aes( x=xs, y=.data$prob ), color = 1+k )
  }
  if( !is.null(cps) ) {
    gp2 <- gp2 + 
      geom_vline( xintercept = cps[which(!is.infinite(cps))], color = "orange" )  
  }
  
  gridExtra::grid.arrange( gp2, gp1, nrow=1 )
}
