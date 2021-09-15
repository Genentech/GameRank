#
# Formula rewrite tools to automatically impute missing values by maximum-likelihood estimate
#

saf <- Vectorize( function( x ) ifelse( !is.na(x), x, 0 ) )
psu <- Vectorize( function( x ) ifelse(  is.na(x), 1, 0 ) )

rewrite_formula <- local( {
  function( fo, dat ) {
    dat_tmp <- model.frame( fo, dat, na.action = na.pass )
    labs <- attr( attr(dat_tmp,'terms'), 'term.labels' )
    
    fo_m <- fo
    for( la in labs ) {
      if( any( is.na(dat_tmp[,la] )  ) ) {
        cat ( ">> ", la,"\n" )
        sup <- sprintf( ". ~ . - %s + saf(%s) + psu(%s)", la, la, la )
        cat( "-- ",sup,"\n")
        fo_m <- update( fo_m, formula( sup ) )
      }
    }
    dat_m <- model.frame( fo_m, dat, na.action=na.pass )
    environment(fo_m) <- new.env()
    return( fo_m )
  }
}, envir = list( saf=saf, psu=psu ) )

