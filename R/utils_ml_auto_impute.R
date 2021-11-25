#
# Formula rewrite tools to automatically impute missing values by maximum-likelihood estimate
#

saf <- Vectorize( function( x ) ifelse( !is.na(x), x, 0 ) )
psu <- Vectorize( function( x ) ifelse(  is.na(x), 1, 0 ) )

#' @title Helper function for Maximum-Likelihood Estimation of Imputation Values
#' 
#' @description Missing values are dummy-encoded such that maximum-likelihood 
#' estimation via lm or glm models will identify a value maximizing the model 
#' likelihood for missing observations.
#' 
#' @param fo Formula specifying lhs ~ rhs of a model
#' @param dat data.frame or tibble that contains data for the model formula
#' 
#' @return Formula object that includes dummy terms and transforms for variables 
#' with missing data.
#' 
#' @examples 
#' 
#' vars <- grep( "the_|rnd", colnames(toy_data), value=TRUE )
#' resp <- "resp"
#' df <- toy_data[,c(resp,vars[1:8])]
#' df$the_cubed[ sample.int(nrow(df), 5) ] <- NA
#' fo <- as.formula(df)
#' rfo <- rewrite_formula( fo, df )
#' rfo
#' glm( rfo, "binomial", df )
#' # saf(the_cubed) denotes the effect for the non-missing variables
#' # psu(the_cubed) denotes the maximum-likelihood imputation value for the missing data
#' 
#' @export
rewrite_formula <- local( {
  function( fo, dat ) {
    dat_tmp <- stats::model.frame( fo, dat, na.action = stats::na.pass )
    labs <- attr( attr(dat_tmp,'terms'), 'term.labels' )
    
    fo_m <- fo
    for( la in labs ) {
      if( any( is.na(dat_tmp[,la] )  ) ) {
        # cat ( ">> ", la,"\n" )
        sup <- sprintf( ". ~ . - %s + saf(%s) + psu(%s)", la, la, la )
        # cat( "-- ",sup,"\n")
        fo_m <- stats::update( fo_m, stats::formula( sup ) )
      }
    }
    dat_m <- stats::model.frame( fo_m, dat, na.action=stats::na.pass )
    environment(fo_m) <- new.env()
    return( fo_m )
  }
}, envir = list( saf=saf, psu=psu ) )

