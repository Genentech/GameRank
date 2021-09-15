#
# Utils and Analysis functions for variable construction
#

# Simple transformations

#
# Add some simple transformations, following https://rcompanion.org/handbook/I_12.html
#
# Add transformed variable if Normality is improved, as found by increase in Shapiro-Wilk W statistic
#
simple_transforms <- function( dat, vars, 
                               transforms = c("sqrt","cubert","log", "zscore" ),
                               name_pattern = function( varname, transform_name ) sprintf( "%s_%s", varname, transform_name )
                               )
{
  stopifnot( is.data.frame(dat) | is_tibble(dat) )
  cat( "Adding simple transformed variables if they are better by Shapiro-Wilk statistic W (larger is better) \n")
  ret <- dat
  lst <- list()
  for( var in vars ) {
    if( is.numeric(dat[[var]]) ) {
      w <- shapiro.test( x = as.numeric(dat[[var]]) )$statistic
      cat( sprintf( "Evaluating %s with W = %1.4f \n", var, w ))
      var_label(ret[,var]) <- sprintf( "original variable: identity (W=%1.4f)", w )
      lst[[length(lst)+1]] <- list( variable = var, transformed_var = "", term = sprintf( "( %s )", var ), W = w, transform = "identity" )
      
      # sqrt transforms
      if( "sqrt" %in% transforms ) {
        nwname <- name_pattern( var, "sqrt" )
        stopifnot( !(nwname %in% colnames(ret) ) )
        ret[,nwname] <- sqrt( ret[,var] )
        wn <- shapiro.test(ret[[nwname]])$statistic
        if( is.nan(wn) || (wn < w) ) {
          ret[,nwname] <- NULL # Delete as it is not better
        } else {
          var_label(ret[,nwname]) <- sprintf( "simple_transforms: sqrt (W=%1.4f)", wn )
          lst[[length(lst)+1]] <- list( variable = var, transformed_var = nwname, term = sprintf( "sqrt( %s )", var ), W = wn, transform = "sqrt" )
        }
      }
      # cube root transforms
      if( "cubert" %in% transforms ) {
        nwname <- name_pattern( var, "cubert" )
        stopifnot( !(nwname %in% colnames(ret) ) )
        ret[,nwname] <- ( ret[,var] )^(1/3)
        wn <- shapiro.test(ret[[nwname]])$statistic
        if( is.nan(wn) || (wn < w) ) {
          ret[,nwname] <- NULL # Delete as it is not better
        } else {
          var_label(ret[,nwname]) <- sprintf( "simple_transforms: cubert (W=%1.4f)", wn )
          lst[[length(lst)+1]] <- list( variable = var, transformed_var = nwname, term = sprintf( "( %s )^(1/3)", var ), W = wn, transform = "cube root" )
        }
      }
      # log transforms
      if( "log" %in% transforms ) {
        nwname <- name_pattern( var, "log" )
        stopifnot( !(nwname %in% colnames(ret) ) )
        ret[,nwname] <- log( ret[,var] )
        wn <- shapiro.test(ret[[nwname]])$statistic
        if( is.nan(wn) || (wn < w) ) {
          ret[,nwname] <- NULL # Delete as it is not better
        } else {
          var_label(ret[,nwname]) <- sprintf( "simple_transforms: log (W=%1.4f)", wn )
          lst[[length(lst)+1]] <- list( variable = var, transformed_var = nwname, term = sprintf( " log( %s ) ", var ), W = wn, transform = "log" )
        }
      }
      # z-score transform
      if( "zscore" %in% transforms ) {
        nwname <- name_pattern( var, "zscore" )
        stopifnot( !(nwname %in% colnames(ret) ) )
        mn <- mean( ret[[var]], na.rm=TRUE )
        st <- sd( ret[[var]], na.rm=TRUE )
        ret[,nwname] <- qnorm( p = ret[[var]], mean = mn, sd = st )
        var_label(ret[,nwname]) <- sprintf( "simple_transforms: zscore (mean=%1.4f, sd=%1.4f)", mn, st )
        lst[[length(lst)+1]] <- list( variable = var, transformed_var = nwname, term = sprintf( "qnorm( %s, mean=%1.4f, sd=%1.4f )", var, mn, st ), W = NA, transform = "zscore" )
      }
    }
  }
  return( list( data = ret, transformations = lst ) )
}

# Box-Cox transformation for regression

# Box-Cox transformation for binomial regression

# Box-Cox transformation for survival regression


source( "R/utils_histogram.R")

# Density estimators

# Mixture models

