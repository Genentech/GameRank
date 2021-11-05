

#' @title Utility function that performs high-level checks for a variable in a dataset
#' 
#' @param dat data.frame or tibble comprising data for model generation and validation.
#' @param var Character indicating the variable to check.
#' @param min_cases Minimal number of cases for which the entropy is estimated or outlier tests are performed.
#' @param c_out Constant for the Robust Outlier Test determining the scaling of IQR to determine non-outlier range from Q1 - c_out x IQR to Q3 + c_out x IQR.
#' 
#' @return A tibble with entries
#' \describe{
#' \item{N}{Overall number of values}
#' \item{n}{Overall number of non-missing values}
#' \item{nmiss}{Number of missing values}
#' \item{p}{Percentage of non-missing values.}
#' \item{check_missing}{Categorization of variable into to 'Drop','Bad', 'Try', 'Good', 'Perfect' depending on the percentage of missingness <70%, <80%, <90%, <99%, or 100%}
#' \item{type}{Type of the variable binary, categorical, integer, or real}
#' \item{entropy}{Entropy estimated, eventually by package entropy function, for the variable. This value describes the information content of the variable distribution. The more the better.}
#' \item{check_entropy}{Variables with entropy close to 0.0 (<0.001) are recommended to be dropped.}
#' \item{rot.nmin}{Number of samples lower than Robust Outlier Test cut-off of Q1 - c_out * IQR.}
#' \item{rot.nmax}{Number of samples lower than Robust Outlier Test cut-off of Q3 + c_out * IQR.}
#' \item{rng_sd}{Ratio of sample range by sample standard deviation.}
# \item{dixon.Q}{Dixons Q-statistics.}
# \item{dixon.pval}{Dixons p-value that there is an outlier in the sample.}
#' }
#' 
#' @export 
check_variable <- function( dat, var, min_cases = 35L, c_out = 1.5 ) {
  cat( sprintf( "Evaluating variable %s \n", var ) )
  ret <- list()
  ret[["variable"]] <- var
  
  xval <- dat[[var]]
  N <- length( xval )
  n <- length( which(!is.na(xval)) )
  p <- n / N

  ret[["N"]] <- N
  ret[["n"]] <- n
  ret[["nmiss"]] <- N - n
  ret[["p"]] <- p * 100.0
  
  cats_nmiss <- c(0.0, 0.7, 0.8, 0.9, 0.99, 1.0)
  names(cats_nmiss) <- c("Drop","Drop","Bad","Try","Good","Perfect")
  lbl <- as.character( cut(p,breaks = cats_nmiss, labels = names(cats_nmiss)[-1], include.lowest = TRUE ) )
  ret[["check_missing"]] <- lbl
  
  x <- xval[ which(!is.na(xval)) ]
  
  # Compute the entropy of each type of variable
  # Entropy H(x) = sum( p(x) * log( p(x) ) )
  epsi <- 1E-3
  if( min_cases < length( x ) ) { # Only check entropy if there are enough non-missing values
    if( is.logical(x) )  {
      
      # x <- c( TRUE, TRUE, FALSE, FALSE, FALSE )
      # x <- c( TRUE, TRUE,TRUE, TRUE )
      tab <- table( x )
      en <- entropy::entropy( y = tab, method = "ML" )
      
      # tab <- table( x )
      # ptab <- prop.table( tab )
      # 
      # en <- -sum( ptab * log( ptab ) )
      # en
      
      ret[["type"]] <- "binary"
      ret[["entropy"]] <- en
      ret[["check_entropy"]] <- ifelse( entropy < epsi, "Entropy too low", "Entropy ok" )
      
    } else if( is.factor(x) | is.character(x) ) {
      
      # x <- c( "TRUE", "TRUE", "FALSE", "FALSE", "FALSE" )
      # x <- c( rep.int( "TRUE", 100 ), "FALSE" )
      # x <- x %>% factor( levels = c("FALSE","TRUE") )
      tab <- table( x )
      en <- entropy::entropy( y = tab, method = "ML" )
      
      # tab <- table( x )
      # ptab <- prop.table( tab )
      # 
      # en <- -sum( ptab * log( ptab ) )
      # en
      
      ret[["type"]] <- "categorical"
      ret[["entropy"]] <- en
      ret[["check_entropy"]] <- ifelse( en < epsi, "Entropy too low", "Entropy ok" )
      
    } else if( is.integer(x) ) {
      
      # x <- as.integer( rnorm(100) * 100 )
      # x <- rep.int( 1L, 100 )
      tab <- table( x )
      en <- entropy::entropy( y = tab, method = "ML" )
      
      # tab <- table( x )
      # ptab <- prop.table( tab )
      # 
      # en <- -sum( ptab * log( ptab ) )
      # en
      
      ret[["type"]] <- "integer"
      ret[["entropy"]] <- en
      ret[["check_entropy"]] <- ifelse( en < epsi, "Entropy too low", "Entropy ok" )
      
      
      # Evaluate for outliers
      iqr <- IQR( x, na.rm=TRUE )
      q1 <- quantile( x, prob = 0.25, na.rm=TRUE )
      q3 <- quantile( x, prob = 0.75, na.rm=TRUE )
      rot.nmin <- length( which( x < q1 - c_out * iqr ) )
      rot.nmax <- length( which( x > q3 + c_out * iqr ) )
      rng_sd <- range( x ) / sd( x )
      # tst <- outliers::dixon.test( x=x, type=22 )
      # tst_stat <- tst$statistic
      # tst_pval <- tst$p.value
      ret[["rot.nmin"]] <- rot.nmin
      ret[["rot.nmax"]] <- rot.nmax
      ret[["rng_sd"]] <- rng_sd
      # ret[["dixon.Q"]] <- tst_stat
      # ret[["dixon.pval"]] <- tst_pval
    } else if( is.double(x) ) {
      
      # x <- rnorm( 100 )
      # x <- rep.int( 1.0, 100 )
      
      # In the continuous variable case, we first need to determine the optimal binwidht by Leave-One-Out cross-validation
      # before we can use it to estimate the entropy.
      bb <- bins_ucv( x )
      tab <- as.integer( table( cut( x, breaks = bb, include.lowest = TRUE ) ) )
      en <- entropy::entropy( y = tab, method = "ML" )
      
      # Attempt 1:
      # n_points <- 512 # Default
      # bw <- bw.ucv( x = x )
      # bw
      # dens <- density( x = x, bw = bw, kernel = "gaussian", n = n_points )
      # fdens <- approxfun( x = dens$x, y = dens$y, yleft = 0.0, yright = 0.0 )
      # # ii_val <- integrate( f = fdens, lower = -Inf, upper = +Inf )
      # # ii_val # We could adjust n_points depending on how close ii_val is to 1.0
      # fentropy <- function( x ) fdens(x) * log( fdens(x) )
      # n_divs <- c(100, 250, 500, 1000, 1500, 2000 )
      # ii <- NULL
      # k <- 1
      # while( is.null(ii) && k <= length(n_divs) ) {
      #   ii <- tryCatch({ integrate( f = fentropy, lower = min(x), upper = max(x), subdivisions = n_divs[k] ) },
      #                  error = function(e) NULL )
      #   k <- k + 1
      # }
      # if( !is.null(ii) ) {
      #   entropy <- -ii$value
      # } else {
      #   entropy <- NA
      # }
      # entropy
      
      ret[["type"]] <- "real"
      ret[["entropy"]] <- en
      ret[["check_entropy"]] <- ifelse( en < epsi, "Entropy too low", "Entropy ok" )
      
      # Evaluate for outliers
      iqr <- IQR( x, na.rm=TRUE )
      q1 <- quantile( x, prob = 0.25, na.rm=TRUE )
      q3 <- quantile( x, prob = 0.75, na.rm=TRUE )
      rot.nmin <- length( which( x < q1 - c_out * iqr ) )
      rot.nmax <- length( which( x > q3 + c_out * iqr ) )
      rng_sd <- range( x ) / sd( x )
      # tst <- outliers::dixon.test( x=x, type=22 )
      # tst_stat <- tst$statistic
      # tst_pval <- tst$p.value
      ret[["rot.nmin"]] <- rot.nmin
      ret[["rot.nmax"]] <- rot.nmax
      ret[["rng_sd"]] <- rng_sd
      # ret[["dixon.Q"]] <- tst_stat
      # ret[["dixon.pval"]] <- tst_pval
    }
  }

  ret <- as_tibble( ret )
  return( ret )
}

#' @title Utility function that applies check_variable to a dataset
#' 
#' @description Applies check_variable to each response and each variable in vars, and returns a merged dataset.
#' 
#' @param dat data.frame or tibble comprising data for model generation and validation.
#' @param resp Character string defining the response (lhs) of the model formula.
#' @param vars Character indicating the variables to check.
#' @param min_cases Minimal number of cases for which the entropy is estimated or outlier tests are performed.
#' @param c_out Constant for the Robust Outlier Test determining the scaling of IQR to determine non-outlier range from Q1 - c_out x IQR to Q3 + c_out x IQR.
#' 
#' @return A tibble with entries
#' \describe{
#' \item{N}{Overall number of values}
#' \item{n}{Overall number of non-missing values}
#' \item{nmiss}{Number of missing values}
#' \item{p}{Percentage of non-missing values.}
#' \item{check_missing}{Categorization of variable into to 'Drop','Bad', 'Try', 'Good', 'Perfect' depending on the percentage of missingness <70%, <80%, <90%, <99%, or 100%}
#' \item{type}{Type of the variable binary, categorical, integer, or real}
#' \item{entropy}{Entropy estimated, eventually by package entropy function, for the variable. This value describes the information content of the variable distribution. The more the better.}
#' \item{check_entropy}{Variables with entropy close to 0.0 (<0.001) are recommended to be dropped.}
#' \item{rot.nmin}{Number of samples lower than Robust Outlier Test cut-off of Q1 - c_out x IQR.}
#' \item{rot.nmax}{Number of samples lower than Robust Outlier Test cut-off of Q3 + c_out x IQR.}
#' \item{rng_sd}{Ratio of sample range by sample standard deviation.}
# \item{dixon.Q}{Dixons Q-statistics.}
# \item{dixon.pval}{Dixons p-value that there is an outlier in the sample.}
#' }
#' 
#' @export 
check_variables <- function( dat, resp, vars, min_cases = 35L, c_out = 1.5 ) {
  
  ret <- list()
  if( !is.null(resp) ) {
    evl <- check_variable( dat, resp, min_cases = min_cases, c_out = c_out )
    evl$is_response <- TRUE
    ret[[resp]] <- evl
  }
  
  for( var in vars ) {
    evl <- check_variable( dat, var, min_cases = min_cases, c_out = c_out )
    evl$is_response <- FALSE
    ret[[var]] <- evl
  }
  
  ret <- Reduce( bind_rows, ret, NULL )
  return( ret )  
}
