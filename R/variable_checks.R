
#' @import survival entropy
#' @importFrom rlang .data

#' @title Utility function that performs high-level checks for a variable in a dataset
#' 
#' @param dat data.frame or tibble comprising data for model generation and 
#' validation.
#' @param var Character indicating the variable to check.
#' @param min_cases Minimal number of cases for which the entropy is estimated or
#' outliers are evaluated.
#' @param c_out Constant for the Robust Outlier Test determining the scaling 
#' of IQR to determine non-outlier range from Q1 - c_out x IQR to Q3 + c_out x IQR.
#' @param resp_cat Categorical array with response variable to estimate 
#' mutation information in bits.
#' 
#' @return A tibble with entries
#' \describe{
#' \item{N}{Overall number of values}
#' \item{n}{Overall number of non-missing values}
#' \item{nmiss}{Number of missing values}
#' \item{p}{Percentage of non-missing values.}
#' \item{check_missing}{Categorization of variable into to 'Drop','Bad', 'Try', 
#' 'Good', 'Perfect' depending on the percentage of 
#' missingness <70\%, <80\%, <90\%, <99\%, or 100\%}
#' \item{type}{Type of the variable binary, categorical, integer, or real}
#' \item{entropy}{Entropy estimated, eventually by package entropy function, for
#' the variable. This value describes the information content in bits of the 
#' variable distribution. The more the better.}
#' \item{mutual_information}{Mutual information between variable and response 
#' (if provided) in bits. For survival endpoints, the event flag (Yes/No) is 
#' used.}
#' \item{check_entropy}{Variables with entropy close to 0.0 (<0.001) are 
#' recommended to be dropped.}
#' \item{rot.nmin}{Number of samples lower than Robust Outlier Test cut-off
#'  of Q1 - c_out * IQR.}
#' \item{rot.nmax}{Number of samples lower than Robust Outlier Test cut-off 
#' of Q3 + c_out * IQR.}
#' \item{rot.p}{Percentage of outliers.}
#' \item{rng_sd}{Ratio of sample range by sample standard deviation.}
#' }
check_variable <- function( dat, var, min_cases = 25L, 
                            c_out = 1.5, resp_cat = NULL ) {
  message( sprintf( "check_variable: Evaluating variable %s.", var ) )
  ret <- list()
  ret[["variable"]] <- var
  
  # xval <- dat[[var]]
  xval <- stats::model.frame(stats::formula(sprintf("%s ~ 1",var)), dat, 
                             na.action = stats::na.pass)[[1]]
  N <- length( xval )
  n <- length( which(!is.na(xval)) )
  p <- n / N

  ret[["N"]] <- N
  ret[["n"]] <- n
  ret[["nmiss"]] <- N - n
  ret[["p"]] <- p * 100.0
  
  cats_nmiss <- c(0.0, 0.7, 0.8, 0.9, 0.99, 1.0)
  names(cats_nmiss) <- c("Drop","Drop","Bad","Try","Good","Perfect")
  lbl <- as.character( cut(p,breaks = cats_nmiss, 
                           labels = names(cats_nmiss)[-1], 
                           include.lowest = TRUE ) )
  ret[["check_missing"]] <- lbl
  
  x <- xval[ which(!(is.na(xval) | is.infinite(xval) | is.nan(xval) ) ) ]
  
  # Compute the entropy of each type of variable
  # Entropy H(x) = sum( p(x) * log( p(x) ) )
  epsi <- 1E-3
  if( min_cases < length( x ) ) { # Only check entropy if there are enough non-missing values
    if( is.Surv(x) ) {
      # Do nothing
    } else if( is.logical(x) )  {
      
      # x <- c( TRUE, TRUE, FALSE, FALSE, FALSE )
      # x <- c( TRUE, TRUE,TRUE, TRUE )
      en <- entropy::entropy( y =  table( x ), method = "ML", unit = "log2" )
      mi <- NA_real_
      if( !is.null(resp_cat) ) {
        mi <- entropy::mi.empirical( y2d = table( resp_cat, xval ), unit = "log2"  )  
      }
      
      # tab <- table( x )
      # ptab <- prop.table( tab )
      # 
      # en <- -sum( ptab * log( ptab ) )
      # en
      
      ret[["type"]] <- "binary"
      ret[["entropy"]] <- en
      ret[["mutual_information"]] <- mi
      ret[["check_entropy"]] <- ifelse( en < epsi, "Entropy too low", "Entropy ok" )
      
    } else if( is.factor(x) | is.character(x) ) {
      
      # x <- c( "TRUE", "TRUE", "FALSE", "FALSE", "FALSE" )
      # x <- c( rep.int( "TRUE", 100 ), "FALSE" )
      # x <- x %>% factor( levels = c("FALSE","TRUE") )
      en <- entropy::entropy( y = table( x ), method = "ML", unit = "log2" )
      mi <- NA_real_
      if( !is.null(resp_cat) ) {
        mi <- entropy::mi.empirical( y2d = table( resp_cat, xval ), unit = "log2"  )  
      }
      
      # tab <- table( x )
      # ptab <- prop.table( tab )
      # 
      # en <- -sum( ptab * log( ptab ) )
      # en
      
      ret[["type"]] <- "categorical"
      ret[["entropy"]] <- en
      ret[["mutual_information"]] <- mi
      ret[["check_entropy"]] <- ifelse( en < epsi, "Entropy too low", "Entropy ok" )
      
    } else if( is.integer(x) ) {
      
      # x <- as.integer( rnorm(100) * 100 )
      # x <- rep.int( 1L, 100 )
      en <- entropy::entropy( y = table( x ), method = "ML", unit = "log2" )
      mi <- NA_real_
      if( !is.null(resp_cat) ) {
        mi <- entropy::mi.empirical( y2d = table( resp_cat, xval ), unit = "log2"  )  
      }
      
      # tab <- table( x )
      # ptab <- prop.table( tab )
      # 
      # en <- -sum( ptab * log( ptab ) )
      # en
      
      ret[["type"]] <- "integer"
      ret[["entropy"]] <- en
      ret[["mutual_information"]] <- mi
      ret[["check_entropy"]] <- ifelse( en < epsi, "Entropy too low", "Entropy ok" )
      
      
      # Evaluate for outliers
      iqr <- stats::IQR( x, na.rm=TRUE )
      q1 <- stats::quantile( x, prob = 0.25, na.rm=TRUE )
      q3 <- stats::quantile( x, prob = 0.75, na.rm=TRUE )
      rot.nmin <- length( which( x < q1 - c_out * iqr ) )
      rot.nmax <- length( which( x > q3 + c_out * iqr ) )
      rot.p <- (rot.nmin+rot.nmax) / n * 100.0
      rng_sd <- abs( max(x) - min(x) ) / stats::sd( x )
      # tst <- outliers::dixon.test( x=x, type=22 )
      # tst_stat <- tst$statistic
      # tst_pval <- tst$p.value
      ret[["rot.nmin"]] <- rot.nmin
      ret[["rot.nmax"]] <- rot.nmax
      ret[["rot.p"]] <- rot.p
      ret[["rng_sd"]] <- rng_sd
      # ret[["dixon.Q"]] <- tst_stat
      # ret[["dixon.pval"]] <- tst_pval
    } else if( is.double(x) ) {
      
      # x <- rnorm( 100 )
      # x <- rep.int( 1.0, 100 )
      
      # In the continuous variable case, we first need to determine the 
      # optimal binwidth by Leave-One-Out cross-validation
      # before we can use it to estimate the entropy.
      en <- NA_real_
      en <- tryCatch({
        bb <- bins_ucv( x )
        bb <- sort(unique(bb))
        if( length(bb)<2 ) bb <- c( bb-1E-4, bb, bb+1E-4 )
        tab <- as.integer( table( cut( x, breaks = bb, include.lowest = TRUE ) ) )
        en <- entropy::entropy( y = tab, method = "ML", unit = "log2" )
        en  
      }, error = function(ee) NA_real_ )
      
      mi <- NA_real_
      if( !is.null(resp_cat) ) {
        mi <- tryCatch({
          # CHE/2023-03-19: Adding tryCatch to make more robust
          entropy::mi.empirical( y2d = table( resp_cat, 
                                              cut( xval, breaks = bb, 
                                                   include.lowest = TRUE ) ),
                                 unit = "log2"  )  
          }, error = function( ee ) NA_real_ )
      }
      
      ret[["type"]] <- "real"
      ret[["entropy"]] <- en
      ret[["mutual_information"]] <- mi
      ret[["check_entropy"]] <- ifelse( en < epsi, "Entropy too low", "Entropy ok" )
      
      # Evaluate for outliers
      iqr <- stats::IQR( x, na.rm=TRUE )
      q1 <- stats::quantile( x, prob = 0.25, na.rm=TRUE )
      q3 <- stats::quantile( x, prob = 0.75, na.rm=TRUE )
      rot.nmin <- length( which( x < q1 - c_out * iqr ) )
      rot.nmax <- length( which( x > q3 + c_out * iqr ) )
      rot.p <- (rot.nmin+rot.nmax) / n * 100.0
      rng_sd <- abs( max(x) - min(x) ) / stats::sd( x )
      # tst <- outliers::dixon.test( x=x, type=22 )
      # tst_stat <- tst$statistic
      # tst_pval <- tst$p.value
      ret[["rot.nmin"]] <- rot.nmin
      ret[["rot.nmax"]] <- rot.nmax
      ret[["rot.p"]] <- rot.p
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
#' @description Applies check_variable to each response and each variable in 
#' vars, and returns a merged dataset.
#' 
#' @param dat data.frame or tibble comprising data for model generation and 
#' validation.
#' @param resp Character string defining the response (lhs) of the model formula.
#' @param vars Character indicating the variables to check.
#' @param min_cases Minimal number of cases for which the entropy is estimated 
#' or outliers are evaluated.
#' @param c_out Constant for the Robust Outlier Test determining the scaling of 
#' IQR to determine non-outlier range from Q1 - c_out x IQR to Q3 + c_out x IQR.
#' 
#' @return A tibble with entries
#' \describe{
#' \item{N}{Overall number of values}
#' \item{n}{Overall number of non-missing values}
#' \item{nmiss}{Number of missing values}
#' \item{p}{Percentage of non-missing values.}
#' \item{check_missing}{Categorization of variable into to 'Drop','Bad', 'Try', 
#' 'Good', 'Perfect' depending on the percentage of 
#' missingness <70\%, <80\%, <90\%, <99\%, or 100\%}
#' \item{type}{Type of the variable binary, categorical, integer, or real}
#' \item{entropy}{Entropy estimated, eventually by package entropy function, 
#' for the variable. This value describes the information content in bits of 
#' the variable distribution. The more the better.}
#' \item{mutual_information}{Mutual information between variable and response 
#' (if provided) in bits. For survival endpoints, the event flag (Yes/No) 
#' is used.}
#' \item{check_entropy}{Variables with entropy close to 0.0 (<0.001) are 
#' recommended to be dropped.}
#' \item{rot.nmin}{Number of samples lower than Robust Outlier Test cut-off 
#' of Q1 - c_out x IQR.}
#' \item{rot.nmax}{Number of samples lower than Robust Outlier Test cut-off 
#' of Q3 + c_out x IQR.}
#' \item{rng_sd}{Ratio of sample range by sample standard deviation.}
#' }
#' 
#' @examples 
#' library( ggplot2 )
#' resp <- "resp"
#' vars <- grep( "the_|rnd", colnames(toy_data), value=TRUE )
#' 
#' vck <- check_variables( toy_data, resp, vars )
#' vck %>% summary
#' vck %>% filter( !is_response ) %>% arrange( desc(entropy) )
#' 
#' vck %>% 
#'   ggplot(aes(x=entropy, y=mutual_information) ) +
#'   geom_point()
#' @export 
check_variables <- function( dat, resp, vars, min_cases = 25L, c_out = 1.5 ) {
  message( sprintf( "check_variables: Obtaining variable screening information for %d variables (min_cases=%d, outlier c=%1.2f).",
                    length( c(resp,vars) ), min_cases, c_out ))
  resp_cat <- NULL
  ret <- list()
  if( !is.null(resp) ) {
    evl <- check_variable( dat, resp, min_cases = min_cases, 
                           c_out = c_out, resp_cat = NULL )
    evl$is_response <- TRUE
    ret[[resp]] <- evl
    
    fo <- stats::formula( sprintf( "%s ~ 1", resp ) )
    mf <- stats::model.frame( fo, dat, na.action = stats::na.pass )
    
    if( is.Surv( mf[[1]] ) ) {
      # TODO For now take Event/No-Event as category
      mf$resp_cat <- as.factor( mf[[1]][,"status"] )
    } else if( is.integer(mf[[1]]) | is.double(mf[[1]]) ) {
      bb <- bins_ucv( mf[[1]] )
      mf$resp_cat = factor( cut( mf[[1]], breaks = bb, include.lowest = TRUE ) ) 
    } else if( is.character(mf[[1]]) | is.factor(mf[[1]]) | is.logical(mf[[1]]) ) {
      mf$resp_cat <- as.factor( mf[[1]] )
    } 
    resp_cat <- mf$resp_cat
    stopifnot( length(resp_cat)==nrow(dat) )
  }
  
  for( var in vars ) {
    evl <- check_variable( dat, var, min_cases = min_cases,
                           c_out = c_out, resp_cat = resp_cat )
    evl$is_response <- FALSE
    ret[[var]] <- evl
  }
  
  ret <- Reduce( bind_rows, ret, NULL )
  ret <- ret %>% 
    mutate( check_missing = factor( .data$check_missing, 
                                    levels = c("Drop","Bad","Try","Good","Perfect") ), 
            check_entropy = factor( .data$check_entropy, 
                                    levels = c("Entropy too low", "Entropy ok") ),
            type = factor( ifelse( is.na(.data$type), "Entropy not done", .data$type ), 
                           levels = c("Entropy not done","real", "integer", "categorical", "binary" ) ) )
  return( ret )  
}
