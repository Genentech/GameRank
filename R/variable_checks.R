


check_variable <- function( dat, var ) {
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
  ret[["nmiss_pct"]] <- p * 100.0
  
  cats_nmiss <- c(0.0, 0.7, 0.8, 0.9, 0.99, 1.0)
  names(cats_nmiss) <- c("Drop","Drop","Bad","Try","Good","Perfect")
  lbl <- as.character( cut(p,breaks = cats_nmiss, labels = names(cats_nmiss)[-1], include.lowest = TRUE ) )
  ret[["check_missing"]] <- lbl
  
  x <- xval[ which(!is.na(xval)) ]
  
  # Compute the entropy of each type of variable
  # Entropy H(x) = sum( p(x) * log( p(x) ) )
  epsi <- 1E-3
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
  }
  
  ret <- as_tibble( ret )
  return( ret )
}

check_variables <- function( dat, resp, vars ) {
  
  ret <- list()
  if( !is.null(resp) ) {
    evl <- check_variable( dat, resp )
    evl$is_response <- TRUE
    ret[[resp]] <- evl
  }
  
  for( var in vars ) {
    evl <- check_variable( dat, var )
    evl$is_response <- FALSE
    ret[[var]] <- evl
  }
  
  ret <- Reduce( bind_rows, ret, NULL )
  return( ret )  
}
