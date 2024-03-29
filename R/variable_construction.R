#
# Utils and Analysis functions for variable construction
#

#' @import purrr
#' @importFrom stats dnorm
#' @importFrom rlang .data

# Simple transformations ----
#' @title Function to evaluate and add simple variable transformations
#' 
#' @description This function adds some simple transformations, as described in 
#' \url{https://rcompanion.org/handbook/I_12.html} to a dataset for a set of 
#' given variables. However, variables are only added if they improve Normality 
#' as measured by an increase in Shapiro-Wilk W statistic.
#' 
#' @param dat A data.frame or tibble with data for each variable in vars
#' @param vars A character vector of variable names to process.
#' @param transforms A list of transformations to try. These can be of 
#' square root 'sqrt', cube root 'cubert', natural logarithm 'log', and 
#' z-Score to Normal distribution 'zscore'
#' @param name_pattern A function that takes two arguments 'varname' and 
#' 'transform_name' and generates a new variable name under with the 
#' transformed data is stored. Default is 'varname + '_' + transform_name.
#' 
#' @return list with two elements
#' \describe{
#' \item{data}{A dataset comprising additional columns with transformed 
#' variables, if they improve Normality.}
#' \item{transformations}{A list of lists with elements for variable, 
#' transformed variable, a formula term, Shapiro-Wilk W statistic, and the 
#' name of the transform.}
#' } 
#' @examples 
#' vars <- grep( "the_|rnd", colnames(toy_data), value=TRUE )
#' smp <- simple_transforms( toy_data, vars = vars )
#' tfs <- smp$transformations %>% Reduce( bind_rows, ., NULL )
#' tfs %>% group_by( variable ) %>% filter( max(W)==W )
#' tfs %>% pull( transform ) %>% table 
#' @export
simple_transforms <- function( dat, vars, 
                               transforms = c("sqrt","cubert","log", "zscore" ),
                               name_pattern = function( varname, transform_name ) sprintf( "%s_%s", varname, transform_name )
                               )
{
  stopifnot( is.data.frame(dat) | tibble::is_tibble(dat) )
  message( sprintf( "simple_transforms: Adding standard transformed for %d variables if they are better by Shapiro-Wilk statistic W (larger is better)", length(c(vars)) ) )
  ret <- dat
  lst <- list()
  for( var in vars ) {
    if( is.numeric(dat[[var]]) & 3 < length(stats::na.omit(dat[[var]])) ) {
      w <- stats::shapiro.test( x = as.numeric(dat[[var]]) )$statistic
      message( sprintf( "Evaluating %s with W = %1.4f.", var, w ))
      labelled::var_label(ret[,var]) <- sprintf( "original variable: identity (W=%1.4f)", w )
      lst[[length(lst)+1]] <- list( variable = var, transformed_var = "", 
                                    term = sprintf( "( %s )", var ), W = w, 
                                    transform = "identity" )
      
      # sqrt transforms
      if( "sqrt" %in% transforms ) {
        nwname <- name_pattern( var, "sqrt" )
        stopifnot( !(nwname %in% colnames(ret) ) )
        ret[,nwname] <- sqrt( ret[,var] )
        wn <- stats::shapiro.test(ret[[nwname]])$statistic
        if( is.nan(wn) || (wn < w) ) {
          ret[,nwname] <- NULL # Delete as it is not better
        } else {
          labelled::var_label(ret[,nwname]) <- sprintf( "simple_transforms: sqrt (W=%1.4f)", wn )
          lst[[length(lst)+1]] <- list( variable = var, transformed_var = nwname,
                                        term = sprintf( "sqrt( %s )", var ), W = wn, 
                                        transform = "sqrt" )
        }
      }
      # cube root transforms
      if( "cubert" %in% transforms ) {
        nwname <- name_pattern( var, "cubert" )
        stopifnot( !(nwname %in% colnames(ret) ) )
        ret[,nwname] <- ( ret[,var] )^(1/3)
        wn <- stats::shapiro.test(ret[[nwname]])$statistic
        if( is.nan(wn) || (wn < w) ) {
          ret[,nwname] <- NULL # Delete as it is not better
        } else {
          labelled::var_label(ret[,nwname]) <- sprintf( "simple_transforms: cubert (W=%1.4f)", wn )
          lst[[length(lst)+1]] <- list( variable = var, transformed_var = nwname, 
                                        term = sprintf( "( %s )^(1/3)", var ), W = wn, 
                                        transform = "cube root" )
        }
      }
      # log transforms
      if( "log" %in% transforms ) {
        nwname <- name_pattern( var, "log" )
        stopifnot( !(nwname %in% colnames(ret) ) )
        ret[,nwname] <- log( ret[,var] )
        wn <- stats::shapiro.test(ret[[nwname]])$statistic
        if( is.nan(wn) || (wn < w) ) {
          ret[,nwname] <- NULL # Delete as it is not better
        } else {
          labelled::var_label(ret[,nwname]) <- sprintf( "simple_transforms: log (W=%1.4f)", wn )
          lst[[length(lst)+1]] <- list( variable = var, transformed_var = nwname, 
                                        term = sprintf( " log( %s ) ", var ), W = wn, 
                                        transform = "log" )
        }
      }
      # z-score transform
      if( "zscore" %in% transforms ) {
        nwname <- name_pattern( var, "zscore" )
        # browser( condition={nwname=="AGE_zscore"})
        stopifnot( !(nwname %in% colnames(ret) ) )
        mn <- mean( ret[[var]], na.rm=TRUE )
        st <- stats::sd( ret[[var]], na.rm=TRUE )
        ret[,nwname] <- stats::pnorm( q = ret[[var]], mean = mn, sd = st )
        wn <- stats::shapiro.test(ret[[nwname]])$statistic
        if( is.nan(wn) || (wn < w) ) {
          ret[,nwname] <- NULL # Delete as it is not better
        } else {
          labelled::var_label(ret[,nwname]) <- sprintf( "simple_transforms: zscore (mean=%1.4f, sd=%1.4f)", mn, st )
          lst[[length(lst)+1]] <- list( variable = var, transformed_var = nwname,
                                        term = sprintf( "pnorm( %s, mean=%1.4f, sd=%1.4f )", var, mn, st ), W = wn, 
                                        transform = "zscore" )
        }
      }
    }
  }
  return( list( data = ret, transformations = lst ) )
}

# Box-Cox Transformations ----
# Box-Cox transformation for regression

#' @title Function to generate Box-Cox transformations for Regression
#' 
#' @description Applies the Box-Cox transformation for set of variables to
#' improve a regression models. Only univariate models are considered.
#' 
#' @param dat A data.frame or tibble as dataset.
#' @param resp Character string for lhs of the formula.
#' @param vars Character vector with variables to transform.
#' @param lambda Real vector of lambda coefficients to evaluate for the 
#' Box-Cox transform.
#' 
#' @return list with two elements
#' \describe{
#' \item{data}{A dataset comprising additional columns with transformed 
#' variables, if they improve Normality.}
#' \item{transformations}{A list of lists with elements for variable, 
#' transformed variable, the Box-Cox results, and the formula terms.}
#' } 
#' @export
box_cox_regression <- function( dat, resp, vars, lambda = seq( -2, +2, 0.1 ) ) {
  stopifnot( is.numeric( dat[[resp]] ) )
  message( sprintf( "box_cox_regression: Generating Box-Cox transformations for %d.", length( c(resp,vars) ) ) )
  trans <- list()
  mod <- dat
  for( var in vars ) {
    if( is.numeric( dat[[var]]) ) {
      mo <- stats::lm( stats::formula( sprintf( "%s ~ %s", resp, var ) ), 
                                       data = dat, y = TRUE, qr = TRUE ) 
      bc <- MASS::boxcox( mo, interp=FALSE, plotit = FALSE )
      py <- bc$x[ which.max( bc$y ) ]
      px <- 1.0/py
      nvar <- sprintf("%s_bc_reg",var)
      mod[[nvar]] <- dat[[var]]^px
      trans[[nvar]] <- list( var = var, new_var = nvar, boxcox_result = bc, 
                             py = py, px = px, 
                             term = sprintf("(%s)^%1.4f",var,px) )
    }
  }
  return( list( transforms = trans, data = mod ) )
  
  # MASS::boxcox()
  # Algorithm in MASS:::boxcox.default:
  # if (is.null(y <- object$y) || is.null(xqr <- object$qr)) 
  #   stop(gettextf("%s does not have both 'qr' and 'y' components", 
  #                 sQuote(deparse(substitute(object)))), domain = NA)
  # if (any(y <= 0)) 
  #   stop("response variable must be positive")
  # n <- length(y)
  # y <- y/exp(mean(log(y)))
  # logy <- log(y)
  # xl <- loglik <- as.vector(lambda)
  # m <- length(xl)
  # for (i in 1L:m) {
  #   if (abs(la <- xl[i]) > eps) 
  #     yt <- (y^la - 1)/la
  #   else yt <- logy * (1 + (la * logy)/2 * (1 + (la * logy)/3 * (1 + (la * logy)/4)))
  #   loglik[i] <- -n/2 * log(sum(qr.resid(xqr, yt)^2))
  # }
  # ...
  # rr <- boxcox(Volume ~ log(Height) + log(Girth), data = trees, lambda = seq(-0.25, 0.25, length = 10))
  # rr
  # plot( rr$x, rr$y )
  
  # xx <- seq( -0.25, 0.25, length.out = 100 )
  # yy <- sapply( xx, FUN=function( ll ) {
  #   # ll <- -0.1
  #   mo <- lm( formula( sprintf( "%s^%1.4f ~ %s", "Volume", ll, "log(Height) + log(Girth)" )), data = trees )
  #   # ret <- logLik(mo)
  #   ret <- sum( log( abs( residuals(mo) ) ) )
  #   # summary(mo)
  #   ret <- as.numeric( ret )
  #   return( ret )
  # } )
  # plot( xx, yy )
  
}

# Box-Cox transformation for binomial regression
#' @title Box-Cox transform for binomial models
#' 
#' @description Determines Box-Cox transformed variables for a binary 
#' outcome modelling.
#' 
#' @param dat A data.frame or tibble comprising the data.
#' @param resp A lhs term for the formula. Must be a binary variable (0/1-coded).
#' @param vars Character vector of variables for which the Box-Cox 
#' transformation should be computed.
#' @param lambda Real numeric vector of lambda values to evaluate.
#' 
#' @return list with two elements
#' \describe{
#' \item{data}{A dataset comprising additional columns with transformed 
#' variables, if they improve Normality.}
#' \item{transformations}{A list of lists with elements for variable, 
#' transformed variable, the Box-Cox results, and the formula terms.}
#' } 
#' @examples 
#' res <- box_cox_binomial( toy_data, "resp", 
#'                          c("the_squared","the_cubed","rnd01"), 
#'                          lambda = seq( -2, +2, 0.1 ) )
#' res$transforms %>% 
#'    map_dfr( function(ee) { ee[["boxcox_result"]] <- NULL; return( ee )}) 
#' res$data[,purrr::map_chr( res$transforms, "new_var")]
#' @export
box_cox_binomial <- function( dat, resp, vars, lambda = seq( -2, +2, 0.1 ) ) {
  stopifnot( is.logical( dat[[resp]] )  )
  message( sprintf( "box_cox_binomial: Generating Box-Cox transformations for %d variables \n", length( c(resp,vars) ) ) )
  # Helper function for probability of observing a
  ll_binomial <- Vectorize( function( x, beta, lambda ) {
    ifelse( 0==lambda,
            1.0 / ( 1.0 + exp( -beta * x ) ),
            1.0 / ( 1.0 + ( 1.0 + lambda * beta * x )^(-1.0/lambda) ) )
  } )
  
  mod <- dat
  trans <- list()
  for( var in vars ) {
    # Note: Formula is lhs var (such that it is evaluated, if containing expressions) and rhs is response (not evaluated)
    mf <- stats::model.frame( stats::formula( sprintf( "%s ~ %s", var, resp )),
                              dat, na.action = stats::na.pass ) %>%
      as_tibble %>% stats::setNames( c("x","y") )
    if( is.numeric( mf[["x"]]) ) {
      # dd <- dat %>% dplyr::select( all_of( c(resp,var ) ) ) %>% setNames( c("y","x") )
      # dd <- dd[which(complete.cases(dd)),]
      dd <- mf[stats::complete.cases(mf),]
      fu_opt <- local( {
        function( px ) {
          beta <- px[1]; lambda <- px[2]
          dd <- dd %>% mutate( pi = ll_binomial( .data$x, beta, lambda ), 
                               pi2 = 2 * pi^.data$y * (1-pi)^(2-.data$y) )
          re <- -sum( dd %>% pull( .data$pi2 ) %>% log )
          re
        }
      } )
      beta0 <- tryCatch({
        or <- stats::optim( par = c(0,0), fn = fu_opt, method="CG" ) # Or: Try BFGS
        or
        beta0 <- round( or$par[2], 4 )
        beta0        
      }, error = function(e) NA )
      if( !is.na(beta0) ) {
        nvar <- sprintf( "%s_bc_bin", var )
        mod[[nvar]] <- mf[["x"]]^beta0
        trans[[nvar]] <- list( var = var, new_var = nvar, boxcox_result = or, 
                               py = 1.0/beta0, px = beta0, 
                               term = sprintf("(%s)^%1.4f",var,beta0) )
      }
    }
  }
  return( list( transforms = trans, data = mod ) )
    
  # From 01_numeric_distributions.R under nhl_safety:
  # pp <- Vectorize( function( x, beta, lambda ) {
  #   ifelse( 0==lambda, 
  #           1.0 / ( 1.0 + exp( -beta * x ) ), 
  #           1.0 / ( 1.0 + ( 1.0 + lambda * beta * x )^(-1.0/lambda) ) )
  # } )
  #
  # dd <- data.frame( y = df.folds[,out_var], x = df.folds[,num_vars[i]], stringsAsFactors = FALSE ) %>% filter( !is.na(x) )
  # fu_opt <- local( {
  #   function( px ) {
  #     beta <- px[1]; lambda <- px[2]
  #     dd <- dd %>% mutate( pi = pp( x, beta, lambda ), pi2 = 2 * pi^y * (1-pi)^(2-y) )
  #     re <- -sum( dd %>% pull( pi2 ) %>% log ) 
  #     re
  #   }
  # } )
  # or <- optim( par = c(0,0), fn = fu_opt, method="CG" )
  # or
  # beta0 <- round( or$par[2], 1 )
  # xvalt <- df.folds %>% pull( num_vars[i] )
  # xvalt <- xvalt^beta0
}

# Box-Cox transformation for survival regression
# TODO

# Mixture models ----

eval_aics <- function( dat, var, n_comp = 5, m_fits = 25, 
                       min_fits_converged = 24 ) {
  ks <- seq_len( n_comp )
  idx <- which(!is.na(dat[,var]))

  models <- list()
  tab <- NULL
  midx <- 1
  for( k in seq_len( n_comp ) ) {
    for( j in seq_len( m_fits ) ) {
      mo <- NULL
      mo <- tryCatch({
        flexmix::flexmix( stats::formula( sprintf( "%s ~ 1", var ) ), 
                          data = dat[idx,], 
                          k = k )
      }, error = function(e) NULL )
      models[[midx]] <- mo
      if( !is.null(mo) ) {
        tab <- bind_rows( tab, tibble( round = j, k = k, 
                                       midx = midx, aic = stats::AIC(mo), 
                                       converged = mo@converged ))
      }
      midx <- midx + 1
    }
  }
  
  agg <- tab %>% 
    filter( .data$converged ) %>%
    group_by(k) %>% 
    summarise( min_aic = min(.data$aic), sum_converged = sum( .data$converged ) )

  best <- agg %>%
    filter( min_fits_converged <= .data$sum_converged ) %>%
    filter( min(.data$min_aic) == .data$min_aic ) %>%
    filter( min(.data$k) == .data$k )
  kmin <- best %>% pull( .data$k )
  maic <- best %>% pull( .data$min_aic )
  
  best_idx <- tab %>% filter( .data$k == kmin & maic == .data$aic & .data$converged ) %>% 
    pull( .data$midx )
  momin <- models[[ best_idx[1] ]]

  ret <- list( var = var, 
               n_comp = n_comp, m_fits = m_fits, 
               min_fits_converged = min_fits_converged,
               aic_tab = tab, aic_aggregate = agg, min_k = kmin, 
               best_model = momin, cut_points = NULL )
  
  if( 1 < kmin ) {
    message( sprintf( "Variable %s is multi-modal with %d Normal components. Determining cut-points.", var, kmin ) )
    prm <- flexmix::parameters(momin)
    ood <- order(prm["coef.(Intercept)",])
    prio <- momin@prior[ood]
    prm <- prm[,ood]
    
    cps <- c()
    i <- 1
    j <- i + 1
    while( j <= kmin ) {
      rt <- NULL
      rt <- tryCatch({
        fu_root <- function( val ) {
          d1 <- stats::dnorm( x = val, 
                              mean = prm["coef.(Intercept)",i], 
                              sd = prm["sigma",i]  ) * prio[i]
          d2 <- stats::dnorm( x = val, 
                              mean = prm["coef.(Intercept)",j], 
                              sd = prm["sigma",j]  ) * prio[j]
          ret <- d1 - d2
          return( ret )
        }
        
        # We are only interested if the cut-point is between the two adjacent modes.
        rt <- stats::uniroot( f = fu_root, 
                              interval = c(prm["coef.(Intercept)",i],prm["coef.(Intercept)",j]),
                              extendInt = "no" ) 
        rt
      }, error = function(e) NULL ) # Catch error, if fu_root is not of opposite signs at both interval ends.
      if( !is.null(rt) ) {
        cp <- rt$root
        cps <- c( cps, cp )
      }
      
      i <- i + 1
      j <- i + 1
    } # while
    
    cps <- c( -Inf, cps, Inf )
    ret[['cut_points']] <- cps
  } # if
  
  return( ret )
} # eval_aics (END)

#' @title A function to evaluate continues variables for multi-modiality using
#' Gaussian-Mixture Models
#' 
#' @description The function receives a list of variables and first performs a 
#' Gaussian-Mixture Model (GMM) model selection using the Akaike Information 
#' Criterion (AIC) to determine the number of components. If this process 
#' indicates more than one componet, it determines cut-points between the 
#' individual mixture distributions and generates a transformed categorical 
#' variable based on those.
#' 
#' @param dat A data.frame or tibble
#' @param resp A response term for the lhs of a formula
#' @param vars A character vector of variables to screen
#' @param n_comp Integer determining the maximum number of components to check
#' @param m_fits Each GMM is fitted multiple times to avoid local minima. m_fits 
#' determines how many GMMs are generated each time.
#' @param min_fits_converged Minimal number of GMMs that need to converge before 
#' considering their AIC for model selection.
#' 
#' @return list with two elements
#' \describe{
#' \item{data}{A dataset comprising additional columns with transformed 
#' variables, if they improve Normality.}
#' \item{transformations}{A list of check_aic results for each variable.}
#' } 
#' @examples 
#' library( ggplot2 )
#' library( flexmix )
#' resp <- "resp"
#' vars <- grep( "the_|rnd", colnames(toy_data), value=TRUE )
#' 
#' toy_data %>%
#'   ggplot( aes( x=the_multi, y=..density.. ) ) +
#'   geom_histogram( bins = 100, alpha=0.5 ) +
#'   geom_density( bw = "ucv" ) 
#' 
#' mumo <- check_multimodality( dat = toy_data, resp = resp, 
#'                              vars = c("the_multi","rnd01","rnd02"),
#'                              n_comp = 3, 
#'                              m_fits = 25, min_fits_converged = 20 )
#' 
#' mumo$transforms$the_multi$aic_aggregate
#' mumo$transforms$the_multi$best_model
#' parameters(mumo$transforms$the_multi$best_model)
#' prior(mumo$transforms$the_multi$best_model)
#' @export
check_multimodality <- function( dat, resp, vars, n_comp = 5,
                                 m_fits = 25, min_fits_converged = 24 ) {
  message( sprintf( "check_mulitmodality: Checking mulit-modality for %d variables (k<=%d, fits=%d, convgd fits=%d).", length(c(resp,vars)), n_comp, m_fits, min_fits_converged ) )
  ret <- list()
  # Check response, if numeric, if it is bimodal
  mod <- dat
  if( ! is.null(resp) && is.numeric( dat[[resp]] ) ) {
    message( sprintf( "Processing response %s.", resp ))
    evl <- eval_aics( dat, resp, n_comp = n_comp, 
                      m_fits = m_fits, 
                      min_fits_converged = min_fits_converged )
    if( !is.null(evl$cut_points) ) {
      nvar <- sprintf( "%s_resp_grp", resp )
      mod[[nvar]] <- as.character( cut( mod[[resp]], 
                                        breaks = evl$cut_points, 
                                        labels = sprintf( "resp_group[%d]", seq_len(length(evl$cut_points)-1) ), 
                                        include.lowest = TRUE ) )
      evl$transformed_var <- nvar
    }
    ret[[resp]] <- evl
  }
  
  # Check all numeric variables if they are bimodal
  for( var in vars ) {
    if( is.numeric(dat[[var]]) ) {
      message( sprintf( "Processing %s.", var ))
      evl <- eval_aics( dat, var, n_comp = n_comp, 
                        m_fits = m_fits, 
                        min_fits_converged = min_fits_converged )
      if( !is.null(evl$cut_points) ) {
        nvar <- sprintf( "%s_grp", var )
        mod[[nvar]] <- as.character( cut( mod[[var]], 
                                          breaks = evl$cut_points, 
                                          labels = sprintf( "group[%d]", seq_len(length(evl$cut_points)-1) ), 
                                          include.lowest = TRUE ) )
        evl$transformed_var <- nvar
      }
      ret[[var]] <- evl
    }
  }
  
  return( list( transforms = ret, data = mod ) )
} # check_modal (END)
