rm( list=ls() )
library( dplyr )
library( tidyr )
library( purrr )

library( ggplot2 )
library( rjags )
library( random )

load.module("mix")

file_tdata <- "inst/data-raw/toy_data.rda"
file_toy <- "inst/data-raw/toy_model.jags"
the_data <- list(
  mu =          c(  1.0,  1.2,  0.7,  0.3,  2.1,  1.3 ),
  cov = matrix( c(  1.5,  0.0,  0.0,  0.0,  0.0,  0.5,
                    0.0,  1.0,  0.0,  0.0,  0.0,  0.0,
                    0.0,  0.0,  1.0,  0.0,  0.0,  0.0,
                    0.0,  0.0,  0.0,  1.0,  0.0,  0.0,
                    0.0,  0.0,  0.0,  0.0,  1.0,  0.0,
                    0.5,  0.0,  0.0,  0.0,  0.0,  2.0
                   ), 
               nrow = 6, ncol = 6, byrow = TRUE ),
  co =          c( -0.5,  0.7, -0.2,  0.4,  1.3,  0.1 ),
  tau = 0.7,
  mimu = c(0.4,2.7),
  mita = 1.0/c(0.3,0.5),
  mipi = c(0.4,0.6)
)

the_model <- '
model {
  
  omega <- inverse( cov )
  ee ~ dmnorm( mu, omega )
  mm ~ dnormmix( mimu, mita, mipi )
  
  vv[1] <- ee[1]             /* Normal distributed variable */
  vv[2] <- ee[2]^2           /* ee[2]^2 squared Normal distributed variable */
  vv[3] <- ee[3]^3           /* ee[3]^3 power of 3 Normal distributed variable */
  vv[4] <- exp( ee[4] )       /* exp( ee[4] ) exponential Normal distributed variable */
  vv[5] <- mm             /* mm Mixture Normal distributed variable */
  vv[6] <- pow( ee[6], 5)             /* Normal distributed variable to the power of +5 for Box-Cox transformations; jags has issues with non-rationa powers like 0.3 */
  
  oo <- co %*% vv  
  
  ll ~ dlogis( oo, tau )
  pp <- plogis( ll, oo, tau  )
}

' # JAGS model (END)
cat( the_model, file = file_toy )

jmod <- jags.model( file = file_toy, data = the_data, n.chains = 2, n.adapt = 5000 )
jsmp <- coda.samples( model = jmod, variable.names = c("ll", "oo", "vv", "ee","mm","pp"), n.iter = 3600, thin = 10 )
df <- jsmp[[1]] %>% as_tibble %>% mutate( resp = (pp > 0.5) )
colnames(df) <- gsub( "\\[|\\]", "", colnames(df) )
df %>% head

ggplot( data = df ) +
  geom_histogram( aes( x=mm, y=..density.. ), bins = 100 ) + 
  geom_density( aes( x=mm, y=..density.. ), bw = "ucv", color=1 )
ggplot( data = df ) +
  geom_density( aes( x=ee1, y=..density.. ), bw = "ucv", color=1 ) +
  geom_density( aes( x=ee2, y=..density.. ), bw = "ucv", color=2 ) +
  geom_density( aes( x=ee3, y=..density.. ), bw = "ucv", color=3 ) +
  geom_density( aes( x=ee4, y=..density.. ), bw = "ucv", color=4 ) +
  geom_density( aes( x=ee5, y=..density.. ), bw = "ucv", color=5 ) +
  geom_density( aes( x=ee6, y=..density.. ), bw = "ucv", color=6 ) +
  geom_density( aes( x=mm, y=..density.. ), bw = "ucv", color=7, lwd=1.5 )
ggplot( data = df ) +
  geom_density( aes( x=vv1, y=..density.. ), bw = "ucv", color=1 ) +
  geom_density( aes( x=vv2, y=..density.. ), bw = "ucv", color=2 ) +
  geom_density( aes( x=vv3, y=..density.. ), bw = "ucv", color=3 ) +
  geom_density( aes( x=vv4, y=..density.. ), bw = "ucv", color=4 ) +
  geom_density( aes( x=vv5, y=..density.. ), bw = "ucv", color=5, lwd=1.5 ) +
  geom_density( aes( x=vv6, y=..density.. ), bw = "ucv", color=6 )
ggplot( data = df ) +
  geom_density( aes( x=oo, y=..density.. ), bw = "ucv", color=1 ) 
ggplot( data = df ) +
  geom_density( aes( x=ll, y=..density.. ), bw = "ucv", color=1 ) 
ggplot( data = df ) +
  geom_density( aes( x=pp, y=..density.. ), bw = "ucv", color=1 ) 
ggplot( data = df ) +
  geom_point( aes( x=ee1, y=ee2 ), color=1 ) +
  geom_smooth(  aes( x=ee1, y=ee2 ), color=2, method= "lm" )

mod <- glm(   oo ~ ee1 + ee2 + ee3 + ee4 + ee5 + ee6, "gaussian", df )
mod

mod <- glm(   pp ~ ee1 + ee2 + ee3 + ee4 + ee5 + ee6, "binomial", df )
mod
summary(mod)

mod <- glm( resp ~ ee1 + ee2 + ee3 + ee4 + ee5 + ee6, "binomial", df )
mod
the_data$co

dat <- df %>% 
  transmute( USUBJID = sprintf( "PAT%0.3d", dplyr::row_number() ),  
             resp, 
             pp, oo, 
             the_normal = vv1, 
             the_squared = vv2, 
             the_cubed = vv3, 
             the_exped = vv4, 
             the_multi = vv5, 
             the_power = vv6  )
dat %>% head

# After sampling from structure generate completely random variables
n_rndvars <- 20
m_rows <- nrow(dat)
m_rows * n_rndvars
rr <- random::randomNumbers( n = m_rows * n_rndvars, min = 0, max = 1000 )
rr <- rr / 1000
mat_rnd <- matrix( rr, ncol = n_rndvars )
colnames(mat_rnd) <- sprintf( "rnd%0.2d", 1:n_rndvars )
# Maybe scale random variables?
dat <- dat %>% bind_cols( as_tibble( mat_rnd )  )
dat %>% head

random::randomQuota()               
# random::randomNumbers()

toy_data <- dat
save( toy_data, file = file_tdata  )
save( toy_data, file = "~/GameRank/data/toy_data.rda"  )

# https://docs.google.com/presentation/d/1bc5ktbty1BOLV6J_Z4P-_WlW6aaQ1xE7hHGKq0pYPPA/edit#slide=id.g6ffbf95dbd_0_111
# https://www.nature.com/articles/s41586-021-03430-5
# http://contributions.bioconductor.org/general.html
# https://bioconductor.org/packages/3.15/bioc/html/BiocCheck.html
# https://r-pkgs.org/vignettes.html
# https://stackoverflow.com/questions/22265837/transfer-git-repositories-from-gitlab-to-github-can-we-how-to-and-pitfalls-i
  

