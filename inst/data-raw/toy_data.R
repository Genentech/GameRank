rm( list=ls() )
library( dplyr )
library( tidyr )
library( purrr )

library( ggplot2 )
library( rjags )
library( random )

library( survival )
library( survminer )

load.module("mix")

file_tdata <- "inst/data-raw/toy_data.rda"
file_data <- "data/toy_data.rda"
file_toy <- "inst/data-raw/toy_model.jags"
the_data <- list(
  mu =          c(  1.00,  0.50,  0.50,  0.10,  0.00,  0.10 ),
  cov = matrix( c(  1.00,  0.00,  0.00,  0.00,  0.00,  0.00,
                    0.00,  0.05,  0.00,  0.00,  0.00,  0.00,
                    0.00,  0.00,  0.05,  0.00,  0.00,  0.00,
                    0.00,  0.00,  0.00,  0.05,  0.00,  0.00,
                    0.00,  0.00,  0.00,  0.00,  1.00,  0.00,
                    0.00,  0.00,  0.00,  0.00,  0.00,  0.02
                   ), 
               nrow = 6, ncol = 6, byrow = TRUE ),
  cof =          c(  2.0, -2.0,  2.0, -2.0,  2.0,  2.0 ),
  mimu = c(0.2,0.7),
  mita = 1.0/c(0.01,0.02),
  mipi = c(0.6,0.4)
  
)
the_data[["b0"]] <- sum( c( the_data[["mu"]], the_data[["mimu"]] ) )

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
  
  reg <- cof %*% vv - b0
  
  logit( pi ) <- reg 
  resp ~ dbern( pi  )
  srv  ~ dpois( pi * 10 ) 
  dys  <- srv * 30.25
}

' # JAGS model (END)
cat( the_model, file = file_toy )

jmod <- jags.model( file = file_toy, data = the_data, n.chains = 2, n.adapt = 5000 )
jsmp <- coda.samples( model = jmod, variable.names = c("vv", "ee","mm", "reg", "pi","resp", "srv", "dys"), n.iter = 3600, thin = 10 )
df <- jsmp[[1]] %>% 
  as_tibble %>% 
  mutate( resp = as.logical(resp) )
colnames(df) <- gsub( "\\[|\\]", "", colnames(df) )
df %>% head
df %>% summary

df %>% ggplot( aes( x=ee1, y=vv1 )) + geom_point()
df %>% ggplot( aes( x=ee2, y=vv2 )) + geom_point() + geom_point( aes( x=ee2, y=sqrt(vv2) ), color = "red" )
df %>% ggplot( aes( x=ee3, y=vv3 )) + geom_point() + geom_point( aes( x=ee3, y=(vv3)^(1/3) ), color = "red" )
df %>% ggplot( aes( x=ee4, y=vv4 )) + geom_point() + geom_point( aes( x=ee4, y=log(vv4) ), color = "red" )
df %>% ggplot( aes( x=vv5, y=vv5 )) + geom_point() 
df %>% ggplot( aes( x=ee6, y=vv6 )) + geom_point() + geom_point( aes( x=ee6, y=(vv6)^(1/5) ), color = "red" )

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
ggplot( data = df ) + geom_density( aes( x=vv1, y=..density.. ), bw = "ucv", color=1 ) 
ggplot( data = df ) + geom_density( aes( x=vv2, y=..density.. ), bw = "ucv", color=2 )
ggplot( data = df ) + geom_density( aes( x=vv3, y=..density.. ), bw = "ucv", color=3 )
ggplot( data = df ) + geom_density( aes( x=vv4, y=..density.. ), bw = "ucv", color=4 )
ggplot( data = df ) + geom_density( aes( x=vv5, y=..density.. ), bw = "ucv", color=5, lwd=1.5 )
ggplot( data = df ) + geom_density( aes( x=vv6, y=..density.. ), bw = "ucv", color=6 )


ggplot( data = df ) +
  geom_density( aes( x=reg, y=..density.. ), bw = "ucv", color=1 ) 

glm( reg ~ vv1 + vv2 + vv3 + vv4 + vv5 + vv6, "gaussian", df )
glm( reg ~ ee1 + I( ee2^2 ) + I( ee3^3 ) + I( exp(ee4) ) + vv5 + I( ee6^5 ), "gaussian", df )

ggplot( data = df, aes( x=reg, y=pi ) ) + geom_point()
ggplot( data = df ) + geom_density( aes( x=pi, y = ..density.. ) )

table( df$resp )
ggplot( data = df, aes( x=reg, y=resp ) ) + geom_point()

ggplot( data = df, aes( x=reg, y=srv ) ) + geom_point()
ggplot( data = df, aes( x=reg, y=dys, color=resp ) ) + geom_point()

ggsurvplot( survfit( Surv( dys, resp ) ~ 1, df ), df, risk.table = TRUE )


dat <- df %>% 
  transmute( USUBJID = sprintf( "PAT%0.3d", dplyr::row_number() ),  
             dys,
             srv,
             resp, 
             pi,
             reg,
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

# create decorrelated variables to reg
mat_rnd <- cbind( df$reg, mat_rnd )
mat_rnd

pc <- princomp( mat_rnd )
cor( cbind( df$reg, pc$scores ) ) 
# 1st PC is regression model, all other PCs are uncorrelated to it
rnds <- pc$scores[,-1]
colnames(rnds) <- sprintf( "rnd%0.2d", seq_len(n_rndvars)  )
rnds %>% head

# Maybe scale random variables?
dat <- dat %>% bind_cols( as_tibble( rnds )  )
dat %>% head

random::randomQuota()               
# random::randomNumbers()

devtools::load_all("~/GameRank/")
vck <- check_variables( dat, "reg", grep( "the_|rnd", colnames(dat), value=TRUE ) )
vck %>% filter( !is_response ) %>% arrange( desc(mutual_information), desc(entropy) ) %>% pull( variable )

fwd <- forward( dat, "resp", vck %>% filter( !is_response ) %>% arrange( desc(mutual_information), desc(entropy) ) %>% pull( variable ), 
                fn_train_binomial, fn_eval_binomial_auroc, 6L, 3L, TRUE )
fwd <- bidirectional( dat, "resp", vck %>% filter( !is_response ) %>% arrange( desc(mutual_information), desc(entropy) ) %>% pull( variable ), 
                fn_train_binomial, fn_eval_binomial_auroc, 6L, 3L, TRUE )
fwd$variable_selections
fwd$agg_results %>% arrange(desc(mean_validation) )

toy_data %>%
  dplyr::select( all_of( setdiff( c("reg", vck$variable ), c("the_multi_grp" ) )) ) %>%
  tidyr::pivot_longer( cols = setdiff( c(vck$variable ), c("reg","the_multi_grp" ) ), names_to = "var", values_to = "value" ) %>%
  mutate( flg = grepl( "the_", var ) ) %>%
  ggplot( aes( x=reg, y=value, group=var, color=var ) ) + 
  facet_grid( . ~ flg ) + 
  geom_point() +
  geom_smooth( method = "lm" )

toy_data <- dat
save( toy_data, file = file_tdata  )
# save( toy_data, file = file_data  )

# https://docs.google.com/presentation/d/1bc5ktbty1BOLV6J_Z4P-_WlW6aaQ1xE7hHGKq0pYPPA/edit#slide=id.g6ffbf95dbd_0_111
# https://www.nature.com/articles/s41586-021-03430-5
# http://contributions.bioconductor.org/general.html
# https://bioconductor.org/packages/3.15/bioc/html/BiocCheck.html
# https://r-pkgs.org/vignettes.html
# https://stackoverflow.com/questions/22265837/transfer-git-repositories-from-gitlab-to-github-can-we-how-to-and-pitfalls-i
  

