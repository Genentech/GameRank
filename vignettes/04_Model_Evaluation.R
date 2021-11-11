#' ---
#' title: "GameRank Algorithm (Vignette)"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{GameRank}
#'   %\VignetteEngine{knitr::rmarkdown}
#'   %\VignetteEncoding{UTF-8}
#' ---
# /* Hit Ctrl + Shift + K in Rstudio to generate html_document for quick view. Note: Vignette YAML header doesn't work then. */
# /* Run knitr::spin("vignettes/01_GameRank.R", format="Rmd") to turn this into a Vignette */

#+ setup, include=FALSE
rm( list=ls() )
library( dplyr )
library( tidyr )
library( purrr )

library( survival )
library( survminer )

devtools::load_all("~/GameRank/")

# setting outputs 
dir_varsel <- "out/vig4/"
file_varchk <- system.file( "data", "luad_nopckyrssmked_varchk.rds", package = "GameRank" )
file_gmr <- file.path( dir_varsel, "luad__gmrsel.rds")

# load data
load( "data/tcga_luad_cna_cnv.Rdata" )
dat <- dat %>% filter( grepl( "01A-11D", dat$SampleID_cna ) ) # Filter only one half of samples 
vck <- readRDS( file = file_varchk )

# Prepare variable selection
re <- vck %>% filter( is_response ) %>% pull( variable )
va <- vck %>% filter( "Perfect"==check_missing ) %>% arrange( desc(mutual_information), desc(entropy), rot.p ) %>% head( 200 ) %>% pull( variable ) %>% unique

