
#
# TCGA Barcode documentation:
# https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
# https://docs.gdc.cancer.gov/Encyclopedia/pages/images/TCGA-TCGAbarcode-080518-1750-4378.pdf
# => e.g. TCGA-2E-A9G8-01A-11D-A402-01 is 01 sample type ( 1-9 are tumor, 10 - 19 are control)
# 


rm( list=ls() )
devtools::load_all("~/GameRank/")
load("~/GameRank/data/tcga_ucec_cna_cnv.Rdata")
dat <- dat %>% filter( Sample %in% c("01","02") ) # Select Tumor samples only
re <- "Surv( days_to_death, vital_status )"
va <- lst_vars
vck <- check_variables( dat, re, va )
comment(vck) <- "Variable check only for tumor samples"
saveRDS( vck, file = "~/GameRank/data/ucec_surv_varchk.rds")

rm( list=ls() )
devtools::load_all("~/GameRank/")
load("~/GameRank/data/tcga_dlbcl_cna_cnv.Rdata")
dat <- dat %>% filter( Sample %in% c("01") ) # Select Tumor samples only
re <- "Surv( days_to_death, vital_status )"
va <- lst_vars
vck <- check_variables( dat, re, va )
comment(vck) <- "Variable check only for tumor samples"
saveRDS( vck, file = "~/GameRank/data/dlbcl_surv_varchk.rds")

rm( list=ls() )
devtools::load_all("~/GameRank/")
load("~/GameRank/data/tcga_coad_cna_cnv.Rdata")
dat <- dat %>% filter( Sample %in% c("01","02","06") ) # Select Tumor samples only
re <- "Surv( days_to_death, vital_status )"
va <- lst_vars
vck <- check_variables( dat, re, va )
comment(vck) <- "Variable check only for tumor samples"
saveRDS( vck, file = "~/GameRank/data/coad_surv_varchk.rds")

rm( list=ls() )
devtools::load_all("~/GameRank/")
load("~/GameRank/data/tcga_luad_cna_cnv.Rdata")
dat <- dat %>% filter( Sample %in% c("01","02") ) # Select Tumor samples only
re <- "Surv( days_to_death, vital_status )"
va <- lst_vars
vck <- check_variables( dat, re, va )
comment(vck) <- "Variable check only for tumor samples"
saveRDS( vck, file = "~/GameRank/data/luad_surv_varchk.rds")

re <- "karnofsky_performance_score"
va <- setdiff( lst_vars, re )
vck <- check_variables( dat, re, va )
comment(vck) <- "Variable check only for tumor samples"
saveRDS( vck, file = "~/GameRank/data/luad_karnofsky_varchk.rds")

re <- "number_pack_years_smoked"
va <- setdiff( lst_vars, re )
vck <- check_variables( dat, re, va )
comment(vck) <- "Variable check only for tumor samples"
saveRDS( vck, file = "~/GameRank/data/luad_nopckyrssmked_varchk.rds")


rm( list=ls() )
library(rmarkdown)
devtools::load_all("~/GameRank/")
vck <- readRDS( file = "~/GameRank/data/coad_surv_varchk.rds")
rmarkdown::render( input = "~/GameRank/inst/templates/rmd_variable_checks.Rmd", output_dir = "~/GameRank/out/coad_varchk/",
                   params = list( vck = vck ), envir = new.env()  )

vck <- readRDS( file = "~/GameRank/data/dlbcl_surv_varchk.rds")
rmarkdown::render( input = "~/GameRank/inst/templates/rmd_variable_checks.Rmd", output_dir = "~/GameRank/out/dlbcl_varchk/",
                   params = list( vck = vck ), envir = new.env()  )

vck <- readRDS( file = "~/GameRank/data/ucec_surv_varchk.rds")
rmarkdown::render( input = "~/GameRank/inst/templates/rmd_variable_checks.Rmd", output_dir = "~/GameRank/out/ucec_varchk/",
                   params = list( vck = vck ), envir = new.env()  )

vck <- readRDS( file = "~/GameRank/data/luad_surv_varchk.rds")
rmarkdown::render( input = "~/GameRank/inst/templates/rmd_variable_checks.Rmd", output_dir = "~/GameRank/out/luad_varchk/",
                   params = list( vck = vck ), envir = new.env()  )

vck <- readRDS( file = "~/GameRank/data/luad_karnofsky_varchk.rds")
rmarkdown::render( input = "~/GameRank/inst/templates/rmd_variable_checks.Rmd", output_dir = "~/GameRank/out/luad_varchk/",
                   output_file = "rmd_variable_checks_karnofsky",
                   params = list( vck = vck ), envir = new.env()  )

vck <- readRDS( file = "~/GameRank/data/luad_nopckyrssmked_varchk.rds")
rmarkdown::render( input = "~/GameRank/inst/templates/rmd_variable_checks.Rmd", output_dir = "~/GameRank/out/luad_varchk/",
                   output_file = "rmd_variable_checks_nopckyrssmked",
                   params = list( vck = vck ), envir = new.env()  )

