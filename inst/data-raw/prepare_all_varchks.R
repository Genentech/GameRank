
rm( list=ls() )
devtools::load_all("~/GameRank/")
load("~/GameRank/data/tcga_ucec_cna_cnv.Rdata")
re <- "Surv( days_to_death, vital_status )"
va <- lst_vars
vck <- check_variables( dat, re, va )
saveRDS( vck, file = "~/GameRank/data/ucec_surv_varchk.rds")

rm( list=ls() )
devtools::load_all("~/GameRank/")
load("~/GameRank/data/tcga_dlbcl_cna_cnv.Rdata")
re <- "Surv( days_to_death, vital_status )"
va <- lst_vars
vck <- check_variables( dat, re, va )
saveRDS( vck, file = "~/GameRank/data/dlbcl_surv_varchk.rds")

rm( list=ls() )
devtools::load_all("~/GameRank/")
load("~/GameRank/data/tcga_coad_cna_cnv.Rdata")
re <- "Surv( days_to_death, vital_status )"
va <- lst_vars
vck <- check_variables( dat, re, va )
saveRDS( vck, file = "~/GameRank/data/coad_surv_varchk.rds")

rm( list=ls() )
devtools::load_all("~/GameRank/")
load("~/GameRank/data/tcga_luad_cna_cnv.Rdata")
re <- "Surv( days_to_death, vital_status )"
va <- lst_vars
vck <- check_variables( dat, re, va )
saveRDS( vck, file = "~/GameRank/data/luad_surv_varchk.rds")

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

