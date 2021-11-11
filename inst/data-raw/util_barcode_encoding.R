
#
# Two functions for decoding the TCGA Barcode into it's components
#
lst_meta_cols <- c("Project","TSS","Participant","Sample_Vial","Portion_Analyte","Plate","Center","Sample","Vial","Portion","Analyte")

process_sample_ids <- local( { function( dat,key_name ) { 
  dat %>%
    tidyr::separate( key_name, c("Project","TSS","Participant","Sample_Vial","Portion_Analyte","Plate","Center"), sep="-", remove=FALSE ) %>%
    mutate( Sample = gsub( "^(..).", "\\1", Sample_Vial ), Vial = gsub( "^..(.)", "\\1", Sample_Vial ) ) %>%
    mutate( Portion = gsub( "^(..).", "\\1", Portion_Analyte ), Analyte = gsub( "^..(.)", "\\1", Portion_Analyte ) ) %>%
    mutate_at( lst_meta_cols, as.factor )
} }, envir = list( lst_meta_cols = lst_meta_cols ) )

