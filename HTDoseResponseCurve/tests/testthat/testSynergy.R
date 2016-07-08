library(HTDoseResponseCurve)
context("Normalize synergy data")

test_that("synergy works", {
   
    pkg = "HTDoseResponseCurve"
    fn_data = system.file("extdata", "sample_synergy_ray.txt", package = pkg)
    syn=read.table( fn_data, header=TRUE, sep='\t', stringsAsFactors=FALSE)
    head(syn)
    ds=create_synergy_dataset( syn$sample_type, syn$treatment, syn$treatment_2, 
         syn$concentration, syn$concentration_2, syn$value, 
         negative_control="DMSO")
    expect_equal( dim(syn)[1], length(ds$sample_type) )
    idx_l1_d = which( syn$sample_type=="line_1" & syn$treatment=="DMSO")
    mu_l1_d = mean(syn$value[idx_l1_d])
    
    idx_l1_6=which( syn$sample_type=="line_1" & syn$treatment=="drug_1" & 
               syn$concentration==6.25 & syn$concentration_2==0)
    expect_equal( round( syn$value[idx_l1_6] / mu_l1_d, 2 ), 
                  round( ds$value_normalized[idx_l1_6], 2))
    
} )