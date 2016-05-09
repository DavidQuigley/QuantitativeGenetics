library(HTDoseResponseCurve)
context("Load data from files")

test_that("loading works", {
    pkg = "HTDoseResponseCurve"
    fn_data = system.file("extdata", "sample_data_96.xlsx", package = pkg)
    
    ds_raw = read_plates_from_Incucyte_export( fn_data, "plate_1", 
                                             number_of_wells=96)
    expect_equal( length(sort(unique((ds_raw$hours)))), 14)
    expect_equal( dim(ds_raw)[1], 112)
    expect_equal( dim(ds_raw)[2], 14)
    
    path_to_file= system.file("extdata", "sample_data_96_platemap.txt",
                              package=pkg )
    plate_map=read_platemap_from_Incucyte_XML( path_to_file )
    expect_equal( length(sort( unique( matrix( t(plate_map$treatment) ) ) )), 5)
    expect_equal( sum(plate_map$treatment=="Vehicle"), 11 )
    expect_equal( plate_map$treatment[1,3], "drug_1")
    
    fn_map= system.file("extdata", "sample_data_96_platemap.xlsx",package=pkg)
    plate_map=read_platemap_from_excel( fn_map, number_of_wells=96)
        
    expect_equal( length(sort( unique( matrix( t(plate_map$treatment) ) ) )), 5)
    expect_equal( sum(plate_map$treatment=="Vehicle"), 11 )
    expect_equal( plate_map$treatment[1,3], "drug_1")
    
} )
    