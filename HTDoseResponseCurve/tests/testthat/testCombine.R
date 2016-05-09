library(HTDoseResponseCurve)
context("Combine data and plate map")

test_that("combine works", {
    pkg = "HTDoseResponseCurve"
    fn_data = system.file("extdata", "sample_data_96.xlsx", package = pkg)
    raw_plate = read_plates_from_Incucyte_export( fn_data, "plate_1", 
                                                 number_of_wells=96)
    path_to_file= system.file("extdata", "sample_data_96_platemap.txt",
                              package=pkg)
    plate_map=read_platemap_from_Incucyte_XML( path_to_file )
    
    # Vehicle is 100
    expect_error( combine_data_and_map( raw_plate, plate_map ), 
      paste("negative_control parameter passed with value 0 so we assume each",
            "treatment has one or more wells with concentration = 0 which is",
            "not the case for Vehicle"))
                  
    # Vehicle is 100
    expect_error( combine_data_and_map( raw_plate, plate_map, -1 ), 
      paste("negative_control parameter passed with value -1 so we assume each",
            "treatment has one or more wells with concentration = -1 which is",
            "not the case for Vehicle"))
    
    # bogus name
    expect_error( combine_data_and_map( raw_plate, plate_map, "vvv" ),
        "negative_control vvv is not among the treatments on this plate")
    
    # correct call
    ds_n = combine_data_and_map( raw_plate, plate_map, "Vehicle" ) 
    expect_equal( sum(!is.na(raw_plate[,1:12])), dim(ds_n)[1])
    expect_equal( round(raw_plate[1,2], 5), round(ds_n$value[1], 5) )
    expect_equal( TRUE, ds_n$is_negative_control[1])
    expect_equal( FALSE, ds_n$is_negative_control[2])
    expect_equal( "Vehicle", ds_n$negative_control[1])
    expect_equal( "Vehicle", ds_n$negative_control[2])
    
    # data frame: missing a drug
    vd = data.frame(
        drug=c("drug_1", "drug_2"), vehicle=c("Vehicle", "Vehicle"),
        stringsAsFactors=FALSE)
    expect_error(combine_data_and_map( raw_plate, plate_map, vd ),
     paste("Not all drugs have been assigned a negative control,",
           "check negative_control parameter"))
    
    # data frame: bogus vehicle
    vd = data.frame(
        drug=c("drug_1", "drug_2"), vehicle=c("vvv", "vvv"),
        stringsAsFactors=FALSE)
    expect_error(combine_data_and_map( raw_plate, plate_map, vd ),
     paste("Not all drugs have been assigned a negative control,",
           "check negative_control parameter"))
    
    # data frame: bogus column named boing instead of drug
    vd = data.frame(
        boing=c("drug_1", "drug_2"), vehicle=c("vvv", "vvv"),
        stringsAsFactors=FALSE)
    expect_error(combine_data_and_map( raw_plate, plate_map, vd ),
     paste("if passed as a data frame, parameter negative_control must have",
           "columns named 'drug' and 'vehicle'"))
    
    # data frame: well-formatted
    vd = data.frame(
        drug=c("drug_1", "drug_2", "drug_3", "Vehicle"), 
        vehicle=rep("Vehicle", 4 ),
        stringsAsFactors=FALSE)
    ds_n = combine_data_and_map( raw_plate, plate_map, vd )
    expect_equal( sum(!is.na(raw_plate[,1:12])), dim(ds_n)[1])
    expect_equal( round(raw_plate[1,2], 5), round(ds_n$value[1], 5) )
    expect_equal( TRUE, ds_n$is_negative_control[1])
    expect_equal( FALSE, ds_n$is_negative_control[2])
    expect_equal( "Vehicle", ds_n$negative_control[1])
    expect_equal( "Vehicle", ds_n$negative_control[2])
    
} )
