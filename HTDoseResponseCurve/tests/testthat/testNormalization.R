library(HTDoseResponseCurve)
context("Normalize loaded data")

test_that("normalization works", {
    
    pkg = "HTDoseResponseCurve"
    fn_data = system.file("extdata", "sample_data_96.xlsx", package = pkg)
    raw_plate = read_plates_from_Incucyte_export( fn_data, "plate_1", 
                                               number_of_wells=96)
    path_to_file= system.file("extdata", "sample_data_96_platemap.txt",
                              package=pkg)
    plate_map=read_platemap_from_Incucyte_XML( path_to_file )
    ds_n = combine_data_and_map( raw_plate, plate_map, "Vehicle" ) 
    ds = normalize_by_vehicle(ds_n, summary_method="mean")

    # CHECK NORMALIZATION CALCULATIONS
    
    V = plyr::ddply( ds, c("sample_type", "concentration", "treatment", 
                           "plate_id", "hours"), 
                       function(x){ data.frame(
                           mu=mean(x$value, na.rm=TRUE),
                           med=median(x$value, na.rm=TRUE) ) } )
    V = V[V$hours==0,]
    
    # Check calculation of mean for vehicle, drug
    mu_v = mean(ds$value[ds$treatment=="Vehicle" & ds$hours==0])
    mu_d10625 = mean(ds$value[ds$treatment=="drug_1" & ds$hours==0 & 
                                  ds$concentration==0.0625]) 
    expect_equal(round(mu_v,5), round( V$mu[V$treatment=="Vehicle"]), 5)
    expect_equal(round(mu_d10625,5),
               round( V$mu[V$treatment=="drug_1" & V$concentration==0.0625]), 5)
    
    # check normalized values
    vals_calc = ds$value_normalized[ds$treatment=="drug_1" & 
                                        ds$hours==0 & ds$concentration==0.0625]
    vals_hand = ds$value[ds$treatment=="drug_1" & ds$hours==0 & 
                             ds$concentration==0.0625] /
        mu_v
    for(i in 1:length(vals_calc)){
        expect_equal( round(vals_calc[i],5), round(vals_hand[i],5))
    }
    
    sample_types = rep( c(rep("line1",3), rep("line2",3)), 5)
    treatments=c(rep("DMSO",6), rep("drug",24))
    concentrations = c( rep(0,6),rep(200,6), rep(500,6),rep(1000,6),rep(5000,6))
    values=c(100,99,100,90,91,92,99,97,99,89,87,88,86,89,88,56,59,58,66,65,67,
            25,23,24,42,43,46,4,5,9)
    D=create_dataset( sample_types=sample_types, 
                      treatments=treatments, 
                      concentrations=concentrations, 
                      values=values)
    expect_equal( dim(D)[1], 30 )
    D = normalize_by_vehicle(D, summary_method = "mean")
    expect_equal( sum( D$value != D$value_normalized), 0 )
    
    D=create_dataset( sample_types, 
                      treatments, 
                      concentrations, 
                      values, 
                      negative_control="DMSO")
    expect_equal( dim(D)[1], 30 )
    D = normalize_by_vehicle(D, summary_method = "mean")
    expect_equal( round( D$value_normalized[29],2 ), 0.05)
} )