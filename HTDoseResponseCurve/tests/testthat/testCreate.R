library(HTDoseResponseCurve)
context("Create plate map and dataset")

test_that("create works", {
    pkg = "HTDoseResponseCurve"

    pm = create_empty_plate_map( 6 )
    expect_equal( sum( names(pm)=="treatment"), 1 )
    expect_equal( sum( names(pm)=="concentration"), 1 )
    expect_equal( sum( names(pm)=="sample_type"), 1 )
    expect_equal( sum( names(pm)=="density"), 1 )
    expect_equal( sum( names(pm)=="passage"), 1 )
    expect_equal( dim(pm$treatment)[1], 2 )
    expect_equal( dim(pm$treatment)[2], 3 )
    
    pm_2 = create_empty_plate_map( 6, max_treatments_per_well = 2 )
    expect_equal( sum( names(pm_2)=="treatment"), 1 )
    expect_equal( sum( names(pm_2)=="concentration"), 1 )
    expect_equal( sum( names(pm_2)=="sample_type"), 1 )
    expect_equal( sum( names(pm_2)=="density"), 1 )
    expect_equal( sum( names(pm_2)=="passage"), 1 )
    expect_equal( sum( names(pm_2)=="treatment_2"), 1 )
    expect_equal( sum( names(pm_2)=="concentration_2"), 1 )
    
    expect_error( create_empty_plate_map( 1 ), 
                  "number_of_wells must be 6, 24, 96, or 384")
    
    pm = create_empty_plate_map( 6 )
    pm$treatment[1,1:3] = c("DMSO","drug1","drug2")
    pm$treatment[2,1:3] = c("DMSO","drug1","drug2")
    ds = create_dataset( 
        sample_types= c("line1","line1","line1","line2","line2","line2"),
        treatments = c("DMSO","drug1","drug2","DMSO","drug1","drug2"),
        concentrations = c(0, 100, 200, 0, 100, 200),
        hours = c(48, 48, 48, 48, 48, 48),
        values = c(100, 90, 20, 100, 89, 87), 
        plate_id = "plate_1",
        negative_control = "DMSO")
    
    ts = sort(unique(ds$treatment))
    expect_equal(ts[1], "DMSO")
    expect_equal(ts[2], "drug1")
    expect_equal(ts[3], "drug2")
    expect_equal(unique(ds$hour), 48)
    expect_equal(min(ds$value), 20)
    expect_equal(ds$value_normalized[1], 1)
    expect_equal(ds$value_normalized[2], 0.9)
    expect_equal(ds$value_normalized[5], 0.89)
    
    ts = get_treatments( ds )
    expect_equal(ts[1], "DMSO")
    expect_equal(ts[2], "drug1")
    expect_equal(ts[3], "drug2")
    
    cs = get_concentrations( ds )
    expect_equal(cs[1], 0)
    expect_equal(cs[2], 100)
    expect_equal(cs[3], 200)
    
    hs = get_hours( ds )
    expect_equal(length(hs), 1)
    expect_equal(hs[1], 48)
})