#' @import utils
plate_dimensions_from_wells = function(number_of_wells){
    if( number_of_wells != 6 & number_of_wells != 24 & 
        number_of_wells != 96 & number_of_wells != 384){
        stop("number_of_wells must be 6, 24, 96, or 384")   
    }
    if( number_of_wells == 384 ){
        list(rows=16, cols=24)
    }else if( number_of_wells == 96 ){
        list(rows=8, cols=12)
    }else if( number_of_wells == 24 ){
        list(rows=4, cols=6)
    }else if( number_of_wells == 6 ){
        list(rows=2, cols=3)
    }
}

#-------------------------------------------------------------------------------

# Function constructor for FIT object
#
HT_fit = function(...){
    FIT=list()
    FIT[[ "unique_conditions" ]] = NA
    FIT[[ "is_fitted" ]] = FALSE
    FIT[[ "sample_types" ]] = NA
    FIT[[ "treatments" ]] = NA
    FIT[[ "fit_stats" ]] = NA
    FIT[[ "input" ]] = NA
    FIT[[ "model" ]] = NA
    FIT[[ "ANOVA_F_test" ]] = NA
    FIT[[ "ANOVA_P_value" ]] = NA
    attr(FIT, "class") = "HT_fit"
    class(FIT) = "HT_fit"
    return( FIT )
}

# If two drugs or lines are on different plates, they should be fit using 
# the vehicle concentrations from their own plate as the 0 concentration.    
# for each unique sample_type + treatment, find the appropriate vehicle
#
subset_treatments = function( D, sample_types, treatments, hour ){
    idx = c()
    conditions_to_fit = c()
    for(i in 1:length(sample_types)){
        for(j in 1:length(treatments)){
            #index of this {treatment, sample, hour} 
            idx_D =  which(D$treatment==treatments[j] & 
                           D$sample_type==sample_types[i] & D$hours==hour )
            if( length(idx_D)==0 ){
                stop( paste("Requested combination of treatment=",
                             treatments[j], "sample_type =", sample_types[j], 
                              "hour =", hour, "not found in D" ))
            }
            cur_plate_ids = unique( D$plate_id[idx_D] )
            vehicle_names = unique( D$negative_control[ idx_D ] )
            idx_V = which(D$treatment %in% vehicle_names & 
                          D$sample_type == sample_types[i] & 
                          D$hours == hour & 
                          D$plate_id %in% cur_plate_ids &
                          D$is_negative_control )
            idx = c(idx, idx_D, idx_V)
            new_condition =rep( paste(sample_types[i], treatments[j],sep="_|_"), 
                                 length(idx_D)+length(idx_V))
            conditions_to_fit = c(conditions_to_fit, new_condition)
        }
    }
    data.frame( 
        value =             D$value_normalized[ idx ], 
        concentration =     D$concentration[idx], 
        sample_type =       D$sample_type[idx], 
        conditions_to_fit = conditions_to_fit, 
        stringsAsFactors=FALSE )
}


# Calculate area under all values in D for a single value of sample_type,
# treatment, and concentration using pracma:trapz. Agnostic to normalization.
#' @import pracma
#' @import plyr
AUC_across_timepoints = function( D, sample_type, treatment, concentration, 
                                  summary_method="mean"){
    if( summary_method != "mean" & summary_method != "median"){
        stop("summary_method parameter must be either mean or median")
    }
    DA = D[D$treatment==treatment & 
               D$concentration==concentration & 
               D$sample_type==sample_type,]
    Dsum = plyr::ddply( DA, c("hours"), 
                        function(x){ data.frame( 
                            mu=mean(x$value, na.rm=TRUE), 
                            med=median(x$value, na.rm=TRUE) ) }
    )
    Dsum = Dsum[order(Dsum$hours),]
    if( summary_method=="mean"){
        values=Dsum$mu
    }
    else{
        values = Dsum$med
    }
    pracma::trapz( Dsum$hours, values  )
}

# Calculate AUC using trapezoids method; works either for curves fit with drc()
# or for point estimates when curve fitting fails.
#
AUC = function(FIT, granularity=0.01){
    if(is.null(attr(FIT, "class"))){ stop("Must be called on a class") }
    else{ UseMethod("AUC") }
}
AUC.HT_fit = function(FIT, granularity=0.01, summary_method="mean"){
    obs_min = rep(NA, length(FIT$unique_conditions))
    obs_max = rep(NA, length(FIT$unique_conditions))
    coef_b = rep(NA, length(FIT$unique_conditions))
    coef_c = rep(NA, length(FIT$unique_conditions))
    coef_d = rep(NA, length(FIT$unique_conditions))
    coef_e = rep(NA, length(FIT$unique_conditions))
    auc = rep(NA, length(FIT$unique_conditions))
    
    # Summarize
    Mstat = plyr::ddply( FIT$input, c("conditions_to_fit", "concentration"), 
                         function(x){ data.frame( 
                             mu=mean(x$value, na.rm=TRUE), 
                             med=median(x$value, na.rm=TRUE) ) } )
    if(summary_method=="mean"){
        Mstat$value = Mstat$mu
    }else{
        Mstat$value = Mstat$med
    }
    
    # Record min/max observed summary data
    for( i in 1:length(FIT$unique_conditions) ){
        observed = Mstat$value[ Mstat$conditions_to_fit==
                                        FIT$unique_conditions[i] ]
        obs_min[i] = min(observed, na.rm=-TRUE)
        obs_max[i] = max(observed, na.rm=-TRUE)
    }
    
    if( !FIT$is_fitted ){    
        # If a curve was not fitted, use the summarized observed values 
        # to calculate AUC but do not try to infer other values and do not 
        # estimate EC50
        for(i in 1:length(FIT$unique_conditions)){
            cur_condition=FIT$unique_conditions[i]
            Mcur=Mstat[ Mstat$conditions_to_fit==cur_condition,]
            Mcur = Mcur[order(Mcur$concentration),]
            auc[i] = pracma::trapz(10^(Mcur$concentration), Mcur$value)
        }
    }else{
        max_concentration = log10( max(FIT$input$concentration) )
        conc_to_fit= 10^seq(from=0, to=max_concentration, by=granularity )
        conc_test = rep(conc_to_fit, length(FIT$unique_conditions) )
        cond_test = c()
        for(i in 1:length(FIT$unique_conditions)){
            cond_test = c( cond_test, 
                           rep(FIT$unique_conditions[i],length(conc_to_fit)))
        }
        model_input_test = data.frame( concentration=conc_test, 
                                       conditions_to_fit=factor(cond_test),
                                       stringsAsFactors=FALSE ) 
        preds = try( predict(FIT$model, model_input_test), silent=TRUE )
        
        if( !"try-error" %in% class( preds ) ){
            for( i in 1:length(FIT$unique_conditions) ){
                cur_cond=FIT$unique_conditions[i]
                observed = FIT$input$value[ FIT$input$conditions_to_fit==
                                                cur_cond ]
                idx_first = min( which(cond_test==cur_cond ) )
                idx_last =  max( which(cond_test==cur_cond ) )
                ys = preds[ idx_first : idx_last ]
                ys[ys<0]=0 # No negative values allowed for AUC
                idx_b = which(names(coefficients(FIT$model))==
                                  paste("b",cur_cond,sep=":"))
                idx_c = which(names(coefficients(FIT$model))==
                                  paste("c",cur_cond,sep=":"))
                idx_d = which(names(coefficients(FIT$model))==
                                  paste("d",cur_cond,sep=":"))
                idx_e = which(names(coefficients(FIT$model))==
                                  paste("e",cur_cond,sep=":"))
                if(length(idx_b)>0){
                    coef_b[i] = as.numeric( coefficients(FIT$model)[idx_b])
                }
                if(length(idx_c)>0){
                    coef_c[i] = as.numeric( coefficients(FIT$model)[idx_c])
                }else{
                    coef_c[i] = 0
                }
                if(length(idx_d)>0){
                    coef_d[i] = as.numeric( coefficients(FIT$model)[idx_d])
                }else{
                    coef_d[i] = 1
                }
                if(length(idx_e)>0){
                    coef_e[i] = as.numeric( coefficients(FIT$model)[idx_e])
                }
                auc[i] = pracma::trapz(conc_to_fit, ys)
                if( min( observed, na.rm=TRUE) > 0.5 ){
                    coef_e[i] = Inf
                }
            }
        }
    }
    
    FIT$fit_stats = data.frame( AUC=auc, 
                                coef_slope = round( coef_b, 3 ), 
                                coef_asymp_low = round( coef_c, 3 ), 
                                coef_asymp_high = round( coef_d, 3 ),
                                coef_EC50 = round( coef_e, 3),
                                obs_min = round( obs_min, 3), 
                                obs_max = round( obs_max, 3 ), 
                                row.names = FIT$unique_conditions )
    FIT
}



#-------------------------------------------------------------------------------

#' Create a plate list with empty elements
#' 
#' 
#' @param number_of_wells plate size, one of {6, 24, 96, 384}
#' @param hour hours since the experiment began
#' @param plate_id string that identifies this plate map
#' @return a data frame 
#' @examples 
#' create_empty_plate( 96, hour=0, plate_id="plate_1")
#' create_empty_plate( 384, hour=12, plate_id="plate_9")
#' @seealso \code{\link{create_empty_plate_map}} 
#' @export
create_empty_plate = function( number_of_wells, hour=0, plate_id ){
    PLATE_ROWS = plate_dimensions_from_wells(number_of_wells)$rows
    PLATE_COLS = plate_dimensions_from_wells(number_of_wells)$cols
    if( !is.character( plate_id ) ){
        stop( "parameter plate_id must be a character string")   
    }
    if( !is.numeric( hour ) ){
        stop( "parameter hour must be a number")
    }
    plate = data.frame(matrix(NA, nrow=PLATE_ROWS, ncol=PLATE_COLS))
    plate = cbind(plate, hours=rep(hour, PLATE_ROWS), 
                         plate_id=rep(plate_id, PLATE_ROWS), 
                         stringsAsFactors=FALSE)
    names(plate)[1:PLATE_COLS] = 1:PLATE_COLS
    dimnames(plate)[[1]] = abc[1:PLATE_ROWS]
    plate
}


#' Create an empty plate map 
#' 
#' The plate map object is a list of identically sized matrixes named 
#' treatment, concentration, sample_type, density, and passage. 
#' 
#' @param number_of_wells plate size, one of {6, 24, 96, 384}
#' @return a plate_map list with matrixes for treatment, concentration, 
#' sample_types, density, and passage
#' @seealso \code{\link{create_empty_plate}} 
#' @examples 
#' create_empty_plate_map( 96 )
#' @export
create_empty_plate_map = function( number_of_wells ){
    PLATE_ROWS = plate_dimensions_from_wells(number_of_wells)$rows
    PLATE_COLS = plate_dimensions_from_wells(number_of_wells)$cols 
    list(
        treatment = matrix("", ncol=PLATE_COLS, nrow=PLATE_ROWS),
        concentration = matrix(NA, ncol=PLATE_COLS, nrow=PLATE_ROWS),
        sample_type = matrix("", ncol=PLATE_COLS, nrow=PLATE_ROWS),
        density = matrix("", ncol=PLATE_COLS, nrow=PLATE_ROWS),
        passage = matrix("", ncol=PLATE_COLS, nrow=PLATE_ROWS)
    )
}


#' Sorted unique list of sample types
#' 
#' Convenience method to extract sample types. get_sample_types() can be called 
#' on either a plate map list or a data frame generated by 
#' \code{combine_data_and_map}.
#' 
#' @param D A plate map list or a data frame generated by 
#' \code{combine_data_and_map}
#' @return sorted vector of sample types
#' @examples 
#' pm = create_empty_plate_map( 6 )
#' pm$sample_type[1,1:3] = "line1"
#' pm$sample_type[2,1:3] = "line2"
#' get_sample_types( pm )
#' ds = create_dataset( 
#'   sample_types= c("line1","line1","line1","line2","line2","line2"),
#'   treatments = c("DMSO","drug1","drug2","DMSO","drug1","drug2"),
#'   concentrations = c(0, 100, 200, 0, 100, 200),
#'   hours = c(48, 48, 48, 48, 48, 48),
#'   values = c(100, 90, 20, 100, 89, 87), 
#'   plate_id = "plate_1",
#'   negative_control = "DMSO")
#' get_sample_types( ds )
#' @seealso  \code{\link{get_treatments}}
#' @seealso  \code{\link{get_concentrations}}
#' @seealso  \code{\link{get_hours}}
#' @export
get_sample_types = function( D ){
    setdiff( sort(unique(as.vector(t(D$sample_type)))), "" )
}

#' Sorted unique list of treatments
#'
#' Convenience method to extract treatments
#' @param D A plate map list or a data frame generated by 
#' \code{combine_data_and_map}
#' @return sorted vector of treatments
#' @examples 
#' pm = create_empty_plate_map( 6 )
#' pm$treatment[1,1:3] = c("DMSO","drug1","drug2")
#' pm$treatment[2,1:3] = c("DMSO","drug1","drug2")
#' get_treatments( pm )
#' ds = create_dataset( 
#'   sample_types= c("line1","line1","line1","line2","line2","line2"),
#'   treatments = c("DMSO","drug1","drug2","DMSO","drug1","drug2"),
#'   concentrations = c(0, 100, 200, 0, 100, 200),
#'   hours = c(48, 48, 48, 48, 48, 48),
#'   values = c(100, 90, 20, 100, 89, 87), 
#'   plate_id = "plate_1",
#'   negative_control = "DMSO")
#' get_treatments( ds )
#' @seealso  \code{\link{get_sample_types}}
#' @seealso  \code{\link{get_concentrations}}
#' @seealso  \code{\link{get_hours}}
#' @export
get_treatments = function( D ){
    setdiff( sort(unique(as.vector(t(D$treatment)))), "")
}

#' Sorted unique list of treatment concentrations
#' 
#' Convenience method to extract treatment concentrations
#' 
#' @param D A plate map list or a data frame generated by 
#' \code{combine_data_and_map}
#' @return sorted vector of concentrations, including 0 if present
#' @examples 
#' pm = create_empty_plate_map( 6 )
#' pm$concentration[1,1:3] = c( 0, 100, 200 )
#' pm$concentration[2,1:3] = c( 0, 100, 300 )
#' get_concentrations( pm )
#' ds = create_dataset( 
#'   sample_types= c("line1","line1","line1","line2","line2","line2"),
#'   treatments = c("DMSO","drug1","drug2","DMSO","drug1","drug2"),
#'   concentrations = c(0, 100, 200, 0, 100, 200),
#'   hours = c(48, 48, 48, 48, 48, 48),
#'   values = c(100, 90, 20, 100, 89, 87), 
#'   plate_id = "plate_1",
#'   negative_control = "DMSO")
#' get_concentrations( ds )
#' @seealso  \code{\link{get_sample_types}}
#' @seealso  \code{\link{get_treatments}}
#' @seealso  \code{\link{get_hours}}
#' @export
get_concentrations = function( D ){
    setdiff( sort(unique(as.vector(t(D$concentration)))), "")
}

#' Sorted unique list of hours since time zero
#' 
#' Convenience method to extract a sorted list of the hours 
#' 
#' @param D A dataframe generated by \code{combine_data_and_map}
#' @return sorted vector of hours
#' @examples 
#' ds = create_dataset( 
#'   sample_types= c("line1","line1","line1","line2","line2","line2"),
#'   treatments = c("DMSO","drug1","drug2","DMSO","drug1","drug2"),
#'   concentrations = c(0, 100, 200, 0, 100, 200),
#'   hours = c(48, 48, 48, 48, 48, 48),
#'   values = c(100, 90, 20, 100, 89, 87), 
#'   plate_id = "plate_1",
#'   negative_control = "DMSO")
#' get_hours( ds )
#' @seealso  \code{\link{get_sample_types}}
#' @seealso  \code{\link{get_treatments}}
#' @seealso  \code{\link{get_hours}}
#' @export
get_hours = function( D ){
    sort(unique(D$hours))
}


# helper function to normalize data in D relative to vehicle
prepare_for_normalization = function( D, negative_control ){
    
    # is_neg_ctl is required for the case where a drug at concentration 0 is 
    # its own negative control; we need a way to remember which of the wells 
    # to summarize for the negative control
    is_neg_ctl = rep( FALSE, dim(D)[1] )
    neg_ctl = rep(NA, dim(D)[1] )
    
    if( is.numeric(negative_control) ){
        neg_ctl = D$treatment
        treatments = unique(D$treatment)
        # Check treatments have a well with concentration==negative_control
        for(i in 1:length(treatments)){
            idx = which( D$concentration[ D$treatment==treatments[i] ]==
                             negative_control )
            if( length( idx ) == 0 ){
                stop(paste("negative_control parameter passed with value",
                           negative_control,"so we assume each",
                           "treatment has one or more wells with concentration",
                           "=",negative_control, "which is not the case for", 
                           treatments[i] ) )
            }
            is_neg_ctl[ D$concentration==negative_control & 
                        D$treatment==treatments[i] ] = TRUE
            neg_ctl[ D$treatment==treatments[i] ] = treatments[i]
        }
        
    }else if( is.character(negative_control) ){
        if( length( intersect( negative_control, D$treatment) ) != 1 ){
            stop(paste("negative_control", negative_control,
                       "is not among the treatments on this plate"))
        }
        neg_ctl = rep( negative_control, dim(D)[1])
        is_neg_ctl = D$treatment==negative_control
        
    }else if( is.data.frame( negative_control ) ){
        VD = negative_control
        if( sum(names(VD)=="drug" ) == 0 | 
            sum(names(VD)=="vehicle") == 0 ){
            stop(paste("if passed as a data frame, parameter negative_control",
                       "must have columns named 'drug' and 'vehicle'"))
        }
        drugs = unique( VD$drug )
        vehicles = unique( VD$vehicle )
        d2v = hsh_from_vectors( VD$drug, VD$vehicle )
        # vehicles are their own vehicles
        for(i in 1:length(vehicles)){
            hsh_set( d2v, VD$vehicle[i], VD$vehicle[i] ) 
        }
        treatments = c(drugs, vehicles)
        for( i in 1:length(treatments)){
            neg_ctl[ which(D$treatment==treatments[i]) ]=
                hsh_get(d2v,treatments[i])
        }
        
        for( i in 1:length(vehicles)){
            is_neg_ctl[ which(D$treatment==vehicles[i]) ] = TRUE
        }
        if( sum(is.na(neg_ctl))>0 ){
            stop(paste( "Not all drugs have been assigned a",
                        "negative control, check negative_control parameter"))
        }
        if( length( setdiff(neg_ctl, D$treatment)>0 ) ){
            stop( paste( "parameter negative_control contain an element in",
                         "the vehicle column that is not found in the",
                         "treatments for this plate") )
        }
    }
    D$is_negative_control = is_neg_ctl
    D$negative_control = neg_ctl
    D
}


#' Create a dataset from raw data without plate  
#' 
#' @param sample_types vector of sample types
#' @param treatments vector of treatments
#' @param concentrations vector of concentrations
#' @param hours vector of timepoints
#' @param values vector of measured response to treatment
#' @param negative_control If an empty string, there are no negative control 
#' wells. If a number such as 0, each treatment will be expected to contain 
#' wells with this concentration as a negative control. If a string, the 
#' name of treatment wells which contain the universal negative control for this 
#' plate. If a data frame, a mapping of each vehicle to each drug individually.
#' @param plate_id experiment identification string, useful if multiple datasets 
#' are later combined.
#' @return a data frame where columns indicate the sample type, treatment, 
#' concentration, observed raw value, normalized value (raw until normalization 
#' is run), name of the negative_control treatment, whether a particular row 
#' is a negative control for at least one other row, hours since the start time,
#' and plate of origin 
#' @examples 
#' # six measurements: DMSO, 100, and 200 nM for two drugs. 
#' # plan to normalize each line against DMSO for that line
#' ds = create_dataset( 
#'   sample_types= c("line1","line1","line1","line2","line2","line2"),
#'   treatments = c("DMSO","drug1","drug2","DMSO","drug1","drug2"),
#'   concentrations = c(0, 100, 200, 0, 100, 200),
#'   hours = c(48, 48, 48, 48, 48, 48),
#'   values = c(98, 90, 20, 99, 89, 87), 
#'   plate_id = "plate_1",
#'   negative_control = "DMSO")
#'   
#' # six measurements; drug1 at 0, 100, 200 nM and drug2 at 0, 100, 200 nM. 
#' # plan to normalize against zero concentration for each line
#' ds = create_dataset( 
#'   sample_types= c("line1","line1","line1","line2","line2","line2"),
#'   treatments = c("drug1","drug1","drug2","drug2","drug2","drug2"),
#'   concentrations = c(0, 100, 200, 0, 100, 200),
#'   hours = c(48, 48, 48, 48, 48, 48),
#'   values = c(98, 90, 20, 99, 89, 87), 
#'   plate_id = "plate_1",
#'   negative_control = 0)
#'   
#' # six measurements; drug1 at 0, 100, 200 nM and drug2 at 0, 100, 200 nM. 
#' # plan to normalize drug1 against DMSO and drug2 against ethanol
#' individual_vehicles = data.frame(
#'   drug=c("drug1", "drug2"), 
#'   vehicle=c("DMSO", "ethanol"),
#'   stringsAsFactors=FALSE)
#' ds = create_dataset( 
#'   sample_types= c("line1","line1","line1","line2","line2","line2"),
#'   treatments = c("DMSO","drug1","drug1","ethanol","drug2","drug2"),
#'   concentrations = c(0, 100, 200, 0, 100, 200),
#'   hours = c(48, 48, 48, 48, 48, 48),
#'   values = c(98, 90, 20, 99, 89, 87), 
#'   plate_id = "plate_1",
#'   negative_control = individual_vehicles)
#' @export
create_dataset = function( sample_types, treatments, concentrations, hours,
                           values, plate_id, negative_control = NA){
    if( length(sample_types) != length(treatments) |
        length(treatments) != length(concentrations) | 
        length(concentrations) != length(hours) | 
        length(hours) != length(values) ){
        stop( paste("parameters sample_types, treatments, concentrations,",
                    "values, and hours must have the same length"))
    }
    if( length(plate_id) != 1 ){
        stop("parameter plate_id should be a single string")   
    }
    if( sum( !is.numeric(hours))>0 ){
        stop( "Hours must be a vector of numbers")
    }
    df = data.frame( sample_type=sample_types, 
                     treatment=treatments, 
                     concentration=concentrations, 
                     hours, 
                     value=values, 
                     value_normalized=values,
                     plate_id = rep(plate_id, length(treatments)),
                     negative_control = rep(NA, length(values) ),
                     is_negative_control = rep(NA, length(values) ),
                     stringsAsFactors=FALSE )
    prepare_for_normalization( df, negative_control )
}


#' Combine raw data and plate map to produce a single tall data frame
#'
#' Combine a single data plate matching the output of 
#' \code{\link{create_empty_plate}} with a plate map matching the output of 
#' \code{\link{create_empty_plate_map}} to create a single tall data frame with 
#' one row per observation. 
#' 
#' This code does not normalize the combined data, but the user specifies which 
#' wells (if any) bear negative controls (i.e. vehicles) that can be normalized 
#' against using the \code{\link{normalize_plates_by_vehicle}} function.
#' Negative controls:
#' 
#' If the experiment lacks negative controls, pass an empty string to the 
#' negative_control parameter.
#' 
#' If each treatment has wells with a concentration in the plate map that 
#' indicates the negative control for that treatment, pass that value to the 
#' negative_control parameter. The default value is 0.
#' 
#' If more all wells share a common negative control (e.g. DMSO), pass the name 
#' of the negative control treatment to the negative_control parameter. 
#' 
#' If there is more than one distinct common negative control (e.g. DMSO for 
#' some treatments and ethanol for others), pass a data frame with column names 
#' "drug" and "vehicle" to indicate the appropriate vehicle for each drug.
#' If a data frame is passed, the drug column should contain an entry for each 
#' treatment on the plate(s). 
#' 
#' @param raw_plate data frame with format matching that produced by a call to 
#' \code{read_incucyte_exported_to_excel} 
#' @param plate_map list where each item specifies concentrations, treatments, 
#' and cell lines for a plate specified in raw_plates, matching the format of a 
#' file created by \code{read_platemap_from_excel} 
#' @param negative_control If an empty string, there are no negative control 
#' wells. If a number such as 0, each treatment will be expected to contain 
#' wells with this concentration as a negative control. If a string, the 
#' name of treatment wells which contain the universal negative control for this 
#' plate. If a data frame, a mapping of each vehicle to each drug individually.
#' @return a data frame where columns indicate the sample type, treatment, 
#' concentration, observed raw value, normalized value (raw until normalization 
#' is run), name of the negative_control treatment, whether a particular row 
#' is a negative control for at least one other row, hours since the start time,
#' and plate of origin 
#' @examples 
#' # Create a six well plate testing two drugs on two cell lines. Each drug has 
#' # two concentrations plus a DMSO well. Plan to normalize against DMSO.
#' # Obviously real data would include replicates, more concentrations, etc.
#' plate = create_empty_plate( 6, hour=0, plate_id="plate_1")
#' plate[1:2,1:3] = c(99,98,90,87,77,20)
#' plate_map = create_empty_plate_map( number_of_wells = 6 )
#' plate_map$concentration[1:2,1:3] = c( 0, 0, 100, 100, 200, 200 )
#' plate_map$treatment[1:2,1:3] = c("DMSO","DMSO","drug1","drug2","drug1",
#'                                  "drug2")
#' plate_map$sample_type[1:2,1:3] = rep( c("line1", "line2"), 3)
#' combine_data_and_map( plate, plate_map, negative_control="DMSO" )
#' 
#' # Create a six well plate testing two drugs on two cell lines. Each drug has 
#' # three concentrations, one of which is vehicle-only and labeled zero.
#' plate = create_empty_plate( 6, hour=0, plate_id="plate_1")
#' plate[1:2,1:3] = c(99,98,90,87,77,20)
#' plate_map = create_empty_plate_map( number_of_wells = 6 )
#' plate_map$concentration[1:2,1:3] = c( 0, 0, 100, 100, 200, 200 )
#' plate_map$treatment[1:2,1:3] = c("drug1","drug2","drug1","drug2","drug1",
#'                                  "drug2")
#' plate_map$sample_type[1:2,1:3] = rep( c("line1", "line2"), 3)
#' combine_data_and_map( plate, plate_map, negative_control=0 )
#' @export
combine_data_and_map = function( raw_plate, plate_map, negative_control=0 ){
    
    COLS = dim(plate_map$concentration)[2]
    ROWS = dim(plate_map$concentration)[1]
    
    if( dim(raw_plate)[2] != COLS + 2 ){
        stop("raw_plate must have two more columns than plate map")    
    }
    N_data_cols = COLS
    
    if( length( which(names(raw_plate)=="hours") )==0 ){
        stop("raw_plate must have hours column")
    }
    if( length( which(names(raw_plate)=="plate_id") )==0 ){
        stop("raw_plate must have plate_id column")
    }
    
    if( length( which(names(plate_map)=="concentration") )==0 ){
        stop("plate_map must have concentration matrix")
    }
    if( length( which(names(plate_map)=="treatment") )==0 ){
        stop("plate_map must have treatment matrix")
    }
    if( length( which(names(plate_map)=="sample_type") )==0 ){
        stop("plate_map must have sample_type matrix")
    }
    if( dim(plate_map$treatment)[1] != ROWS | 
        dim(plate_map$treatment)[2] != N_data_cols ){
        stop("plate_map treatment matrix dimensions don't match raw_data")
    }
    if( dim(plate_map$sample_type)[1] != ROWS | 
        dim(plate_map$sample_type)[2] != N_data_cols ){
        stop("plate_map sample_type matrix dimensions don't match raw_data")
    }    
    N_obs = sum( !is.na( raw_plate[, 1:N_data_cols] ) )
    sample_type = rep("", N_obs)
    treatment = rep("", N_obs)
    concentration = rep(0, N_obs)
    value = rep("", N_obs)
    hours = rep(0, N_obs)
    plate_ids = rep("", N_obs)
    idx=1
    rr_p = 0
    # dim(raw_plate)[1] is number of rows per plate * number of timepoints
    for(rr in 1:dim(raw_plate)[1]){
        rr_p = rr_p + 1
        if(rr_p==ROWS+1){
            rr_p = 1 # roll back to first row of plate
        }
        for(cc in 1:N_data_cols){
            if( !is.na(raw_plate[rr,cc]) ){
               sample_type[idx] =   plate_map[["sample_type"]][rr_p, cc]
               treatment[idx] =     plate_map[["treatment"]][rr_p, cc]
               concentration[idx] = plate_map[["concentration"]][rr_p, cc]
               value[idx] =      raw_plate[rr,cc]
               hours[idx] =      raw_plate$hours[rr]
               plate_ids[idx] =  raw_plate$plate_id[rr]
               idx = idx+1
            }
        }
    }
    concentration = as.numeric( concentration )
    
    D = data.frame( 
        sample_type, 
        treatment, 
        concentration, 
        value = as.numeric( value ),
        value_normalized = as.numeric( value ),
        negative_control = rep(NA, length(value) ),
        is_negative_control = rep(NA, length(value) ),
        hours = as.numeric(hours),
        plate_id = plate_ids, 
        stringsAsFactors=FALSE
    )
    prepare_for_normalization( D, negative_control )
}


#' Normalize raw data on each plate against the vehicle control wells
#'
#' Given a data frame of measurements generated by 
#' \code{\link{combine_data_and_map}}
#' or matching that format, divide each raw measurement on the specified plate 
#' by the summarized value for the appropriate negative control measurement on 
#' the same plate. Summary can be either mean or median. If the plate_id 
#' parameter is an empty string (default), all plates will be normalized. 
#'
#' If no normalization has been specified for D, the value_normalized column 
#' will be identical to the value column.
#' 
#' @param D experiment dataset with columns matching the output of 
#' \code{\link{combine_data_and_map}}
#' @param summary_method string indicating how to summarize vehicle replicates, 
#' one of "mean" or "median"
#' @return a data frame with the same columns and dimenson as D
#' @examples 
#' # Create a six well plate testing two drugs on two cell lines. Each drug has 
#' # three concentrations, one of which is vehicle-only and labeled zero.
#' plate = create_empty_plate( 6, hour=0, plate_id="plate_1")
#' plate[1:2,1:3] = c(99,98,90,87,77,20)
#' plate_map = create_empty_plate_map( number_of_wells = 6 )
#' plate_map$concentration[1:2,1:3] = c( 0, 0, 100, 100, 200, 200 )
#' plate_map$treatment[1:2,1:3] = c("drug1","drug2","drug1","drug2","drug1",
#'                                  "drug2")
#' plate_map$sample_type[1:2,1:3] = rep( c("line1", "line2"), 3)
#' ds = combine_data_and_map( plate, plate_map, negative_control=0 )
#' normalize_plates_by_vehicle( ds, summary_method="mean")
#' @export
normalize_plates_by_vehicle = function( D, summary_method ){
    if( summary_method != "mean" & summary_method != "median"){
        stop("summary_method parameter must be either mean or median")
    }
    
    D$value_normalized = D$value
    
    plates = sort(unique(D$plate_id))
    
    for( pp in 1:length(plates)){
        plate_cur = plates[pp]
        hours = unique( D$hours[ D$plate_id == plate_cur ] ) 
        for( tt in 1:length( hours ) ){
            hour_cur = hours[tt]
            vehicles=unique( D$negative_control[ D$hours==hour_cur & 
                                                 D$plate_id == plate_cur] )
            for(vv in 1:length(vehicles)){
                vehicle_cur = vehicles[vv]
                D_v = D[D$hours==hour_cur &
                        D$plate_id==plate_cur &
                        D$treatment==vehicle_cur & 
                        D$is_negative_control,]
                D_v = plyr::ddply( D_v, c("sample_type", "hours"), 
                                   function(x){ data.frame(
                                       mu=mean(x$value, na.rm=TRUE),
                                       med=median(x$value, na.rm=TRUE) ) } 
                )
                if( summary_method=="mean" ){
                    D_v$value = D_v$mu
                }else{
                    D_v$value = D_v$med
                }
                for(i in 1:length( D_v$sample_type ) ){
                    sample_type_cur = D_v$sample_type[i]
                    idx = which( D$hours==hour_cur & 
                                 D$plate==plate_cur & 
                                 D$negative_control==vehicle_cur &
                                 D$sample_type==sample_type_cur )
                    D$value_normalized[idx] = D$value[idx] / D_v$value[i]
                } 
            }
        }
    }
    D
}


    
#' reshape a set of values to a matrix with axes defined by X any Y vectors
#' 
#' Useful for extracting two columns out of a data frame that will serve as 
#' column and row names, with a third column that will represent the value at 
#' that (row, column). The pairs identified by {x_axis_labels,y_axis_labels} 
#' should be unique; if not, the last value seen will be used.
#' 
#' @param x_axis_labels vector of values that define the X axis
#' @param y_axis_labels vector of values that define the Y axis
#' @param values value to position at (X,Y) defined by the labels
#' @examples 
#' xlab = c("a","a","a","b","b","b","c","c","c")
#' ylab = c("x","y","z","x","y","z","x","y","z")
#' vals = c(1, 2 , 3 , 4 , 5 , 6 , 7 , 8, 9)
#' metric_to_grid(xlab, ylab, vals)
#' @return Matrix with columns named for x_axis_labels, rows named for 
#' y_axis_labels, and elements corresponding to the values parameter. X and Y 
#' labels will be sorted.
#' @export
metric_to_grid = function( x_axis_labels, y_axis_labels, values ){
    
    if( length(x_axis_labels) != length(y_axis_labels) ){
        stop("axis labels must have the same length")   
    }
    if( length(x_axis_labels) != length(values) ){
        stop("values must have the same length as axis labels")   
    }
    # sort first, then convert to strings. If convert first, labels that are 
    # numbers will be sorted incorrectly.
    xl = sort(unique(x_axis_labels))
    yl = sort(unique(y_axis_labels))
    
    x_axis_labels=as.character(x_axis_labels)
    y_axis_labels=as.character(y_axis_labels)
    xl=as.character(xl)
    yl=as.character(yl)
    xl2idx = hsh_from_vectors( xl, 1:length(xl))
    yl2idx = hsh_from_vectors( yl, 1:length(yl))
    G = matrix( NA, nrow=length(yl), ncol=length(xl), dimnames=
                    list(yl, xl))
    for(i in 1:length(x_axis_labels)){
        G[ hsh_get(yl2idx, y_axis_labels[i] ), 
           hsh_get(xl2idx, x_axis_labels[i] )] = values[i] 
    }
    G
}


#' Attempt to fit a dose-response curve from the dataset and calculate AUC.
#'
#' Given a data frame of measurements generated by 
#' \code{\link{combine_data_and_map}} or matching the format generated by that 
#' function, fit a dose-response curve for each unique sample_type/treatment 
#' specified by the sample_types and treatments parameters. A curve will 
#' be fit for each sample in sample_types, at each treatment in treatments
#' 
#' @section Non-linear function for curve fitting:
#' 
#' Curve fitting is performed by the \code{drm()} function in the \code{drc} 
#' library. To fit the curve, you need to select a non-linear function. To 
#' estimate the slope, upper asymptote, lower asymptote, and EC50, pass 
#' drc::LL.4(). To fix the lower asymptote at 1 and estimate the other 
#' parameters, pass drc::LL.3(). To fix the upper asympotote at 1 and the lower 
#' asymptote at 0, pass dcr::LL.2. For a list of available functions, see 
#' \code{drc::getMeanFunctions()}. 
#' 
#' To call this function you must load the drc package in your R session.
#' @param D experiment dataset with columns matching the output of 
#' \code{combine_data_and_map}
#' @param sample_types sample types (e.g. distinct cell lines) to fit
#' @param treatments treatments to fit
#' @param hour hour at which to fit
#' @param fct Non-linear function to fit, e.g. drc::LL.3(). See summary. 
#' @examples 
#' sample_types = rep( c(rep("line1",3), rep("line2",3)), 5)
#' treatments = c(rep("DMSO",6), rep("drug",24))
#' concentrations = c( rep(0,6),rep(200,6), rep(500,6),rep(1000,6),rep(5000,6))
#' values=c(100,99,100,90,91,92,99,97,99,89,87,88,86,89,88,56,59,58,66,65,67,
#'          25,23,24,42,43,46,4,5,9)
#' hours = rep(48, length(values))
#' plate_id = "plate_1"
#' ds = create_dataset( sample_types, treatments, concentrations, 
#'                       hours, values, plate_id, negative_control = "DMSO")
#' ds = normalize_plates_by_vehicle(ds, summary_method = "mean")
#' library(drc)
#' # Fit model using three-parameter log-logistic function
#' fit_DRC(ds, sample_types=c("line1", "line2"), treatments=c("drug"), 
#'         hour = 48, fct=drc::LL.3() )
#'  
#' # Fit model using four-parameter log-logistic function
#' fit_DRC(ds, sample_types=c("line1", "line2"), treatments=c("drug"), 
#'         hour = 48, fct=drc::LL.4() )
#' @return A HT_fit object
#' @import drc
#' @export
fit_DRC = function(D, sample_types, treatments, hour, fct ){
    
    for(i in 1:length(sample_types)){
        if( sum(D$sample_type==sample_types[i])==0){
            stop(paste("sample_types parameter",sample_types[i],
                       "not present in D"))
        }
    }
    for(i in 1:length(treatments)){
        if( sum(D$treatment==treatments[i])==0 ) 
            stop(paste("treatments parameter",treatments[i],"not present in D"))
    }
    if( sum( D$hours==hour )==0 ){
        stop(paste("hour parameter", hour, "not present in data frame D"))
    }
    FIT = HT_fit()
    FIT$input = subset_treatments( D, sample_types, treatments, hour )
    FIT$unique_conditions = sort( unique( FIT$input$conditions_to_fit ) ) 
    FIT$sample_types = sample_types
    FIT$treatments = treatments
    drm_settings = drc::drmc(errorm = FALSE)
    model_input = FIT$input
    model_input$conditions_to_fit = factor(model_input$conditions_to_fit)
    if( length( unique( model_input$conditions_to_fit ) == 1 ) ){
        FIT$model = try( drc::drm( value~concentration, 
                                   model_input$conditions_to_fit, 
                                   data=model_input, 
                                   fct = fct) )
    }
    else{
        FIT$model = try( 
            drc::drm( value~concentration, 
                      model_input$conditions_to_fit, 
                      data=model_input, 
                      fct = fct), 
            silent=TRUE )
    }
    if( !"try-error" %in% class( FIT$model ) ){
        FIT$is_fitted=TRUE
        # Difference between fits F statistic and P value
        M = try( drc::drm(value~concentration, data=model_input, fct = fct), 
                 silent=TRUE )
        if( !"try-error" %in% class( M ) ){
            M_anova = anova( FIT$model, M, details=FALSE )
            FIT$ANOVA_F_test = M_anova$`F value`[2]
            FIT$ANOVA_P_value = M_anova$`p value`[2]
            FIT$ANOVA=M_anova
        }
    }
    FIT=AUC( FIT )
    drm_settings = drc::drmc( errorm = TRUE )
    FIT
}

#' Summarize the results of a dose response fit
#' 
#' @param object fit object
#' @param ... other parameters, currently ignored
#' @examples 
#' sample_types = rep( c(rep("line1",3), rep("line2",3)), 5)
#' treatments = c(rep("DMSO",6), rep("drug",24))
#' concentrations = c( rep(0,6),rep(200,6), rep(500,6),rep(1000,6),rep(5000,6))
#' values=c(100,99,100,90,91,92,99,97,99,89,87,88,86,89,88,56,59,58,66,65,67,
#'          25,23,24,42,43,46,4,5,9)
#' hours = rep(48, length(values))
#' plate_id = "plate_1"
#' ds = create_dataset( sample_types, treatments, concentrations, 
#'                       hours, values, plate_id, negative_control = "DMSO")
#' ds = normalize_plates_by_vehicle(ds, summary_method = "mean")
#' library(drc)
#' # Fit model using three-parameter log-logistic function
#' fit_1 = fit_DRC(ds, sample_types=c("line1", "line2"), treatments=c("drug"), 
#'                 hour = 48, fct=drc::LL.3() )
#' summary(fit_1)
#' @return none
#' @export
summary.HT_fit = function(object, ...){ 
    cat("Sample types: ")
    st = sort(unique(object$sample_types))
    tt = sort(unique(object$treatments))
    cat( paste( paste(st, collapse=' '), "\n", sep=""))
    cat("Treatments: ")
    cat( paste( paste(tt, collapse=' '), "\n", sep=""))
    cat("Unique Conditions: ")
    cat( paste( paste(object$unique_conditions, collapse=' '), "\n", sep=""))
    
    if( object$is_fitted ){
        cat("\nFit converged.\n")
        print( round(object$fit_stats,3) )
        cat( paste("F-statistic: ",round(object$ANOVA_F_test,2), 
                     " p-value: ", signif(object$ANOVA_P_value,3),"\n", 
                   collapse=""))
    }else{
        cat("\nFit did not converge\n")
    }
}


#' Fit all observations and return summary information from fit_stats data frame
#' 
#' @section Non-linear function for curve fitting:
#' 
#' Curve fitting is performed by the \code{drm()} function in the \code{drc} 
#' library. To fit the curve, you need to select a non-linear function. To 
#' estimate the slope, upper asymptote, lower asymptote, and EC50, pass 
#' drc::LL.4(). To fix the lower asymptote at 1 and estimate the other 
#' parameters, pass drc::LL.3(). To fix the upper asympotote at 1 and the lower 
#' asymptote at 0, pass dcr::LL.2. For a list of available functions, see 
#' \code{drc::getMeanFunctions()}. 
#' 
#' @param D dataset
#' @param fct Non-linear function to fit, e.g. drc::LL.3(). See summary.
#' @return data frame with fit summary data (is_fitted, hour, treatment, 
#' sample_type, fit_stats, ANOVA_F_test, ANOVA_P_value)
#' @examples 
#' sample_types = rep( c(rep("line1",3), rep("line2",3)), 5)
#' treatments = c(rep("DMSO",6), rep("drug",24))
#' concentrations = c( rep(0,6),rep(200,6), rep(500,6),rep(1000,6),rep(5000,6))
#' values=c(100,99,100,90,91,92,99,97,99,89,87,88,86,89,88,56,59,58,66,65,67,
#'          25,23,24,42,43,46,4,5,9)
#' hours = rep(48, length(values))
#' plate_id = "plate_1"
#' ds = create_dataset( sample_types, treatments, concentrations, 
#'                       hours, values, plate_id, negative_control = "DMSO")
#' ds = normalize_plates_by_vehicle(ds, summary_method = "mean")
#' library(drc)
#' fit_statistics( ds, fct=LL.3() )
#' @seealso \code{\link{fit_DRC}} 
#' @export
fit_statistics = function(D, fct){
    sample_types = get_sample_types(D)
    treatments = sort(unique(D$treatment[!D$is_negative_control]))
    res = c()
    for( tt in 1:length(treatments ) ){
        cur_treat = treatments[tt]
        hours = get_hours(D[D$treatment==cur_treat,])
        for(hh in 1:length(hours) ){
            cur_hour = hours[hh]
            FIT=fit_DRC(D, sample_types, cur_treat, cur_hour, fct=fct)
            rn=rownames(FIT$fit_stats)
            treatment = get.split.col(rn, "_|_", last=TRUE)
            sample_type = get.split.col(rn, "_|_", first=TRUE)
            summ = cbind( is_fitted = rep(FIT$is_fitted, length(treatment)),
                          hour = rep(cur_hour, length(treatment)),
                          treatment, 
                          sample_type, 
                          FIT$fit_stats, 
                        ANOVA_F_test = rep(FIT$ANOVA_F_test, length(treatment)),
                        ANOVA_P_value=rep(FIT$ANOVA_P_value, length(treatment)),
                          stringsAsFactors=FALSE)
            if(is.null(dim(res)) ){
                res = summ
            }else{
                res = rbind(res, summ)
            }
        }
    }
    res
}


#' Convert fit statistics data frame to matrixes for further analysis
#' 
#' This function is convenient for producing heat maps or other plots of 
#' fit summary data across a time series.
#' 
#' @param fits data frame produced by call to \code{\link{fit_statistics}}
#' @return list of lists, indexed first by sample type and then by 
#' ANOVA_P_value, AUC, and EC50. These elements are identically sized matrixes 
#' with the dimensions (treatment x hour) 
#' @examples 
#' sample_types = rep( c(rep("line1",3), rep("line2",3)), 5)
#' treatments = c(rep("DMSO",6), rep("drug",24))
#' concentrations = c( rep(0,6),rep(200,6), rep(500,6),rep(1000,6),rep(5000,6))
#' values=c(100,99,100,90,91,92,99,97,99,89,87,88,86,89,88,56,59,58,66,65,67,
#'          25,23,24,42,43,46,4,5,9)
#' hours = rep(48, length(values))
#' plate_id = "plate_1"
#' ds = create_dataset( sample_types, treatments, concentrations, 
#'                       hours, values, plate_id, negative_control = "DMSO")
#' ds = normalize_plates_by_vehicle(ds, summary_method = "mean")
#' library(drc)
#' fits = fit_statistics( ds, fct=LL.3() )
#' fit_statistics_matrixes( fits )
#' @export
fit_statistics_matrixes = function( fits ){
    expected_names = c("hour","treatment","sample_type","AUC","coef_slope",
               "coef_asymp_low","coef_asymp_high","coef_EC50","obs_min",
               "obs_max","ANOVA_F_test","ANOVA_P_value")
    if( !is.data.frame( fits ) ){
        stop("parameter fits must be a data frame")
    }
    for(i in 1:length(expected_names) ){
        if( !expected_names[i] %in% names(fits) ){
            stop(paste("fits data frame must have column",expected_names[i]))
        }
    }
    sample_types = sort(unique(fits$sample_type))
    MS <- list()
    for(i in 1:length(sample_types)){
        cur_s = sample_types[i]
        idx_cur_s = which( fits$sample_type==cur_s )
        M = list(
             ANOVA_P_value =  metric_to_grid( fits$hour[idx_cur_s], 
                                              fits$treatment[idx_cur_s],
                                              fits$ANOVA_P_value[idx_cur_s] ),
             AUC =  metric_to_grid( fits$hour[idx_cur_s], 
                                              fits$treatment[idx_cur_s],
                                              fits$AUC[idx_cur_s] ),
             EC50 =  metric_to_grid( fits$hour[idx_cur_s], 
                                    fits$treatment[idx_cur_s],
                                    fits$coef_EC50[idx_cur_s] )
        )
        MS[[ length(MS)+1 ]] = M
    }
    names(MS) = sample_types
    MS
}

#' Convert values in M1, M2 to fold-change M1 vs M2
#' @param M1 matrix for numerator of fold-change
#' @param M2 matrix for denominator of fold-change
#' @return matrix with fold-change for M1/M2. Values that are NA in either M1 
#' or M2 are NA in M
#' @examples 
#' M1 = matrix( c(1, 1, 4, 1), nrow=2, ncol=2 )
#' M2 = matrix( 1:4, nrow=2, ncol=2 )
#' convert_to_foldchange(M1, M2)
#' @export
convert_to_foldchange = function(M1, M2){
    if( !is.matrix( M1 ) | !is.matrix(M2) ){
        stop("M1 and M2 must be matrixes")   
    }
    if( dim(M1)[1] != dim(M2)[1] | dim(M1)[2] != dim(M2)[2] ){
        stop("M1 and M2 must have the same dimensions")   
    }
    M = log2( M2/M1 )
    is_na = is.na(M) | is.nan(M)
    is_lt0 = !is_na & M<0 
    M[is_lt0] = M[is_lt0]*-1
    M=2^M
    M[is_lt0] = M[is_lt0] * -1
    M
}


#' Calculate fold-change in Dose Response Curve AUC and EC50 at all time points 
#' for two samples 
#' 
#' This is a convenience function that reproduces analysis that can also be 
#' performed using the \code{\link{fit_statistics}} function. The rationale 
#' for using this method is it contains some logic for screening out results 
#' where:
#'  
#' \enumerate{
#'  \item{the P value for ANOVA testing whether observations are  
#' significantly different from each other is below a given threshold}
#'  \item{at least one sample type was observed with a normalized value 
#' below a given threshold (e.g. achieved a SF50)}
#'  \item{the log10(difference in AUC) exceeds a threshold}
#' }
#' 
#' @param D experiment dataset with columns matching the output of 
#' \code{combine_data_and_map}
#' @param sample_types vector of two sample types to compare
#' @param fct Non-linear function to fit, defaults to four variable LL.4 model 
#' that estimates slope, upper asymptote, lower asymptote, and EC50. To fix 
#' lower asymptote at 1, pass drc::LL.3(). To fix lower asympotote ant 1 and 
#' upper asymptote at 0, pass dcr::LL.2(). For a list of available functions, 
#' see drc::getMeanFunctions(). To pass a function you must load the drc 
#' package in your R session.
#' @param max_Pval mark comparisons where fit P_value exceeds this as NA
#' @param min_obs mark comparisons where neither sample type has an observation 
#' below this value to be NA
#' @param min_log10_EC50_diff mark comparisons where the log10( difference in 
#' AUC ) at that timepoint is less than this value to be NA
#' @return list of matrixes for fold-change in AUC, fold-change in EC50, 
#' P value for ANOVA, sample AUC, and sample EC50. Matrix X-axis is hours and 
#' Y-axis is treatments.
#' @examples 
#' # set up two timepoints; this function is more useful when there are many 
#' # drugs and multiple timepoints
#' sample_types = rep( c(rep("line1",3), rep("line2",3)), 5)
#' treatments = c(rep("DMSO",6), rep("drug",24))
#' concentrations = c( rep(0,6),rep(200,6), rep(500,6),rep(1000,6),rep(5000,6))
#' values=c(100,99,100,90,91,92,99,97,99,77,76,79,92,93,91,66,65,67,77,76,74,
#'          36,35,35,57,56,55,22,25,24,100,99,100,90,91,92,99,97,99,89,87,88,
#'          86,89,88,56,59,58,66,65,67,25,23,24,42,43,46,4,5,9)
#' hours = c( rep(24, length(values)/2), rep(48, length(values)/2))
#' treatments = rep(treatments, 2)
#' concentrations = rep(concentrations, 2)
#' sample_types = rep(sample_types, 2)
#' plate_id = "plate_1"
#' ds = create_dataset( sample_types, treatments, concentrations, 
#'                      hours, values, plate_id, negative_control = "DMSO")
#' ds = normalize_plates_by_vehicle(ds, summary_method = "mean")
#' # calculate grid
#' timecourse_relative_grid(ds, c("line1", "line2"), drc::LL.3() )
#' @export
timecourse_relative_grid = function( D, sample_types, fct, 
                             max_Pval=1, 
                             min_obs=Inf,
                             min_log10_EC50_diff=0 ){
    if(length(sample_types)!=2){
        stop( "Must pass exactly two sample_types" )
    }
    if( length(D$sample_type==sample_types[1])==0 ){
        stop(paste("sample_type[1]",sample_types[1],"not found in D"))
    }
    if( length(D$sample_type==sample_types[2])==0 ){
        stop(paste("sample_type[1]",sample_types[2],"not found in D"))
    }
    fits = fit_statistics(D[D$sample_type %in% sample_types,], fct = fct )
    odd = which( fits$sample_type==sample_types[1] )
    even = which( fits$sample_type==sample_types[2] )
    lbl_hours = fits$hour[odd]
    lbl_treat = fits$treatment[odd]
    fits_P = metric_to_grid( lbl_hours, lbl_treat, fits$ANOVA_P_value[odd] )
    
    
    fits_minobs_o = metric_to_grid( lbl_hours, lbl_treat, fits$obs_min[odd] )
    fits_minobs_e = metric_to_grid( lbl_hours, lbl_treat, fits$obs_min[even] )
    fits_EC50_o = metric_to_grid( lbl_hours, lbl_treat, fits$coef_EC50[odd] )
    fits_EC50_e = metric_to_grid( lbl_hours, lbl_treat, fits$coef_EC50[even] )
    
    HAS_OBS_LT50 = (!is.na( fits_minobs_o ) & fits_minobs_o < min_obs) | 
                   (!is.na( fits_minobs_e ) & fits_minobs_e < min_obs) 
    
    HAS_EC50_o = ( !is.na( fits_EC50_o ) & !is.infinite( fits_EC50_o )  )
    HAS_EC50_e = ( !is.na( fits_EC50_e ) & !is.infinite( fits_EC50_e )  )
    HAS_EC50 = HAS_EC50_e | HAS_EC50_o
    HAS_EC50_DIFF = ( is.infinite( fits_EC50_o ) & HAS_EC50_e ) | 
        ( is.infinite( fits_EC50_e ) & HAS_EC50_o ) |
        ( HAS_EC50_o & HAS_EC50_e & abs(log10(fits_EC50_o)-log10(fits_EC50_e)) 
          > min_log10_EC50_diff )
    HAS_PVAL = !is.na( fits_P ) & fits_P <= max_Pval
    
    fits_AUC_even = metric_to_grid( lbl_hours, lbl_treat, fits$AUC[even] )
    fits_AUC_odd = metric_to_grid( lbl_hours, lbl_treat, fits$AUC[odd] )
    fits_AUC_even[  !HAS_PVAL | !HAS_OBS_LT50 | !HAS_EC50_DIFF ] = NA
    fits_AUC_odd[  !HAS_PVAL | !HAS_OBS_LT50 | !HAS_EC50_DIFF ] = NA
    fits_AUC_ratio=convert_to_foldchange(fits_AUC_even, fits_AUC_odd)
    
    fits_EC50_odd = metric_to_grid( lbl_hours, lbl_treat, fits$coef_EC50[odd] )
    fits_EC50_even = metric_to_grid( lbl_hours, lbl_treat, fits$coef_EC50[even] )
    fits_EC50_even[  !HAS_PVAL | !HAS_OBS_LT50 | !HAS_EC50_DIFF ] = NA
    fits_EC50_odd[  !HAS_PVAL | !HAS_OBS_LT50 | !HAS_EC50_DIFF ] = NA
    fits_EC50_ratio=convert_to_foldchange( fits_EC50_even, fits_EC50_odd )
    list( AUC_foldchange=fits_AUC_ratio, 
          EC50_foldchange=fits_EC50_ratio,
          P_value = fits_P,
          AUC_1 = metric_to_grid( lbl_hours, lbl_treat, fits$AUC[odd]),
          AUC_2 = metric_to_grid( lbl_hours, lbl_treat, fits$AUC[even]),
          EC50_1 = metric_to_grid( lbl_hours, lbl_treat, fits$coef_EC50[odd]),
          EC50_2 = metric_to_grid( lbl_hours, lbl_treat, fits$coef_EC50[even]))
}


#' Calculate AUC for sample/treatment at a concentration across multiple 
#' timepoints
#'
#' This function calculates the AUC at a single concentration 
#' across all timepoints. This is very distinct from a typical dose response 
#' curve, which can be used to calculate the AUC at a single timepoint across 
#' multiple concentrations. For that purpose, use \code{\link{fit_DRC}}. 
#' 
#' AUC for a single sample type/treatment condition across the 
#' whole timecourse is calculated using raw data normalized by the vehicle 
#' control specified by the user when the dataset was created. AUC is 
#' calculated by connecting the normalized data points and then using the 
#' trapezoids method for AUC implemented in pracma::trapz.
#'
#' For each treatment, in each concentration,
#' \itemize{
#'  \item{ AUC_v= AUC( raw vehicle matched plate_id )} 
#'  \item{ AUC_sample1 = AUC( raw sample1 after treatment ) / AUC_v }
#'  \item{ AUC_sample2 = AUC( raw sample2 after treatment ) / AUC_v }   
#'  \item{ AUC ratio = AUC_sample1 / AUC_sample2}
#' }
#' Return value is a list. The matrix in the list contains one column for each 
#' AUC value followed by one column for each AUC ratios of sample_types[1] / 
#' sample_types[2]. The concentrations vector in the list indicates the 
#' treatment concentration in each column. The sample_types vector in the list 
#' indicates the sample type if the column is a single sample AUC, or the word 
#' "ratio" if the value in the column is an AUC ratio.
#' 
#' @param D experiment dataset with columns matching the output of 
#' \code{read_incucyte_exported_to_excel}
#' @param treatments which treatments to calculate
#' @param sample_types samples on which to calculate 
#' @param concentrations which concentrations to draw
#' @param summary_method mean or median
#' @return list( AUC matrix, concentrations, sample_types )
#' @examples 
#' # Create a dataset with two drugs normalized against DMSO, tested at a single
#' # concentration. Timepoints are 0, 12, 24, 36, 48 hours. drug1 affected 
#' # line1 cells at this concentration, while drug2 affected line2 cells.
#' ds = create_dataset( 
#'      sample_types= rep( c("line1", "line2"), 15),
#'      treatments = c( rep("DMSO", 10), rep("drug1", 10), rep("drug2", 10)),
#'      concentrations = c( rep(0, 10), rep(200, 20) ),
#'      hours = rep( c(0,0,12,12,24,24,36,36,48,48), 3),
#'      values = c(10,10,15,15,30,25,60,50,80,65,
#'                 10,10,12,15,22,28,26,48,30,60,
#'                 10,10,15,12,30,11,60,13,80,9), 
#'      plate_id = "plate_1",
#'      negative_control = "DMSO")
#' ds=normalize_plates_by_vehicle(ds, summary_method="mean")
#' tc=timecourse_AUC_ratio(ds, 
#'                       treatments=c("drug1", "drug2"), 
#'                       sample_types=c("line1", "line2"), 
#'                       concentrations=c(0, 200), 
#'                       summary_method="mean")
#' # create a data frame from individual list components
#' data.frame( sample_type = tc$sample_types,
#'             concentration = tc$concentrations,
#'             t(tc$AUC), 
#'             row.names=1:length(tc$sample_types),
#'             stringsAsFactors=FALSE)
#' @export
timecourse_AUC_ratio = function( D, treatments, sample_types, concentrations,
                           summary_method){
    
    if( summary_method != "mean" & summary_method != "median"){
        stop("summary_method parameter must be either mean or median")
    }
    if( length(sample_types) != 2 ){
        stop("Must pass exactly two sample types in parameter sample_types")
    }
    if( sum(D$sample_type==sample_types[1]) == 0  ){
        stop(paste("sample_types[1], '",sample_types[1], 
                   "', not found in D",sep="" ) )
    }
    if( sum(D$sample_type==sample_types[2]) == 0  ){
        stop(paste("sample_types[2], '",sample_types[2], 
                   "', not found in D",sep="" ) )
    }
    for(i in 1:length(treatments)){
        if(sum(D$sample_type==sample_types[1] & D$treatment==treatments[i])==0){
            stop(paste("treatment",treatments[i],"not found for sample type",
                       sample_types[1]))
        }
        if(sum(D$sample_type==sample_types[2] & D$treatment==treatments[i])==0){
            stop(paste("treatment",treatments[i],"not found for sample type",
                       sample_types[2]))
        }
    }
    for(i in 1:length(concentrations)){
        if(sum(D$sample_type==sample_types[1] & 
               D$concentration==concentrations[i])==0){
            stop(paste("concentration",concentrations[i],
                       "not found for sample type", sample_types[1]))
        }
        if(sum(D$sample_type==sample_types[2] & 
               D$concentration==concentrations[i])==0){
            stop(paste("concentration",concentrations[i],
                       "not found for sample type",sample_types[2]))
        }
    }
    Mauc = matrix(0, nrow=length(treatments), 
                  ncol=length(concentrations)*length(sample_types))
    col_concentrations = rep(0, dim(Mauc)[2] +  length(concentrations) )
    col_sample_types = rep(0, dim(Mauc)[2] + length(concentrations) )
    col_identities = rep(0, dim(Mauc)[2] + length(concentrations) )
    for(idx_t in 1:length(treatments) ){
        cur_t = treatments[idx_t]
        plate_id = unique(D$plate_id[D$treatment==cur_t])
        # AUC for vector in each line at the appropriate plate
        AUC_v = c()
        for( idx_s in 1:length(sample_types)){
            cur_s = sample_types[idx_s]
            neg_ctl = unique( D$negative_control[ 
               D$plate==plate_id & D$sample_type==cur_s & D$treatment == cur_t])
            neg_ctl_conc = unique(D$concentration[ D$plate==plate_id &
                                            D$sample_type==cur_s &
                                            D$treatment == neg_ctl &
                                            D$is_negative_control ])
            AUC_v = c(AUC_v,
                      AUC_across_timepoints(D, cur_s, neg_ctl, 
                                            neg_ctl_conc, summary_method) )
        }
        ctr=1
        for(idx_c in 1:length(concentrations)){
            for(idx_s in 1:length(sample_types)){
                cur_s = sample_types[idx_s]
                cur_c = concentrations[idx_c]
                col_concentrations[ctr] = cur_c
                col_identities[ctr] = paste("AUC_", cur_s, collapse="", sep="")
                col_sample_types[ctr] = cur_s
                if( sum( D$sample_type==cur_s & D$concentration==cur_c & 
                         D$treatment==cur_t)==0 ){
                    ratio=NA
                }
                else{
                    ratio = AUC_across_timepoints( D, cur_s, cur_t, cur_c,
                                               summary_method) / AUC_v[idx_s]
                }
                Mauc[idx_t, ctr] = ratio
                ctr = ctr+1
            }
        }
    }
    Mauc_ratio = matrix(0, nrow=length(treatments), 
                        ncol = length(concentrations) )
    for( idx_c in 1:length(concentrations)){
        cur_c = concentrations[idx_c]
        idx_conc_s1 = which(col_concentrations==cur_c &
                            col_sample_types == sample_types[1])
        idx_conc_s2 = which(col_concentrations==cur_c &
                            col_sample_types == sample_types[2])
        Mauc_ratio[,idx_c] = Mauc[,idx_conc_s1] / Mauc[,idx_conc_s2 ]
        col_concentrations[ctr] = cur_c
        col_sample_types[ctr] = "ratio"
        col_identities[ctr] = paste("AUCratio_", sample_types[1],":",
                                    sample_types[2], collapse="", sep="")
        ctr = ctr + 1
    }
    Mauc = cbind(Mauc, Mauc_ratio)
    dimnames(Mauc)[[2]] = col_identities
    dimnames(Mauc)[[1]] = treatments
    
    list( AUC = Mauc, 
          concentrations=col_concentrations,
          sample_types=col_sample_types)
}
