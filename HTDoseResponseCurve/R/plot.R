#' convert a matrix to colors and plot a heatmap grid
#' 
#' If color_bounds
#' @param M matrix of values to plot
#' @param block.height element height, default 20
#' @param block.width element width, default 10
#' @param space.X space between elements on horizontal axis, default 3
#' @param space.Y space between elements on vertical axis, default 10
#' @param cex.x cex value for x axis
#' @param cex.y cex value for y axis
#' @param border boolean, draw a border around grid elements
#' @param color_palatte vector of colors c(cold,middle,hot), defaults to 
#' c("blue", "white", "red") 
#' @param color_bounds two element vector with lower and upper limit for color 
#' palatte, defaults to c( min(M), max(M) )
#' @param num_colors number of color gradations in palatte, default 50
#' @return M, including possible truncation from color bounds
#' @examples 
#' M = matrix(1:12, nrow=3, ncol=4,byrow=TRUE)
#' plot_color_grid(M)
#' dimnames(M)[[1]] = c("row_a","row_b","row_c")
#' dimnames(M)[[2]] = c("col1", "col2", "col3","col4")
#' plot_color_grid(M, color_palatte=c("blue","white","yellow" ) )
#' @export
plot_color_grid=function(M, block.height=20, block.width=10, space.X=3, 
                         space.Y=10, cex.x=1, cex.y=1, border=TRUE, 
                         color_palatte=c("blue","white","red"), 
                         color_bounds=NA, num_colors=50 ) {
    
    cmap = colorRampPalette( colors=color_palatte )( num_colors )
    if( is.na(color_bounds)[1] ){
        color_bounds = c( min(M, na.rm=TRUE), max(M, na.rm=TRUE) )
    }
    if( sum(M>color_bounds[2], na.rm=TRUE)>0 ){
        warning("value(s) in M exceed maximum color bound, truncating value(s)")
    }
    if( sum(M<color_bounds[1], na.rm=TRUE)>0 ){
        warning("value(s) in M exceed minimum color bound, truncating value(s)")
    }
    M[M>color_bounds[2]] = color_bounds[2]
    M[M<color_bounds[1]] = color_bounds[1]
    MC = color_scale( M, cmap, color_bounds = color_bounds )
    
    n.rows = dim(MC)[1]
    n.cols = dim(MC)[2]
    total.width =  ( block.width*n.cols) + ( (n.cols-1) * space.X)
    total.height = ( block.height * n.rows ) + ( (n.rows-1) * space.Y)
    plot(0,0,col="white", xlim=c(0,total.width), ylim=c(0,total.height), 
         axes=FALSE, xlab="", ylab="", bg="azure2")
    xlab_locs = rep(0, n.cols)
    ylab_locs = rep(0, n.rows)
    cur.y = total.height
    for(rr in 1:n.rows){
        for(cc in 1:n.cols){  
            this.x.left = (cc-1)*block.width + (cc-1)*space.X
            this.x.right = this.x.left+block.width
            xlab_locs[cc] = this.x.right - (block.width)
            if( border )
                rect( this.x.left , cur.y - block.height, this.x.right, cur.y, 
                      col=MC[rr,cc], border="azure2")
            else
                rect( this.x.left , cur.y - block.height, this.x.right, cur.y, 
                      col=MC[rr,cc], border=NA)
        }
        ylab_locs[rr] = cur.y - (0.5*block.height)
        cur.y = cur.y - block.height - space.Y
    }
    axis(1, at=xlab_locs, labels=dimnames(MC)[[2]], las=2, cex.axis=cex.x, 
         tick=FALSE, padj=1, line=-1.5 )
    axis(2, at=ylab_locs, labels=dimnames(MC)[[1]], las=2, cex.axis=cex.y, 
         tick=FALSE, hadj=1, line=-1.5)
    M
}

#' Plot a dose response curve
#' 
#' If the curve cannot be fit (meaning that the optimization method used by 
#' the drm() method from the drc package failes to converge) then the summarized 
#' points will be plotted without a curve connecting them.
#' 
#' @param x HTfit object
#' @param ... standard parameters for \code{plot} function
#' @param log10_xmin x minimum for plot on log10 scale, default 0
#' @param log10_xmax x maximum for plot on log10 scale, default 5
#' @param bar_multiple multiplier for standard error bars, default 2
#' @param summary_method summary method for points to plot in timecourse, one 
#' of ("mean", "median"), defaults "mean"
#' @param show_EC50 show dotted vertical lines at EC50 values if they are fit,
#' default TRUE
#' @param show_legend show a legend, default TRUE
#' @param show_x_log_tics if axes==TRUE and this parameter is TRUE, draw tics 
#' between each power of 10 scaled appropriately. If axes==TRUE and this 
#' parameter is FALSE, omit ticks between powers of 10. Default TRUE.
#' @param show_x_exponent If axes==TRUE and this parameter is TRUE, show labels 
#' on X axis as 10^2, 10^3, ... . If axes==TRUE and this parameter is FALSE, 
#' show labels on X axis as an exponent (2, 3, ... ). Default TRUE.
#' @return none
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
#' fit_1=fit_DRC(ds, sample_types=c("line1", "line2"), treatments=c("drug"), 
#'         hour = 48, fct=drc::LL.3() )
#' plot(fit_1)
#' plot(fit_1, show_EC50=FALSE, show_legend=FALSE, lwd=3,col=c("black", "gold"),
#'      log10_xmin=0, log10_xmax=4, xlab="concentration nM", 
#'      ylab="surviving fraction")
#' legend(1.5, 0.3, c("Line 1", "Line 2"), col=c("black", "gold"), pch=15)
#' @export
plot.HT_fit = function( x, ..., log10_xmin=1, log10_xmax=5, bar_multiple=2, 
                        summary_method="mean", show_EC50=TRUE, 
                        show_legend=TRUE, show_x_log_tics=TRUE, 
                        show_x_exponent=TRUE){
    if( summary_method!="mean" & summary_method !="median"){
        stop("parameter summary_method must be one of {mean, median}")   
    }
    
    N = length(x$unique_conditions)
    plot_parameters = list(...)
    if( length( which( names(plot_parameters)== "pch"  ) )==0 )
        plot_parameters[["pch"]] = rep(15, N)
    if( length( which( names(plot_parameters)== "col"  ) )==0 ){
        plot_parameters[["col"]] = INC_colors_DRC[1:N]
    }else{
        if(length(plot_parameters[["col"]])==1)
            plot_parameters[["col"]] = rep(plot_parameters[["col"]], N)
    }
    if( length( which( names(plot_parameters)== "lty"  ) )==0 )
        plot_parameters[["lty"]] = rep(1, N)
    if( length( which( names(plot_parameters)== "ylim"  ) )==0 )
        plot_parameters[["ylim"]] = c(0, 1.2)
    if( length( which( names(plot_parameters)== "xlab"  ) )==0 )
        plot_parameters[["xlab"]] = ""
    if( length( which( names(plot_parameters)== "ylab"  ) )==0 )
        plot_parameters[["ylab"]] = ""
    if( length( which( names(plot_parameters)== "cex"  ) )==0 )
        plot_parameters[["cex"]] = 1
    if( length( which( names(plot_parameters) == "cex.axis"))==0 )
        plot_parameters[["cex.axis"]] = 2
    if( length( which( names(plot_parameters) == "cex.lab"))==0 )
        plot_parameters[["cex.lab"]] = 1.5
    plot_parameters[["xlim"]] = c(10^log10_xmin, 10^log10_xmax)
    plot_parameters[["yaxs"]] = "i"
    plot_parameters[["xaxs"]] = "i"
    uc = as.character(unique(x$input$conditions_to_fit))
    
    cond2color = hsh_from_vectors( uc, plot_parameters[["col"]][1:length(uc)] )
    Mstat = plyr::ddply( x$input, c("conditions_to_fit", "concentration"), 
                         function(po){ data.frame( 
                             mu=mean(po$value, na.rm=TRUE), 
                             med=median(po$value, na.rm=TRUE),
                             sterr = se(po$value), 
                             bar_width=(po$concentration/10), 
                             stringsAsFactors=FALSE ) }
    )
    if( summary_method=="mean"){
        Mstat$value = Mstat$mu
    }
    else{
        Mstat$value = Mstat$med
    }
    
    # if user wants axes (either by passing axes=TRUE or not specifying), we're 
    # going to draw our own axes, so set the plot parameter axes to FALSE. If 
    # the user passed axes=FALSE and therefore DOESN'T want axes, respect that.
    if( length( which( names(plot_parameters)== "axes"  ) )==0 ){
        plot_parameters[["axes"]] = FALSE
        draw_axes = TRUE
    }else{
        draw_axes = plot_parameters[["axes"]]
        if( plot_parameters[["axes"]] ){
            plot_parameters[["axes"]] = FALSE
        }
    }
    
    if( x$is_fitted ){
        plot_parameters[["type"]] = "none"
        plot_parameters[["legend"]] = FALSE # for plot.drc()
        plot_parameters[["bp"]] = 1 # for plot.drc()
        
        plot_parameters = c( list(x$model), plot_parameters)
        do.call( plot, plot_parameters )
        box(col="white", lwd=4)
        if( show_EC50 ){
            EC50 = x$fit_stats$coef_EC50 
            for( i in 1:length(EC50) ){
                if( !is.na(EC50[i]) & !is.infinite(EC50[i]) ){
                    lines( c( EC50[i], EC50[i] ), c(0, 0.5), 
                           col=hsh_get( cond2color,
                                        rownames(x$fit_stats)[i] ),
                           lwd=plot_parameters[["lwd"]], lty=3 )
                }
            }
        }
        bandwidth_factor=10
        x_legend=10
    }else{
        Mstat$concentration=log10(Mstat$concentration)
        Mstat$concentration[is.infinite( Mstat$concentration)] = 0
        plot_parameters[["x"]] = -1
        do.call( plot, plot_parameters)
        Mstat$bar_width = 0.2
        bandwidth_factor=50
        x_legend=0.2
    }
    if( draw_axes ){
        axis( 2, at=c(0, 0.5, 1), labels=c("0", "0.5", "1"), las=1, 
              cex.axis=plot_parameters[["cex.axis"]],
              font=plot_parameters[["font"]])
        tics = logtics( log10_xmin, log10_xmax, show_x_log_tics,show_x_exponent)
        axis(1, at=tics$values, labels=tics$labels, las=1, 
             font=plot_parameters[["font"]],
             cex.axis=plot_parameters[["cex.axis"]], lwd.ticks=1)   
        axis(1, at=tics$values[tics$lwd==2],
             labels=rep("", sum(tics$lwd==2)), las=1, 
             font=plot_parameters[["font"]],
             cex.axis=plot_parameters[["cex.axis"]], lwd.ticks=2)    
    }
        
    Mstat$concentration[Mstat$concentration==0] = 1
    points( Mstat$concentration, Mstat$value, pch = 19, 
            col=hsh_get( cond2color, as.character(Mstat$conditions_to_fit)), 
            cex=plot_parameters[["cex"]])
    if(show_legend){
        legend( x_legend, 0.45, x$unique_conditions, pch=19,
                col=hsh_get( cond2color, uc), bty="n", cex=0.75 )    
    }
    
    for(i in 1:dim(Mstat)[1]){
        px = Mstat$concentration[i]
        py = Mstat$value[i]
        pse = Mstat$sterr[i] * bar_multiple
        pcol = hsh_get( cond2color, as.character(Mstat$conditions_to_fit[i]) )
        bwdth = px/bandwidth_factor # varies if fit or not
        lines(c( px, px ), c( py-pse, py+pse ), lwd=1, col=pcol )
        lines(c( px-bwdth, px+bwdth ), c( py-pse, py-pse ), lwd=1, col=pcol )
        lines(c( px-bwdth, px+bwdth ), c( py+pse, py+pse ), lwd=1, col=pcol )        
    }
}

#' Plot an image of the raw intensities for a plate
#' 
#' @param plate data frame in format matching that produced by 
#' \code{read_incucyte_from_excel}.
#' @param hour hour of the experiment to plot
#' @param color_bounds values at which to draw coldest and hottest colors, 
#' defaults to c(0, 100)
#' @param color_palatte colors to use for cold-middle-hot, defaults to 
#' c("white", "#fdbb84","#e34a33")
#' @param main title to draw above plot, defaults to ""
#' @param cex point size, defaults to 1.5
#' @return none
#' @examples 
#' # small toy example
#' plate_1 = create_empty_plate( 6, hour=0, plate_id="plate_1")
#' plate_1[1,1:3] = c(100, 40, 10 )
#' plate_1[2,1:3] = c(90, 70, 30 )
#' plot_values_by_plate( plate_1, hour=0 )
#' # real-world example
#' pkg = "HTDoseResponseCurve"
#' fn_data = system.file("extdata", "sample_data_384.xlsx", package = pkg)
#' plate_data = read_plates_from_Incucyte_export( fn_data, "p1", 
#'                                                number_of_wells=384)
#' plot_values_by_plate(plate_data, hour=96)
#' @export
plot_values_by_plate = function( plate, hour, color_bounds=c(0,100), 
                                 color_palatte=c("white", "#fdbb84","#e34a33"),
                                 main="", cex=1.5 ){
    if( sum( plate$hours==hour )==0 ){
        stop(paste("hour parameter", hour, "not present in plate"))
    }
    # plot raw value for a plate at the specified hour
    n_row = sum(plate$hours==hour)
    idx_max = dim(plate)[2] - 2
    n_col = idx_max
    plot( -1, -1, xlim=c(0, n_col+1), ylim=c(0, n_row+1), 
          axes=FALSE, xlab="", ylab="",yaxs="i", xaxs="i", main=main )
    abc = abc[1:n_row]
    abc_rev = abc[ length(abc):1]
    axis(1, 1:n_col, cex=0.5, las=2)
    axis(2, at=1:n_row, labels=abc_rev, las=1, cex=0.5)
    cmap = colorRampPalette( colors=color_palatte )(color_bounds[2]-
                                                        color_bounds[1])
    vals = ceiling( data.matrix(plate[plate$hours==hour,1:idx_max] ) )
    colors = color_scale( vals, cmap, color_bounds = color_bounds )
    pches = matrix(19, ncol=n_col, nrow=n_row)
    pches[ is.na(colors) ] = 1
    colors[is.na(colors)] = "grey"
    xs = matrix( rep( 1:n_col, n_row ), ncol=n_col, nrow=n_row, byrow=TRUE)
    ys = matrix( rep( 1:n_row, n_col ), ncol=n_col, nrow=n_row)
    ys = n_row-ys+1
    points( as.numeric(xs), as.numeric(ys), col=colors, 
            pch=as.numeric(pches), cex=cex )
    box()
}

#' Plot an time course of the raw intensities 
#' 
#' @param D experiment dataset with columns matching the output of 
#' \code{combine_data_and_maps}
#' @param sample_types which sample types to draw
#' @param treatments which treatments to draw
#' @param concentrations which concentrations to draw
#' @param ... standard parameters for \code{plot} function
#' @param cex.yaxis size coefficient for Y axis, default 2
#' @param axis.font font value for axis, default 2 (bold)
#' @param summary_method summary method for points to plot in timecourse, one 
#' of (mean, median), defaults to mean
#' @return data frame reporting points plotted, sample types, and concentrations
#' @examples 
#' pkg = "HTDoseResponseCurve"
#' fn_map = system.file("extdata", "sample_data_384_platemap.txt",package=pkg)
#' fn_data = system.file("extdata", "sample_data_384.xlsx", package = pkg)
#' plate_map = read_platemap_from_Incucyte_XML( fn_map )
#' plate_data = read_plates_from_Incucyte_export( fn_data, "p1", 
#'                                                number_of_wells=384)
#' plate_data$hours = round(plate_data$hours)
#' ds = combine_data_and_map( plate_data, plate_map, negative_control = "DMSO" )
#' ds = normalize_plates_by_vehicle( ds, summary_method="mean")
#' ds = ds[ds$treatment=="drug13" | ds$treatment=="DMSO",]
#' plot_timecourse_raw( ds, sample_types=c("line_1", "line_2"), 
#'                      treatments="DMSO", concentrations=0)
#' @export
plot_timecourse_raw = function( D, sample_types, treatments, 
                                concentrations, ..., cex.yaxis=2, axis.font=2,
                                summary_method="mean"){
    if( summary_method != "mean" & summary_method != "median"){
        stop("summary_method parameter must be either mean or median")
    }
    for(i in 1:length(sample_types)){
        if( sum(D$sample_type==sample_types[i])==0 ){
            stop( paste("sample type", sample_types[i],
                        "in parameter sample_types not found in D") )
        }
    }
    for(i in 1:length(treatments)){
        if( sum(D$treatment==treatments[i])==0 ){
            stop( paste("treatment", treatments[i],
                        "in parameter treatments not found in D") )
        }
    }
    for(i in 1:length(concentrations)){
        if( sum(D$concentration==concentrations[i])==0 ){
            stop( paste("concentration",concentrations[i],
                        "in parameter concentrations not found in D") )
        }
    }
    ds_cur = D[ which(D$sample_type %in% sample_types &
                          D$treatment %in% treatments & 
                          D$concentration %in% concentrations),]
    if(dim(ds_cur)[1]==0){
        stop(paste("No matches for this combination of sample_types,",
                   "treatments, and concentrations"))
    }
    hours = sort( unique( ds_cur$hours ) )
    unique_concs = sort( unique(concentrations) )
    unique_lines = sort( unique(sample_types) )
    unique_treatments = sort( unique (treatments ) )
    conc2color = hsh_from_vectors( as.character(unique_concs) , 
                                   INC_colors_DRC[ 1:length(unique_concs) ] )
    
    line2pch = hsh_from_vectors( as.character(unique_lines) , 
                                 INC_pches[1:length(unique_lines)] )
    Mstat = plyr::ddply( ds_cur, c("sample_type", "hours", "concentration"), 
                         function(x){ 
        data.frame(
            mu=mean(x$value, na.rm=TRUE), med=median(x$value, na.rm=TRUE), 
            color=hsh_get( conc2color, as.character(x$concentration) ), 
            pch = hsh_get( line2pch, x$sample_type),
            sample_type=x$sample_type,
            concentration=x$concentration, stringsAsFactors=FALSE)}
    )
    if( summary_method=="mean"){
        values = Mstat$mu
    }else{
        values = Mstat$med
    }
    plot_parameters = list(...)
    plot_parameters[["xlim"]] = c(0, max(hours))
    plot_parameters[["axes"]] = FALSE
    plot_parameters[["pch"]] = Mstat$pch
    plot_parameters[["col"]] = Mstat$color
    plot_parameters[["x"]] = Mstat$hours
    plot_parameters[["y"]] = values
    if( length( which( names(plot_parameters)== "ylim"  ) )==0 ){
        plot_parameters[["ylim"]] = c(0, 100)
    }
    if( length( which( names(plot_parameters)== "ylab"  ) )==0 ){
        plot_parameters[["ylab"]] = ""
    }
    if( length( which( names(plot_parameters)== "xlab"  ) )==0 ){
        plot_parameters[["xlab"]] = ""
    }
    ylim = plot_parameters[["ylim"]]
    do.call( plot, plot_parameters )
    axis(2, c(ylim[1], (ylim[2])/2, ylim[2]), las=1, cex.axis=cex.yaxis,
         font= axis.font )
    axis(1, round(hours,1), las=2, cex.axis=1, font=axis.font )
    box()
    unique(Mstat)
}



#' Plot a dose response curve summary values across one or more time points
#' 
#' Used to generate a line and dot plot of one or more AUC or EC50 values 
#' across a series of time points. All of the values in the fit_stats data 
#' frame will be plotted. Observations that do not meet the criteria specified 
#' in the alpha parameter are plotted using open values (e.g. pch=1). 
#' Observations where the minimum EC50 exceeds the obs_min parameter can be 
#' plotted using a user-specified value for the pch parameter; by default this 
#' is pch=13.
#' 
#' @param fit_stats data frame returned by call to \code{fit_statistics}
#' @param statistic statistic to plot, must be one of AUC, EC50
#' @param ... standard parameters for \code{plot} function
#' @param alpha cut-off to distinguish results with an ANOVA_P_value that is 
#' considered statistically significant, plotted as closed circles. Defaults to 
#' 1. 
#' @param obs_min cut-off to distinguish results where the minimum observed 
#' response across any dose is no greater than obs_min. Useful for restricting 
#' plots to those values that acheived a result such as a Surviving Fraction 
#' below 50 percent. Defaults to NA, which plots everything.
#' @return data frames with plotted points, colors. Useful for drawing a legend.
#' @examples 
#' pkg = "HTDoseResponseCurve"
#' fn_map = system.file("extdata","sample_data_384_platemap.txt",package=pkg)
#' fn_data = system.file("extdata", "sample_data_384.xlsx", package = pkg)
#' plate_map = read_platemap_from_Incucyte_XML( fn_map )
#' plate_data = read_plates_from_Incucyte_export( fn_data, "p1", 
#'                                                number_of_wells=384)
#' plate_data$hours = round(plate_data$hours)
#' ds = combine_data_and_map( plate_data, plate_map, negative_control = "DMSO" )
#' ds = normalize_plates_by_vehicle( ds, summary_method="mean")
#' ds = ds[ds$treatment=="drug13",]
#' fits = fit_statistics(ds, fct = drc::LL.3() )
#' res=plot_fit_statistic( fits, "AUC", ylim=c(0, 6) ) 
#' @export
plot_fit_statistic = function( fit_stats, statistic, ..., alpha = 1, 
                               obs_min=NA){
    
    if( !( statistic == "AUC" | statistic=="EC50" ) ){
        stop("parameter statistic must be one of AUC, EC50" )
    }
    
    if( is.na(obs_min) ){
        obs_min = max(fit_stats$obs_min, na.rm=TRUE) + 1   
    }
    plot_par = list(...)
    if( length( which( names(plot_par)== "ylim"  ) )==0 ){
        plot_par[["ylim"]] = c(0, 10000)
    }
    if( length( which( names(plot_par)== "col"  ) )==0 ){
        colors = INC_colors_DRC    
    }else{
        colors = plot_par[["col"]]
    }
    if( length( which( names(plot_par)== "xlab"  ) )==0 ){
        plot_par[["xlab"]] = ""    
    }
    if( length( which( names(plot_par)== "cex"  ) )==0 ){
        plot_par[["cex"]] = 2
    }
    if( length( which( names(plot_par)== "ylab"  ) )==0 ){
        plot_par[["ylab"]] = ""    
    }
    plot_par[["x"]] = -100
    plot_par[["y"]] = -100
    if( length( which( names(plot_par)== "xlim"  ) )==0 ){
        plot_par[["xlim"]] = c(0, max(fit_stats$hour))
    }
    fit_stats = fit_stats[fit_stats$hour <= plot_par[["xlim"]][2],]
    plot_par[["axes"]] = FALSE
    ylim = plot_par[["ylim"]]
    
    if( length( which( names(plot_par)== "xlim"  ) )==0 ){
        plot_par[["xlim"]] = c(0, max(fit_stats$hour, na.rm=TRUE) )
    }
    
    do.call( plot, plot_par )
    ylims = 1:ylim[2]
    line_y = seq(from=ylim[1], to=ylim[2], by=(ylim[2]-ylim[1])/5 )
    for(i in 1:length(line_y)){
        abline( line_y[i], 0, col="lightgray")   
    }
    
    if( length( which( names(plot_par)== "las"  ) )==0 ){
        plot_par[["las"]] = 1
    }
    
    axis(2, line_y, las=plot_par[["las"]])
    axis(1, sort(unique(fit_stats$hour)), las=plot_par[["las"]])
    ctr=1
    sample_types = sort(unique(fit_stats$sample_type))
    treatments = sort(unique(fit_stats$treatment))
    for( i in 1:length( sample_types ) ){
        cur_samp = sample_types[i]
        for( j in 1:length(treatments) ){
            cur_treat = treatments[j]
            fits_cur = fit_stats$sample_type==cur_samp & 
                fit_stats$treatment==cur_treat
            if( statistic=="AUC" ){
                y=fit_stats$AUC[ fits_cur ]
            }else if( statistic=="EC50" ){
                y=fit_stats$coef_EC50[ fits_cur ]
            }else{
                stop(paste("cannot plot statistic",statistic))
            }
            xs = unique(fit_stats$hour[ fits_cur ] )
            pvals = fit_stats$ANOVA_P_value[ fits_cur ]
            min_observations = fit_stats$obs_min[ fits_cur ]
            if( length( y ) >1 ){
                for(k in 2:length(y)){
                    if( !is.na( pvals[k-1] ) & !is.na( pvals[k] ) & 
                                pvals[k-1] <= alpha & pvals[k] <= alpha ){
                        lines( c(xs[k], xs[k-1] ), 
                               c( y[k],  y[k-1] ), col=colors[ctr] )
                    }
                }
            }
            
            pches = rep(19, length(pvals))
            pches[ pvals > alpha ] = 13    
            pches[ pvals <= alpha & min_observations > obs_min ] = 18 
            
            points( xs, y, col=colors[ctr] , pch=pches, 
                    cex=plot_par[["cex"]])
            for(i in 1:length(y)){
                if( is.infinite( y[i] ) ){
                    text( xs[i], ylim[2], "Inf", cex=1, font=2, col=colors[ctr])
                }
            }
            cur_res = data.frame( sample_type=rep(cur_samp, length(y)),
                                  treatment=rep(cur_treat, length(y)),
                                  AUC=y, hour=xs, 
                                  color=rep(colors[ctr], length(y) ), 
                                  stringsAsFactors=FALSE)
            if( ctr==1 ){
                res = cur_res    
            }
            else{
                res = rbind(res, cur_res)   
            }
            ctr = ctr+1
        }
    }
    res
}



#' Standard boxplot with labeled outliers
#' 
#' Any point outside of the whiskers will be labeled. Passes all standard 
#' \code{boxplot()} parameters through to the R \code{boxplot()} function.
#' 
#' @param M one-dimensional matrix of values with row names for labels
#' @param ... standard commands for \code{boxplot()} function
#' @return standard output of \code{boxplot()} function
#' @examples
#' M = matrix( c(3,4,7,5,6,10), nrow=6, ncol=1 )
#' dimnames(M)[[1]] = paste("drug",1:6)
#' boxplot_label_outliers(M[,1])
#' @export
boxplot_label_outliers = function( M, ... ){
    plot_parameters = list(...)
    if( length( which( names(plot_parameters)== "ylim"  ) )==0 )
        plot_parameters[["ylim"]] = c(0, ceiling(max(M) ) )
    if( length( which( names(plot_parameters)== "xlim" ) )== 0 )
        plot_parameters[["xlim"]] = c(0.25, 1.25)
    y_max = plot_parameters[["ylim"]][2]
    plot_parameters[["x"]] = M 
    b=do.call( boxplot, plot_parameters )
    if( length(b$out)>0 ){
        outs = b$out
        outs = sort(outs, decreasing=TRUE)
        xs = rep(1, length(M))
        not_out = which( !( round(M,3) %in% round(outs,3)) )
        xs[not_out] = jitter( xs[not_out], 2 )
        points( xs, M, pch=19, cex=0.5 )
        for( i in 1:length(b$out)){
            text( 0.7, (y_max*0.75) - (i*0.1), names(outs)[i], adj=1 )
            lines( c(0.7,1), c((y_max*0.75) - (i*0.1), outs[i]), col="#00000033" ) 
        }
    }
    b
}