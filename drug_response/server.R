library(shiny)
library(drc)
library(HTDoseResponseCurve)

colors = c("black", "red", "cornflowerblue", "gold", "darkgreen", "orange", "pink", "gray")

get.split.col = function(v, string, col=0, last=F, first=F){
    if( last & first )
        stop("Cannot request both last and first column")
    if( col==0 & !last & !first)
        stop("Must request either a column by index, first, or last")
    
    for(i in 1:length(v)){
        x = strsplit( v[i], string, fixed=T)[[1]]
        if(last){
            v[i] = x[length(x)]
        }
        else if(first){
            v[i] = x[1]
        }
        else{
            v[i] = x[col]
        }
    }
    v
}


load.matrix=function( fn ){
    read.table(file=fn, sep="\t",row.names=1, check.names=F, header=T, 
               na.strings=c('NA', '-99999'), stringsAsFactors=F)
}
se=function(x){
    sd(x, na.rm=T)/sqrt(sum(!is.na(x))) 
}

make_fit = function(dd, input){
    treatments = rep("drug", dim(dd)[1] * dim(dd)[2])
    concentrations = rep( as.numeric(rownames(dd)), dim(dd)[2])
    values = as.numeric( data.matrix( dd ) )
    plate_id="plate"
    negative_control=0
    st = c()
    for(i in 1:dim(dd)[2]){
        st = c(st, rep(names(dd)[i], dim(dd)[1] ) )    
    }
    D=create_dataset( sample_types=st, 
                      treatments=treatments, 
                      concentrations=concentrations, 
                      values=values, 
                      plate_id=plate_id, 
                      negative_control=negative_control)
    unique_samples = unique(names(dd))
    
    if( input$curve_fit ==1 ){
        fct=drc::LL.4()   
    }else if( input$curve_fit==2){
        fct=drc::LL.3()   
    }else{
        fct=drc::LL.2()   
    }
    fit_DRC(D, sample_types = unique_samples, treatments = "drug", 
            fct=fct )
}

shinyServer(function(input, output, session) {
    
    output$drugresponsePlot <- renderPlot({
        if( !is.null( input$file$datapath ) ){
            dd = load.matrix(input$file$datapath[1])
            fitted = make_fit(dd, input)
            par(mar=c(5,5,3,1))
            if( as.numeric(input$axis_pointsize) > 1.5 ){
                par(mar=c(8,8,3,1))
                par(mgp=c(5,1,0))
            }
            plot( fitted,
                  cex.axis=as.numeric(input$axis_pointsize),
                xlab=input$xlab,
                ylab=input$ylab,
                show_x_log_tics=as.logical(input$show_x_logtic),
                show_x_exponent=as.logical(input$show_x_exponent),
                show_EC50 = as.logical(input$sf50),
                bar_multiple=as.numeric(input$barmultiple),
                xlim=c(as.numeric(input$min_x), as.numeric(input$max_x) ),
                cex.lab=as.numeric(input$axis_labelsize), 
                ylim=c(0, as.numeric(input$max_y) ),
                show_legend=FALSE)
            
            legend( as.numeric(input$legend_x), 
                    as.numeric(input$legend_y), 
                    fitted$sample_types, 
                    pch=15, cex=1.75, 
                    col=colors[1:length(fitted$sample_types) ],
                    bty="n")
            session$sendCustomMessage("download_PDF", list())
        }
    })
    
    output$table <- renderDataTable({
        
        if( !is.null( input$file$datapath) ){
            dd = load.matrix( input$file$datapath[1] )
            fitted = make_fit(dd, input)
            sn=get.split.col( rownames(fitted$fit_stats), "_|_", first=TRUE)
            data = data.frame( sample=sn,
                               AUC=fitted$fit_stats$AUC,
                               EC50=fitted$fit_stats$coef_EC50,
                               obs_min=fitted$fit_stats$obs_min,
                               obs_max=fitted$fit_stats$obs_max,
                               slope = fitted$fit_stats$coef_slope,
                          stringsAsFactors=FALSE)
        }
    })

    output$download_PDF <- downloadHandler(
        filename = function() {
            "dose_response.pdf"
        },
        
        content = function(file) {
            pdf( file, width=10, height=8 )
            dd = load.matrix( input$file$datapath[1] )
            
            par(mar=c(5,5,3,1))
            if( as.numeric(input$axis_pointsize) > 1.5 ){
                par(mar=c(8,8,3,1))
                par(mgp=c(5,1,0))
            }
            fitted = make_fit(dd, input)
            plot( fitted,
                  xlab=input$xlab,
                  ylab=input$ylab,
                  show_x_log_tics=as.logical(input$show_x_logtic),
                  show_x_exponent=as.logical(input$show_x_exponent),
                  cex.axis=as.numeric(input$axis_pointsize),
                  bar_multiple=as.numeric(input$barmultiple), 
                  show_EC50 = as.logical(input$sf50),
                  xlim=c(as.numeric(input$min_x), as.numeric(input$max_x) ),
                  cex.lab=as.numeric(input$axis_labelsize), 
                  ylim=c(0, as.numeric(input$max_y) ),
                  show_legend = FALSE)
            legend( as.numeric(input$legend_x), 
                    as.numeric(input$legend_y), 
                    fitted$sample_types, 
                    pch=15, cex=1.75, 
                    col=colors[1:length(fitted$sample_types) ],
                    bty="n")
            
            dev.off()
        }
    )
})
