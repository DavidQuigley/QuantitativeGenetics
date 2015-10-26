library(shiny)
library(drc)

colors = c("black", "red", "cornflowerblue", "gold", "darkgreen")

fit_model = function( dd, concentrations, conditions ){
    conc = rep(concentrations, dim(dd)[2] )
    response=as.numeric(matrix(data.matrix(dd), nrow=dim(dd)[1] * dim(dd)[2], ncol=1))
    dd.long = data.frame( response=response, conc=conc  )
    curveids = c()
    for(i in 1:length(conditions)){
        curveids = c( curveids, rep( conditions[i], dim(dd)[1] ) )
    }
    curveids = factor(curveids)
    drm( dd.long, curveid=curveids, fct = LL2.2() )
}

extract_sf = function( model, y_val ){
    y_val = 0.5
    b_fitted = coef(m1)[ 1 : (length(coef(m1))/2) ]
    e_fitted = coef(m1)[ ((length(coef(m1))/2)+1) : length(coef(m1)) ]
    y_fitted = rep(0, length(b_fitted))
    for( i in 1:length(b_fitted)){
        y_fitted[i] = ((1/y_val)-1)/exp(b_fitted[i]) + exp(e_fitted[i])
    }
    y_fitted
}

plot.dose.response = function( dd, concentrations, conditions, xlab, ylab, show50=TRUE ){
    
    m1 = fit_model( dd, concentrations, conditions )
    plot(m1, type="none", lty=c(1,1), pch=c(15,15),  xlim=c(0,1e5), ylim=c(0,1), 
         axes=FALSE, legend=FALSE, cex.lab=1.25, col=c("black", "red"),
         xlab=xlab, ylab=ylab, lwd=2 )
    axis(2, at=c(0, 0.5, 1), labels=c("0", "", "1"), las=1)
    axis(1, at=c(1e1, 1e2, 1e3, 1e4, 1e5), labels=c(1,2,3,4,5), las=1)
    unique.conditions = unique(conditions)
    for(i in 1:length(unique.conditions) ){
        points( concentrations, rowMeans(dd[,conditions==unique.conditions[i]]), pch=15, col=colors[i], cex=1.2)
    }
    for(i in 1:dim(dd)[1]){
        for(j in 1:length(unique.conditions)){
            v = as.numeric( dd[i,conditions==unique.conditions[j]] )
            x = concentrations[i]
            y = mean(v)
            bar_width = x/10
            lines(c( x, x ), c( y-(3*se(v)), y+(3*se(v)) ), lwd=1, col=colors[j] )
            lines(c( x-bar_width, x+bar_width ), c( y-(3*se(v)), y-(3*se(v)) ), lwd=1, col=colors[j] )
            lines(c( x-bar_width, x+bar_width ), c( y+(3*se(v)), y+(3*se(v)) ), lwd=1, col=colors[j] )        
        }
    }
    if( show50 ){
        y_fitted = extract_sf( m1, 0.5 )
        for(i in 1:length(y_fitted)){
            lines( c( y_fitted[i], y_fitted[i] ), c(0.5, 0), col=colors[i], lwd=2, lty=2 )   
        }
    }
}

load.matrix=function( fn ){
    read.table(file=fn, sep="\t",row.names=1, check.names=F, header=T, na.strings=c('NA', '-99999'), stringsAsFactors=F)
}
se=function(x){
    sd(x, na.rm=T)/sqrt(sum(!is.na(x))) 
}

shinyServer(function(input, output) {
    
    output$drugresponsePlot <- renderPlot({
        if( !is.null( input$file$datapath) ){
            dd = load.matrix(input$file$datapath[1])
            concentrations =  as.numeric(rownames(dd))
            conditions = factor(names(dd))
            legends = sort(unique(names(dd)))
            plot.dose.response( dd, concentrations, conditions, input$xlab, input$ylab, input$sf50 )
            legend( 10, 0.1, legends, pch=15, cex=1.25, col=colors[1:length(legends)], bty="n")
        }
    })
    
    output$table <- renderDataTable({
        
        if( !is.null( input$file$datapath) ){
            dd = load.matrix(input$file$datapath[1])
            concentrations =  as.numeric(rownames(dd))
            conditions = factor(names(dd))
            legends = sort(unique(names(dd)))
            y_fitted = extract_sf( fit_model( dd, concentrations, conditions), 0.5 )
            data = data.frame( legends, SF50=y_fitted, stringsAsFactors=FALSE )
        }
    })
})
