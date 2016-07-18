# Utility functions that are not specialized to DRC calculations

#-------------------------------------------------------------------------------
# local hardcoded variables
#-------------------------------------------------------------------------------

# Hat tip to http://colorbrewer2.org/ for color choices
INC_colors = list()
INC_colors[["colors1"]] = c("#fcae91") 
INC_colors[["colors2"]] = c("#fcae91", "#cb181d") 
INC_colors[["colors3"]] = c("#fcae91", "#fb6a4a", "#cb181d")
INC_colors[["colors4"]] = c("#fcae91", "#fb6a4a", "#de2d26", "#a50f15")
INC_colors[["colors5"]] = c("#fcbba1", "#fc9272", "#fb6a4a", "#de2d26", 
                            "#a50f15")
INC_colors[["colors6"]] = c("#fcbba1", "#fc9272", "#fb6a4a", "#ef3b2c", 
                            "#cb181d", "#99000d")
INC_colors[["colors7"]] = c("#fee0d2", "#fcbba1", "#fc9272", "#fb6a4a", 
                            "#ef3b2c", "#cb181d", "#99000d")
INC_colors[["colors8"]] = c("#fee0d2", "#fcbba1", "#fc9272", "#fb6a4a", 
                            "#ef3b2c", "#cb181d", "#99000d", "black")

INC_colors_DRC = c("black", "red", "cornflowerblue", "gold", "darkgreen", 
                   "orange", "pink", "gray", "springgreen","indianred1",
                   "yellow")

INC_pches = c(19, 22, 23, 24, 1, 4, 5)

# to label plates
abc = c('A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P')

# Calculate standard error for error bars
#
se = function(x){
    stats::sd(x, na.rm=TRUE)/sqrt(sum(!is.na(x))) 
}


logtics = function( log10_min, log10_max, show_x_log_tics, show_x_exponent ){
    if( log10_min < 0 ){
        stop("log10_min must be at least 0")
    }
    if( log10_min >= log10_max ){
        stop("log10_min must be less than log10_max")   
    }
    if( show_x_log_tics ){
        tics = c( 10^log10_min )
        if( show_x_exponent ){
            labels = c( 10^log10_min )
        }else{
            labels = c( log10_min )
        }
        lwd_tics = c(2) # allows for thicker ticks at whole exponent numbers
        for(n in (log10_min+1):log10_max){
            tics = c(tics, seq( from=10^(n-1)+10^(n-1), to=10^n, by=10^(n-1) ) )
            if( show_x_exponent ){
                labels = c(labels, rep("", 8), 
                           parse( text=paste("10^",n,sep="") ) )
            }else{
                labels = c(labels, rep("", 8), parse( text=n ) )
            }
            lwd_tics = c(lwd_tics, rep(1,8), 2)
        }
        
    }
    else{
        tics = 10^(log10_min:log10_max)
        if( show_x_exponent ){
            labels = parse( text=paste("10^",log10_min:log10_max,sep="") ) 
        }else{
            labels = log10_min:log10_max
        }
        lwd_tics = rep(2, length(tics))
    }
    list(values=tics, labels=labels, lwd=lwd_tics)
}

# given a vector of real number V and a vector of colors color_map, scale 
# values in V to the appropriate value in color_map. If color_bounds is not 
# passed, bounds are set to c( min(V), max(V)). If color_NA is passed, it is 
# used for indexes which( is.na(V) )
# 
color_scale = function( V, color_map, color_bounds=NA, color_NA=NA ){
    if( is.na( color_bounds[1] ) ){
        Vmax = max(V, na.rm=TRUE)
        Vmin = min(V, na.rm=TRUE)
    }else{
        Vmin = color_bounds[1]
        Vmax = color_bounds[2]
    }
    increment = (Vmax-Vmin) / (length(color_map)-1)
    lookup = seq(from=Vmin, to=Vmax, by=increment )
    if( is.vector(V) ){
        out = rep("", length(V))
        for(i in 1:length(V)){
            out[i] = color_map[ which( lookup>V[i] )[1] ]
        }
    }else if( is.matrix(V) ){
        out = matrix("", nrow=dim(V)[1], ncol=dim(V)[2], 
                     dimnames=list( dimnames(V)[[1]], dimnames(V)[[2]]))
        for(rr in 1:dim(V)[1]){
            for(cc in 1:dim(V)[2]){
                col=color_map[ which( lookup>=V[rr,cc] )[1]]
                out[rr,cc] = col 
            }
        }
    }
    if( !is.na( color_NA) )
        out[is.na(out)] = color_NA
    out
}


# given vector V, of strings, split elements using string and col/first/last
get.split.col = function(v, string, col=0, last=FALSE, first=FALSE){
    if( last & first )
        stop("Cannot request both last and first column")
    if( col==0 & !last & !first)
        stop("Must request either a column by index, first, or last")
    
    for(i in 1:length(v)){
        x = strsplit( v[i], string, fixed=TRUE)[[1]]
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


#-------------------------------------------------------------------------------
# Simple hash functions, implemented using an environment
#-------------------------------------------------------------------------------

# Create 
#
hsh_new = function(){
    new.env(hash=TRUE, parent=emptyenv()) 
}

# Does key exist in H
#
hsh_in = function(H, key){
    exists(key, H)
}

# retrieve value from hash
#
hsh_get = function( H, key, na.if.not.found=FALSE ){
    if( length(key)==1 ){
        if( na.if.not.found ){
            if( exists(key, H) )
                get(key, H)
            else
                NA
        }else{
            get(key, H)
        }
    }
    else{
        results = rep(0, length(key) )
        if( !na.if.not.found ){
            for(i in 1:length(key) ){
                if( exists(key[i], H) ){
                    results[i] = get(key[i], H )
                }
                else{
                    results[i] = NA
                }
            }
        }else{
            for(i in 1:length(key) ){
                results[i] = get(key[i], H )
            }
        }
        results
    }
}

# Set value in a hash
#
hsh_set = function( H, key, value ){
    assign(key, value, envir=H)
}

# Create a hash from vectors v1, v2 with keys from v1 and values from v2
# if v2 is null, set it to 1:length(v1)
#
hsh_from_vectors = function( v1, v2=NULL ){
    if( is.null(v2) )
        v2 = 1:length(v1)
    if( length(v1) != length(v2) ){
        stop("Length of v1 != length of v2")
    }
    H = hsh_new()
    for( i in 1:length(v1) ){
        hsh_set(H, v1[i], v2[i] )
    }
    H
}
