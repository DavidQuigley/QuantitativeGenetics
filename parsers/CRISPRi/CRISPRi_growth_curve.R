cc = read.table('/datasets/human_lines_CAPAN1/CRISPRi/2015_11_10_CAPAN_cell_count_worksheet.txt', header=TRUE, stringsAsFactors=FALSE)


staggered_slopes = function( T, counts, seeds ){
    # intuition here is that if data were full we would want the ratio 
    # this measurement's count / last measurement's seed, but the data contain gaps and 
    # we want to use the last data point where there is a measurement
    slopes = c()
    intervals = c()

    i = 2
    while( i < length(T) ){
        if(!is.na( counts[i] )){
            j=i-1
            while( is.na( seeds[j] ) ){
                j=j-1
            }
            slopes = c(slopes, counts[i]/seeds[j])
            intervals = c(intervals, T[i] - T[j] )
        }
        i=i+1
    }
    list(slopes=round(slopes,3), intervals=intervals)
}
res_H1v =staggered_slopes( cc$T, cc$H1_veh_cnt, cc$H1_veh_seed )
res_H1o = staggered_slopes( cc$T, cc$H1_olap_cnt, cc$H1_olap_seed )
res_H2v = staggered_slopes( cc$T, cc$H2_veh_cnt, cc$H2_veh_seed )
res_H2o = staggered_slopes( cc$T, cc$H2_olap_cnt, cc$H2_olap_seed )

# Plot on actual time scale

layout(matrix(1:2,1,2))
plot(0,0, col="white", xlim=c(0,35), ylim=c(0,4), main="H1 count multiples", 
     xlab="days", ylab="count / previous seed", axes=FALSE)
axis(1, seq(from=0, to=35, by=5))
axis(2, 0:4, las=1)

xs = cumsum(res_H1v$intervals)
rect(0,0,xs[1],res_H1v$slopes[1], col="#33333333")
for(i in 2:length(xs)){
    rect(xs[i-1],0,xs[i],res_H1v$slopes[i], col="#33333333")
}

xs = cumsum(res_H1o$intervals)
rect(0,0,xs[1],res_H1o$slopes[1], col="#0000ff33")
for(i in 2:length(xs)){
    rect(xs[i-1],0,xs[i],res_H1o$slopes[i], col="#0000ff33")
}

plot(0,0, col="white", xlim=c(0,35), ylim=c(0,4), main="H2 count multiples", 
     ylab="count / previous seed", xlab="days", axes=FALSE)
axis(1, seq(from=0, to=35, by=5))
axis(2, 0:4, las=1)

xs = cumsum(res_H2v$intervals)
rect(0,0,xs[1],res_H2v$slopes[1], col="#33333333")
for(i in 2:length(xs)){
    rect(xs[i-1],0,xs[i],res_H2v$slopes[i], col="#33333333")
}

xs = cumsum(res_H2o$intervals)
rect(0,0,xs[1],res_H2o$slopes[1], col="#0000ff33")
for(i in 2:length(xs)){
    rect(xs[i-1],0,xs[i],res_H2o$slopes[i], col="#0000ff33")
}

layout(matrix(1:2,1,2))
barplot( res_H1v$slopes / res_H1v$intervals, col="#33333333", main="Growth per day, H1", las=1)
barplot( res_H1o$slopes / res_H1o$intervals, col="#0000ff33", add=TRUE, las=1)

barplot( res_H2v$slopes / res_H2v$intervals, col="#33333333", main="Growth per day, H2", las=1)
barplot( res_H2o$slopes / res_H2o$intervals, col="#0000ff33", add=TRUE, las=1)