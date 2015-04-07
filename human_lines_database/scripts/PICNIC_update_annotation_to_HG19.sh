python

id2pos = {}

f = open('/notebook/annotations/manufacturers/affymetrix/SNP6_ProbePositions_na34.txt')
for line in f:
    if line[0] != "#":
        a = line.rstrip('\r\n').split('\t')
        id2pos[a[0]] = a[2]

f.close()
f = open('/notebook/annotations/manufacturers/affymetrix/SNP6_CopyNumberPositions_na34.txt')
for line in f:
    if line[0] != "#":
        a = line.rstrip('\r\n').split('\t')
        id2pos[a[0]] = a[2]

f.close()

f = open('/datasets/human_lines_CCLE/CN/PICNIC/annotation/old_ProbeRef.csv')
fo = open('/datasets/human_lines_CCLE/CN/PICNIC/annotation/ProbeRef.csv', 'w')
fo.write( f.readline() )
for line in f:
    a = line.rstrip('\r\n').split(',')
    if a[0] in id2pos:
        if id2pos[ a[0] ] !=  "---":
            a[3] = id2pos[a[0]]
    
    fo.write( ','.join( a ) + '\n' )

fo.close()
f.close()

f = open('/datasets/human_lines_CCLE/CN/PICNIC/annotation/old_Snp6FeatureMappings.csv')
fo = open('/datasets/human_lines_CCLE/CN/PICNIC/annotation/Snp6FeatureMappings.csv', 'w')
for line in f:
    a = line.rstrip('\r\n').split(',')
    if a[0] in id2pos:
        if id2pos[ a[0] ] !=  "---":
            a[2] = id2pos[a[0]]
    
    fo.write( ','.join( a ) + '\n' )

fo.close()
f.close()



f = open('/datasets/human_lines_CCLE/CN/PICNIC/annotation/old_probe_id.csv')
fo = open('/datasets/human_lines_CCLE/CN/PICNIC/annotation/probe_id.csv', 'w')
fo.write( f.readline() )
probe_locs = []
for line in f:
    a = line.rstrip('\r\n').split(',')
    if a[0] in id2pos:
        if id2pos[a[0]] !=  "---":    
            a[3] = id2pos[a[0]]
    
    probe_locs.append( a[3] )
    fo.write( ','.join( a ) + '\n' )

fo.close()
f.close()

quit()

#
# Reorder so everything is in ascending order
# Rewrite SNP_pos.csv
#
R 
ref = read.table('/datasets/human_lines_CCLE/CN/PICNIC/annotation/ProbeRef.csv', 
                 sep=',', header=TRUE, stringsAsFactors=FALSE)
probe_id = read.table('/datasets/human_lines_CCLE/CN/PICNIC/annotation/probe_id.csv',  
                 sep=',', header=TRUE, stringsAsFactors=FALSE)
sum(ref$probe_set != probe_id$probe_set)

A_coeffs = read.table('/datasets/human_lines_CCLE/CN/PICNIC/annotation/old_A_coeffs.csv',  
                      sep=',', stringsAsFactors=FALSE)
B_coeffs = read.table('/datasets/human_lines_CCLE/CN/PICNIC/annotation/old_B_coeffs.csv',  
                      sep=',', stringsAsFactors=FALSE)
T_coeffs = read.table('/datasets/human_lines_CCLE/CN/PICNIC/annotation/old_T_coeffs.csv',  
                      sep=',', stringsAsFactors=FALSE)                      
no_features = read.table('/datasets/human_lines_CCLE/CN/PICNIC/annotation/old_no_features.csv',  
                      sep=',', stringsAsFactors=FALSE)
p_coeffs = read.table('/datasets/human_lines_CCLE/CN/PICNIC/annotation/old_p_coeffs.csv',  
                      sep=',', stringsAsFactors=FALSE)
mappings = read.table('/datasets/human_lines_CCLE/CN/PICNIC/annotation/Snp6FeatureMappings.csv',  
                      sep=',', stringsAsFactors=FALSE)

new_order = order(probe_id$chr, probe_id$pos)
ref = ref[new_order,]
probe_id = probe_id[new_order,]

A_coeffs=A_coeffs[new_order,]
B_coeffs=B_coeffs[new_order,]
T_coeffs=T_coeffs[new_order,]
no_features=no_features[new_order,]
p_coeffs=p_coeffs[new_order,]
mappings = mappings[new_order,]
write.table( probe_id, file="/datasets/human_lines_CCLE/CN/PICNIC/annotation/probe_id.csv", 
             quote=F, sep=",", row.names=FALSE, col.names=TRUE)
write.table( ref, file="/datasets/human_lines_CCLE/CN/PICNIC/annotation/ProbeRef.csv", 
             quote=F, sep=",", row.names=FALSE, col.names=TRUE)
write.table( probe_id$pos, file="/datasets/human_lines_CCLE/CN/PICNIC/annotation/SNP_pos.csv", 
             quote=F, row.names=FALSE, col.names=FALSE)
write.table( ref$ids, file="/datasets/human_lines_CCLE/CN/PICNIC/annotation/SNP_ind.csv", 
             quote=F, row.names=FALSE, col.names=FALSE)
write.table( mappings, file="/datasets/human_lines_CCLE/CN/PICNIC/annotation/Snp6FeatureMappings.csv", 
             quote=F, sep=',', row.names=FALSE, col.names=FALSE)             
# Update chr_info.csv
mins = rep(0, 24)
maxs = rep(0, 24)
for(i in 1:24){
    mins[i] = min( which(ref$chr==i) )
    maxs[i] = max( which(ref$chr==i) )
}
write.table( df, file="/datasets/human_lines_CCLE/CN/PICNIC/annotation/chr_info.csv", 
             quote=F, sep=",", row.names=FALSE, col.names=FALSE)

write.table( A_coeffs, file="/datasets/human_lines_CCLE/CN/PICNIC/annotation/A_coeffs.csv", 
             quote=F, sep=",", row.names=FALSE, col.names=FALSE)
write.table( B_coeffs, file="/datasets/human_lines_CCLE/CN/PICNIC/annotation/B_coeffs.csv", 
             quote=F, sep=",", row.names=FALSE, col.names=FALSE)
write.table( T_coeffs, file="/datasets/human_lines_CCLE/CN/PICNIC/annotation/T_coeffs.csv", 
             quote=F, sep=",", row.names=FALSE, col.names=FALSE)
write.table( no_features, file="/datasets/human_lines_CCLE/CN/PICNIC/annotation/no_features.csv", 
             quote=F, sep=",", row.names=FALSE, col.names=FALSE)
write.table( p_coeffs, file="/datasets/human_lines_CCLE/CN/PICNIC/annotation/p_coeffs.csv", 
             quote=F, sep=",", row.names=FALSE, col.names=FALSE)

quit()
