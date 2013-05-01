fn_GO_h = '/notebook/annotations/gene_association.goa_human'
fn_UCSC_genes_h = '/notebook/annotations/UCSC_all_known_human_genes_Feb_2009_simplified.txt'
fn_out_h = '/notebook/code/release/human_annotation.txt'
fn_GO_m = '/notebook/annotations/gene_association.mgi'
fn_UCSC_genes_m = '/notebook/annotations/UCSC_all_known_mouse_genes_July_2007.txt'
fn_out_m = '/notebook/code/release/mouse_annotation.txt'
fn_GO_obo = '/notebook/annotations/gene_ontology_ext_may_2011.obo'
fn_GO_obo_compressed = '/notebook/code/release/gene_ontology.1_2.obo.compressed'

def generate_species_file(fn_GO, fn_UCSC_genes, fn_out):
    g2GO = {}
    f = open(fn_GO)
    for line in f:
        if line[0] != '!':
            a = line.rstrip('\r\n').split('\t')
            if g2GO.has_key( a[2] ):
                g2GO[a[2]].append( a[4] )
            else:
                g2GO[a[2]] = [ a[4] ]
    
    f = open(fn_UCSC_genes)
    h = f.readline().rstrip('\r\n').split('\t')
    i_gn = h.index('Gene Name')
    i_chr = h.index('Chr')
    i_uchr = h.index('UCSC_Chr')
    i_ls = h.index('cds_start')
    i_fullname = h.index('description')
    i_rs = h.index('refseq')
    g2annot = {}
    for line in f:
        a = line.rstrip('\r\n').split('\t')
        blurb = a[i_chr] + '\t' + a[i_uchr] + '\t' + a[i_ls] + '\t' + a[i_fullname]
        gene = a[i_gn]
        if g2annot.has_key( gene ):
            if a[i_rs] != 'NA':
                g2annot[ a[i_gn] ] = blurb # prefer annotation with refseq, as gene name is more accurate
        else:
            g2annot[ a[i_gn] ] = blurb             
    f = open(fn_out, 'w')
    f.write('symbol\tchr\tucsc_chr\tloc_start\tfullname\tGO\n')
    for gene in sorted(g2annot):
        if gene in g2GO:
            go = ','.join(g2GO[gene])
        else:
            go = 'NA'
        f.write( gene + '\t' + g2annot[gene] + '\t' + go + '\n')
    f.close()


def generate_species_file_from_simplified(fn_GO, fn_UCSC_genes, fn_out):
    # PREFERRED METHOD
    g2GO = {}
    f = open(fn_GO)
    for line in f:
        if line[0] != '!':
            a = line.rstrip('\r\n').split('\t')
            if g2GO.has_key( a[2] ):
                g2GO[a[2]].append( a[4] )
            else:
                g2GO[a[2]] = [ a[4] ]

    f = open(fn_UCSC_genes)
    h = f.readline().rstrip('\r\n').split('\t')
    # symbol	chr	loc_start	loc_end	refseq	description
    i_gn = h.index('symbol')
    i_chr = h.index('chr')
    i_uchr = i_chr
    i_ls = h.index('loc_start')
    i_fullname = h.index('description')
    i_rs = h.index('refseq')
    g2annot = {}
    for line in f:
        a = line.rstrip('\r\n').split('\t')
        blurb = a[i_chr] + '\t' + a[i_uchr] + '\t' + a[i_ls] + '\t' + a[i_fullname]
        gene = a[i_gn]
        if g2annot.has_key( gene ):
            if a[i_rs] != 'NA':
                g2annot[ a[i_gn] ] = blurb # prefer annotation with refseq, as gene name is more accurate
        else:
            g2annot[ a[i_gn] ] = blurb             
    f = open(fn_out, 'w')
    f.write('symbol\tchr\tucsc_chr\tloc_start\tfullname\tGO\n')
    for gene in sorted(g2annot):
        if gene in g2GO:
            go = ','.join(g2GO[gene])
        else:
            go = 'NA'
        f.write( gene + '\t' + g2annot[gene] + '\t' + go + '\n')
    f.close()

def generate_GO_file(fn_obo, fn_out):
    f = open( fn_obo )
    fo = open(fn_out, 'w')
    fo.write('ID\tcategory\tname\n')
    for line in f:
        if line[0:3] == "id:":
            goid = line.rstrip('\r\n').split(' ')[1]
        elif line[0:5] == "name:":
            name = line.rstrip('\r\n')[6:]
        elif line[0:10] == "namespace:":
            ns = line.rstrip('\r\n').split(' ')[1]
            if ns=="biological_process":
                ns = "BP"
            elif ns=="cellular_component":
                ns = "CC"
            elif ns=="molecular_function":
                ns = "MF"
            fo.write( goid + '\t' + ns + '\t' + name + '\n')
        elif line[0:7] == "alt_id:":
            goid = line.rstrip('\r\n').split(' ')[1]
            fo.write( goid + '\t' + ns + '\t' + name + '\n')
generate_species_file_from_simplified(fn_GO_h, fn_UCSC_genes_h, fn_out_h)
#generate_species_file(fn_GO_m, fn_UCSC_genes_m, fn_out_m)
#generate_GO_file(fn_GO_obo, fn_GO_obo_compressed)
