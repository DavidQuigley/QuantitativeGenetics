#cut -f 3,5 /notebook/annotations/2015_gene_association.goa_human > /notebook/annotations/gene_ontology/2015_gene_association

fn_ontology = '/notebook/annotations/2015_go-basic.obo'
fn_assn = '/notebook/annotations/gene_ontology/2015_gene_association'
fn_onto_compressed = '/notebook/annotations/gene_ontology/2015_gene_association_lookup'
fn_output_R = '/notebook/annotations/gene_ontology/load_GO.R'

f = open(fn_ontology)
fo = open(fn_onto_compressed, 'w')
fo.write( "GO\tname\n")
for line in f:
    if line[0:3]=="id:":
        cur_id=line.rstrip('\r\n').split(": ")[1]
    
    elif line[0:5]=="name:":
        cur_name=line.rstrip('\r\n').split(": ")[1].replace("'", "prime")
        fo.write( cur_id  + '\t' + cur_name + '\n')
    
    elif line[0:5]=="[Type":
        break

fo.close()

f = open(fn_assn)
go2s = {}
s2go = {}
for line in f:
    if line[0] == "!":
        continue
    a = line.rstrip('\r\n').split('\t')
    try:
        go2s[a[1]].append(a[0])
    
    except KeyError:
        go2s[a[1]] = [ a[0] ]
    
    try:
        s2go[a[0]].append(a[1])
    
    except KeyError:
        s2go[a[0]] = [ a[1] ]

f.close()

fo = open(fn_output_R, 'w')
fo.write("go2name = read.table('" + fn_onto_compressed + "', sep='\\t', header=TRUE, stringsAsFactors=FALSE)\n")
fo.write("go2name = hsh_from_vectors( go2name$GO, go2name$name)\n")

fo.write( 'go2s = hsh_new()\n')
for go in go2s:
    fo.write( "hsh_set( go2s, '" + go + "', c('" + "','".join( list(set(go2s[go])) ) + "') )\n" )

fo.write( 's2go = hsh_new()\n')
for s in s2go:
    fo.write( "hsh_set( s2go, '" + s + "', c('" + "','".join( list(set(s2go[s])) ) + "') )\n" )

fo.close()

