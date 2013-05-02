# Two annotation files are used by CARMEN, one for mice and one for human data
# These combine symbol, official name, a reasonable location (as many genes have 
# multiple locations for distinct isoforms) and a list of GO accessions associated
# with the symbol.
#  
# NCBI Gene is used for symbol lists and gene names
# UCSC refFlat is used for symbol locations
# Ensembl/MGI are used for GO identifiers

wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz
wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
gunzip Mus_musculus.gene_info
gunzip Homo_sapiens.gene_info
cut -f3,12 Mus_musculus.gene_info | sort > symbol2name_MM
cut -f3,12 Homo_sapiens.gene_info | sort > symbol2name_HS
python
f = open("symbol2name_MM")
fo = open("symbol2name_MM_clean", "w")
for line in f:
    fo.write(line.replace("'", "prime"))

fo.close()
f = open("symbol2name_HS")
fo = open("symbol2name_HS_clean", "w")
for line in f:
    fo.write(line.replace("'", "prime"))

fo.close()
quit()
mv symbol2name_HS_clean symbol2name_HS
mv symbol2name_MM_clean symbol2name_MM

wget http://viewvc.geneontology.org/viewvc/GO-SVN/trunk/gene-associations/gene_association.goa_human.gz?rev=HEAD
wget http://viewvc.geneontology.org/viewvc/GO-SVN/trunk/gene-associations/gene_association.mgi.gz?rev=HEAD
gunzip gene_association.goa_human.gz
gunzip gene_association.mgi.gz

cut -f3,5 gene_association.goa_human | sort > symbol2GO_HS
cut -f3,5 gene_association.mgi | sort > symbol2GO_MM

wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/refFlat.txt.gz
gunzip refFlat.txt
mv refFlat.txt refFlat_MM.txt
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz
gunzip refFlat.txt
mv refFlat.txt refFlat_HS.txt
cut -f1,3,5 refFlat_MM.txt | sort > symbol2loc_MM
cut -f1,3,5 refFlat_HS.txt | sort > symbol2loc_HS

# JOIN THE FILES IN R
R
source('http://www.davidquigley.com/software/quantitative_genetics.R')
s2go.MM = read.table('symbol2GO_MM', sep='\t', stringsAsFactors=F, comment.char="!")
names(s2go.MM) = c("symbol", "GO")
s2go.HS = read.table('symbol2GO_HS', sep='\t', stringsAsFactors=F, comment.char="!")
names(s2go.HS) = c("symbol", "GO")
s2loc.MM = read.table('symbol2loc_MM', sep='\t', stringsAsFactors=F)
s2loc.HS = read.table('symbol2loc_HS', sep='\t', stringsAsFactors=F)
names(s2loc.MM) = c("symbol", "chrom", "loc.start")
names(s2loc.HS) = c("symbol", "chrom", "loc.start")

s2n.MM = read.table('symbol2name_MM', sep='\t', stringsAsFactors=F)
s2n.HS = read.table('symbol2name_HS', sep='\t', stringsAsFactors=F)
names(s2n.MM) = c("symbol", "name")
names(s2n.HS) = c("symbol", "name")
s2n.MM$symbol[is.na(s2n.MM$symbol)]="NA"
s2n.HS$symbol[is.na(s2n.HS$symbol)]="NA"

# There are numerous symbols that are not annotated on the reference assembly
symbols.MM = sort(unique( s2n.MM$symbol ) )
symbols.HS = sort(unique( s2n.HS$symbol ) )
chrom.MM = rep("NA", length(symbols.MM))
loc.start.MM = rep("NA", length(symbols.MM))
fullname.MM = symbols.MM
go.MM = rep("NA", length(symbols.MM))
chrom.HS = rep("NA", length(symbols.HS))
loc.start.HS = rep("NA", length(symbols.HS))
fullname.HS = symbols.HS
go.HS = rep("NA", length(symbols.HS))

m = match.idx(symbols.MM, s2loc.MM$symbol)
chrom.MM[ m$idx.A ] = s2loc.MM$chrom[m$idx.B]
loc.start.MM[ m$idx.A ] = s2loc.MM$loc.start[m$idx.B]
m = match.idx(symbols.MM, s2n.MM$symbol)
fullname.MM[ m$idx.A ] = s2n.MM$name[m$idx.B]


m = match.idx(symbols.HS, s2loc.HS$symbol)
chrom.HS[ m$idx.A ] = s2loc.HS$chrom[m$idx.B]
loc.start.HS[ m$idx.A ] = s2loc.HS$loc.start[m$idx.B]
m = match.idx(symbols.HS, s2n.HS$symbol)
fullname.HS[ m$idx.A ] = s2n.HS$name[m$idx.B]

for(i in 1:length(symbols.MM)){
    go = paste(s2go.MM$GO[s2go.MM$symbol==symbols.MM[i]], collapse=",")
    if( go != "" ){
        go.MM[i] = go
    }
}

for(i in 1:length(symbols.HS)){
    go = paste(s2go.HS$GO[s2go.HS$symbol==symbols.HS[i]], collapse=",")
    if( go != "" ){
        go.HS[i] = go
    }
}

mm = data.frame(symbol=symbols.MM, chr=get.split.col(chrom.MM,"chr",last=T), ucsc_chr=chrom.MM, loc_start=loc.start.MM, fullname=fullname.MM)
hs = data.frame(symbol=symbols.HS, chr=get.split.col(chrom.HS,"chr",last=T), ucsc_chr=chrom.HS, loc_start=loc.start.HS, fullname=fullname.HS)
write.matrix(mm, 'mouse_annotation.txt')
write.matrix(hs, 'human_annotation.txt')
quit()

# BACK IN THE SHELL. CLEAN UP
rm Mus_musculus.gene_info
rm Homo_sapiens.gene_info
rm gene_association.goa_human
rm gene_association.mgi
rm refFlat_MM.txt
rm refFlat_HS.txt
