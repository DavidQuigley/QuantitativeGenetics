# Example call for Affymetrix Gene ST:
# python equalizer.py \
# -f gene \
# -p MoGene-1_1_noSprFVBSnps \
# -v /datasets/mouse_sequence_reference/SNPS/snps-all.SPR.vcf,/datasets/mouse_sequence_reference/SNPS/2012-0612-snps+indels_FVBNJ_annotated.vcf \
# -a data/annotation/probe_bed/MoGene_probes_MM9.bed \
# -b /notebook/code/BEDTools-Version-2.15.0/bin \
# -s na32 \
# -g mm9 \
# -d 37 \
# -c MoGene-1_1-st-v1 \
# -i 1.1 \
# -r "Mus musculus" \
# -w mouse \
# -u "David Quigley" \
# -l dquigley@cc.ucsf.edu \
# -y data/annotation/MoGene-1_1/MoGene-1_1-st-v1.na32.mm9.probeset.csv \
# -t data/annotation/MoGene-1_1/MoGene-1_1-st-v1.na32.mm9.transcript.csv \
# -q data/annotation/MoGene-1_1/MoGene-1_1-st-v1.r4.pgf \
# -m data/annotation/MoGene-1_1/MoGene-1_1-st-v1.r4.mps \
# -k data/annotation/MoGene-1_1/MoGene-1_1-st-v1.r4.clf \
# -o results/MoGene_FVB_SPR
 
# Example call for Affymetrix M430 v2
# python equalizer.py \
# -f IVT \
# -p M430_noSprFVBSnps \
# -v /datasets/mouse_sequence_reference/SNPS/snps-all.SPR.vcf,/datasets/mouse_sequence_reference/SNPS/2012-0612-snps+indels_FVBNJ_annotated.vcf \
# -a data/annotation/probe_bed/Mouse430_2.mm9.bed \
# -b /notebook/code/BEDTools-Version-2.15.0/bin \
# -s na32 \
# -g mm9 \
# -d 37 \
# -c Mouse430_2 \
# -i 1.0 \
# -r "Mus musculus" \
# -w mouse \
# -u "David Quigley" \
# -l dquigley@cc.ucsf.edu \
# -n data/annotation/M430/Mouse430_2.cdf \
# -x data/annotation/M430/Mouse430_2.probe_tab \
# -z data/annotation/M430/Mouse430_2.CEL \
# -o results/M430_FVB_SPR

import subprocess
import sys
import os
import datetime
from optparse import OptionParser
import uuid

    
def require_file(fn, what_is_it):
    if len(fn)==0:
        sys.exit("ERROR: did not pass a value for " + what_is_it )
    
    if not os.path.exists(fn):
        sys.exit("ERROR: " + what_is_it + " " + fn + " does not exist, exiting")    

def require_not_zero(str, msg):    
    if len(str)==0:
        sys.exit(msg)
 
def validate_input(options):
    """ check input parameters, ensuring required files exist """
    
    new_package_name = options.newpackage
    if len(new_package_name)==0:
        sys.exit("ERROR: must provide a new package name")
    
    chip_format = options.chip_format
    if chip_format not in ["gene", "IVT"]:
        sys.exit("ERROR: must specify chip format as one of {gene,IVT}")
    
    fn_vcf = options.fn_list_vcf.split(",")
    require_not_zero( fn_vcf, "no VCF files passed" )
    
    for fn in fn_vcf:
        require_file( fn, "VCF file" )
    
    require_file( options.fn_affy_bed, "affymetrix probe BED file")
    require_file( options.dir_bedtools, "BEDtools folder" )
    require_file( options.dir_bedtools + "/intersectBed", "bedtools folder missing intersectBed" )
    require_file( options.dir_bedtools + "/sortBed", "bedtools folder missing sortBed" )
    require_file( options.fn_ps, "probeset.csv" )
    
    if chip_format=="gene":
        require_file( options.fn_pgf, "affymetrix PGF" )
        require_file( options.fn_clf, "affymetrix CLF" )
        require_file( options.fn_mps, "affymetrix MPS" )
        require_file( options.fn_ts, "transcript.csv" )
    
    if chip_format=="IVT":
        require_file( options.fn_CDF, "affymetrix CDF" )
        require_file( options.fn_probe_tab, "affymetrix probe_tab" )
        require_file( options.fn_CEL, "affymetrix CEL file (text format, not binary)" )
    
    require_file( options.dir_out, "output directory" )
    
    require_not_zero( options.annot_version, "must indicate annotation version (e.g. na32)")
    require_not_zero( options.annot_revision, "must indicate annotation verison (e.g. 4)")
    require_not_zero( options.genome_version_ucsc, "must indicate annotation version (e.g. mm9)")
    require_not_zero( options.genome_version_ncbi, "must indicate annotation version (e.g. 37)")
    require_not_zero( options.chip, "must indicate chip (e.g. MoGene-1_1-st-v1)")
    require_not_zero( options.species, "must pass an species (e.g. Mus musculus)")
    
    require_not_zero( options.author, "must pass an Author (e.g. Gregor Mendel)")
    require_not_zero( options.email, "must pass an email address (e.g. gregor@wrinkledpea.org)")
    require_not_zero( options.organism, "must pass an organism (e.g. mouse)")        
    
    if options.dir_out == ".":
        print "MESSAGE: no output directory passed, writing files to current directory"
    
    P ={}
    P["new_package_name"] = options.newpackage
    P["fn_vcf"] = options.fn_list_vcf.split(",")
    P["fn_affy_bed"] = options.fn_affy_bed
    P["dir_bedtools"] = options.dir_bedtools
    P["annot_version"] = options.annot_version # na32
    P["annot_revision"] = options.annot_revision # 4
    P["genome_version_ucsc"] = options.genome_version_ucsc.lower() # mm9
    P["genome_version_ncbi"] = options.genome_version_ncbi.lower() # 37
    P["chip"] = options.chip # MoGene-1_1-st-v1
    P["chip_version"] = options.chip_version
    P["chip_format"] = options.chip_format
    P["organism"] = options.organism # Mouse
    P["species"] = options.species # Mus musculus
    P["author"] = options.author # David Quigley
    P["email"] = options.email # dquigley@cc.ucsf.edu
    P["fn_ps_orig"] = options.fn_ps
    P["fn_ts_orig"] = options.fn_ts
    P["fn_pgf_orig"] = options.fn_pgf
    P["fn_mps_orig"] = options.fn_mps
    P["fn_clf_orig"] = options.fn_clf
    P["fn_cdf_orig"] = options.fn_CDF
    P["fn_CEL"] = options.fn_CEL    
    P["fn_probetab_orig"] = options.fn_probe_tab
    P["dir_out"] = options.dir_out
    if P["dir_out"]==".":
            P["dir_out"] = os.getcwd()

    P["fn_probe_count_changes"] = P["dir_out"] + "/probe_count_changes.txt"
    P["fn_R_script"] = P["dir_out"] + '/create_package.R'
    P["fn_probesets_with_snps"] = P["dir_out"] + '/combined.bed'
    
    P["fn_cdf_new"] = P["dir_out"] + "/" + P["new_package_name"] + ".cdf"    
    P["fn_mps_new"] = P["dir_out"] + "/" + P["new_package_name"] + ".mps"
    P["fn_pgf_new"] = P["dir_out"] + "/" + P["new_package_name"] + ".pgf"
    P["fn_ps_new"] =  P["dir_out"] + "/" + P["new_package_name"] + ".probeset.csv"
    P["fn_ts_new"] =  P["dir_out"] + "/" + P["new_package_name"] + ".transcript.csv"
    P["fn_probetab_new"] =  P["dir_out"] + "/" + P["new_package_name"] + ".probe_tab"

    P["fn_combined_bed"] = P["dir_out"] + "/" + "combined.bed"
    
    return P

def find_probe_intersections( P ):
    """ Call intersectBed to build a single large BED for all intersections between SNPs and affy """
    
    if os.path.exists( P["fn_combined_bed"] ):
        print "MESSAGE: Combined BED file indicating SNP intersections already exists: " + P["fn_combined_bed"]
        print "MESSAGE: Skipping BEDtools step. To re-create this file, delete it and re-run the script."
        return
    
    cmd = ""
    for i,fn in enumerate(P["fn_vcf"]):
        cmd += dir_bedtools + "/intersectBed -u -a " + P["fn_affy_bed"] + " -b " + fn + " > " + P["dir_out"] + "/snps_from_vcf_" + str(i+1) + ".bed\n"
    
    cmd += "cat" 
    for i,fn in enumerate(P["fn_vcf"]):
        cmd += " " + P["dir_out"] + "/snps_from_vcf_" + str(i+1) + ".bed"
    cmd += " | sortBed > " + P["fn_combined_bed"] + "\n"

    fo = open( P["dir_out"] + "/combine_vcf.sh", 'w')
    fo.write(cmd)    
    fo.close()

    print "MESSAGE: Executing BEDtools commands to combine VCF files."
    print "MESSAGE: Please be patient, as this can take several minutes...\n"
    subprocess.call(["sh", P["dir_out"] + "/combine_vcf.sh"])


class Probeset():
    def __init__(self, probe_id, probe_type, probe_name):
        self.probe_id = probe_id
        self.probe_type = probe_type
        self.probe_name = probe_name        
        self.probe_atoms = []
        self.at_least_one_valid_probe = True
    
    def __str__( self ):
        s = self.probe_id + '\t' + self.probe_type + '\t' + self.probe_name + '\n'
        for atom in self.probe_atoms:
            s += str(atom)
            for probe in atom.probes:
                s += str(probe)
        
        return s
    
    def remove_probes_by_idx(self, idx):
        # Having problems with probesets containing no probes.
        # to work around, always keeping one probe in probeset, but marking invalid if 
        # that probeset should have been removed.
        tmp = []
        for i,atom in enumerate(self.probe_atoms):
            if not i in idx:
                tmp.append( atom )
        if len(tmp)==0:
            tmp.append(atom)
            self.at_least_one_valid_probe=False
        
        self.probe_atoms = tmp
        return self.at_least_one_valid_probe
        
    def n_probes(self):
        n = 0
        for atom in self.probe_atoms:
            n += len(atom.probes)
        
        return n

class Probe_atom():
    # one atom can have multiple probes in the case of controls probes, e.g. 10338002
    def __init__(self, probe_atom ):
        self.probe_atom = probe_atom
        self.probes = []
    
    def __str__(self):
        if len(self.probes)>0: # Modified to only write atom if there are probes
            return '\t' + self.probe_atom + '\n'

class Probe():
    def __init__(self, probe_id, type,gc_count,probe_length,interrogation_position,probe_sequence):
        self.probe_id = probe_id
        self.type = type
        self.gc_count = gc_count
        self.probe_length = probe_length
        self.interrogation_position = interrogation_position
        self.probe_sequence = probe_sequence
    
    def __str__( self ):
        s='\t\t' + self.probe_id + '\t' + self.type + '\t' + self.gc_count + '\t' + self.probe_length
        s += '\t' + self.interrogation_position + '\t' + self.probe_sequence + '\n'
        return s

def parse_PGF(fn_pgf):
    # Read in the PGF file and store it as a vector of probeset objects
    PS = None
    f = open(fn_pgf)
    probesets = []
    for line in f:
        if line[0]=="#":
            continue
        linestrip = line.rstrip('\r\n')
        a = linestrip.split('\t')
        if len(a[0])>0:
            # this is a new probeset
            if not PS is None:
                probesets.append(PS)
            
            PS = Probeset( a[0], a[1], a[2] )
        
        else:
            if len(a)==2:
                ATOM = Probe_atom(a[1])
                PS.probe_atoms.append(ATOM)
            else:
                x,y, probe_id, type,gc_count,probe_length,interrogation_position,probe_sequence = a
                P = Probe( probe_id, type,gc_count,probe_length,interrogation_position,probe_sequence )
                ATOM.probes.append(P)
    
    probesets.append(PS)
    return probesets

def parse_pre_post(fn):
    # transcript.id	orig.n.probes	probeset.id	orig.probeset.n.probes	final.probeset.n.probes
    f = open(fn)
    h = f.readline()
    ps2n = {}
    ts2n = {}
    ps2n_orig = {}
    for line in f:
        a = line.rstrip('\r\n').split('\t')
        try:
            ts2n[a[0]] += int(a[4])
        except:
            ts2n[a[0]] = int(a[4])
        ps2n[a[2]] = int(a[4])
        ps2n_orig[a[2]] = int(a[3])
    
    f.close()
    
    return ps2n, ts2n, ps2n_orig


def check_sequence_chr_consistency( fn_vcf, fn_affy_bed ):
    """ Confirm that all VCF files and the affymetrix BED file agree about whether column one
    starts with the string 'chr' or not"""
    
    bed_has_chr = False
    for line in open(fn_affy_bed):
        if line[0] != "#":
            bed_has_chr = "chr" in line.split('\t')[0]
            break
    
    print "MESSAGE: Does Affymetrix bed file uses 'chr' format? " + str(bed_has_chr)
    agrees = True
    for fn in fn_vcf:
        f = open(fn)
        for line in f:
            if line[0] != "#":
                if ("chr" in line.split('\t')[0]) != bed_has_chr:
                    agrees = False
                break
    
    return agrees


def process_gene_format(P):
    """ Generate new files for Affymetrix Gene ST format 
        I am well aware of how comically large this function is; however, it works, and 
        there is only one path through the code for all of these statements.
    """

    #---------------------------------------------------
    # identify probesets containing a SNP (BED file)
    # load these into a remove hash. We need both the probeset
    # identifier and the probe start location to uniquely identify a probe, 
    # since multiple probes are identified by the same probeset ID in the 
    # Affymetrix probe BED file.
    #---------------------------------------------------
    print "MESSAGE: Reading probesets with SNPs file:"
    print "         " + P["fn_probesets_with_snps"]
    f = open(P["fn_probesets_with_snps"])
    ps_start_remove = {}
    for line in f:
        a = line.rstrip('\r\n').split('\t')
        chrom, start, stop, probeset = a[0:4]
        ps_start_remove[ probeset + ' ' + start ] = 1

    f.close()
    print "MESSAGE: VCF file analysis identified " + str(len(ps_start_remove.keys())) + " probes containing SNPs."

    #---------------------------------------------------
    # spin through probe BED, identifying the index of probes within a 
    # probeset that are within the remove hash. The probes in the PGF file
    # within a probeset can be identified only by index which seems to match 
    # the index in which they are present into the Affymetrix BED file.
    #
    # probe_indexes_to_remove is a hash of probeset->list of indexes to remove
    #---------------------------------------------------
    f = open(P["fn_affy_bed"])
    current_probe = ""
    remove_list = []
    probe_indexes_to_remove = {}
    marked_probes_found_in_BED = 0
    for line in f:
        a = line.rstrip('\r\n').split('\t')
        chrom, start, stop, probeset = a[0:4]
        if probeset != current_probe:
            probe_idx = 0
            current_probe = probeset
    
        if probeset + " " + start in ps_start_remove:
            transcript_id, probeset_id = probeset.split("_")
            remove_list.append( [transcript_id, probeset_id, probe_idx ] )
            try:
                probe_indexes_to_remove[probeset_id].append(probe_idx)
            except KeyError:
                probe_indexes_to_remove[probeset_id] = [probe_idx]
            marked_probes_found_in_BED += 1
        probe_idx += 1

    f.close()

    print "MESSAGE: Found " + str(marked_probes_found_in_BED) + " of these probes in the affymetrix BED file"

    #---------------------------------------------------
    # rewrite the PGF file
    # Representing this nested file format as objects makes it easier to parse
    #-------------------------------------------------
    print "MESSAGE: Parsing affymetrix PGF file..."
    probesets = parse_PGF(P["fn_pgf_orig"])

    # find the index of each probeset
    ids = [x.probe_id for x in probesets]
    id2idx = {}
    for idx, id in enumerate(ids):
        id2idx[id] = idx

    # How many probes does each probeset have before we touch it?
    ps2before = {}
    n_probes = 0
    for id in ids:
        ps2before[id] = probesets[id2idx[id]].n_probes()
        n_probes += ps2before[id]

    print "MESSAGE: Counted " + str(len(probesets)) + " probesets containing " + str(n_probes) + " probes in the original PGF file"

    # Remove probes from probesets. Count remaining probes in ps2after.
    ps2after = {}
    n_probes_removed = 0 
    probesets_no_valid_probes = {}
    for probeset in probe_indexes_to_remove.keys():
        indexes = sorted( probe_indexes_to_remove[probeset], reverse=True)
        n_probes_removed += len(indexes)
        if not probesets[ id2idx[probeset] ].remove_probes_by_idx(indexes):
            probesets_no_valid_probes[ probeset ] = 1
        ps2after[probeset] = probesets[ id2idx[probeset] ].n_probes() 

    n_probes_after = 0
    for id in ids:
        n_probes_after += probesets[id2idx[id]].n_probes()

    print "MESSAGE: Removed " + str(n_probes_removed) + " probes for SNP adjustement."
    print "MESSAGE: Counted " + str(n_probes_after) + " probes after SNP adjustement."

    # Write out a summary table of transcript_id probeset_id before after
    # Use the stock MPS file to map probesets to transcript IDs and get probe counts
    print "MESSAGE: Parsing Affymetrix MPS file. This file assigns probesets to transcripts."
    fo = open(P["fn_probe_count_changes"], 'w')
    fo.write( "transcript.id\torig.n.probes\tprobeset.id\torig.probeset.n.probes\tfinal.probeset.n.probes\tinclude.in.MPS\n" )
    f = open(P["fn_mps_orig"])
    n_empty_probesets = 0
    n_mps_probesets_at_least_one_probe = 0
    n_mps_probesets_total = 0
    for line in f:
        if line[0]!="#":
            ps, ts, probes, length = line.rstrip('\r\n').split('\t')
            if ts=="":
                ts = ps
        
            if length=="":
                length="0"
                fo.write( ps + "\t" + length + "\t" + ps + "\t" + length + '\t' + length + '\tno\n')
                n_empty_probesets += 1
        
            elif ps != "probeset_id":
                for probe in probes.split(" "):
                    n_mps_probesets_total += 1
                    n_before = str(ps2before[probe])
                    if probe in ps2after:
                        n_after = str(ps2after[probe])
                
                    else:
                        n_after = n_before
                    if probe in probesets_no_valid_probes:
                        include_in_MPS = "no"
                    else:
                        include_in_MPS = "yes"
                        n_mps_probesets_at_least_one_probe += 1
                    fo.write( ps + "\t" + length + "\t" + probe + "\t" + n_before + '\t' + n_after + '\t' + include_in_MPS + '\n')

    fo.close()

    print "MESSAGE: The Affymetrix MPS file contained " + str(n_mps_probesets_total) + " probesets"
    print "MESSAGE: The Affymetrix MPS file contained " + str(n_mps_probesets_at_least_one_probe) + " probesets with >= 1 probe after correction"
    print "MESSAGE: The Affymetrix MPS file contained " + str(n_empty_probesets) + " probesets that did not contain a probe assignment"
    print "MESSAGE: NOTE: Not every probeset in the PGF file will be assigned to a transcript,"
    print "MESSAGE:       even in the stock Affymetrix assignments. This accounts for the difference"
    print "MESSAGE:       in probeset count between the PGF and MPS files and for the difference "
    print "MESSAGE:       in probe counts as well."
    print "MESSAGE: Wrote a table of changes to probeset composition to " + P["fn_probe_count_changes"]

    ps2n, ts2n, ps2n_orig = parse_pre_post(P["fn_probe_count_changes"])
    n_probesets_orig=0
    for probeset in ps2n_orig.keys():
        n_probesets_orig += ps2n_orig[probeset]
    
    n_probes_after=0
    for probeset in ps2n.keys():
        n_probes_after += ps2n[probeset]

    print "MESSAGE: Summary of changes:"
    print "MESSAGE: Before accounting for SNPs, PGF contains " + str(len(ts2n)) + " transcripts " + str(len(ps2n_orig)) + " probesets " + str(n_probesets_orig) + " probes"
    print "MESSAGE: After accounting for SNPs, PGF contains  " + str(len(ts2n)) + " transcripts "  + str(len(ps2n)) + " probesets " + str(n_probes_after) + " probes"

    timestamp = datetime.datetime.today().strftime("%a %b %d %H:%M:%S PST %Y")
    guid = str(uuid.uuid1())

    # Write new PGF file
    fo = open(P["fn_pgf_new"], 'w')
    fo.write("""#%chip_type=""" + new_package_name + """
    #%lib_set_name=""" + new_package_name + """
    #%lib_set_version=""" + annot_revision + """
    #%create_date=""" + timestamp + """
    #%guid=""" + guid + """
    #%pgf_format_version=1.0
    #%rows=1190
    #%cols=990
    #%probesets=""" + str(len(probesets)) + """
    #%datalines=2023948
    #%sequential=1
    #%order=col_major
    #%header0=probeset_id	type	probeset_name
    #%header1=	atom_id
    #%header2=		probe_id	type	gc_count	probe_length	interrogation_position	probe_sequence""" + '\n')
    for probeset in probesets:
        fo.write( str(probeset) )

    fo.close()

    print "MESSAGE: wrote new PGF file to " + P["fn_pgf_new"]
 
    #------------------------------------------------------------------
    # Rewrite the MPS file
    # In the MPS format the first two columns are identical.
    # The third is space-delimited list of probesets.
    # The fourth is count of probes.
    # If a probeset contains zero probes, we will remove it.
    # If the number of probes has changed, we'll change sixth column in transcript.csv 
    # to update with new number of probes
    #
    # Read probes for exclusion out of probe_count_changes
    # For each transcript we want 
    #   the number of probes remaining
    #   the list of probesets with at least one probe
    # 
    # Rewrite MPS file
    # MPS links individual probesets to a transcript
    # Apparently the case of no probesets for a transcript is
    # ps {empty}  ps {empty}
    #----------------------------------------
    f = open(P["fn_mps_orig"])
    fo = open(P["fn_mps_new"], 'w')
    fo.write("#%chip_type=" + P["new_package_name"] + "\n" )
    fo.write("#%lib_set_name=" + P["new_package_name"] + "\n" )
    fo.write("#%lib_set_version=r" + P["annot_revision"] + "\n" )
    fo.write("#%guid=" + guid + "\n")
    fo.write("#%create_date=" + timestamp + "\n" )
    fo.write("#%genome-species=" + P["organism"] + "\n" )
    fo.write("#%genome-version=" + P["genome_version_ucsc"] + "\n" )
    fo.write("#%genome-version-ucsc=" + P["genome_version_ucsc"] + "\n" )
    fo.write("#%genome-version-ncbi=" + P["genome_version_ncbi"] + "\n" )
    fo.write("#%genome-version-create_date=2007 July\n")
    for line in f:
        if line[0] != "#":
            ps, ts, probes, length = line.rstrip('\r\n').split('\t')
            if ps == "probeset_id":
                fo.write(line.rstrip('\r\n') + '\r')
        
            else:
                length = str(ts2n[ps])
                probes_not_zero = []
                for probe in probes.split(" "):
                    if ps2n[probe]>0 and not probe in probesets_no_valid_probes:
                        probes_not_zero.append(probe)
                if len(probes_not_zero) == 0:
                    fo.write( ps + '\t\t' + ps + '\t\n') # matches MPS entries
                else:
                    fo.write( ps + '\t' + ps + '\t' + ' '.join(probes_not_zero) + '\t' + length + '\n' )

    fo.close()
    f.close()

    print "MESSAGE: wrote new MPS file to " + P["fn_mps_new"]

    #------------------------------------------------------------------------
    # Modify probeset and transcript files 
    # For probeset, do not write if no probes are left.
    # For transcript, update number of probes. 
    # TOTAL PROBESETS: 212,639
    # TOTAL TRANSCRIPTS: 35,556
    #------------------------------------------------------------------------

    n_probes, n_transcripts = 0,0
    f = open(P["fn_ps_orig"])
    fo = open(P["fn_ps_new"], "w")
    for line in f:
        if line[0] =="#":
            fo.write(line.rstrip('\r\n') + '\n' )
    
        else:
            a = line.rstrip('\r\n')[1:-1].split('","')
            if a[0]=="probeset_id":
                fo.write(line.rstrip('\r\n') + '\n' )
            elif int( ps2n[a[0]]) > 0:
                fo.write(line.rstrip('\r\n') + '\n' )
                n_probes += 1

    fo.close()
    f.close()

    print "MESSAGE: wrote new probeset description file to " + P["fn_ps_new"]

    f = open(P["fn_ts_orig"])
    fo = open(P["fn_ts_new"], "w")
    for line in f:
        if line[0] =="#":
            fo.write(line.rstrip('\r\n') + '\n' )
    
        else:
            a = line.rstrip('\r\n')[1:-1].split('","')
            if a[0]=="transcript_cluster_id":
                fo.write(line.rstrip('\r\n') + '\n' )
        
            else:
                a[6] = str( ts2n[a[0]] )
                out = '"' + '","'.join(a) + '"'
                fo.write(out + '\n' )
                n_transcripts += 1

    fo.close()
    f.close()

    print "MESSAGE: wrote new transcript description file to " + P["fn_ts_new"]

def process_block(block_lines, probeset2idx_remove, unithead2details):
    """ Input is an array of lines read from the cdf
        we want to catch blocks which contain _Block1 in the first (header) line
        for these, we need to read past the boilerplate header and then check if
        any of the probe definition lines are in probeset2idx_remove
        If they are, we need to remove them from block_header (which will be rewritten)
        and also capture their X,Y values.
    
    """
    if len(block_lines)==0:
        return [], []
    
    standard_header = "CellHeader=X\tY\tPROBE\tFEAT\tQUAL\tEXPOS\tPOS\tCBASE\tPBASE\tTBASE\tATOM\tINDEX\tCODONIND\tCODON\tREGIONTYPE\tREGION"
    
    block_header = block_lines[0][1:-1]
    XY_removed = []
    lines_kept = 11 # by default keep all
    if "_Block" in block_header:
        # This is a probeset block, which may need modification
        probeset = block_lines[1].split("=")[1]
        if probeset in probeset2idx_remove:
            idx_to_be_removed = probeset2idx_remove[probeset]
            new_block_lines = []
            lines_kept = 0
            for line in block_lines[1:]:
                line_type, line_value=line.split("=")
                if "Cell" in line_type and not "Header" in line_type and not "NumCells" in line_type:
                    a = line.split('\t')
                    X, Y, line_idx = a[0], a[1], a[5]
                    if line_idx in idx_to_be_removed:
                        XY_removed.append( [ probeset, X,Y ] )
                    else:
                        new_block_lines.append( line )
                        lines_kept += 1
            
            # rebuild block_lines, first boilerplate header than the remaining probe lines
            block_lines = ["[" + block_header + "]" ]
            block_lines.append( "Name=" + probeset )
            block_lines.append( "BlockNumber=1" )
            block_lines.append( "NumAtoms=" + str(lines_kept/2) )
            block_lines.append( "NumCells=" + str(lines_kept) )
            block_lines.append( "StartPosition=0" )
            block_lines.append( "StopPosition=10" ) # not clear if this should change or not
            block_lines.append( standard_header )
            for nbl in new_block_lines:
                block_lines.append( nbl )
    
        unit_head = block_header.split("_")[0]
        unit_direction, unit_type = unithead2details[ unit_head ] 
        unit_number = unit_head.replace("Unit","")
        header = [ "[" + unit_head + "]" ]
        header.append( "Name=NONE" )
        header.append( unit_direction )
        header.append( "NumAtoms="+str(lines_kept) )
        header.append( "NumCells=" + str(2 * lines_kept) )
        header.append("UnitNumber=" + unit_number )
        header.append(unit_type)
        header.append("NumberBlocks=1")
        header.append("")

        header.extend(block_lines)
        block_lines = header
        block_lines.append("")
    
    block_lines.append("")    
    return block_lines, XY_removed

def remove_probes_from_CDF(targets, fn_orig, fn_final):
    f = open(fn_orig)
    fo = open(fn_final, "w")
    XY_removed = {}
    block = []
    unithead2details = {}
    for line in f:
        line = line.rstrip('\r\n')
        if len(line)==0:
            continue
        
        if line[0] == "[":
            if len(block)>0:
                block_header = block[0][1:-1]
                if "Unit" in block_header and not "_Block" in block_header:
                    unithead2details[block_header] = [ block[2], block[6] ] # store direction, UnitType
                else:    
                    new_block, XY_removed_from_block = process_block(block, targets, unithead2details)
                    for new_line in new_block:
                        fo.write( new_line + '\n' )
                    for X, Y, probeset in XY_removed_from_block:
                        XY_removed[ probeset + "_" + X + "_" + Y ] = 1
                block = []
        
        block.append(line)
    
    
    new_block, XY_removed_from_block = process_block(block, targets,unithead2details)
    
    for new_line in new_block:
        fo.write( new_line + '\n' )
    for ps_x_y in XY_removed:
        XY_removed[ ps_x_y ] = 1
    
    return XY_removed

def update_probetab(XY_removed, fn_probetab_orig, fn_probetab_new):
    f = open(fn_probetab_orig)
    fo = open(fn_probetab_new, "w")
    fo.write(f.readline()) # header
    for line in f:
        a = line.rstrip('\r\n').split('\t')
        id = '_'.join( [a[0], a[1], a[2]] )
        if not id in XY_removed:
            fo.write(line)
    
    
def process_IVT_format(P):
    
    #---------------------------------------------------
    # identify probesets containing a SNP (BED file)
    # load these into a remove hash. 
    # for IVT, probesets in the hash will have the form 1415670_at:1
    #---------------------------------------------------
    print "MESSAGE: Reading probesets with SNPs file:"
    print "         " + P["fn_probesets_with_snps"]
    f = open(P["fn_probesets_with_snps"])
    ps_remove = {}
    for line in f:
        a = line.rstrip('\r\n').split('\t')
        chrom, start, stop, probeset = a[0:4]
        probeset, probe_index = probeset.split(":")
        if probeset in ps_remove:
            ps_remove[ probeset ].append(probe_index)
        else:
            ps_remove[ probeset ] = [probe_index]

    f.close()
    print "MESSAGE: VCF file analysis identified " + str(len(ps_remove.keys())) + " probes containing SNPs."
    print "MESSAGE: Reading CDF file:"
    print "         " + P["fn_cdf_orig"]
    
    # spin through cdf, removing atoms which are overlapped. Store probeset_id\tX\tY
    # Then spin through Mouse430_2.probe_tab, removing overlaps
    XY_removed=remove_probes_from_CDF(ps_remove, P["fn_cdf_orig"], P["fn_cdf_new"])
    update_probetab(XY_removed, P["fn_probetab_orig"], P["fn_probetab_new"])


#------------------------------------------------------------------------------------------
# Begin main 
#------------------------------------------------------------------------------------------

parser = OptionParser()
parser.add_option("-f", "--format", dest="chip_format", help="Chip format, one of {gene,IVT}", default="")
parser.add_option("-p", "--packagename", dest="newpackage", help="New package name, default MoGenenoSnps-1_1-st-v1", default="MoGenenoSnps-1_1-st-v1")
parser.add_option("-v", "--vcf", dest="fn_list_vcf", help="Path to VCF files, comma-delimited", default="")
parser.add_option("-a", "--affy", dest="fn_affy_bed", help="Path to affymetrix probe BED file", default="")
parser.add_option("-b", "--bedtools", dest="dir_bedtools", help="Path to bedTools bin folder", default="")

parser.add_option("-s", "--annot_version", dest="annot_version", help="Annotation version, default na32", default="na32")
parser.add_option("-e", "--annot_revision", dest="annot_revision", help="Affy chip library revision, default 4", default="4")
parser.add_option("-g", "--genome_version_ucsc", dest="genome_version_ucsc", help="UCSC Genome version, default mm9", default="mm9")
parser.add_option("-d", "--genome_version_ncbi", dest="genome_version_ncbi", help="NCI Genome version, default 37", default="37")
parser.add_option("-c", "--chip", dest="chip", help="Affy chip, default MoGene-1_1-st-v1", default="MoGene-1_1-st-v1")
parser.add_option("-i", "--chip_version", dest="chip_version", help="chip version string, default 1.1", default="1.1")
parser.add_option("-r", "--species", dest="species", help="species name, default Mus musculus", default="Mus musculus")

parser.add_option("-u", "--author", dest="author", help="Author Name for package, default David Quigley", default="David Quigley")
parser.add_option("-l", "--email", dest="email", help="Author email for package, default dquigley@cc.ucsf.edu", default="dquigley@cc.ucsf.edu")
parser.add_option("-w", "--organism", dest="organism", help="species name, default Mouse", default="Mouse")

parser.add_option("-n", "--file_CDF", dest="fn_CDF", help="Path to affymetrix CDF file (IVT format only)", default="")
parser.add_option("-x", "--file_probe_tab", dest="fn_probe_tab", help="Path to affymetrix probe_tab (IVT format only)", default="")
parser.add_option("-z", "--file_CEL", dest="fn_CEL", help="Path to sample Affymetrix CEL file (IVT format only)", default="")

parser.add_option("-y", "--file_probeset_csv", dest="fn_ps", help="Path to affymetrix probeset.csv", default="")
parser.add_option("-t", "--file_transcript_csv", dest="fn_ts", help="Path to affymetrix transcript.csv (gene format only)", default="")
parser.add_option("-q", "--file_pgf", dest="fn_pgf", help="Path to chip PGF file from Affymetrix library (gene format only)", default="")
parser.add_option("-m", "--file_mps", dest="fn_mps", help="Path to chip MPS file from Affymetrix library (gene format only)", default="")
parser.add_option("-k", "--file_clf", dest="fn_clf", help="Path to chip CLF file from Affymetrix library (gene format only)", default="")

parser.add_option("-o", "--output", dest="dir_out", help="Path to write package files, default current directory", default=".")

(options, args) = parser.parse_args()
P = validate_input( options ) # Check to make sure requested files exist and commands are coherent

print "MESSAGE: New package will be named: " + P["new_package_name"]
print "MESSAGE: New files will be written to: " + P["dir_out"]
print "MESSAGE: Affymetrix probe sequence file: " + P["fn_affy_bed"]
print "MESSAGE: User passed the following VCF files to match against probes:"
for fn in P["fn_vcf"]:
    print "         " + fn

print "MESSAGE: User indicated chip was " + P["chip"] + ", version " + P["annot_revision"] + " annotation version " + P["annot_version"]

if not check_sequence_chr_consistency( P["fn_vcf"], P["fn_affy_bed"] ):
    print "ERROR: Affymetrix probe BED file does not have the same chromosome string encoding as one or"
    print "       or more of the VCF files. Chromosomes can be encoded as 'chr1' or '1'; the probe_bed "
    print "       VCF files must all use the same encoding."
    sys.exit(0)

find_probe_intersections( P )

if P["chip_format"]=="gene":
    process_gene_format( P )
else:
    process_IVT_format( P )

r_script = """
if(!require("pdInfoBuilder")){
    print("Attempting to install pdInfoBuilder package")
    source("http://bioconductor.org/biocLite.R")
    biocLite("pdInfoBuilder")
}
library(pdInfoBuilder)
"""
if P["chip_format"]=="gene":
    r_script += "fn.probes = '" + P["fn_ps_new"] + "'\n"
    r_script += "fn.pgf = '" + P["fn_pgf_new"] + "'\n"
    r_script += "fn.clf = '" + P["fn_clf_new"] + "'\n"
    r_script += "fn.mps = '" + P["fn_mps_new"] + "'\n"
    r_script += "fn.trans = '" + P["fn_trans_new"] + "'\n\n"
    r_script += "seed <- new('AffyExpressionPDInfoPkgSeed', "
else:
    r_script += "fn.tab = '" + P["fn_probetab_new"] + "'\n"
    r_script += "fn.cel = '" + P["fn_CEL"] + "'\n"
    r_script += "fn.cdf = '" + P["fn_cdf_new"] + "'\n\n"
    r_script += "seed <- new('AffyExpressionPDInfoPkgSeed', "

r_script += "chipName = '" + P["new_package_name"] + "', "
r_script += "version = '" + P["chip_version"] + "', "

if P["chip_format"]=="gene":
    r_script += "pgfFile = fn.pgf, "
    r_script += "clfFile = fn.clf, "
    r_script += "coreMps = fn.mps, "
    r_script += "transFile = fn.trans, "
else:
    r_script += "cdfFile = fn.cdf, "
    r_script += "celFile = fn.cel, "
    r_script += "tabSeqFile = fn.tab, "

r_script += "author = '" + P["author"] + "', email = '" + P["email"] + "', biocViews = 'AnnotationData', "
r_script += "genomebuild = '" + P["genome_version_ucsc"].upper() + "', organism = '" + P["organism"] + "', "
r_script += "species = '" + P["species"] + "', url ='')\n"
r_script += "makePdInfoPackage(seed, destDir = '" + P["dir_out"] + "', unlink=T)\n"
pd_name = "pd." + P["new_package_name"].lower().replace("_",".")
r_script += "#install.packages('" + P["dir_out"] + "/" + pd_name + "', repos = NULL, type='source')"


fo = open(P["fn_R_script"], 'w')
fo.write(r_script)
fo.close()

print "MESSAGE: Wrote R script to create package: " + P["fn_R_script"]
print "MESSAGE: There are several ways to run this script:"
print "         1) Open the file, run R, and paste the contents into R"
print "         2) From the command line, type: "
print "            R " + P["fn_R_script"]
print "         The script will attempt to install the bioconductor 'pdInfoBuilder' package if "
print "         it is not already present on your machine. If you do not have the ability to "
print "         install packages on your machine, contact your local system administrator."
print "         To install the new package, uncomment and run the last line of the R script."
print ""
print "         IMPORTANT FOR OS X users running R v2.14.0 or newer!"
print "         If you have installed the binary build of affxparser v1.26.2, it will crash"
print "         when you attempt to build the package. See: "
print "              https://groups.google.com/forum/#!topic/aroma-affymetrix/lEfDanThLEA/discussion"
print "              https://stat.ethz.ch/pipermail/bioc-devel/2011-November/002969.html"
print "         You can install a working binary build of affxparser v1.26.2 and avoid the crash with: "
print "         remove.packages('affxparser')"
print "         source('http://www.braju.com/R/hbLite.R)"
print "         hbBiocLite('affxparser');"

