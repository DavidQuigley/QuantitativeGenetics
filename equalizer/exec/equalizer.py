# Copyright 2015 David Quigley
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import subprocess
import sys
import os
import datetime
from optparse import OptionParser
import uuid

VERSION = "0.3"
RELEASE_DATE = "June 05 2015"
    
def require_file(fn, what_is_it):
    if len(fn)==0:
        print( "ERROR: did not pass a value for " + what_is_it + "\n" )
        sys.exit("For usage instructions, call equalizer.py --help or visit http://http://github.com/DavidQuigley/QuantitativeGenetics/equalizer")
    
    if not os.path.exists(fn):
        print( "ERROR: " + what_is_it + " " + fn + " does not exist, exiting\n" )
        sys.exit("For usage instructions, call equalizer.py --help or visit http://http://github.com/DavidQuigley/QuantitativeGenetics/equalizer")
        
def require_not_zero(str, msg):    
    if len(str)==0:
        print( msg + '\n' )
        sys.exit("For usage instructions, call equalizer.py --help or visit http://http://github.com/DavidQuigley/QuantitativeGenetics/equalizer")


def validate_input(options):
    """ check input parameters, ensuring required files exist """
    
    new_package_name = options.newpackage
    if len(new_package_name)==0:
        print( "ERROR: must provide a new package name\n" )
        sys.exit("For usage instructions, call equalizer.py --help or visit http://http://github.com/DavidQuigley/QuantitativeGenetics/equalizer")
    
    chip_format = options.chip_format
    if chip_format not in ["gene", "IVT"]:
        print( "ERROR: must specify chip format as one of {gene,IVT}\n" )
        sys.exit("For usage instructions, call equalizer.py --help or visit http://http://github.com/DavidQuigley/QuantitativeGenetics/equalizer")
    
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
    require_not_zero( options.annot_revision, "must indicate annotation revision (e.g. 4)")
    require_not_zero( options.genome_version_ucsc, "must indicate annotation version (e.g. mm9)")
    require_not_zero( options.genome_version_ncbi, "must indicate annotation version (e.g. 37)")
    require_not_zero( options.chip, "must indicate chip (e.g. MoGene-1_1-st-v1)")
    require_not_zero( options.species, "must pass an species (e.g. Mus musculus)")
    
    require_not_zero( options.author, "must pass an Author (e.g. Gregor Mendel)")
    require_not_zero( options.email, "must pass an email address (e.g. gregor@wrinkledpea.org)")
    require_not_zero( options.organism, "must pass an organism (e.g. mouse)")        
    
    if options.dir_out == ".":
        print( "MESSAGE: no output directory passed, writing files to current directory" )
    
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
    P["fn_clf_new"] = P["fn_clf_orig"] # do not rewrite the clf file
    P["fn_pgf_new"] = P["dir_out"] + "/" + P["new_package_name"] + ".pgf"
    P["fn_ps_new"] =  P["dir_out"] + "/" + P["new_package_name"] + ".probeset.csv"
    P["fn_ts_new"] =  P["dir_out"] + "/" + P["new_package_name"] + ".transcript.csv"
    P["fn_probetab_new"] =  P["dir_out"] + "/" + P["new_package_name"] + ".probe_tab"

    P["fn_combined_bed"] = P["dir_out"] + "/" + "combined.bed"
    
    return P

def find_probe_intersections( P ):
    """ Call intersectBed to build a single large BED for all intersections between SNPs and affy """
    
    if os.path.exists( P["fn_combined_bed"] ):
        print( "MESSAGE: Combined BED file indicating SNP intersections already exists: " + P["fn_combined_bed"] )
        print( "MESSAGE: Skipping BEDtools step. To re-create this file, delete it and re-run the script." )
        return
    
    cmd = ""
    for i,fn in enumerate(P["fn_vcf"]):
        cmd += 'echo "MESSAGE: processing ' + fn + '"\n'
        cmd += P["dir_bedtools"] + "/intersectBed -u -a " + P["fn_affy_bed"] + " -b " + fn + " > " + P["dir_out"] + "/snps_from_vcf_" + str(i+1) + ".bed\n"
    
    cmd += "cat" 
    for i,fn in enumerate(P["fn_vcf"]):
        cmd += " " + P["dir_out"] + "/snps_from_vcf_" + str(i+1) + ".bed"
    cmd += " | " + P["dir_bedtools"] + "/sortBed > " + P["fn_combined_bed"] + "\n"

    fo = open( P["dir_out"] + "/combine_vcf.sh", 'w')
    fo.write(cmd)    
    fo.close()

    print( "MESSAGE: Executing BEDtools commands to combine VCF files (combine_vcf.sh)" )
    print( "MESSAGE: Please be patient, as this can take several minutes...\n" )
    subprocess.call(["sh", P["dir_out"] + "/combine_vcf.sh"])


class Probeset():
    def __init__(self, probe_id, probe_type, probe_name):
        self.probe_id = probe_id
        self.probe_type = probe_type
        self.probe_name = probe_name        
        self.probe_atoms = []
        self.at_least_one_valid_probe = True
    
    def __str__( self ):
        if self.probe_type.count("control")==0 and len(self.probe_atoms)==0:
            s = ""
        else:
            s = self.probe_id + '\t' + self.probe_type + '\t' + self.probe_name + '\n'
            for atom in self.probe_atoms:
                s += str(atom)
                for probe in atom.probes:
                    s += str(probe)
        
        return s
    
    def remove_probes_by_idx(self, idx):
        tmp = []
        for i,atom in enumerate(self.probe_atoms):
            if not i in idx:
                tmp.append( atom )
        if len(tmp)==0:
            self.at_least_one_valid_probe=False
        
        self.probe_atoms = tmp
        return self.at_least_one_valid_probe
    
        
    def n_probes(self):
        n = 0
        for atom in self.probe_atoms:
            n += len(atom.probes)
        
        return n

class Probe_atom():
    """ one atom can have multiple probes in the case of controls probes, e.g. 10338002 """
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
    """ Read in the PGF file and store it as a vector of probeset objects """
    
    print( "MESSAGE: Parsing affymetrix PGF file..." )
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


def check_sequence_chr_consistency( fn_vcf, fn_affy_bed ):
    """ Confirm that all VCF files and the affymetrix BED file agree about 
        whether column one starts with the string 'chr' or not
    """
    
    bed_has_chr = False
    for line in open(fn_affy_bed):
        if line[0] != "#":
            bed_has_chr = "chr" in line.split('\t')[0]
            break
    
    print( "MESSAGE: Does Affymetrix bed file use 'chr' format? " + str(bed_has_chr) )
    agrees = True
    for fn in fn_vcf:
        f = open(fn)
        for line in f:
            if line[0] != "#":
                if ("chr" in line.split('\t')[0]) != bed_has_chr:
                    agrees = False
                break
    
    return agrees


def read_probes_with_SNPs(P):
    """ identify probesets containing a SNP (BED file)
        load these into a remove hash. We need both the probeset
        identifier and the probe start location to uniquely identify a probe, 
        since multiple probes are identified by the same probeset ID in the 
        Affymetrix probe BED file.
        
        RETURNS:
        ps_start_remove, hash of probeset<space>genomic_start_location
    """
    print( "MESSAGE: Reading probesets with SNPs file:" )
    print( "         " + P["fn_probesets_with_snps"] )
    f = open(P["fn_probesets_with_snps"])
    ps_start_remove = {}
    for line in f:
        a = line.rstrip('\r\n').split('\t')
        chrom, start, stop, probeset = a[0:4]
        ps_start_remove[ probeset + ' ' + start ] = 1

    f.close()
    print( "MESSAGE: VCF file analysis identified " + str(len(ps_start_remove.keys())) + " probes containing SNPs." )
    return ps_start_remove


def parse_probe_BED(P, ps_start_remove):
    """ spin through probe BED, identifying the index of probes within a 
        probeset that are within the remove hash. The probes in the PGF file
        within a probeset can be identified only by index which seems to match 
        the index in which they are present into the Affymetrix BED file.
        
        RETURNS:
        probe_indexes_to_remove is a hash of probeset->list of indexes to remove
    """
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
    print( "MESSAGE: Found " + str(marked_probes_found_in_BED) + " of these probes in the affymetrix BED file" )

    return probe_indexes_to_remove
    
def process_gene_format(P):
    """ Generate new files for Affymetrix Gene ST format 
    """
    
    ps_start_remove = read_probes_with_SNPs(P)
    probe_indexes_to_remove = parse_probe_BED(P, ps_start_remove)
    probesets = parse_PGF(P["fn_pgf_orig"])

    # find the index of each probeset
    ids = [x.probe_id for x in probesets]
    id2idx = {}
    for idx, id in enumerate(ids):
        id2idx[id] = idx

    # How many probes does each probeset have before we touch it?
    ps2before = {}
    n_probes_total = 0
    for id in ids:
        ps2before[id] = probesets[id2idx[id]].n_probes()
        n_probes_total += ps2before[id]
    print( "MESSAGE: Counted " + str(len(probesets)) + " probesets containing " + str(n_probes_total) + " probes in the original PGF file")

    # Remove marked probes from probesets. update the new number of probes in ps2after.
    # ps2after is only set for probesets that we modify
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

    print( "MESSAGE: Removed " + str(n_probes_removed) + " probes for SNP adjustement." )
    print( "MESSAGE: Counted " + str(n_probes_after) + " probes after SNP adjustement." )
    
    # the changes we make are tracked in the pre_post file
    summaries = parse_MPS(P, ps2before, ps2after, probesets_no_valid_probes)
    summaries.write_summary_file( P["fn_probe_count_changes"], "ST" )
        
    timestamp = datetime.datetime.today().strftime("%a %b %d %H:%M:%S PST %Y")
    guid = str(uuid.uuid1())
    
    write_new_PGF_file(P, timestamp, guid, probesets)
    
    # TODO: MPS write needs to know which values changed
    write_new_MPS_file(P, guid, timestamp, summaries)
    write_new_probeset_file(P, summaries)
    write_new_transcript_file(P, summaries)

class Probeset_summary():
    """ Track original and final number of probes for an individual probeset """
    def __init__(self, probeset_id, probes_orig, probes_final):
        self.probeset_id = probeset_id
        self.probes_orig = probes_orig
        self.probes_final = probes_final
    
    def n_probes_final(self):
        return self.probes_final

class Transcript_summary():
    """ Track original and final number or probes for all probesets in transcript cluster """
    def __init__(self, transcript_id, transcript_type):
        self.transcript_id = transcript_id
        self.transcript_type = transcript_type
        self.probeset_summaries = {}
        self.probeset_ids = []
    
    def add_probeset_summary(self, probeset_id, n_probes_orig, n_probes_final):
        self.probeset_ids.append( probeset_id )
        self.probeset_summaries[probeset_id] = Probeset_summary(probeset_id, n_probes_orig, n_probes_final)
                        
    def is_control(self):
        return self.transcript_type=="control"
    
    def n_probes_final(self):
        n=0
        for probeset_id in self.probeset_ids:
            n += self.probeset_summaries[probeset_id].probes_final
        return n
    
    def include_in_MPS(self):
        if self.is_control() or self.n_probes_final()>0:
            return "yes"
        else:
            return "no"
    
    def n_probesets(self):
        return len(self.probeset_ids)

    def n_probesets_with_probe_final(self):
        # note we're returning sum(probesets with >0 probes), not sum(probes)
        n=0
        for probeset_id in self.probeset_ids:
            if self.probeset_summaries[probeset_id].probes_final>0:
                n+=1
        return n
            
    def n_probes_original(self):
        n=0
        for probeset_id in self.probeset_ids:
            n += self.probeset_summaries[probeset_id].probes_orig
        return n

 
class Transcript_summaries():
    """ Wrapper class to manage transcript summaries """
    def __init__(self):
        self.ts2summary = {}
        self.transcript_ids = []
    
    def __len__(self):
        return len(self.transcript_ids)
    
    def add_control_summary( self, transcript_id ):
        if not transcript_id in self.ts2summary:
            self.transcript_ids.append( transcript_id )
            self.ts2summary[transcript_id] = Transcript_summary(transcript_id, "control")
            
    def add_probeset_summary( self, transcript_id, probeset_id, n_probes_orig, n_probes_final ):
        if not transcript_id in self.ts2summary:
            self.transcript_ids.append( transcript_id )
            self.ts2summary[transcript_id] = Transcript_summary(transcript_id, "main")
        self.ts2summary[transcript_id].add_probeset_summary( probeset_id, n_probes_orig, n_probes_final ) 
    
    def get_probeset_summary(self, transcript_id, probeset_id):
        if transcript_id in self.ts2summary:
            ts = self.ts2summary[transcript_id]
            if probeset_id in ts.probeset_summaries:
                return ts.probeset_summaries[probeset_id]
        return None
        
    def write_summary_file(self, fn, chiptype):
        if chiptype != "ST" and chiptype != "IVT":
            raise RuntimeError("chiptype must be ST or IVT")
        
        fo = open(fn, "w")
        if chiptype=="ST":
            fo.write( "transcript.id\torig.n.probes\tprobeset.id\torig.probeset.n.probes\tfinal.probeset.n.probes\tinclude.in.MPS\n" )
            for transcript_id in self.transcript_ids:
                ts = self.ts2summary[transcript_id]
                if ts.is_control():
                    fo.write( transcript_id + "\t0\t" + transcript_id + "\t0\t0\tyes\n")
                else:
                    for probeset_id in ts.probeset_ids:
                        n_orig = str(ts.probeset_summaries[probeset_id].probes_orig)
                        n_final = str(ts.probeset_summaries[probeset_id].probes_final)
                        include_in_MPS = ts.include_in_MPS()
                        fo.write( transcript_id + "\t" + str(ts.n_probes_original()) + "\t" + probeset_id + "\t" + n_orig + '\t' + n_final + '\t' + include_in_MPS + '\n')
        elif chiptype=="IVT":
            fo.write( "probeset.id\torig.probeset.n.probes\tfinal.probeset.n.probes\n" )
            for transcript_id in self.transcript_ids:
                ts = self.ts2summary[transcript_id]
                if ts.is_control():
                    fo.write( transcript_id + "\t0\t0\n")
                else:
                    for probeset_id in ts.probeset_ids:
                        n_orig = str(ts.probeset_summaries[probeset_id].probes_orig)
                        n_final = str(ts.probeset_summaries[probeset_id].probes_final)
                        fo.write( probeset_id + "\t" + n_orig + '\t' + n_final + '\n')
        
        fo.close()
    
    def n_control_transcripts(self):
        n=0
        for transcript_id in self.transcript_ids:
            if self.ts2summary[transcript_id].is_control():
                n += 1
        return n
    
    def n_probesets(self):
        n=0
        for transcript_id in self.transcript_ids:
            n += self.ts2summary[transcript_id].n_probesets()  
        return n
    
    def n_probesets_with_a_probe_final(self):
        n=0
        for transcript_id in self.transcript_ids:
            n += self.ts2summary[transcript_id].n_probesets_with_probe_final()  
        return n

def parse_MPS(P, ps2before, ps2after, probesets_no_valid_probes):
    """ Use the stock MPS file to map probesets to transcript IDs and get probe counts
        RETURN:
        probe_statistics object
    """
    print( "MESSAGE: Parsing Affymetrix MPS file. This file assigns probesets to transcripts." )
    f = open(P["fn_mps_orig"])
    n_empty_probesets = 0
    n_mps_probesets_at_least_one_probe = 0
    n_mps_probesets_total = 0
    n_control_probesets = 0
    
    sums = Transcript_summaries()
    
    for line in f:
        if line[0]!="#":
            ps, ts, probes, length = line.rstrip('\r\n').split('\t')            
            if ts=="":
                ts = ps    
            if length=="":
                # these are control probes
                sums.add_control_summary(ps)
        
            elif ps != "probeset_id":
                for probe in probes.split(" "):
                    n_before = ps2before[probe]
                    if probe in ps2after:
                        n_after = ps2after[probe]
                    else:
                        n_after = n_before
                    sums.add_probeset_summary( ps, probe, n_before, n_after )
    f.close()
    
    print( "MESSAGE: The Affymetrix MPS file contained " + str(len(sums)) + " transcripts." )
    print( "MESSAGE: The Affymetrix MPS file contained " + str(sums.n_control_transcripts()) + " control transcripts" )
    print( "MESSAGE: The Affymetrix MPS file contained " + str( sums.n_probesets() ) + " probesets" )
    print( "MESSAGE: The Affymetrix MPS file contained " + str( sums.n_probesets_with_a_probe_final() ) + " probesets with >= 1 probe after correction" )
    print( "MESSAGE: Wrote a table of changes to probeset composition to " + P["fn_probe_count_changes"] )
    
    return sums


def write_new_PGF_file(P, timestamp, guid, probesets):
    # Write new PGF file
    # probesets that contain no probes will NOT be written because the overridden str()
    # command knows to return an empty string 
    fo = open(P["fn_pgf_new"], 'w')
    fo.write("#%chip_type=" + P["new_package_name"] + "\n")
    fo.write("#%lib_set_name=" + P["new_package_name"] + "\n")
    fo.write("#%lib_set_version=" + P["annot_revision"] + "\n")
    fo.write("#%create_date=" + timestamp + "\n")
    fo.write("#%guid=" + guid + "\n")
    fo.write("#%pgf_format_version=1.0\n")
    fo.write("#%rows=1190\n")
    fo.write("#%cols=990\n")
    fo.write("#%probesets=" + str(len(probesets)) + "\n")
    fo.write("#%datalines=2023948\n") # TODO: update this correctly
    fo.write("#%sequential=1\n")
    fo.write("#%order=col_major\n")
    fo.write("#%header0=probeset_id	type	probeset_name\n")
    fo.write("#%header1=	atom_id\n")
    fo.write("#%header2=		probe_id	type	gc_count	probe_length	interrogation_position	probe_sequence" + '\n')
    for probeset in probesets:
        fo.write( str(probeset) )
    fo.close()
    print( "MESSAGE: wrote new PGF file to " + P["fn_pgf_new"] )


def write_new_MPS_file(P, guid, timestamp, summaries ):
    """ Rewrite the MPS file
        The MPS links individual probesets to a transcript
    
        In the MPS format the first two columns are identical.
        The third is space-delimited list of probesets.
        The fourth is count of probes.
        If a probeset contains zero probes, we will remove it.
        If the number of probes has changed, we'll change sixth column in transcript.csv 
        to update with new number of probes
    
        Read probes for exclusion out of probe_count_changes
        For each transcript we want 
          the number of probes remaining
          the list of probesets with at least one probe
    
        Apparently the case of no probesets for a transcript is
        ps<tab>{empty}<tab>ps<tab>{empty}<newline>
        We will only write this for cases where the transcript contains no probesets but 
        not because we didn't remove them
    """
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
                fo.write(line.rstrip('\r\n') + '\r') # header
            else:
                transcript_summary = summaries.ts2summary[ps]
                if transcript_summary.include_in_MPS:
                    if transcript_summary.is_control():
                        fo.write( ps + '\t\t' + ps + '\t\n') # control probe
                    else:
                        
                        if transcript_summary.n_probesets_with_probe_final()>0:
                            length = str( transcript_summary.n_probes_final() )
                            probes_not_zero = []
                            for probeset_id in transcript_summary.probeset_ids:
                                if transcript_summary.probeset_summaries[probeset_id].n_probes_final() > 0:
                                    probes_not_zero.append(probeset_id)
                            length = str(transcript_summary.n_probes_final())
                            fo.write( ps + '\t' + ps + '\t' + ' '.join(probes_not_zero) + '\t' + length + '\n' )
    fo.close()
    f.close()

    print( "MESSAGE: wrote new MPS file to " + P["fn_mps_new"] )


def write_new_probeset_file(P, summaries):
    # If all probes have been remove from the probeset, do not write it
    f = open(P["fn_ps_orig"])
    fo = open(P["fn_ps_new"], "w")
    for line in f:
        if line[0] =="#":
            fo.write(line.rstrip('\r\n') + '\n' ) # header
        else:
            a = line.rstrip('\r\n')[1:-1].split('","')
            if a[0]=="probeset_id":
                fo.write(line.rstrip('\r\n') + '\n' ) # header
            else:
                probeset_id = a[0]
                transcript_id = a[6]
                probeset_summary = summaries.get_probeset_summary(transcript_id, probeset_id)
                if not probeset_summary is None and probeset_summary.probes_final>0:
                    a[5] = str( probeset_summary.probes_final ) # rewrite number of probes
                    fo.write( '"' + '","'.join(a) + '"\n' )
    fo.close()
    f.close()
    print( "MESSAGE: wrote new probeset description file to " + P["fn_ps_new"] )


def write_new_transcript_file(P, summaries):
    # If all probesets have been remove from the transcript cluster, do not write it
    f = open(P["fn_ts_orig"])
    fo = open(P["fn_ts_new"], "w")
    for line in f:
        if line[0] =="#": 
            fo.write(line.rstrip('\r\n') + '\n' ) # header
        else:
            a = line.rstrip('\r\n')[1:-1].split('","')
            if a[0]=="transcript_cluster_id":
                fo.write(line.rstrip('\r\n') + '\n' ) # header
            else:
                transcript_id = a[0]
                if summaries.ts2summary[ transcript_id ].n_probesets_with_probe_final()>0:
                    a[6] = str( summaries.ts2summary[ transcript_id ].n_probes_final() )
                    out = '"' + '","'.join(a) + '"'
                    fo.write(out + '\n' )    
    fo.close()
    f.close()
    print( "MESSAGE: wrote new transcript description file to " + P["fn_ts_new"] )


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
    print( "MESSAGE: Reading probesets with SNPs file:" )
    print( "         " + P["fn_probesets_with_snps"] )
    f = open(P["fn_probesets_with_snps"])
    probeset_index_bearing_snp = {}
    ps_remove = {}
    for line in f:
        a = line.rstrip('\r\n').split('\t')
        chrom, start, stop, probeset_plus_index = a[0:4]
        probeset_index_bearing_snp[probeset_plus_index]=1
        try:
            probeset, probe_index = probeset_plus_index.split(":")
            if probeset in ps_remove:
                ps_remove[ probeset ].append(probe_index)
            else:
                ps_remove[ probeset ] = [probe_index]
        except ValueError:
            # certain probe_bed tracks may intermingle individual probe and probe 
            # consensus lines; identify individual probes by presence of colon
            pass
    
    f.close()
    # restrict ps_remove to unique entries; multiple SNPs result in duplicates
    for probeset in ps_remove.keys():
        ps_remove[probeset] = sorted(list(set(ps_remove[probeset])))
    
    print( "MESSAGE: VCF file analysis identified " + str(len(ps_remove.keys())) + " probes containing SNPs." )
    print( "MESSAGE: Reading CDF file:" )
    print( "         " + P["fn_cdf_orig"] )
    
    # Count the number of probes in each probeset
    f = open(P["fn_probetab_orig"])
    h = f.readline()
    probeset2NXY = {}
    current_probeset = ""
    current_N = 0
    for line in f:
        a = line.rstrip('\r\n').split('\t')
        probeset,x,y = a[0], a[1], a[2]
        if probeset != current_probeset:
            current_N=0
            current_probeset = probeset
        else:
            current_N = current_N + 1
        try:
            probeset2NXY[ probeset ].append( [str(current_N), x, y] )
        except:
            probeset2NXY[ probeset ] = [ [str(current_N), x, y] ]
    f.close()
    
    # Write probe count changes
    f = open( P["fn_probe_count_changes"], 'w' )
    f.write( 'probeset\tprobe_index\tX\tY\tprobe_removed\n')
    for probeset in probeset2NXY.keys():
        for n,x,y in probeset2NXY[probeset]:
            probe_removed = "FALSE"
            if probeset + ":" + n in probeset_index_bearing_snp:
                probe_removed = "TRUE"
            f.write( '\t'.join([probeset, n, x, y, probe_removed] ) + '\n')
    f.close()
    
    # spin through cdf, removing atoms which are overlapped. Store probeset_id\tX\tY
    # Then spin through Mouse430_2.probe_tab, removing overlaps
    XY_removed=remove_probes_from_CDF(ps_remove, P["fn_cdf_orig"], P["fn_cdf_new"])
    update_probetab(XY_removed, P["fn_probetab_orig"], P["fn_probetab_new"])
    

#------------------------------------------------------------------------------------------
# Begin main 
#------------------------------------------------------------------------------------------

parser = OptionParser()
parser.add_option("-f", "--format", dest="chip_format", help="Chip format, one of {gene,IVT}", default="")
parser.add_option("-p", "--packagename", dest="newpackage", help="New package name", default="")
parser.add_option("-v", "--vcf", dest="fn_list_vcf", help="Path to VCF files, comma-delimited", default="")
parser.add_option("-a", "--affy", dest="fn_affy_bed", help="Path to affymetrix probe BED file", default="")
parser.add_option("-b", "--bedtools", dest="dir_bedtools", help="Path to bedTools bin folder", default="")

parser.add_option("-s", "--annot_version", dest="annot_version", help="Annotation version", default="")
parser.add_option("-e", "--annot_revision", dest="annot_revision", help="Affy chip library revision", default="")
parser.add_option("-g", "--genome_version_ucsc", dest="genome_version_ucsc", help="UCSC Genome version", default="")
parser.add_option("-d", "--genome_version_ncbi", dest="genome_version_ncbi", help="NCI Genome version", default="")
parser.add_option("-c", "--chip", dest="chip", help="Affy chip", default="")
parser.add_option("-i", "--chip_version", dest="chip_version", help="chip version string", default="")
parser.add_option("-r", "--species", dest="species", help="species name", default="")

parser.add_option("-u", "--author", dest="author", help="Author Name for package", default="")
parser.add_option("-l", "--email", dest="email", help="Author email for package", default="")
parser.add_option("-w", "--organism", dest="organism", help="species name", default="")

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

print( "*** EQUALIZER version " + VERSION )
print( "*** Copyright 2015 David Quigley" )
print( "*** Contact David Quigley (dquigley@cc.ucsf.edu, davidquigley.com)" )
print( "*** Released " + RELEASE_DATE  )

P = validate_input( options ) # Check to make sure requested files exist and commands are coherent

print( "MESSAGE: New package will be named: " + P["new_package_name"] )
print( "MESSAGE: New files will be written to: " + P["dir_out"] )
print( "MESSAGE: Affymetrix probe sequence file: " + P["fn_affy_bed"] )
print( "MESSAGE: User passed the following VCF files to match against probes:" )
for fn in P["fn_vcf"]:
    print( "         " + fn )

print( "MESSAGE: User indicated chip was " + P["chip"] + ", version " + P["annot_revision"] + " annotation version " + P["annot_version"] )

if not check_sequence_chr_consistency( P["fn_vcf"], P["fn_affy_bed"] ):
    print( "ERROR: Affymetrix probe BED file does not have the same chromosome string encoding as one or" )
    print( "       or more of the VCF files. Chromosomes can be encoded as 'chr1' or '1'; the probe_bed " )
    print( "       VCF files must all use the same encoding." )
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
    r_script += "fn.trans = '" + P["fn_ts_new"] + "'\n\n"
    r_script += "seed <- new('AffyGenePDInfoPkgSeed', "
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
    r_script += "probeFile = fn.probes, "
else:
    r_script += "cdfFile = fn.cdf, "
    r_script += "celFile = fn.cel, "
    r_script += "tabSeqFile = fn.tab, "

r_script += "author = '" + P["author"] + "', email = '" + P["email"] + "', biocViews = 'AnnotationData', "
r_script += "genomebuild = '" + P["genome_version_ucsc"].upper() + "', organism = '" + P["organism"] + "', "
r_script += "species = '" + P["species"] + "', url ='')\n"
r_script += "makePdInfoPackage(seed, destDir = '" + P["dir_out"] + "', unlink=T)\n"
pd_name = "pd." + P["new_package_name"].lower().replace("_",".").replace("-",".")
r_script += "#install.packages('" + P["dir_out"] + "/" + pd_name + "', repos = NULL, type='source')"

fo = open(P["fn_R_script"], 'w')
fo.write(r_script)
fo.close()

print( "MESSAGE: Wrote R script to create package: " + P["fn_R_script"] )
print( "MESSAGE: There are several ways to run this script:" )
print( "         1) Open the file, run R, and paste the contents into R" )
print( "         2) From the command line, type: " )
print( "            R < " + P["fn_R_script"] + " --no-save" )
print( " " )
print( "         The script will attempt to install the bioconductor 'pdInfoBuilder' package if " )
print( "         it is not already present on your machine. If you do not have the ability to " )
print( "         install packages on your machine, contact your local system administrator." )
print( "         To install the new package, uncomment and run the last line of the R script." )
print( "         This line reads:" )
print( "         #install.packages('" + P["dir_out"] + "/" + pd_name + "', repos = NULL, type='source')" )
