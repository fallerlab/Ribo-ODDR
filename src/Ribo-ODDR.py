#!/usr/bin/env python3

################################################################################
#
#  ---  Ribo-ODDR - Ribo-seq focused Oligo Design pipeline for ---
#  ---     experiment-specific Depletion of Ribosomal RNAs     ---
#
#  Oligo-design script
#  Version 1.0 (See the ChangeLog.md file for changes.)
#
#  Copyright 2019   Faller Lab <fallerlab@gmail.com>
#                   Ferhat Alkan <f.alkan@nki.nl>
#                   William Faller <w.faller@nki.nl>
#  (See the AUTHORS file for other contributors.)
#
#  This file is part of the Ribo-ODDR pipeline.
#
#  Ribo-ODDR is a free software: you can
#  redistribute it and/or modify it under the terms of the
#  GNU General Public License as published by the
#  Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Ribo-ODDR is distributed in the hope that it will
#  be useful, but WITHOUT ANY WARRANTY; without even
#  the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with the Ribo-ODDR Pipeline, see file LICENSE.
#  If not, see <http://www.gnu.org/licenses/>.
################################################################################

import sys, gzip, pysam, os, argparse, itertools, subprocess, time
from Bio import SeqIO, Seq, SeqRecord
from Bio.Alphabet import generic_dna, generic_rna
from shutil import copyfile

__version__ = "1.0"
__author__ = "Ferhat Alkan <f.alkan@nki.nl>"

## GLOBAL DECLARATIONS ##
# add intermediate files here that should be cleared/deleted at the exit
INT_FILES = set()


# store reads in rRNA specific matrices
ALL_READS_START_END = {}
ALL_rRNA_READ_COUNT = {}
# rRNA sequences dictionary
TX_dic = {}

# Exiting function
def EXIT(ex_type=0,message=""):
    for f in INT_FILES:
        os.remove(f)
    if ex_type==1:
        print("#### Pipeline successfully completed. ####")
        sys.exit(1)
    else:
        sys.stderr.write("ERROR: "+message+"\n")
        sys.exit(ex_type)

# File opener function
def f_open(filename, mode='r'):
    if filename.endswith(".gz"):
        if mode=="r":
            return(gzip.open(filename,"rt"))
        else:
            return(gzip.open(filename, mode))
    else:
        return(open(filename,mode))

# Build the argument parser
def get_parser():
    parser = argparse.ArgumentParser(description='Ribo-ODDR: Oligo Design for Depleting rRNAs in Ribo-seq experiments.')
    parser.add_argument('-r','--rRNA_seqs', metavar='<filepath>', help='<Required> Fasta file for rRNA sequences (Same file used when aligning).', required=True)
    parser.add_argument('-o','--output', metavar='<directory_path>', help='<Required> FULLPATH to output DIRECTORY for designed oligos. EMPTY DIRECTORY MUST BE PASSED TO AVOID OVERWRITING!', required=True)
    parser.add_argument('--RIsearch2_exe', metavar='<exe>', help='<Required> Path to RIsearch2 executable. (Pass "--RIsearch2_exe d" if available by default). It is required for accurate depleting potential computation, cross-species optimization and off-target search.', required=True)
    parser.add_argument('--bowtie2-build_exe', metavar='<exe>', help='<Optional> Path to bowtie2-build executable. (Pass "--bowtie2-build_exe d" if available by default)', required=False)
    parser.add_argument('--tophat_exe', metavar='<exe>', help='<Optional> Path to tophat executable. (Pass "--tophat_exe d"if available by default)', required=False)
    parser.add_argument('--samtools_exe', metavar='<exe>', help='<Optional> Path to samtools executable. (Pass "--samtools_exe" with no argument if available by default)', required=False)
    parser.add_argument('--RNAfold_exe', metavar='<exe>', help='<Optional> Path to RNAfold executable. (Pass "--RNAfold_exe d" if available by default)', required=False)
    parser.add_argument('--OFFtargets', metavar='<filepath>', help='Path to protein-coding transcripts fasta file for potential off-target search.', required=False)
    parser.add_argument('--OFFtargetsSUF', metavar='<filepath>', help='Alternatively, RIsearch2 index file for pc_transcript sequences.', required=False)
    design = parser.add_argument_group("Novel oligo design mode","Design new oligos for rRNA depletion based on given Ribo-seq data.")
    design.add_argument('-s','--samples', metavar='<multiple labels>', nargs='+', help='<Required> Sample names. These labels will also be used to access fastq or bam files. (Prefix+Label+Suffix)')
    design.add_argument('-fp','--fastq_prefix', metavar='<prefix>', help='FULL prefix for trimmed Ribo-seq reads (fastq files). Do not use "./" or "../".')
    design.add_argument('-fs','--fastq_suffix', metavar='<suffix>', help='Suffix for trimmed Ribo-seq reads (fastq files).')
    design.add_argument('-ap','--alignment_prefix', metavar='<prefix>', help='FULL prefix for alignment files (Ribo-seq reads aligned to rRNAs). Do not use "./" or "../".')
    design.add_argument('-as','--alignment_suffix', metavar='<suffix>', help='Suffix for alignment files (Ribo-seq reads aligned to rRNAs)')
    design.add_argument('-l','--range', metavar='<int-int>', default='25-30', help='Length range of oligos given in "min:max" format. Default:"25-30"')
    design.add_argument('-f','--force_design', type=int, metavar='<int>', default=20, help='For every rRNA, report at least given number of oligos at every valid length, even if they don'+"'"+'t pass the filter. Default:"20"')
    optimize = parser.add_argument_group("Cross-species optimization mode","Optimize source oligos for depleting given rRNAs.")
    optimize.add_argument('-op','--oligos2optimize', metavar='<filepath>', help="Fasta file for oligos that you want to optimize. (Overrides the design mode)")
    advanced = parser.add_argument_group("Advanced","Advanced options.")
    advanced.add_argument('-pt','--rrna_perc_threshold', type=float, metavar='<float>', default=5.0, help="Filtering threshold (rRNA depletion percentage), applied within oligo score calculation. (default=5)")
    advanced.add_argument('-st','--score_threshold', type=float, metavar='<float>',  default=0.25, help="Score threshold. Minimum ratio of samples that the oligo can deplete at least the selected percentage of rRNAs within. (default:0.25)")
    # advanced.add_argument('--no_OFF_search', action='store_true', help='Do not look for offtargets!')
    # parser.add_argument('--DNA_DNA_params', help='<Required> Path to DNA-DNA energy parameters. (misc/dna_mathews2004.par file from ViennaRNA package installation folder)', required=True)
    return parser

# Parse fasta, return dic, make sure it's dna sequence
def get_seq(fa_file):
    txdic = {}
    with f_open(fa_file) as handle:
        for record in SeqIO.parse(handle, "fasta") :
            txdic[record.id] = str(record.seq).replace("U","T").replace("u","t")
    return(txdic)

# Parse BAM alignment file, return positional depths (number of reads for every position)
def get_reads(bam_file):
    all_TX_reads = {}
    read_count = 0
    parsed_count = 0

    samfile = pysam.AlignmentFile(bam_file, "rb")
    for TX_id, TX_seq in TX_dic.items():
        TX_l = len(TX_seq)
        TX_reads = [[0]*(TX_l+1) for i in range(TX_l)]

        n=0
        for read in samfile.fetch(TX_id):
            if parsed_count%1000000 == 0:
                sys.stdout.write(str(parsed_count)+" reads read \r")
            parsed_count += 1
            # only count primary aligned reads that are on plus strand
            # Not counting secondary aligned reads let us calculate depleting potential accurately (avoiding multiple counts of same read)
            if read.is_reverse==False and read.is_secondary==False:
                read_count += 1
                TX_reads[read.reference_start][read.reference_end] += 1
                n += 1
        print('#',n,' primary-aligned reads on',TX_id,'for sample',bam_file)
        all_TX_reads[TX_id] = TX_reads
    return(all_TX_reads,read_count)

# oligo class, hold information on designed oligos
class Oligo:
    def __init__(self, id, TX_id, i, oligo_target, score, avg_target_depth, dep_pot):
        self.id = id
        self.targetTX = TX_id
        self.startPos = i
        self.l = len(oligo_target)
        self.seq =  str(Seq.Seq(oligo_target, generic_dna).reverse_complement()).replace('T','U')
        self.score = score
        self.avg_target_depth = avg_target_depth
        self.avg_pot = sum([x[1] for x in dep_pot.values()])
        self.dep_pot = dep_pot
        self.GC = sum([1 for c in self.seq if c in ['G','C','g','c']]) / self.l
        self.offs = set()
        self.MFE = "NA"
        self.structure = "NA"
        self.BPper = "NA"
        self.Emin = "NA"
    def add_fold_info(self, strc, en):
        self.MFE = en
        self.structure = strc
        self.BPper = sum([c=='(' for c in strc]) / self.l
    def add_self_inter(self, str_en):
        self.Emin = float(str_en)
    def add_offs(self, off_list):
        for line in off_list:
            cols = line.decode('utf-8').rstrip().split('\t')
            # Filter negative strand and ofss with high energy
            if cols[6]=='+' and float(cols[7]) < (self.Emin*0.5):
                self.offs.add(line.decode('utf-8'))


# Iterate over samples and find oligos that can deplete rRNAs in all
def get_oligos_with_reads_start_end(sample_cutoff=0.25, rrna_cutoff=5, oligo_range=range(25,30), forceALL=20, risearch="risearch2.x", rRNAs_suf="./rRNAs.suf", working_dir="./"):
    sample_size = len(ALL_READS_START_END.keys())
    TXspec_oligos = {}
    o_l = 0
    print("# Oligo design starts now")
    print("# First locating suboptimal binding regions")
    try:
        os.mkdir(working_dir+"/RIsearch2_suboptimal_oligobind_results")
    except FileExistsError:
        print("# WARNING: ",working_dir+"/RIsearch2_suboptimal_oligobind_results","already exists")
    except:
        print("# Unexpected error:", sys.exc_info()[0])
        EXIT()
    RIs2_dir = os.path.abspath(working_dir+"/RIsearch2_suboptimal_oligobind_results/")
    for l in oligo_range:
        print("# Running RIsearch2 for length="+str(l)+" oligos")
        done_seqs=set()
        min_aimed_pair = min(10, (l-int(l/3)))
        with open(RIs2_dir+"/all_length_"+str(l)+"_oligos_for_RIsearch2.fa","w") as outf:
            for TX_id,TX_seq in TX_dic.items():
                TX_l = len(TX_seq)
                # iterating over positions, not oligos
                for i in range(0,TX_l+1-l):
                    oligo_seq = str(Seq.Seq(TX_seq[i:i+l], generic_dna).reverse_complement()).replace('T','U')
                    if oligo_seq not in done_seqs:
                        outf.write(">"+oligo_seq+"\n"+oligo_seq+"\n")
                        done_seqs.add(oligo_seq)
        # Run RIsearch allowing GUs
        out,err = subprocess.Popen(risearch+' -q '+os.path.abspath(RIs2_dir+"/all_length_"+str(l)+"_oligos_for_RIsearch2.fa")+" -i "+os.path.abspath(rRNAs_suf)+" -l 0 -e -10 -s "+str(min_aimed_pair)+" -p2",
                                                shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=RIs2_dir).communicate()
        if(err.decode('utf-8')!=""):
            EXIT("# There is a potential error in RIsearch2 run:", err.decode('utf-8'))

    print("# Calculating oligo potentials")
    for l in oligo_range:
        print("## Designing for size", l,"                                    ")
        for TX_id,TX_seq in TX_dic.items():
            print("## Designing for TX", TX_id,"                              ")
            if TX_id not in TXspec_oligos:
                TXspec_oligos[TX_id] = {ol:[] for ol in oligo_range}
            TX_l = len(TX_seq)

            # iterating over positions, not oligos
            for i in range(0,TX_l+1-l):
                if i%100 == 1:
                    print("## Designing for position", i, "("+str(int(100*(i/TX_l)))+"% finished)    ",end='\r')
                deplete_potential_summary = {}
                # calculate potential in each sample
                for sample in ALL_READS_START_END.keys():
                    dep_pot = 0

                    # necessary to avoid multiple counting
                    added_pos_pairs = set()

                    # Heuristic part to compute deplete potential
                    # Depleteable Reads must have at least min(10,l-int(l/3)) nt pairing with designed oligo
                    # Depleteable Reads can have max(10,int(l/3)) unpaired nts with designed oligo
                    max_allowed_unpair = max(10, int(l/3))
                    min_aimed_pair = min(10, (l-int(l/3)))

                    # Oligo is designed perfect-complementary to the position starting at i for a length of l
                    # thats why we allow the read-to-deplete to start leftmost at (i-max_allowed_unpair)
                    # and rightmost at (i+l+1-min_aimed_pair)
                    for r in range(i-max_allowed_unpair, i+l+1-min_aimed_pair): # iterating over positions where reads might start
                        if r>=0: # obvious
                            read_basin = ALL_READS_START_END[sample][TX_id][r]
                            maxdist = len(read_basin)  # it's the TXlength+1 probably
                            prior_unpaired = -min(r-i,0) # to keep track of unpaired nts in reads starting before i
                            if r < maxdist:
                                end_lmp=max((r+min_aimed_pair), (i+min_aimed_pair))-1 # depletable read can end here leftmost
                                end_rmp=min((i+l+max_allowed_unpair-prior_unpaired), maxdist) # depletable read can end here rightmost

                                if end_rmp > end_lmp: # it's the
                                    for pep in range(end_lmp, end_rmp):
                                        if (TX_id, r, pep) not in added_pos_pairs:
                                            added_pos_pairs.add((TX_id, r, pep)) # add the (TX_id, read_SP,read_EP) pair to the set
                                            dep_pot += read_basin[pep] # Add the depletable reads to potential

                    # Add secondary potentials to oligo potentials based on predicted interaction sites
                    # Use the same heuristic approach
                    oligo_seq = str(Seq.Seq(TX_seq[i:i+l], generic_dna).reverse_complement()).replace('T','U')
                    if not os.path.isfile(RIs2_dir+"/risearch_"+oligo_seq+".out.gz"):
                        EXIT("No risearch prediction with previously reported potential error:")
                    with(f_open(RIs2_dir+"/risearch_"+oligo_seq+".out.gz")) as inf:
                        for line in inf: # iterate over prediction lines
                            rt_os, rt_osp, rt_oep, rt_tarTX, rt_tar_sp, rt_tar_ep, rt_tar_str, rt_energy, rt_pair = line.rstrip().split("\t")
                            if rt_tar_str=="+": # for every "+" strand interaction
                                out_unp = int(rt_osp)-1
                                out_unp += int(rt_oep)-1 # number of oligo unpaired bases outside the predicted interaction
                                in_p = (int(rt_oep)-int(rt_osp))+1 # number of oligo paired bases inside the predicted interaction

                                # iterating over positions where depletable reads start
                                for r in range((int(rt_tar_sp)-1-max_allowed_unpair)+out_unp, int(rt_tar_sp)+(in_p-min_aimed_pair)):
                                    if r>=0: # obvious
                                        if rt_tarTX not in ALL_READS_START_END[sample]:
                                            EXIT("WEIRD that risearch target is not in parsed reads dic")
                                        read_basin = ALL_READS_START_END[sample][rt_tarTX][r]
                                        maxdist = len(read_basin)  # It's the TXlength+1 probably
                                        prior_unpaired = (-min(r-(int(rt_tar_sp)-1),0)) # to keep track of leftunpaired nts in reads starting before i
                                        if r < maxdist:
                                            end_lmp=(int(rt_tar_sp)-1)+min_aimed_pair # depletable read can end here leftmost
                                            end_rmp=min((int(rt_tar_ep)+1+max_allowed_unpair)-prior_unpaired, maxdist) # depletable read can end here rightmost
                                            if end_rmp>end_lmp:
                                                for pep in range(end_lmp, end_rmp):
                                                    if (rt_tarTX,r,pep) not in added_pos_pairs:
                                                        added_pos_pairs.add((rt_tarTX, r, pep)) # add the (rt_tarTX, read_SP,read_EP) pair to the set
                                                        dep_pot += read_basin[pep] # Add the depletable reads to potential

                    deplete_potential_summary[sample] = (dep_pot, round(100*dep_pot/ALL_rRNA_READ_COUNT[sample],2))
                    #os.remove(working_dir+"/risearch_"+oligo_seq+".out.gz")

                total_score = sum([v[1]>=rrna_cutoff for k,v in deplete_potential_summary.items()])*(1/sample_size)
                o_l+=1
                oligo = Oligo("oligo-"+str(o_l+1), TX_id, i, TX_seq[i:i+l], total_score, None, deplete_potential_summary)
                TXspec_oligos[TX_id][l].append(oligo)

    # Get final oligos based on being able to deplete min 5(rrna_cutoff) percent rRNA in at least 0.25(sample_cutoff) of samples
    final_oligos = {}
    identical_oligos = {}
    for TX_id, l_spec_oligos in TXspec_oligos.items():
        add_l = 0
        size_spec_add_l =  {ol:0 for ol in oligo_range}
        for ol, oligos in l_spec_oligos.items():
            for oligo in oligos:
                if oligo.score >= sample_cutoff:
                    final_oligos[oligo.id] = oligo
                    add_l += 1
                    size_spec_add_l[ol] += 1
                    # Add it into identical oligos list
                    if oligo.seq not in identical_oligos:
                        identical_oligos[oligo.seq] = set()
                    identical_oligos[oligo.seq].add(oligo.id)

        # if none passes filters, force some designs
        if forceALL>0:
            for ol, oligos in l_spec_oligos.items():
                if size_spec_add_l[ol]<forceALL:
                    for oligo in sorted(oligos, key=lambda x: x.avg_pot, reverse=True):
                        final_oligos[oligo.id] = oligo
                        if oligo.seq not in identical_oligos:
                            identical_oligos[oligo.seq] = set()
                            size_spec_add_l[ol] += 1
                            add_l += 1

                        identical_oligos[oligo.seq].add(oligo.id)

                        if size_spec_add_l[ol] == forceALL:
                            break

    print("## Designed", len(final_oligos), "oligos                           ", )
    return(final_oligos,identical_oligos)

# Main function for the pipeline, Use functional programming when necessary to simplify the code
def main():
    global TX_dic
    global ALL_READS_START_END
    global ALL_rRNA_READ_COUNT

    print("#### Ribo-ODDR: Ribo-seq Oligo Design for Depleting rRNAs ####")
    print("#### PIPELINE STARTED ####\n#### These are the given args ####")
    parser = get_parser()
    args = parser.parse_args()
    for arg in vars(args):
        print ("#", arg, getattr(args, arg))
    print("#### END of the given args ####")
    print("#### RUN started ####")

    output_path = os.path.abspath(args.output)
    rrna_seq_file = os.path.abspath(args.rRNA_seqs)
    o2o_file = os.path.abspath(args.oligos2optimize) if args.oligos2optimize not in [None] else args.oligos2optimize
    OFFtargetsSUF = os.path.abspath(args.OFFtargetsSUF) if args.OFFtargetsSUF not in [None] else args.OFFtargetsSUF
    OFFtargets = os.path.abspath(args.OFFtargets) if args.OFFtargets not in [None] else args.OFFtargets

    risearch_exe = os.path.abspath(args.RIsearch2_exe) if args.RIsearch2_exe not in [None,'d'] else args.RIsearch2_exe
    risearch_exe = "risearch2.x" if str(risearch_exe)=="d" else risearch_exe
    samtools_exe = os.path.abspath(args.samtools_exe) if args.samtools_exe not in [None,'d'] else args.samtools_exe
    samtools_exe = "samtools" if str(samtools_exe)=="d" else samtools_exe
    bowtie2_build_exe = os.path.abspath(args.bowtie2_build_exe) if args.bowtie2_build_exe not in [None,'d'] else args.bowtie2_build_exe
    bowtie2_build_exe = "bowtie2-build" if str(bowtie2_build_exe)=="d" else bowtie2_build_exe
    tophat_exe = os.path.abspath(args.tophat_exe) if args.tophat_exe not in [None,'d'] else args.tophat_exe
    tophat_exe = "tophat" if str(tophat_exe)=="d" else tophat_exe
    RNAfold_exe = os.path.abspath(args.RNAfold_exe) if args.RNAfold_exe not in [None,'d'] else args.RNAfold_exe
    RNAfold_exe = "RNAfold" if str(RNAfold_exe)=="d" else RNAfold_exe


    # Work in output directory
    if not os.path.isdir(output_path):
        EXIT(message=output_path+" is not a directory.")
    os.chdir(output_path)
    if not os.path.isdir(output_path):
        EXIT(message=output_path+" is not a directory.")
    print("#### WORKING IN "+output_path+" directory ####")


    # Get rRNA sequences
    TX_dic = get_seq(rrna_seq_file)
    for k,v in TX_dic.items():
        print('#rRNA sequence read:',k,len(v),'nts')

    # Create suf file for rRNA seqs. It's necessary for determining the additional secondary potential of oligos
    rRNAs_suf = output_path+"/rRNA_seqs.suf"
    out,err = subprocess.Popen(risearch_exe+' -c '+rrna_seq_file+" -o "+rRNAs_suf, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if not os.path.isfile(rRNAs_suf):
        EXIT(message=err.decode("utf-8"))

    # Get oligo length range
    oligo_range = None
    if len(args.range.split('-')) != 2:
        sys.stderr.write("ERROR: Check given oligo length range !!!\n")
        EXIT()
    else:
        oligo_min, oligo_max = [int(c) for c in args.range.split('-')]
        if oligo_min > oligo_max:
            sys.stderr.write("ERROR: Check given oligo length range. MIN larger than MAX. !!!\n")
            EXIT()
        else:
            oligo_range = range(oligo_min, oligo_max+1)
    print('#Allowed lengths for oligo designs: ',str(oligo_range),' nts')

    # STORE ALL OLIGOS HERE, KEEP TRACK OF IDENTICAL ONES
    oligos = {} # key is oligo id
    identical_oligos = {} # key is oligo sequence

    # GET run mode
    # if o2o file is given, run cross-species optimization
    if o2o_file != None:
        ols2opt = get_seq(o2o_file)
        ol_l_min = min([len(ol) for ol in ols2opt.values()])
        out,err = subprocess.Popen(risearch_exe+' -q '+o2o_file+" -i "+output_path+"/rRNA_seqs.suf -t 1 -s "+str(ol_l_min/2),
                                    shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
        for ol_id, ol_seq in ols2opt.items():
            if not os.path.isfile("risearch_"+ol_id+".out.gz"):
                EXIT(message=err.decode("utf-8"))
            with gzip.open("risearch_"+ol_id+".out.gz") as self_int:
                if not os.path.isfile("risearch_"+ol_id+".out.gz"):
                    EXIT(message=err.decode("utf-8"))
                min_eng = 0
                min_cols = None
                for line in self_int:
                    cols = line.decode('utf-8').rstrip().split('\t')
                    if float(cols[7])<min_eng:
                        min_eng = float(cols[7])
                        min_cols = cols

                spos = (int(min_cols[4])-1)+(int(min_cols[2])-len(ol_seq))
                epos = int(min_cols[5])+int(min_cols[1])-1
                oligo = Oligo("optimized-"+ol_id, min_cols[3], spos, TX_dic[min_cols[3]][spos:epos], 0, None, {"NA":(0,0,0)})
                oligos[oligo.id] = oligo
            INT_FILES.add("risearch_"+ol_id+".out.gz")
            #oligo-12	2	24	ENSMUST00000099123	1228	1250 +	-20.66	PPUPPPUWWWWPPPWPPPPWPPP

    else: # otherwise run DESIGN mode
        deplete_potentials = {}

        if len(args.samples)<1:
            EXIT(message="ERROR: Do not forget to pass sample labels!!!\n")

        alignment_prefix = args.alignment_prefix
        alignment_suffix = args.alignment_suffix
        # Perform alignments if not given
        if alignment_prefix==None and alignment_suffix==None:
            print("### Alignment suffix&prefix is not provided. Will check for fastq files instead.")
            if args.fastq_prefix==None and args.fastq_suffix==None:
                EXIT(message="ERROR: NO fastq suffix and prefix provided.!!!\n")
            else:
                print("### Pilot data will be aligned to rRNAs.")
                print("### Creating rRNA index with bowtie2-build.")
                if bowtie2_build_exe==None or tophat_exe==None or samtools_exe==None:
                    EXIT(message="ERROR: 'bowtie2-build' and 'tophat' and 'samtools' executables are needed.!!!\n")
                copyfile(rrna_seq_file, output_path+"/rRNA_sequences.fa")
                out,err = subprocess.Popen((bowtie2_build_exe+" "+output_path+"/rRNA_sequences.fa "+output_path+"/rRNA_sequences"),
                            shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
                if not os.path.isfile(output_path+"/rRNA_sequences.1.bt2"):
                    EXIT(message=err.decode("utf-8"))
                for sample_label in args.samples:
                    print("### Aligning "+sample_label)
                    os.mkdir(output_path+"/"+sample_label+"_vs_rRNAsequences")
                    command = [tophat_exe, "-p 8 --no-novel-juncs","--no-novel-indels","--no-coverage-search", "--segment-length", "25",
                                    "-o",output_path+"/"+sample_label+"_vs_rRNAsequences/", output_path+"/rRNA_sequences",
                                        ("" if args.fastq_prefix==None else args.fastq_prefix)+sample_label+("" if args.fastq_suffix==None else args.fastq_suffix)]
                    out,err=subprocess.Popen((" ".join(command)), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
                    if not os.path.isfile(output_path+"/"+sample_label+"_vs_rRNAsequences/accepted_hits.bam"):
                        EXIT(message=err.decode("utf-8"))
                    out,err=subprocess.Popen(samtools_exe+" index "+output_path+"/"+sample_label+"_vs_rRNAsequences/accepted_hits.bam",
                                        shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
                    if not os.path.isfile(output_path+"/"+sample_label+"_vs_rRNAsequences/accepted_hits.bam.bai"):
                        EXIT(message=err.decode("utf-8"))
                alignment_prefix = output_path+"/"
                alignment_suffix = "_vs_rRNAsequences/accepted_hits.bam"

        # Read alignment files to get rRNA fragments
        for sample_label in args.samples:
            alignment_file = alignment_prefix + sample_label + alignment_suffix
            print("### Reading "+sample_label+" alignment file and getting depth profile.")

            reads_start_end, read_count = get_reads(alignment_file)
            ALL_READS_START_END[sample_label] = reads_start_end

            # deplete_potentials[sample_label] = dep_pot
            ALL_rRNA_READ_COUNT[sample_label] = read_count
            print(str(ALL_rRNA_READ_COUNT[sample_label])+ ' rRNA mapping reads')

        # oligos fasta
        oligos, identical_oligos = get_oligos_with_reads_start_end(sample_cutoff=args.score_threshold, rrna_cutoff=args.rrna_perc_threshold,
                                                                    oligo_range=oligo_range, forceALL=args.force_design,
                                                                    risearch=risearch_exe, rRNAs_suf=rRNAs_suf, working_dir=output_path)


    if oligos != None:
        SeqIO.write([SeqRecord.SeqRecord(Seq.Seq(oligo.seq,generic_rna), id=oligo.id, description='|'.join([str(round(oligo.score,2)),oligo.targetTX,str(oligo.startPos)+"-"+str(oligo.startPos+oligo.l)])) for oligo in oligos.values()], output_path+"/oligos.fa", "fasta")

        # oligos order.fa
        with open(output_path+"/oligos_order.fa",'w') as out_f:
            for oligo in oligos.values():
                out_f.write('>'+oligo.id+'\n/5Biosg/'+''.join([('r'+c) for c in oligo.seq])+'\n')

        # RNAfold results
        oligo_folds = {}
        if RNAfold_exe == None:
            sys.stderr.write("# WARNING: RNAfold_exe is not set. Self-folding information will not be reported.\n")
            print("# WARNING: RNAfold_exe is not set. Self-folding information will not be reported.")
        else:
            out, err = subprocess.Popen(RNAfold_exe +' --noPS < '+output_path+"/oligos.fa",
                                        shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
            if("not found" in out.decode('utf-8') or "not found" in err.decode('utf-8')):
                EXIT(message=out.decode('utf-8')+"|"+err.decode('utf-8')+" within RNAfold subprocess")
            fold_res = [i.decode('utf-8') for i in out.strip().splitlines()]
            # if len(fold_res)%3 != 0:
            #    print("ERROR",len(fold_res))
            for i in range(int(len(fold_res)/3)):
                oligos[fold_res[3*i][1:].split()[0]].add_fold_info(fold_res[(3*i)+2].split(None,1)[0].strip(), float(fold_res[(3*i)+2].rsplit("(",1)[1].strip("()").strip()))

        # Adding the Binding ENergy and off-target data/score on oligos
        if risearch_exe == None:
            sys.stderr.write("# WARNING: RIsearch_exe is not set. Binding energy information and OFF-TARGETS will not be reported.\n")
            print("# WARNING: RIsearch_exe is not set. Binding energy information and OFF-TARGETS will not be reported.")
        else:
            self_eng = {}
            out,err = subprocess.Popen(risearch_exe+' -c '+output_path+"/oligos.fa -o "+output_path+"/oligos.suf", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
            out,err = subprocess.Popen(risearch_exe+' -q '+output_path+"/oligos.fa -i "+output_path+"/oligos.suf -e 0 -t 1 -s "+str(min(oligo_range)),
                                shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
            print("# Analyzing Emins")
            min_Emin = 0
            for oligo_id in oligos.keys():
                if not os.path.isfile("risearch_"+oligo_id+".out.gz"):
                    EXIT(message=err.decode("utf-8"))
                # Get Emin first. Energy of perfect complementary ineraction.
                with gzip.open("risearch_"+oligo_id+".out.gz") as self_int:
                    for line in self_int:
                        cols = line.decode('utf-8').rstrip().split('\t')
                        if cols[0]==cols[3] and cols[1]==cols[4] and cols[2]==cols[5] and cols[6]=='-':
                            oligos[oligo_id].add_self_inter(cols[7])
                            # update min Emin
                            min_Emin = oligos[oligo_id].Emin if oligos[oligo_id].Emin<min_Emin else min_Emin
                            INT_FILES.add("risearch_"+oligo_id+".out.gz")
                            break

            if (OFFtargets==None and OFFtargetsSUF==None):
                sys.stderr.write("# WARNING: Off-target transcripts are not given. OFF-TARGETS will not be reported.\n")
                print("# WARNING: Off-target transcripts are not given. OFF-TARGETS will not be reported.")
            else:
                print("# RUNning off-targets")

                suffile = OFFtargetsSUF
                if suffile==None:
                    if not os.path.isfile(OFFtargets):
                        EXIT(message=OFFtargets+" is not reachable!")
                    elif OFFtargets.endswith("suf"):
                        EXIT(message=OFFtargets+" is not a .fa, .fa.gz, .FASTA, .FASTA.gz file.")
                    out,err = subprocess.Popen(risearch_exe+' -c '+OFFtargets+" -o "+output_path+"/OFFtargets.suf",
                                shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
                    suffile = output_path+"/OFFtargets.suf"
                    if not os.path.isfile(output_path+"/OFFtargets.suf"):
                        EXIT(message=err.decode('utf-8')+"\n Above problem in RIsearch2 off-target index generation")

                # max energy -20, but if no oligo has BE > -40 there is no need to predict with -e -20. that's why there is max option here
                out,err = subprocess.Popen(risearch_exe+' -q '+output_path+"/oligos.fa -t 1 -i "+suffile+
                                            " -s "+str(min(round(0.75*min(oligo_range)), 10))+" -e "+str(max(min_Emin*0.5, -20)),
                            shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
                t_c = 0
                for oligo_id in oligos.keys():
                    if t_c%10 == 0:
                        sys.stdout.write(str(round((100.0*t_c)/len(oligos)))+" percent complete  \r")
                    t_c+=1

                    if not os.path.isfile(output_path+"/risearch_"+oligo_id+".out.gz"):
                        EXIT(message=err.decode('utf-8')+"\n Above problem in RIsearch2 off-target prediction")
                    with gzip.open(output_path+"/risearch_"+oligo_id+".out.gz") as self_int:
                        oligos[oligo_id].add_offs(self_int.readlines())

        # Report oligos in bed file
        with open(output_path+"/oligos.bed",'w') as out_f:
            for oligo in oligos.values():
                out_f.write('\t'.join([oligo.targetTX,str(oligo.startPos),str(oligo.startPos+oligo.l),oligo.id])+'\n')


        # Report oligos in gff file
        with open(output_path+"/oligos.gff3",'w') as out_f:
            out_f.write("##gff-version 3\n")
            for oligo in oligos.values():
                out_f.write('\t'.join([oligo.targetTX,'CUSTOM','oligo',str(oligo.startPos+1),str(oligo.startPos+oligo.l),'.','-','.',
                            ';'.join(['ID='+oligo.id,'score='+str(oligo.score),'seq='+oligo.seq,
                            'gc_content='+str(round(oligo.GC,2)),'mfe='+str(oligo.MFE),'structure='+oligo.structure,
                            'bp_per='+("NA" if oligo.BPper=="NA" else str(round(oligo.BPper,2))),'eng_min='+("NA" if oligo.Emin=="NA" else str(round(oligo.Emin,1))),
                            'no_of_OFFS='+("NA" if (OFFtargets==None and OFFtargetsSUF==None) else str(len(oligo.offs)))]+[(k+'='+str(v[0])+"-reads_"+str(v[1])+"-pc") for k,v in oligo.dep_pot.items()])])+'\n')

        # Report oligos in CSV file
        with open(output_path+"/oligos.csv",'w') as out_f:
            out_f.write("# Ribo-ODDR "+__version__+" oligo output file\n")
            samples = [k for k,v in oligo.dep_pot.items()]

            out_f.write(",".join([ 'oligoID','length','sequence','target','avg_dep_per','score','gc_content','mfe','structure','bp_per','eng_min','no_of_OFFs', 'OFFs']+samples)+"\n")
            for oligo in oligos.values():
                avg_dep_per = round((sum([oligo.dep_pot[k][1] for k in samples]) / float(len(samples))) , 2)
                out_f.write(','.join([oligo.id, str(oligo.l), oligo.seq,
                                        ('"'+oligo.targetTX+':'+str(oligo.startPos+1)+'-'+str(oligo.startPos+oligo.l)+'"'),
                                        str(avg_dep_per), str(oligo.score), str(round(oligo.GC,2)),
                                        str(oligo.MFE), oligo.structure,
                                        ("NA" if oligo.BPper=="NA" else str(round(oligo.BPper,2))),
                                        ("NA" if oligo.Emin=="NA" else str(round(oligo.Emin,1))),
                                        ("None" if (OFFtargets==None and OFFtargetsSUF==None) else str(len(oligo.offs))),
                                        ('"'+';'.join(list(set([off.split("\t")[3].replace('"','_') for off in oligo.offs])))+'"') if len(oligo.offs)<10 else '"More than 10. See oligos_offtargets.txt file for details."']+
                                     [('"'+str(oligo.dep_pot[k][0])+"-reads_"+str(oligo.dep_pot[k][1])+"-pc"+'"') for k in samples])+'\n')

        # Report oligos off-targets in txt file if OFFtargets target seq is given
        if (OFFtargets!=None or OFFtargetsSUF!=None):
            with open(output_path+"/oligos_offtargets.txt",'w') as out_f:
                for oligo in oligos.values():
                    if len(oligo.offs)>0:
                        out_f.write(''.join([off for off in oligo.offs])+'\n')
    else:
        print("## NO OLIGOS FOUND ##")




# Run the main when executed
if __name__ == "__main__":
    main()
    EXIT(1)
