#!/usr/bin/env python3

import sys, gzip, pysam, os, argparse, itertools
from Bio import SeqIO, Seq

## GLOBAL DECLARATIONS ##
# intermediate files that should be cleared at the exit
INT_FILES = []
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
    parser.add_argument('-a','--alignments', nargs='+', help='<Required> Alignment files (Ribo-seq reads vs rRNAs). Each alignment file should be given as label:folder:file_basename', required=True)
    parser.add_argument('-r','--rRNA_seqs', help='<Required> Fasta file for rRNA sequences (Same file used when aligning).', required=True)
    parser.add_argument('-o','--output', help='<Required> Output file for designed oligos.', required=True)
    return parser

# Parse fasta, return dic
def get_seq(fa_file):
    TX_dic = {}
    with f_open(fa_file) as handle:
        for record in SeqIO.parse(handle, "fasta") :
            TX_dic[record.id] = str(record.seq)
    return(TX_dic)

# Parse bam, return positional depths
def get_positional_depths(bam_file, TX_dic, total_mapped):
    TX_depths = {}
    read_count = 0
    parsed_count = 0
    samfile = pysam.AlignmentFile(bam_file, "rb")
    for TX_id, TX_seq in TX_dic.items():
        TX_depths[TX_id] = [0]*len(TX_seq)
        for read in samfile.fetch(TX_id):
            if parsed_count%1000000 == 0:
                sys.stdout.write(str(round((100.0*parsed_count)/total_mapped))+" percent complete                  \r")
            parsed_count += 1
            if not read.is_reverse:
                read_count += 1
                for i in range (read.reference_start, read.reference_end):
                    TX_depths[TX_id][i] += 1
    return(TX_depths,read_count)

# Iterate over samples and find oligos that can deplete rRNAs in all
def get_oligos(all_rRNA_depths, all_rRNA_read_count, all_total_read_count, TX_dic, total_cutoff=0.05, rrna_cutoff=0.1, length_range=range(20,35)):
    over_threshold = {sample:{TX[0]:[False]*len(TX[1]) for TX in TX_dic.items()} for sample in all_rRNA_depths.keys()}
    for sample,rRNA_depths in all_rRNA_depths.items():
        for TX_id, rRNA_depth in rRNA_depths.items():
            per_profile = [((float(x)/float(all_total_read_count[sample])),((float(x)/float(all_rRNA_read_count[sample])))) for x in rRNA_depth]
            for i in range(len(per_profile)):
                totper,rrnaper = per_profile[i]
                if totper > total_cutoff or rrnaper > rrna_cutoff:
                    over_threshold[sample][TX_id][i] = True
                    print(sample,TX_id,i)


# Main function for the pipeline, Use functional programming when necessary to simplify the code
def main():
    print("#### Ribo-ODDR: Ribo-seq Oligo Design for Depleting rRNAs ####")
    print("#### PIPELINE STARTED ####\n#### These are the given args ####")
    parser = get_parser()
    args = parser.parse_args()
    for arg in vars(args):
        print ("#", arg, getattr(args, arg))
    print("#### END of the given args ####")
    print("#### RUN started ####")

    TX_dic = get_seq(args.rRNA_seqs)

    all_rRNA_depths = {}
    all_rRNA_read_count = {}
    all_total_read_count = {}
    for sample_info in args.alignments:
        sample_label, alignment_folder = sample_info.rstrip().split(':')
        print("### Reading "+sample_label+" alignment file and getting depth profile.")

        total_mapped = 0
        with open(alignment_folder+'/align_summary.txt') as inf:
            lines = inf.readlines()
            all_total_read_count[sample_label] = int(lines[1].rstrip().replace(' ','').split(':')[1])
            total_mapped = float(lines[2].rstrip().replace(' ','').split(':')[1].split('(')[0])

        rRNA_depths, read_count = get_positional_depths('/'.join([alignment_folder,'accepted_hits.bam']), TX_dic, total_mapped)
        all_rRNA_depths[sample_label] = rRNA_depths
        all_rRNA_read_count[sample_label] = read_count
        print(str(all_total_read_count[sample_label])+' reads; '+str(all_rRNA_read_count[sample_label])+' (%'+str(round((100.0*all_rRNA_read_count[sample_label])/float(all_total_read_count[sample_label])))+') rRNA reads')

    oligos = get_oligos(all_rRNA_depths,all_rRNA_read_count, all_total_read_count, TX_dic)


# Run the main when executed
if __name__ == "__main__":
    main()
    EXIT(1)
