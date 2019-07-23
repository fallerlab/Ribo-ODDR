#!/usr/bin/env python3

##################################################
#
#  --- Ribo-ODDR - Oligo Design for Depleting rRNAs ---
#
#  Version 0.1 (See the ChangeLog.md file for changes.)
#
#  Copyright 2019   Ferhat Alkan <f.alkan@nki.nl>
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
#################################################

import sys, gzip, pysam, os, argparse, itertools, subprocess
from Bio import SeqIO, Seq, SeqRecord
from Bio.Alphabet import generic_dna, generic_rna

__version__ = "0.1"
__author__ = "Ferhat Alkan <f.alkan@nki.nl>"

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
    parser.add_argument('-l','--range', default='25-34', help='Length range of oligos given in min:max format. Default:"25-35"')
    parser.add_argument('--RNAfold_exe', help='<Required> Path to RNAfold executable. (Default:RNAfold)', required=True)
    parser.add_argument('--RIsearch_exe', help='<Required> Path to RIsearch2 executable. (Default:risearch2.x)', required=True)
    parser.add_argument('--OFFtargets', help='<Required> Path to suffix file for potential off-targets.', required=True)
    #parser.add_argument('--DNA_DNA_params', help='<Required> Path to DNA-DNA energy parameters. (misc/dna_mathews2004.par file from ViennaRNA package installation folder)', required=True)
    return parser

# Parse fasta, return dic
def get_seq(fa_file):
    TX_dic = {}
    with f_open(fa_file) as handle:
        for record in SeqIO.parse(handle, "fasta") :
            TX_dic[record.id] = str(record.seq)
    return(TX_dic)

# Parse bam, return positional depths
def get_positional_depths(bam_file, TX_dic, total_mapped, oligorange=range(25,35)):
    TX_depths = {}
    deplete_potential = {}
    read_count = 0
    parsed_count = 0

    samfile = pysam.AlignmentFile(bam_file, "rb")
    for TX_id, TX_seq in TX_dic.items():
        TX_l = len(TX_seq)
        TX_depths[TX_id] = [0]*TX_l
        deplete_potential[TX_id] = [[0]*TX_l]*len(oligorange)
        n=0
        for read in samfile.fetch(TX_id):
            if parsed_count%1000000 == 0:
                sys.stdout.write(str(round((100.0*parsed_count)/total_mapped))+" percent complete  \r")
            parsed_count += 1
            if not read.is_reverse:
                read_count += 1
                for i in range (read.reference_start, read.reference_end):
                    TX_depths[TX_id][i] += 1

                for o_l in oligorange:
                   for i in range(read.reference_start-round(0.1*o_l), read.reference_end-round(0.9*o_l)):
                       if i>=0:
                           deplete_potential[TX_id][o_l-oligorange[0]][i] += 1

                n += 1
        print('#',n,round(n/total_mapped,2),'reads on',TX_id,'for sample',bam_file)
    return(TX_depths,deplete_potential,read_count)

# Parse bam, return positional depths
def get_reads(bam_file, TX_dic, total_mapped):
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
                sys.stdout.write(str(round((100.0*parsed_count)/total_mapped))+" percent complete  \r")
            parsed_count += 1
            if not read.is_reverse:
                read_count += 1
                TX_reads[read.reference_start][read.reference_end] += 1
                n += 1
        print('#',n,round(n/total_mapped,2),'reads on',TX_id,'for sample',bam_file)
        all_TX_reads[TX_id] = TX_reads
    return(all_TX_reads,read_count)

# oligo class
class Oligo:
    def __init__(self, id, TX_id, i, oligo_target, score, avg_target_depth, dep_pot):
        self.id = id
        self.targetTX = TX_id
        self.startPos = i
        self.l = len(oligo_target)
        self.seq =  str(Seq.Seq(oligo_target, generic_dna).reverse_complement()).replace('T','U')
        self.score = score
        self.avg_target_depth = avg_target_depth
        self.dep_pot = dep_pot
        self.GC = sum([1 for c in self.seq if c in ['G','C','g','c']]) / self.l
        self.offs = set()
    def add_fold_info(self, strc, en):
        self.MFE = en
        self.structure = strc
        self.BPper = sum([c=='(' for c in strc]) / self.l
    def add_self_inter(self, str_en):
        self.Emin = float(str_en)
    def add_offs(self, off_list):
        for line in off_list:
            cols = line.decode('utf-8').rstrip().split('\t')
            if cols[6]=='+':
                self.offs.add(line.decode('utf-8'))

# Iterate over samples and find oligos that can deplete rRNAs in all
def get_oligos(all_rRNA_depths, dep_pots, all_rRNA_read_count, all_total_read_count, TX_dic, total_cutoff=5, rrna_cutoff=10, oligo_range=range(25,35)):

    sample_size = len(all_rRNA_depths.keys())

    over_threshold = {sample:{TX[0]:[False]*len(TX[1]) for TX in TX_dic.items()} for sample in all_rRNA_depths.keys()}
    pass_score = {TX[0]:[0]*len(TX[1]) for TX in TX_dic.items()}
    for sample,rRNA_depths in all_rRNA_depths.items():
        for TX_id, rRNA_depth in rRNA_depths.items():
            per_profile = [((float(x)/float(all_total_read_count[sample])), ((float(x)/float(all_rRNA_read_count[sample])))) for x in rRNA_depth]
            for i in range(len(per_profile)):
                totper,rrnaper = per_profile[i]
                if totper > total_cutoff:
                    pass_score[TX_id][i] += (0.5 / sample_size)
                    over_threshold[sample][TX_id][i] = True
                if rrnaper > rrna_cutoff:
                    pass_score[TX_id][i] += (0.5 / sample_size)
                    over_threshold[sample][TX_id][i] = True

    oligos = {}
    for l in oligo_range:
        for TX_id,TX_seq in TX_dic.items():
            TX_l = len(TX_seq)
            for i in range(0,TX_l+1-l):
                # total_score = sum(pass_score[TX_id][i:i+l])/l
                # if total_score > 0.5:
                #     avg_target_depths = {}
                #     deplete_potential_summary = {}
                #     for sample in all_rRNA_depths.keys():
                #         average = float(sum(all_rRNA_depths[sample][TX_id][i:i+l]))/l
                #         avg_target_depths[sample] = (round(average), round(100*average/all_rRNA_read_count[sample]), round(100*average/all_total_read_count[sample]))
                #         dep_pot = dep_pots[sample][TX_id][l-oligo_range[0]][i]
                #         deplete_potential_summary[sample] = (dep_pot,round(100*dep_pot/all_rRNA_read_count[sample]),round(100*dep_pot/all_total_read_count[sample]))
                #     oligo = Oligo("oligo-"+str(len(oligos)+1), TX_id, i, TX_seq[i:i+l], total_score, avg_target_depths, deplete_potential_summary)
                #     oligos["oligo-"+str(len(oligos)+1)] = oligo

                avg_target_depths = {}
                deplete_potential_summary = {}
                for sample in dep_pots.keys():
                    average = float(sum(all_rRNA_depths[sample][TX_id][i:i+l]))/l
                    avg_target_depths[sample] = (round(average), round(100*average/all_rRNA_read_count[sample]), round(100*average/all_total_read_count[sample]))
                    dep_pot = dep_pots[sample][TX_id][l-oligo_range[0]][i]
                    deplete_potential_summary[sample] = (dep_pot,round(100*dep_pot/all_rRNA_read_count[sample]),round(100*dep_pot/all_total_read_count[sample]))

                total_score = sum([v[1]>=rrna_cutoff for k,v in deplete_potential_summary.items()])*(0.5/sample_size)
                total_score+= sum([v[2]>=total_cutoff for k,v in deplete_potential_summary.items()])*(0.5/sample_size)

                if total_score > 0.5:
                    oligo = Oligo("oligo-"+str(len(oligos)+1), TX_id, i, TX_seq[i:i+l], total_score, avg_target_depths, deplete_potential_summary)
                    oligos["oligo-"+str(len(oligos)+1)] = oligo



# Iterate over samples and find oligos that can deplete rRNAs in all
def get_oligos_with_reads_start_end(all_reads_start_end, all_rRNA_read_count, all_total_read_count, TX_dic, total_cutoff=5, rrna_cutoff=10, oligo_range=range(25,35)):
    sample_size = len(all_reads_start_end.keys())

    oligos = {}
    for l in oligo_range:
        for TX_id,TX_seq in TX_dic.items():
            TX_l = len(TX_seq)
            for i in range(0,TX_l+1-l):
                deplete_potential_summary = {}
                for sample in all_reads_start_end.keys():
                    dep_pot = 0
                    for r in range(i-round(0.1*l), i+round(0.1*l)):
                        if r>=0:
                            dep_pot += sum(all_reads_start_end[sample][TX_id][r][r:])

                    deplete_potential_summary[sample] = (dep_pot,round(100*dep_pot/all_rRNA_read_count[sample]),round(100*dep_pot/all_total_read_count[sample]))

                total_score = sum([v[1]>=rrna_cutoff for k,v in deplete_potential_summary.items()])*(0.5/sample_size)
                total_score+= sum([v[2]>=total_cutoff for k,v in deplete_potential_summary.items()])*(0.5/sample_size)

                if total_score > 0.25:
                    oligo = Oligo("oligo-"+str(len(oligos)+1), TX_id, i, TX_seq[i:i+l], total_score, None, deplete_potential_summary)
                    oligos["oligo-"+str(len(oligos)+1)] = oligo


    return(oligos)

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

    # Get oligo length range
    oligo_range = None
    if len(args.range.split('-')) != 2:
        sys.stderr.write("ERROR: Check given oligo length range !!!\n")
        EXIT()
    else:
        oligo_min, oligo_max = [int(c) for c in args.range.split('-')]
        if oligo_min >= oligo_max:
            sys.stderr.write("ERROR: Check given oligo length range. MIN larger than MAX. !!!\n")
            EXIT()
        else:
            oligo_range = range(oligo_min, oligo_max+1)
    print('#Allowed lengths for oligo designs: ',str(oligo_range),' nts')

    # Get rRNA sequences
    TX_dic = get_seq(args.rRNA_seqs)
    for k,v in TX_dic.items():
        print('#rRNA sequence read:',k,len(v),'nts')

    all_rRNA_depths = {}
    all_reads_start_end = {}

    deplete_potentials = {}
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

        # rRNA_depths, dep_pot, read_count = get_positional_depths('/'.join([alignment_folder,'accepted_hits.bam']), TX_dic, total_mapped)
        # all_rRNA_depths[sample_label] = rRNA_depths

        reads_start_end, read_count = get_reads('/'.join([alignment_folder,'accepted_hits.bam']), TX_dic, total_mapped)
        all_reads_start_end[sample_label] = reads_start_end

        # deplete_potentials[sample_label] = dep_pot
        all_rRNA_read_count[sample_label] = read_count
        print(str(all_total_read_count[sample_label])+' reads; '+str(all_rRNA_read_count[sample_label])+' (%'+str(round((100.0*all_rRNA_read_count[sample_label])/float(all_total_read_count[sample_label])))+') rRNA reads')

    # oligos fasta
    # oligos = get_oligos(all_rRNA_depths, deplete_potentials, all_rRNA_read_count, all_total_read_count, TX_dic)
    oligos = get_oligos_with_reads_start_end(all_reads_start_end, all_rRNA_read_count, all_total_read_count, TX_dic)

    SeqIO.write([SeqRecord.SeqRecord(Seq.Seq(oligo.seq,generic_rna), id=oligo.id, description='|'.join([str(round(oligo.score,2)),oligo.targetTX,str(oligo.startPos)])) for oligo in oligos.values()], args.output+"/oligos.fa", "fasta")

    # oligos order.fa
    with open(args.output+"/oligos_order.fa",'w') as out_f:
        for oligo in oligos.values():
            out_f.write('>'+oligo.id+'\n/5Biosg/'+''.join([('r'+c) for c in oligo.seq])+'\n')

    # RNAfold results
    oligo_folds = {}
    fold_res = [i.decode('utf-8') for i in subprocess.Popen(args.RNAfold_exe+' --noPS < '+args.output+"/oligos.fa", shell=True, stdout=subprocess.PIPE).communicate()[0].strip().splitlines()]
    if len(fold_res)%3 != 0:
        print("ERROR",len(fold_res))
    for i in range(int(len(fold_res)/3)):
        oligos[fold_res[3*i][1:].split()[0]].add_fold_info(fold_res[(3*i)+2].split(None,1)[0].strip(), float(fold_res[(3*i)+2].rsplit("(",1)[1].strip("()").strip()))

    # Adding the off-taget data/score on oligos
    self_eng = {}
    sub = subprocess.Popen(args.RIsearch_exe+' -c '+args.output+"/oligos.fa -o "+args.output+"/oligos.suf", shell=True, stdout=subprocess.PIPE).communicate()
    sub = subprocess.Popen(args.RIsearch_exe+' -q '+args.output+"/oligos.fa -i "+args.output+"/oligos.suf -s "+str(min(oligo_range)), shell=True, stdout=subprocess.PIPE).communicate()
    for oligo_id in oligos.keys():
        # Get Emin first. Energy of perfect complementary ineraction.
        with gzip.open("risearch_"+oligo_id+".out.gz") as self_int:
            for line in self_int:
                cols = line.decode('utf-8').rstrip().split('\t')
                if cols[0]==cols[3] and cols[1]==cols[4] and cols[2]==cols[5] and cols[6]=='-':
                    oligos[oligo_id].add_self_inter(cols[7])
                    INT_FILES.append("risearch_"+oligo_id+".out.gz")
                    print("#Emin for "+oligo_id,cols[7])
                    break
        oligo = oligos[oligo_id]
        SeqIO.write(SeqRecord.SeqRecord(Seq.Seq(oligo.seq,generic_rna), id=oligo.id), args.output+"/"+oligo_id+".fa", "fasta")
        INT_FILES.append(args.output+"/"+oligo_id+".fa")
        sub = subprocess.Popen(args.RIsearch_exe+' -q '+args.output+"/"+oligo_id+".fa -i "+args.OFFtargets+" -s "+str(round(0.75*oligo.l))+" -e "+str(oligo.Emin*0.5), shell=True, stdout=subprocess.PIPE).communicate()
        with gzip.open("risearch_"+oligo_id+".out.gz") as self_int:
            oligos[oligo_id].add_offs(self_int.readlines())


    # oligos bed
    with open(args.output+"/oligos.bed",'w') as out_f:
        for oligo in oligos.values():
            out_f.write('\t'.join([oligo.targetTX,str(oligo.startPos),str(oligo.startPos+oligo.l),oligo.id])+'\n')


    # oligos gff
    with open(args.output+"/oligos.gff3",'w') as out_f:
        for oligo in oligos.values():
            out_f.write('\t'.join([oligo.targetTX,'CUSTOM','oligo',str(oligo.startPos+1),str(oligo.startPos+oligo.l),'.','-','.',
                        ';'.join(['ID='+oligo.id,'score='+str(oligo.score),'seq='+oligo.seq,
                        'GC_content='+str(round(oligo.GC,2)),'MFE='+str(oligo.MFE),'structure='+oligo.structure,
                        'BPper='+str(round(oligo.BPper,2)),'Emin='+str(round(oligo.Emin,1)),'NOofOFFS='+str(len(oligo.offs))]+[(k+'='+str(v[0])+"-reads,"+str(v[1])+"pc-rRNA-mapping,"+str(v[2])+"pc-total") for k,v in oligo.dep_pot.items()])])+'\n')

    # oligos off
    with open(args.output+"/oligos.off",'w') as out_f:
        for oligo in oligos.values():
            out_f.write(''.join([off for off in oligo.offs])+'\n')





# Run the main when executed
if __name__ == "__main__":
    main()
    EXIT(1)
