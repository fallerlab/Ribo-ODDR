#!/usr/bin/env python3

##################################################
#
#  --- Ribo-ODDR - Oligo Design for Depleting rRNAs ---
#  ---      in Ribosome profiling experiments       ---
#  Version 1.0 (See the ChangeLog.md file for changes.)
#
#  Copyright 2019   Ferhat Alkan <f.alkan@nki.nl>
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
    parser.add_argument('-r','--rRNA_seqs', help='<Required> Fasta file for rRNA sequences (Same file used when aligning).', required=True)
    parser.add_argument('-o','--output', help='<Required> Output directory for designed oligos.', required=True)
    parser.add_argument('--RNAfold_exe', help='<Required> Path to RNAfold executable. (Default:RNAfold)', required=True)
    parser.add_argument('--RIsearch_exe', help='<Required> Path to RIsearch2 executable. (Default:risearch2.x)', required=True)
    parser.add_argument('--OFFtargets', help='<Required> Path to RIsearch2 suffix file for potential off-target search.', required=True)
    design = parser.add_argument_group("Design-mode","Design new oligos for rRNA depletion based on given samples.")
    design.add_argument('-a','--alignments', nargs='+', help='<Required> Alignment files (Ribo-seq reads vs rRNAs). Each alignment file should be given as "label1:folder1 label2:folder2" or just as "label1 label2" with --prefix and --suffix.')
    design.add_argument('-pre','--prefix', help='Prefix for sample rRNA-alignment folder')
    design.add_argument('-suf','--suffix', help='Suffix for sample rRNA-alignment folder')
    design.add_argument('-l','--range', default='25-30', help='Length range of oligos given in "min:max" format. Default:"25-30"')
    design.add_argument('-f','--force_design', type=int, default=20, help='For every given rRNA, force designing given number of oligos at every length in the pre-given range. Default:"10"')
    optimize = parser.add_argument_group("Optimize-mode","Optimize given oligos for depleting given rRNAs.")
    optimize.add_argument('-op','--oligos2optimize', help="Fasta file for oligos that you want to optimize. (Overrides the design mode)")
    advanced = parser.add_argument_group("Advanced","Advanced options.")
    advanced.add_argument('-pt','--rrna_perc_threshold', type=float, default=5.0, help="Filtering threshold (rRNA depletion percentage), applied within oligo score calculation. (default=5)")
    advanced.add_argument('-st','--score_threshold', type=float,  default=0.25, help="Score threshold. Minimum ratio of samples that the oligo can deplete at least the selected percentage of rRNAs within. (default:0.25)")
    advanced.add_argument('--no_OFF_search', action='store_true', help='Do not look for offtargets!')

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
        self.avg_pot = sum([x[1] for x in dep_pot.values()])
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
            if cols[6]=='+' and float(cols[7]) < min(self.Emin*0.5, -25):
                self.offs.add(line.decode('utf-8'))

# Iterate over samples and find oligos that can deplete rRNAs in all
def get_oligos_with_reads_start_end(all_reads_start_end, all_rRNA_read_count, all_total_read_count, TX_dic, sample_cutoff=0.25, rrna_cutoff=5, oligo_range=range(25,30), forceALL=10):
    sample_size = len(all_reads_start_end.keys())

    TXspec_oligos = {}
    o_l = 0
    print("# Oligo design starts now")
    for l in oligo_range:
        print("## Designing for size", l)
        for TX_id,TX_seq in TX_dic.items():
            if TX_id not in TXspec_oligos:
                TXspec_oligos[TX_id] = {ol:[] for ol in oligo_range}
            TX_l = len(TX_seq)
            for i in range(0,TX_l+1-l):
                deplete_potential_summary = {}
                for sample in all_reads_start_end.keys():
                    dep_pot = 0
                    for r in range(i-10, i+10): # Heuristic part
                        if r>=0:
                            dep_pot += sum(all_reads_start_end[sample][TX_id][r][r:])

                    deplete_potential_summary[sample] = (dep_pot, round(100*dep_pot/all_rRNA_read_count[sample],2), round(100*dep_pot/all_total_read_count[sample],2))

                total_score = sum([v[1]>=rrna_cutoff for k,v in deplete_potential_summary.items()])*(1/sample_size)
                o_l+=1
                oligo = Oligo("oligo-"+str(o_l+1), TX_id, i, TX_seq[i:i+l], total_score, None, deplete_potential_summary)
                TXspec_oligos[TX_id][l].append(oligo)

    # Get final oligos based on being able to deplete min 5(rrna_cutoff) percent rRNA in at least 0.25(sample_cutoff) of samples
    final_oligos = {}
    for TX_id, l_spec_oligos in TXspec_oligos.items():
        add_l = 0
        size_spec_add_l =  {ol:0 for ol in oligo_range}
        for ol, oligos in l_spec_oligos.items():
            for oligo in oligos:
                if oligo.score >= sample_cutoff:
                    final_oligos[oligo.id] = oligo
                    add_l += 1
                    size_spec_add_l[ol] += 1

        # if none passes filters, force some designs
        if forceALL>0:
            for ol, oligos in l_spec_oligos.items():
                if size_spec_add_l[ol]<forceALL:
                    for oligo in sorted(oligos, key=lambda x: x.avg_pot, reverse=True):
                        final_oligos[oligo.id] = oligo
                        size_spec_add_l[ol] += 1
                        add_l += 1
                        if size_spec_add_l[ol] == forceALL:
                            break

    print("## Designed", len(final_oligos), "oligos", )
    return(final_oligos)

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


    # Get rRNA sequences
    TX_dic = get_seq(args.rRNA_seqs)
    for k,v in TX_dic.items():
        print('#rRNA sequence read:',k,len(v),'nts')

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

    oligos = {}
    # GET run mode
    if args.oligos2optimize != None:
        ols2opt = get_seq(args.oligos2optimize)
        ol_l_min = min([len(ol) for ol in ols2opt.values()])
        sub = subprocess.Popen(args.RIsearch_exe+' -c '+args.rRNA_seqs+" -o "+args.output+"/rRNA_seqs.suf", shell=True, stdout=subprocess.PIPE).communicate()
        sub = subprocess.Popen(args.RIsearch_exe+' -q '+args.oligos2optimize+" -i "+args.output+"/rRNA_seqs.suf -t 1 -s "+str(ol_l_min/2), shell=True, stdout=subprocess.PIPE).communicate()
        for ol_id, ol_seq in ols2opt.items():
            with gzip.open("risearch_"+ol_id+".out.gz") as self_int:
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
            INT_FILES.append("risearch_"+ol_id+".out.gz")
                #oligo-12	2	24	ENSMUST00000099123	1228	1250 +	-20.66	PPUPPPUWWWWPPPWPPPPWPPP


    else:
        all_rRNA_depths = {}
        all_reads_start_end = {}

        deplete_potentials = {}
        all_rRNA_read_count = {}
        all_total_read_count = {}
        for sample_info in args.alignments:
            sample_label = None
            alignment_folder = None
            if args.prefix!=None and args.suffix!=None:
                sample_label = sample_info
                alignment_folder = args.prefix+sample_info+args.suffix
            else:
                sample_label, alignment_folder = sample_info.rstrip().split(':')
            print("### Reading "+sample_label+" alignment file and getting depth profile.")

            total_mapped = 0
            with open(alignment_folder+'/align_summary.txt') as inf:
                lines = inf.readlines()
                all_total_read_count[sample_label] = int(lines[1].rstrip().replace(' ','').split(':')[1])
                total_mapped = float(lines[2].rstrip().replace(' ','').split(':')[1].split('(')[0])


            reads_start_end, read_count = get_reads('/'.join([alignment_folder,'accepted_hits.bam']), TX_dic, total_mapped)
            all_reads_start_end[sample_label] = reads_start_end

            # deplete_potentials[sample_label] = dep_pot
            all_rRNA_read_count[sample_label] = read_count
            print(str(all_total_read_count[sample_label])+' reads; '+str(all_rRNA_read_count[sample_label])+' (%'+str(round((100.0*all_rRNA_read_count[sample_label])/float(all_total_read_count[sample_label])))+') rRNA reads')

        # oligos fasta
        oligos = get_oligos_with_reads_start_end(all_reads_start_end, all_rRNA_read_count, all_total_read_count, TX_dic,
                                    sample_cutoff=args.score_threshold, rrna_cutoff=args.rrna_perc_threshold, oligo_range=oligo_range, forceALL=args.force_design)

    if oligos != None:
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
        sub = subprocess.Popen(args.RIsearch_exe+' -q '+args.output+"/oligos.fa -i "+args.output+"/oligos.suf -t 1 -s "+str(min(oligo_range)), shell=True, stdout=subprocess.PIPE).communicate()

        print("# Analyzing Emins")
        min_Emin = 0
        for oligo_id in oligos.keys():
            # Get Emin first. Energy of perfect complementary ineraction.
            with gzip.open("risearch_"+oligo_id+".out.gz") as self_int:
                for line in self_int:
                    cols = line.decode('utf-8').rstrip().split('\t')
                    if cols[0]==cols[3] and cols[1]==cols[4] and cols[2]==cols[5] and cols[6]=='-':
                        oligos[oligo_id].add_self_inter(cols[7])
                        # update min Emin
                        min_Emin = oligos[oligo_id].Emin if oligos[oligo_id].Emin<min_Emin else min_Emin
                        INT_FILES.append("risearch_"+oligo_id+".out.gz")
                        break


        if args.no_OFF_search==False:
            print("# RUNning off-targets")
            sub = subprocess.Popen(args.RIsearch_exe+' -q '+args.output+"/oligos.fa -t 1 -i "+args.OFFtargets+" -s "+str(round(0.75*min(oligo_range)))+" -e "+str(min(min_Emin*0.5, -25)),
                shell=True, stdout=subprocess.PIPE).communicate()

        t_c = 0
        for oligo_id in oligos.keys():
            if t_c%10 == 0:
                sys.stdout.write(str(round((100.0*t_c)/len(oligos)))+" percent complete  \r")
            t_c+=1

            oligo = oligos[oligo_id]
            with gzip.open("risearch_"+oligo_id+".out.gz") as self_int:
                oligos[oligo_id].add_offs(self_int.readlines())

        # oligos bed
        with open(args.output+"/oligos.bed",'w') as out_f:
            for oligo in oligos.values():
                out_f.write('\t'.join([oligo.targetTX,str(oligo.startPos),str(oligo.startPos+oligo.l),oligo.id])+'\n')


        # oligos gff
        with open(args.output+"/oligos.gff3",'w') as out_f:
            out_f.write("##gff-version 3\n")
            for oligo in oligos.values():
                out_f.write('\t'.join([oligo.targetTX,'CUSTOM','oligo',str(oligo.startPos+1),str(oligo.startPos+oligo.l),'.','-','.',
                            ';'.join(['ID='+oligo.id,'score='+str(oligo.score),'seq='+oligo.seq,
                            'gc_content='+str(round(oligo.GC,2)),'mfe='+str(oligo.MFE),'structure='+oligo.structure,
                            'bp_per='+str(round(oligo.BPper,2)),'eng_min='+str(round(oligo.Emin,1)),
                            'no_of_OFFS='+"None" if args.no_OFF_search else str(len(oligo.offs))]+[(k+'='+str(v[0])+"_reads-"+str(v[1])+"_pc_rRNA_mapping-"+str(v[2])+'_pc_total') for k,v in oligo.dep_pot.items()])])+'\n')

        # oligos off
        if args.no_OFF_search==False:
            with open(args.output+"/oligos.off",'w') as out_f:
                for oligo in oligos.values():
                    if len(oligo.offs)>0:
                        out_f.write(''.join([off for off in oligo.offs])+'\n')
    else:
        print("## NO OLIGOS FOUND ##")




# Run the main when executed
if __name__ == "__main__":
    main()
    EXIT(1)
