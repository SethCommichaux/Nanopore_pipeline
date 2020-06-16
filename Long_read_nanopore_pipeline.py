# This pipeline is a Python3 wrapper for long read assembly
#
# Long read assembly pipeline
# 1) Collect raw read statistics --> NanoPlot
# 2) Filter/trim raw reads --> NanoFilt
# 3) Assembly --> Flye, Miniasm+Racon
# 4) Collect assembly statistics --> Python3
# 5) Gene prediction --> Prokka

##########################################
###### Required Python3 libraries ########
##########################################

import argparse
import os
from Bio import SeqIO

##########################################
##### Required bioinformatic tools #######
##########################################

# NanoPlot
# NanoFilt
# Flye
# Miniasm
# Minimap2
# Racon
# Prokka
miniasm="/lustre/projects/SethCommichaux/Nanopore_pipeline_dev/mark_data/miniasm/miniasm"

#########################################
###### Functions for running tools ######
#########################################

# Generate reads statistics and plots
def read_stat_plots(reads,minqual,minLength,threads):
        print("Generating statistics and plots for your raw reads with NanoPlot...")
        os.system('NanoPlot -t %s --raw --o %s.stats_plots --minqual %s --minlength %s --fastq %s' % (threads,reads,minqual,minLength,reads))
        print("Statistics gathered and plots created!")

# filter/trim the input long reads
def filter_trim_reads(reads,filtered_reads,minqual,minLength,threads):
        print("Trimming/filtering the long reads using NanoFilt...")
        os.system('NanoFilt -l %s -q %s %s > %s' % (minLength,minqual,reads,filtered_reads))
        print("Long reads trimmed and filtered!")

# assemble long reads with flye
def assembly_flye(longReads_filtered,genome_size,threads):
	print("Beginning to assemble reads with Flye...")
	os.system('flye --nano-raw %s -g %s -o %s.flye -t %s' % (longReads_filtered,genome_size,longReads_filtered,threads))
	print("Flye assembly finished!\nCollecting assembly statistics...")
	assembly_stats('%s.flye/assembly.fasta' % (longReads_filtered))
	print("Assembly statistics collected!\nPredicting genes in assembly...")
	predict_genes('%s.flye/assembly.fasta' % (longReads_filtered),threads)
	print("Genes predicted!\nAnalysis finished for flye assembly!!!")

# assemble long reads with Miniasm
def assembly_miniasm(longReads_filtered,threads):
	print("Beginning to assemble reads with Miniasm...")
	# map reads to themselves to determine overlaps
	os.system('minimap2 -x ava-ont -t12 %s %s | gzip -1 > %s.miniasm.paf.gz' % (longReads_filtered,longReads_filtered,longReads_filtered))
	# run miniasm assembly
	os.system('%s -f %s %s.miniasm.paf.gz > %s.miniasm.gfa' % (miniasm,longReads_filtered,longReads_filtered,longReads_filtered))
	# convert GFA file to FASTA
	gfa2fasta(longReads_filtered+'.miniasm.gfa')
	# map reads back to assembly for polishing
	os.system('minimap2 -a -t 12 %s.miniasm.gfa.fasta %s > %s.miniasm.sam' % (longReads_filtered,longReads_filtered,longReads_filtered))
	# polishing with Racon
	os.system('racon -t 12 %s %s.miniasm.sam %s.miniasm.gfa.fasta > %s.miniasm.racon' % (longReads_filtered,longReads_filtered,longReads_filtered,longReads_filtered))
	print("Miniasm assembly finished!\nCollecting assembly statistics...")
	assembly_stats('%s.miniasm.racon' % (longReads_filtered))
	print("Assembly statistics collected!\nPredicting genes in assembly...")
	predict_genes('%s.miniasm.racon' % (longReads_filtered),threads)
	print("Genes predicted!\nAnalysis finished for miniasm assembly!!!")

# convert GFA to fasta file
def gfa2fasta(GFA):
	with open(GFA+'.fasta','w') as out:
		for i in open(GFA):
			tmp = i.strip().split('\t')
			if 'S' == tmp[0]:
				out.write(">"+tmp[1]+"\n"+tmp[2]+"\n")

def assembly_stats(assembly):
        print('Collecting assembly statistics...')
        contig_lens = [len(i.seq) for i in SeqIO.parse(assembly,'fasta')]
        total_assembly_length = sum(contig_lens)
        number_contigs = len(contig_lens)
        longest_contig = max(contig_lens)
        n = 0
        for i in sorted(contig_lens,reverse=True):
                n += i
                if n >= total_assembly_length/2.:
                        n50 = i
                        break
        with open(assembly+'.stats','w') as out:
                out.write('''Total assembly length (bp): %s
Number of contigs: %s
N50: %s
Longest contig: %s
''' % (total_assembly_length,number_contigs,n50,longest_contig))
        print('Assembly statistics collected.')

# predict genes with Prokka
def predict_genes(assembly,threads):
        print('Predicting genes in assembly...')
        os.system('prokka --cpus %s --outdir %s.prokka %s' % (threads,assembly,assembly))
        print('Genes predicted in assembly.')

#########################################
###### User input parameters ############
#########################################

# parse user inputs
parser = argparse.ArgumentParser()
parser.add_argument("-L", help="Long reads fastq file",required=True)
parser.add_argument("-g", help="Estimated genome size of species in sample being assembled [k/m/g] e.g. for 2 million bp use 2m",required=True)
parser.add_argument("-t", help="Number of threads to use",default="12")
parser.add_argument("-min_qual", help="Minimum average qscore for a read to be kept for further analysis. Default is 10.",default="10")
parser.add_argument("-min_len", help="Minimum length for a read to be kept for further analysis. Default is 500.",default="500")
args = parser.parse_args()

longReads = args.L
longReads_filtered = longReads.replace('fastq','trimmed.fastq')
threads = args.t
genome_size = args.g
min_read_qual = args.min_qual
min_read_len = args.min_len

#########################################
###### Tool command calls ###############
#########################################

# collect statistics about the long reads such as their length and base quality distribution before and after trimming/filtering
read_stat_plots(longReads,min_read_qual,min_read_len,threads)

# filter/trim long reads
filter_trim_reads(longReads,longReads_filtered,min_read_qual,min_read_len,threads)

# assemble long reads with Flye
assembly_flye(longReads_filtered,genome_size,threads)

# assemble long reads with Miniasm
assembly_miniasm(longReads_filtered,threads)

