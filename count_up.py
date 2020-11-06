#!/usr/bin/env python3

# this is a python script template
# this next line will download the file using curl

gff="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz"
fasta="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz"

import os,gzip,itertools,csv,re,sys

# this is code which will parse FASTA files
# define what a header looks like in FASTA format
def isheader(line):
    return line[0] == '>'

def aspairs(f):
    seq_id = ''
    sequence = ''
    for header,group in itertools.groupby(f, isheader):
        if header:
            line = next(group)
            seq_id = line[1:].split()[0]
        else:
            sequence = ''.join(line.strip() for line in group)
            yield seq_id, sequence



if not os.path.exists(gff):
    os.system("curl -O ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/gff3/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz")

if not os.path.exists(fasta):
    os.system("curl -O ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz")
    
with gzip.open(gff,"rt") as fh:
    genecount = 0
    total = 0
    gff = csv.reader(fh,delimiter="\t")
    for row in gff:
        if row[0].startswith("#"):
            continue
        if row[2].startswith("gene"):
            genecount += 1       
            intarray =[int(row[3]),int(row[4])]
            diff = intarray[1] - intarray[0]
            total += diff

with gzip.open(fasta, "rt") as f: 
    pairs= aspairs(f)    
    seqs= dict(pairs)	
    for seqid in seqs:
        bp = len(seqs[seqid])
              
print("2) There are {} genes in E.coli".format(genecount))        
print("3) The total length of all the E.coli genes is {}.".format(total))
print("4) There are {} basepairs in the E.coli genome.".format(bp))
print("5) %.2f percent of the E.coli genome is coding."%(100 * total/bp))      
   		
