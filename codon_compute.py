#!/usr/bin/env python3

import os, gzip, itertools, re

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

url1="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2/cds/Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
url2="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/mycobacterium_tuberculosis_h37rv/cds/Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"
file1="Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
file2="Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"

if not os.path.exists(file1):
    os.system("curl -O %s"%(url1))

if not os.path.exists(file2):
    os.system("curl -O %s"%(url2))

with gzip.open(file1,"rt") as fh:
    pair = aspairs(fh)
    seqs = dict(pair)
    seqnumS = 0
    tlS = 0
    totGCS = 0
    for seq in seqs:
        seqnumS += 1
        seqlenS = len(seqs[seq])
        tlS += seqlenS
        GCS = 0
        for bp in seqs[seq]:
            if bp == "G":
                GCS += 1
            if bp == "C":
                GCS += 1
            else:
                continue
        totGCS += GCS
    codonS = {}
    for first in {'A','T','C','G'}:
        for second in {'A','T','C','G'}:
            for third in {'A','T','C','G'}:
                new_codon = first+second+third
                codonS[new_codon] = 0 
    for seq in seqs:     
        for n in range(0, len(seqs[seq]), 3):
            codonS[seqs[seq][n:n+3]] += 1 
    total_codonsS = 0
    for value in codonS: 
        total_codonsS += codonS[value]
    

with gzip.open(file2,"rt") as fh:
    pair = aspairs(fh)
    seqs = dict(pair)
    seqnum = 0
    tl = 0 
    totGC = 0
    for seq in seqs:
        seqnum += 1
        seqlen = len(seqs[seq])
        tl += seqlen
        GC = 0
        for bp in seqs[seq]:
            if bp == "G": 
                GC += 1
            if bp == "C": 
                GC += 1 
            else: 
                continue
        totGC += GC
    codon = {}
    for first in {'A','T','C','G'}:
        for second in {'A','T','C','G'}:
            for third in {'A','T','C','G'}:
                new_codon = first+second+third
                codon[new_codon] = 0
    for seq in seqs:
        for n in range(0, len(seqs[seq]), 3):
            codon[seqs[seq][n:n+3]] += 1
    total_codons = 0
    for value in codon:
        total_codons += codon[value]


print("1) There are {} genes in Salmonella".format(seqnumS))
print("There are {} genes in Mycobacterium".format(seqnum))
print("2) The total length of genes in Salmonella is {} bp.".format(tlS))
print("The total length of genes in Mycobacterium is {} bp.".format(tl))
print("3) The perc GC for Salmonella is %.2f."%(100 * totGCS/tlS))
print("The perc GC for Mycobacterium is %.2f."%(100 * totGC/tl))
print("4) The total num codons in Salmonella is %d."%(total_codonsS))
print("The total num codons in Mycobacterium is %d."%(total_codons))
print("Codon" "\t" "Salmonella" "\t""\t" "Mycobacterium")  
for value in codonS:
    print("\t".join([value,str(100 * codonS[value]/total_codonsS), str(100 * codon[value]/total_codons)]))



