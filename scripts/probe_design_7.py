#Python script to pick probes from the user provided consensus region of a DNA seq
#Inatall primer3, Biopython, re libraries
#input filename should be test.fasta
#to add decremental from 5' end
import primer3
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio import SeqIO
import re

min_len = 15
value = 0
prefix_nt = "T"
poly_nt = 15

for seq_record in SeqIO.parse("test.fasta", "fasta"):
    probe_bind = str(seq_record.seq.reverse_complement())
    probe_confirm = re.search(r"[^ATGC]", probe_bind)
    if probe_confirm:
        print('probe is not a nucleotide sequence')
        non_nt = probe_confirm.group()
        print("The non nt is " + non_nt)
    else:
        probe_list=[]
        probe_bind = str(seq_record.seq.reverse_complement())
        probe_bind_len = len(probe_bind)
        #probe_bind = str(seq_record.seq)
        probe_rc = probe_bind
        while probe_bind_len >= min_len:
            probe_list.append(probe_rc)
            probe_rc = probe_rc[:-1]
            probe_bind_len = len(probe_rc)
        for i in range(0, len(probe_bind)-15, 1):
            probe_list.append(probe_bind[i:i+15])
        probe_bind_len = len(probe_bind)
        probe_bind = str(seq_record.seq)
        probe_rc = probe_bind
        while probe_bind_len >= min_len:
            probe_list.append(probe_rc)
            probe_rc = probe_rc[1:]
            probe_bind_len = len(probe_rc)
print("\n".join(probe_list))
