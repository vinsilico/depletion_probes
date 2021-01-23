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

for seq_record in SeqIO.parse("test2.fasta", "fasta"):
     probe_list=[]
     probe_bind = str(seq_record.seq.reverse_complement())
     #print(probe_bind_rc)
     probe_bind_len = len(probe_bind)
     #probe_bind = str(seq_record.seq)
     probe_rc = probe_bind
     while probe_bind_len >= min_len:
        probe_list.append(probe_rc)
        probe_rc = probe_rc[:-1]
        probe_bind_len = len(probe_rc)
        print(probe_rc)
     print(probe_list[5])

for seq_record in SeqIO.parse("test2.fasta", "fasta"):
     probes_win=[]
     probe_bind = str(seq_record.seq.reverse_complement())
     for i in range(0, len(probe_bind)-15, 1):
         #print(seq_record.id + "\n")
         probes_win.append(probe_bind[i:i+15])
         print(probe_bind[i:i+15])
     print(probes_win[5])


for seq_record in SeqIO.parse("test2.fasta", "fasta"):
     probe_list_f=[]
     probe_bind = str(seq_record.seq.reverse_complement())
     #print(probe_bind_rc)
     probe_bind_len = len(probe_bind)
     #probe_bind = str(seq_record.seq)
     probe_rc = probe_bind
     while probe_bind_len >= min_len:
        probe_list_f.append(probe_rc)
        probe_rc = probe_rc[1:]
        probe_bind_len = len(probe_rc)
        print(probe_rc)
     print(probe_list_f[5])
