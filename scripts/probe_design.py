#Python script to pick probes from the user provided consensus region of a DNA seq
#Inatall primer3, Biopython, re libraries
#input filename should be test.fasta
#to add decremental from 5' end
import primer3
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio import SeqIO
import re

min_len = 30
value = 0
prefix_nt = "T"
poly_nt = 15

for seq_record in SeqIO.parse("test.fasta", "fasta"):
     probe_list=[]
     probe_bind = str(seq_record.seq)
     probe_bind_len = len(probe_bind)
     difference = probe_bind_len - min_len
     probe_confirm = re.search(r"[^ATGC]", probe_bind)
     if probe_confirm:
        print('probe is not a nucleotide sequence')
        non_nt = probe_confirm.group()
        print("The non nt is " + non_nt)
     else:
        probe  = (prefix_nt * poly_nt) + str(seq_record.seq)
        length = len(probe)    
        if (length <= 60):
            len_probe = len(probe)
            idx = 0
            while len_probe >= min_len:
                Tm = primer3.calcTm(probe)
                Tm2 = ("%.2f" % Tm)
                Hairpin  =  primer3.calcHairpin(probe)
                Homodimer = primer3.calcHomodimer(probe)
                GC_percent = GC(probe)
                GC_per = ("%.2f" % GC_percent)
                idx = idx + 1
                print(seq_record.id + "_" + str(idx) +  "\t" + probe + "\t" + str(GC_per) + "\t" + Tm2 + "\t" + str(Hairpin.structure_found) + "\t" + str(Homodimer.structure_found))
                probe = probe[:-1]
                len_probe = len(probe)
        else:
            print(seq_record.id + 'longer than 60 nt')
