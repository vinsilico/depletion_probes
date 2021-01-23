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
        probe_rc = probe_bind
        probe_bind_len = len(probe_rc)
        if(probe_bind_len <= 60):
            #len_probe = len(probe_bind)
            idx = 0
            while probe_bind_len >= min_len:
                probe_list.append(probe_rc)
                Tm = primer3.calcTm(probe_rc)
                Tm2 = ("%.2f" % Tm)
                Hairpin  =  primer3.calcHairpin(probe_rc)
                Homodimer = primer3.calcHomodimer(probe_rc)
                GC_percent = GC(probe_rc)
                GC_per = ("%.2f" % GC_percent)
                idx = idx + 1
                print(seq_record.id + "_" + str(idx) +  "\t" + probe_rc + "\t" + str(GC_per) + "\t" + Tm2 + "\t" + str(Hairpin.structure_found) + "\t" + str(Homodimer.structure_found))
                probe_rc = probe_rc[:-1]
                probe_bind_len = len(probe_rc)
                #for i in range(0, len(probe_bind)-15, 1):
                 #   probe_list.append(probe_bind[i:i+15])
                  #  probe_bind_len = len(probe_bind)
                   # probe_bind = str(seq_record.seq)
                    #probe_rc = probe_bind
                     #   while probe_bind_len >= min_len:
                      #  probe_list.append(probe_rc)
                      #  probe_rc = probe_rc[1:]
                      #  probe_bind_len = len(probe_rc)
        else:
            print(seq_record.id + ' longer than 60 nt')




   # print("\n".join(probe_list))
