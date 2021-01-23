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
    probe_non_nt = re.search(r"[^ATGC]", probe_bind)
    if probe_non_nt:
        print('probe is not a nucleotide sequence')
        non_nt = probe_confirm.group()
        print("The non nt is " + non_nt)
    else:
        probe_list=[]
        long_seq_list=[]
        #probe_bind = str(seq_record.seq.reverse_complement())
        probe_rc = probe_bind
        probe_fc = probe_bind
        probe_bind_len = len(probe_rc)
        probe_bind_len2 = len(probe_rc)
        probe_list.append(probe_bind)
        if(probe_bind_len <= 60):
           idx_sl = 0
           idx_rc = 0
           idx_fc = 0
           for i in range(0, len(probe_bind)-min_len, 1):
               probe_list.append(probe_bind[i:i+15])
               idx_sl = idx_sl + 1
           while probe_bind_len >= min_len:
               probe_rc = probe_rc[:-1]
               probe_list.append(probe_rc)
               probe_bind_len = len(probe_rc)
               idx_rc = idx_rc + 1
           while probe_bind_len2 >= min_len:
               probe_fc = probe_fc[1:]
               probe_list.append(probe_fc)
               probe_bind_len2 = len(probe_fc)
               idx_fc = idx_fc + 1
                                        
 #              probe_bind_len = len(probe_bind)
 #              probe_bind = str(seq_record.seq)
 #              probe_rc = probe_bind
#
 ##           
 #           
  #          while probe_bind_len >= min_len:
         #       probe_list.append(probe_rc)
#
#
#

 #               probe_rc = probe_rc[:-1]
  #              probe_bind_len = len(probe_rc)
           #     idx = idx+1
   #         for i in range(0, len(probe_bind)-min_len, 1):
    #            probe_list.append(probe_bind[i:i+15])
     #           probe_bind_len = len(probe_bind)
            #    probe_bind = str(seq_record.seq)
      #          probe_rc = probe_bind
       #     while probe_bind_len >= min_len:
   ####             probe_list.append(probe_rc)
       #         probe_rc = probe_rc[1:]
        #        probe_bind_len = len(probe_rc)
        else:
            long_seq_list.append(seq_record.id)
            #print(seq_record.id + ' longer than 60 nt')
                
                
              #  Tm2 = ("%.2f" % Tm)
              #  Hairpin  =  primer3.calcHairpin(probe_rc)
              #  Homodimer = primer3.calcHomodimer(probe_rc)
               # GC_percent = GC(probe_rc)
               # GC_per = ("%.2f" % GC_percent)
                #idx = idx + 1
                #print(seq_record.id + "_" + str(idx) +  "\t" + probe_rc + "\t" + str(GC_per) + "\t" + Tm2 + "\t" + str(Hairpin.structure_found) + "\t" + str(Homodimer.structure_found))
                #probe_rc = probe_rc[:-1]
                #probe_bind_len = len(probe_rc)
                #for i in range(0, len(probe_bind)-15, 1):
                 #   probe_list.append(probe_bind[i:i+15])
                  #  probe_bind_len = len(probe_bind)
                   # probe_bind = str(seq_record.seq)
                    #probe_rc = probe_bind
                     #   while probe_bind_len >= min_len:
                      #  probe_list.append(probe_rc)
                      #  probe_rc = probe_rc[1:]
                      #  probe_bind_len = len(probe_rc)
       # else:
       #     print(seq_record.id + ' longer than 60 nt')



    print(idx_sl + idx_rc + idx_fc)
    print("\n".join(probe_list))

    print("\n".join(long_seq_list))

