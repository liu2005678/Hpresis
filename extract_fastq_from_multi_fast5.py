#!/usr/bin/env python
# coding: utf-8


import h5py
import os
import sys

fast5_dir = sys.argv[1] ; print("\nfast5_dir is: ", fast5_dir)
save_dir = sys.argv[2] ; print("file to be writting is: ", save_dir)


# fast5_dir="/media/logen/sdc/h.pylori/demux_fast5/link_by_strains/Hpfe073"
# save_dir="/media/logen/sdc/h.pylori/demux_fast5/fastq_by_strains/Hpfe073.fq"


fs = os.listdir(fast5_dir)
inmemfq=""
for file in fs:
    path= fast5_dir + "/" + file
    print("processing ", file)
    f5= h5py.File(path,"r",driver='core')
    keys=list(f5.keys())
    for key in keys:
        read = key + "/Analyses/Basecall_1D_000/BaseCalled_template/Fastq"
        fq=f5[read]
        stred = fq[()].split()
        for i in range(len(stred)):
            stred[i] = stred[i].decode("utf-8")
        joinstring = " "
        stred[6]=joinstring.join(stred[0:7])
        joinfq = "\n"
        fqseq = joinfq.join(stred[6:10])
        inmemfq += fqseq + "\n"
    f5.close()
f = open(save_dir, "w")
f.write(inmemfq)
f.close()




