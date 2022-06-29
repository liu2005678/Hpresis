#!/usr/bin/env python
# coding: utf-8

import h5py
import pandas as pd
import numpy as np
import os
import math
import sys


fast5_dir= sys.argv[1] ; print("\nfast5_dir is: ", fast5_dir)
sum_dir = sys.argv[2] ; print("directory for sequencing summary file is: ", sum_dir)
batches = sys.argv[3]; print("subfolders for batches is: ", batches)
head = sys.argv[4]; print("file headers for splitted file is: ", head)
savedir = sys.argv[5] + batches + "/" ; print("saving directory is: ", savedir)  
chunksize=sys.argv[6]; print("chunksize is ", chunksize, "\n")
# fast5_dir="/media/logen/sdc/h.pylori/called2/workspace/hpgenomes73-84_20220107/trial2/20220107_0944_MN31740_FAS01060_7f01268a/fast5/"
# sum_dir="/media/logen/sdc/h.pylori/called2/"
# batches="hpgenomes73-84_20220107"
# head="FAS01060_7f01268a"
# savedir="/home/logen/fast5/" + batches + "/" 
# chunksize=20
chunksize=int(chunksize)


fs = os.listdir(fast5_dir)
if not os.path.exists(savedir):
    os.makedirs(savedir)


summary_dir=sum_dir + "sequencing_summary.txt"
summary=pd.read_table(summary_dir)
subset_sum=summary[['read_id','barcode_arrangement']]
subset_sum=subset_sum.assign(readid = lambda dataframe: dataframe['read_id'].map(lambda read_id: "read_" + read_id))

subset_sum['read_id'] = subset_sum['readid']
subset_sum.head()
barcodes=pd.unique(subset_sum["barcode_arrangement"])
for bari in barcodes:
    savedir2=savedir + bari
    if not os.path.exists(savedir2):
        os.makedirs(savedir2)


for i in range(math.ceil(len(fs)/chunksize)):
    chunk = fs[chunksize*i:min(len(fs), chunksize*(i+1) )]
    tmp_fast5=[ {} for _ in range(len(barcodes)) ]
    tmp_path=[[] for _ in range(len(barcodes)) ]
    for bari in range(len(barcodes)):
        tmp_path[bari]=savedir + barcodes[bari] + "/" + head + "_" + str(i) + ".fast5"
        tmp_fast5[bari] = h5py.File(tmp_path[bari],"a",driver='core')
    for fast5 in chunk:    
        files = fast5_dir + fast5
        print("reading ", files)
        f = h5py.File(files, "r",driver='core')
        fast5_ids=list(f.keys())
        demux_ids=[ [] for _ in range(len(barcodes)) ]
        for bari in range(len(barcodes)):
            barcode_temp=subset_sum[subset_sum['barcode_arrangement'] == barcodes[bari]]
            barcode_temp = list(barcode_temp['read_id'])
            demux_ids[bari].extend(np.intersect1d(barcode_temp, fast5_ids))
            for reads in demux_ids[bari]:
                f.copy(reads,tmp_fast5[bari])
        f.close()
    for bari in range(len(barcodes)):
        print("writting ", tmp_path[bari])
        tmp_fast5[bari].close()


