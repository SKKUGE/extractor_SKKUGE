import pandas as pd
import numpy as np
from tqdm import tqdm

import multiprocessing as mp    # Multiprocessing
import gc



FRAG_LENGTH = [278,268,194]
GENEROSITY = 1
FILES ="""
20220314_1.extendedFrags.fastq
20220314_2.extendedFrags.fastq
20220314_3.extendedFrags.fastq
20220314_4.extendedFrags.fastq
20220314_5.extendedFrags.fastq
20220314_6.extendedFrags.fastq
20220314_7.extendedFrags.fastq
20220314_8.extendedFrags.fastq
20220314_9.extendedFrags.fastq
20220314_10.extendedFrags.fastq
20220517_37.extendedFrags.fastq
20220629_6.extendedFrags.fastq
20220629_7.extendedFrags.fastq
20220629_8.extendedFrags.fastq
20220629_9.extendedFrags.fastq
20220629_10.extendedFrags.fastq
20220629_11.extendedFrags.fastq
20220629_12.extendedFrags.fastq
20220629_13.extendedFrags.fastq
20220727_4.extendedFrags.fastq
20220727_12.extendedFrags.fastq
20220727_20.extendedFrags.fastq
20220819_4.extendedFrags.fastq
20220819_12.extendedFrags.fastq
20220819_20.extendedFrags.fastq
20220824_4.extendedFrags.fastq
20220824_12.extendedFrags.fastq
20220824_20.extendedFrags.fastq
20220830_1.extendedFrags.fastq
20220830_2.extendedFrags.fastq
20220830_3.extendedFrags.fastq
20220830_4.extendedFrags.fastq
20220830_5.extendedFrags.fastq
20220830_6.extendedFrags.fastq
20220830_7.extendedFrags.fastq
20220830_8.extendedFrags.fastq
20220902_1.extendedFrags.fastq
20220902_9.extendedFrags.fastq
20220902_17.extendedFrags.fastq
20220915_1.extendedFrags.fastq
20220915_9.extendedFrags.fastq
20220915_17.extendedFrags.fastq
"""

def sep(file):
    with open(file, 'r')as f:
        df = pd.DataFrame(f.readlines())
        df['length'] = df[0].str.len()
        
        for idx, f_len in enumerate(FRAG_LENGTH):
            frag = df[(f_len -GENEROSITY <= df['length']) & (df['length'] <= f_len + GENEROSITY)].copy()
            frag.drop('length', axis=1, inplace=True)
            np.savetxt(f"{file.split('.')[0]}_F{idx+1}.{file.split('.')[1]}.fastq",frag.values, fmt='%s')

            print(f"{file.split('.')[0]}_F{idx+1}.{file.split('.')[1]}.fastq saved")

    del df
    gc.collect()

    return 0

def split(list_a, chunk_size):

  for i in range(0, len(list_a), chunk_size):
    yield list_a[i:i + chunk_size]

def multi_process(files, n_jobs=4):
    """
    Multiprocessing
    """
    CHUNK_SIZE = n_jobs
    for chunk in tqdm(list(split(files,CHUNK_SIZE))):
        pool = mp.Pool(n_jobs)
        pool.map_async(sep, chunk)
        pool.close()
        pool.join()
    

if __name__ == '__main__':

    files = FILES.split('\n')[1:-1] 
    multi_process(files)
    # for f in files:
    #     sep(f)
            

    print()