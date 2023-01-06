import pandas as pd
import numpy as np
from tqdm import tqdm

import multiprocessing as mp  # Multiprocessing
import gc

N_JOBS = 2
FRAG_LENGTH = [278, 268, 194]
GENEROSITY = 1
FILES = """
20220930_1.extendedFrags.fastq
"""

# TODO: seaparating removes the quality metrics of the reads
FASTQ_FORMAT = ["id", "sequence", "spacer", "quality"]


def sep(file):
    with open(file, "r") as f:
        data = f.readlines()

        iterables = [
            [
                i
                for i in range(int(len(data) / len(FASTQ_FORMAT)))
            ],
            FASTQ_FORMAT,
        ]

        idx = pd.MultiIndex.from_product(iterables, names=["first", "properties"])
        df = pd.DataFrame(
            data, index=idx).reset_index(level=1)

        df["length"] = df[lambda x: x["properties"] == "sequence"][0].str.len()

        for idx, f_len in tqdm(enumerate(FRAG_LENGTH)):
            frag = df[

                (f_len - GENEROSITY <= df["length"])
                & (df["length"] <= f_len + GENEROSITY)
                ].copy()
            frag.drop("length", axis=1, inplace=True)
            np.savetxt(
                f"{file.split('.')[0]}_F{idx + 1}.{file.split('.')[1]}.fastq",
                frag[0].values,
                fmt="%s",
            )

            print(f"{file.split('.')[0]}_F{idx + 1}.{file.split('.')[1]}.fastq saved")

    del df
    gc.collect()

    return 0


def split(list_a, chunk_size):
    for i in range(0, len(list_a), chunk_size):
        yield list_a[i: i + chunk_size]


def multi_process(files, n_jobs=N_JOBS):
    """
    Multiprocessing
    """
    CHUNK_SIZE = n_jobs
    for chunk in tqdm(list(split(files, CHUNK_SIZE))):
        pool = mp.Pool(n_jobs)
        pool.map_async(sep, chunk)
        pool.close()
        pool.join()


if __name__ == "__main__":

    files = FILES.split("\n")[1:-1]
    # multi_process(files)

    for f in files:
        sep(f)
