import pandas as pd

import numpy as np
from tqdm import tqdm

import multiprocessing as mp  # Multiprocessing
import gc
import pathlib

N_JOBS = 1
FRAG_LENGTH = [278, 268, 194]
GENEROSITY = 1
FILES = """
./separator/test
"""
# ./separator/20230109_1.extendedFrags.fastq
# ./separator/20230109_2.extendedFrags.fastq
# ./separator/20230109_3.extendedFrags.fastq
# ./separator/20230109_4.extendedFrags.fastq
# ./separator/20230109_5.extendedFrags.fastq

# TODO: seaparating removes the quality metrics of the reads
FASTQ_FORMAT = ["id", "sequence", "spacer", "quality"]


def sep(file):
    from torch import cuda

    df = pd.DataFrame(pd.read_csv(file, header=None).values.reshape(-1, 4), columns=FASTQ_FORMAT)
    cuda_available = cuda.is_available()
    # cuda_available = False    # debug
    if cuda_available:
        # TODO: cudf_dask
        import cudf
        import dask_cudf as dc
        print("Nvidia GPU detected!")
        df = cudf.from_pandas(df)
        df = dc.from_cudf(df, npartitions=mp.cpu_count())

        df["length"] = (
            df["sequence"].compute().str.len()
        )

        for idx, f_len in tqdm(enumerate(FRAG_LENGTH)):
            frag = df[
                (f_len - GENEROSITY <= df["length"])
                & (df["length"] <= f_len + GENEROSITY)
                ].compute()
            frag.drop("length", axis=1, inplace=True)
            np.savetxt(
                f"{pathlib.Path(file).absolute()}_F{idx + 1}.fastq",
                frag.to_pandas().values.reshape(-1, 1),
                fmt="%s",
                newline="",
            )

            print(
                f"{pathlib.Path(file).absolute()}_F{idx + 1}.fastq saved"
            )
    else:
        import dask.dataframe as dd
        print("No Nvidia GPU in system!")

        df = dd.from_pandas(df, npartitions=mp.cpu_count())

        df["length"] = (
            df["sequence"].compute().str.len()
        )

        for idx, f_len in tqdm(enumerate(FRAG_LENGTH)):
            frag = df[
                (f_len - GENEROSITY <= df["length"])
                & (df["length"] <= f_len + GENEROSITY)
                ].compute()
            frag.drop("length", axis=1, inplace=True)
            np.savetxt(
                f"{pathlib.Path(file).absolute()}_F{idx + 1}.fastq",
                frag.values.reshape(-1, 1),
                fmt="%s",
                newline="",
            )

            print(
                f"{pathlib.Path(file).absolute()}_F{idx + 1}.fastq saved"
            )

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
