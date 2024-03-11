#!/usr/bin/python
# -*- coding: utf-8 -*-

# Original code is here, https://github.com/MyungjaeSong/Paired-Library.git
# This is the modified version of the program for academic uses in SKKUGE Lab

__author__ = "forestkeep21@naver.com"
__editor__ = "poowooho3@g.skku.edu"

import gc
import pathlib

import numpy as np
import pandas as pd
import skbio
from tqdm import tqdm


def extract_read_cnts(
    sequence_file: pathlib.Path, barcode_file: pathlib.Path, result_dir, sep=":"
):
    # df index == barcode, column == read count
    # Positional meaning of barcode file: Name, Barcode, 1st barcode-paired substring, 2nd barcode-paired substring, ...

    from torch import cuda
    from dask.dataframe import from_pandas as dask_from_pandas
    import swifter

    tqdm.pandas()
    # Load barcode file
    result_df = pd.read_csv(
        barcode_file, sep=sep, header=None, names=["Gene", "Barcode"]
    ).iloc[:, [0, 1]]  # Use only Gene and Barcode columns
    if not result_df["Barcode"].is_unique:
        # Barcode used as a PK in the database, so duplication is not allowed
        print("Barcode duplication detected!")
        print("Remove duplicated Barcodes... only the first one will be kept.")
        result_df.drop_duplicates(subset=["Barcode"], keep="first", inplace=True)

    result_df["Barcode"] = result_df["Barcode"].str.upper()
    result_df["Barcode_copy"] = result_df["Barcode"]

    result_df = result_df.set_index("Barcode")  # TODO: tentative design

    # Load a split sequencing result using high-level I/O; validating fastq format
    result_df["Read_counts"] = 0
    result_df["ID"] = ""
    seqs = skbio.io.read(
        sequence_file, format="fastq", verify=True, variant="illumina1.8"
    )  # FASTQ format verification using skbio

    seq_df = pd.DataFrame(
        [(seq.metadata["id"], seq._string.decode()) for seq in seqs],
        columns=["ID", "Sequence"],
    )  # Sequence pool
    seq_detection_array = np.zeros(seq_df.shape[0], dtype=bool)  # TODO: remove it?

    def find_occurrences_for_barcode(x: pd.Series):
        query_result = seq_df["Sequence"].str.contains(row["Barcode_copy"])
        # boolean indexing for fast processing
        result_df.at[idx, "ID"] = (
            seq_df.loc[query_result[query_result].index]["ID"].to_numpy().tolist()
        )
        result_df.loc[idx, "Read_counts"] = (query_result.sum(),)

        seq_detection_array[query_result[query_result].index] = True

        # TODO: Sample with replacement option
        # Without replacement from the sequence pool
        if not sample_replacement_mode:
            seq_df.drop(query_result[query_result].index, inplace=True, axis=0)

    result_df.swifter.apply(find_occurrences_for_barcode, axis=1)

    for idx, row in tqdm(result_df.iterrows()):
        # query_result format
        #

        query_result = seq_df["Sequence"].str.contains(row["Barcode_copy"])
        # boolean indexing for fast processing
        result_df.at[idx, "ID"] = (
            seq_df.loc[query_result[query_result].index]["ID"].to_numpy().tolist()
        )
        result_df.loc[idx, "Read_counts"] = (query_result.sum(),)

        seq_detection_array[query_result[query_result].index] = True

        # TODO: Sample with replacement option
        # Without replacement from the sequence pool
        seq_df.drop(query_result[query_result].index, inplace=True, axis=0)

    del seq_df
    gc.collect()

    result_df.drop("Barcode_copy", axis=1, inplace=True)
    result_df.reset_index(inplace=True, drop=False)
    result_df.iloc[1], result_df.iloc[-1] = result_df.iloc[-1], result_df.iloc[1]

    def name():
        from datetime import datetime

        dt_string = datetime.now().strftime("%Y-%m-%d;%H:%M:%S")
        return str(f"{dt_string}")

    result_df.to_parquet(
        f"{result_dir}/{name()}+{pathlib.Path(sequence_file).name}.parquet"
    )

    del result_df
    gc.collect()

    return seq_detection_array


def main(*args) -> pd.DataFrame:
    (sequence, barcode, logger, result_dir, sep) = args[0]

    # start = time.time()
    rval = extract_read_cnts(sequence, barcode, result_dir, sep)
    # end = time.time()

    # logger.info(f"Extraction is done. {end - start}s elapsed.")

    return rval
