#!/usr/bin/python
# -*- coding: utf-8 -*-

# Original code is here, https://github.com/MyungjaeSong/Paired-Library.git
# This is the modified version of the program for academic uses in SKKUGE Lab

__author__ = "forestkeep21@naver.com"
__editor__ = "poowooho3@g.skku.edu"


import time
import traceback

import dask.dataframe as dd
import pandas as pd
from icecream import ic
from tqdm import tqdm


def extractor_main(
    sequence_frame: dd.DataFrame,
    barcode_df: pd.DataFrame,
    result_dir,
    sep=",",
    logger=None,
    chunk_number=None,
):
    start_time = time.time()
    try:
        if chunk_number is None:
            raise ValueError("chunk_number is not defined")

        pbar = tqdm(total=barcode_df.shape[0])
        pbar.set_description(f"Barcode extraction...{chunk_number}")
        for ntp in barcode_df.itertuples():
            gene = ntp.Gene
            barcode = ntp.Query
            sequence_frame[gene] = sequence_frame["Sequence"].str.contains(
                barcode, regex=True
            )
            pbar.update(1)

        pbar.close()

        # Drop reads that do not contain any barcode
        sequence_frame["Detection_count"] = sequence_frame.sum(
            axis=1, numeric_only=True
        ).astype(int)
        sequence_frame = sequence_frame[sequence_frame["Detection_count"] > 0]

        sequence_frame.to_parquet(
            f"{result_dir}/parquets/{chunk_number}",
            compression="snappy",
            engine="pyarrow",
            compute=True,
            write_index=True,
            write_metadata_file=True,
        )
        end_time = time.time()
        ic(f"Barcode extraction finished...{gene} in {end_time-start_time} seconds")
        return 0

    except Exception as e:
        ic(e)
        ic(traceback.format_exc())

        return -1
