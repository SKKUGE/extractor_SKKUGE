#!/usr/bin/python
# -*- coding: utf-8 -*-

# Original code is here, https://github.com/MyungjaeSong/Paired-Library.git
# This is the modified version of the program for academic uses in SKKUGE Lab

__author__ = "forestkeep21@naver.com"
__editor__ = "poowooho3@g.skku.edu"


import pathlib
import traceback

import dask.dataframe as dd
import pandas as pd
from icecream import ic

# def load_test():
#     import time  # debug

#     if chunk_number is not None:
#         ic(chunk_number)

#     start_time = time.time()
#     cnt = 0
#     while True:
#         cnt += 1
#         ic(f"{chunk_number}:{cnt}")
#         # time.sleep(0.1)
#         if time.time() - start_time >= 5:
#             break


def extract_read_cnts(
    sequence_frame: dd.DataFrame,
    barcode_df: pd.DataFrame,
    result_dir,
    sep=",",
    logger=None,
    chunk_number=None,
):
    try:
        if chunk_number is None:
            raise ValueError("chunk_number is not defined")

        if not barcode_df["Gene"].is_unique or not barcode_df["Barcode"].is_unique:
            # Barcode used as a PK in the database, so duplication is not allowed
            ic(
                f"Barcode duplication detected! Check your program run design {chunk_number}"
            )
            ic(
                f"Remove duplicated Barcodes... only the first one will be kept. {chunk_number}"
            )
            barcode_df.drop_duplicates(subset=["Barcode"], keep="first", inplace=True)

        barcode_df["Barcode"] = barcode_df["Barcode"].str.upper()

        ic(f"Barcode extraction initiated...{chunk_number}")

        for i, (gene, barcode) in barcode_df.iterrows():
            # ic(gene, barcode) # DEBUG
            sequence_frame[gene] = sequence_frame["Sequence"].str.contains(
                barcode, regex=True
            )

            if i % (2**12) == 0:  # DEBUG
                # ic(i)
                sequence_frame = sequence_frame.persist()

        # Drop heavy columns
        sequence_frame = sequence_frame.drop(
            columns=["Sequence"],
        )

        # OPTION 1 : Save as parquet
        pathlib.Path(f"{result_dir}/visualize").mkdir(parents=True, exist_ok=True)
        sequence_frame.visualize(filename=f"{result_dir}/visualize/{chunk_number}.png")
        sequence_frame.to_parquet(
            f"{result_dir}/parquets/{chunk_number}",
            compression="snappy",
            engine="pyarrow",
            write_index=True,
            write_metadata_file=True,
            compute=True,
        )
        del sequence_frame, barcode_df
        return f"{result_dir}/parquets/{chunk_number}"

    except Exception as e:
        ic(e)
        ic(traceback.format_exc())
        logger.error(e)
        return -1


def extractor_main(sequence, barcode, logger, result_dir, sep, chunk_number):

    rval = extract_read_cnts(
        sequence, barcode, result_dir, sep, logger, chunk_number=chunk_number
    )  # return 0 upon successful completion

    if rval == -1:
        ic(rval)
        logger.error("Barcode extraction failed")
        return rval

    logger.info("Barcode extraction completed")
    logger.info("Merging parquet files...")
    return rval  # OPTION 1: parquet path
