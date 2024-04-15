#!/usr/bin/python
# -*- coding: utf-8 -*-

# Original code is here, https://github.com/MyungjaeSong/Paired-Library.git
# This is the modified version of the program for academic uses in SKKUGE Lab

__author__ = "forestkeep21@naver.com"
__editor__ = "poowooho3@g.skku.edu"


import time
import traceback

from icecream import ic
from tqdm import tqdm


def extractor_main(sequence_frame, barcode_df, logger, result_dir):
    start_time = time.time()
    try:
        if not barcode_df["Gene"].is_unique:
            # Barcode used as a PK in the database, so duplication is not allowed
            ic("Barcode duplication detected! Check your program run design")
            ic("Remove duplicated Barcodes... only the first one will be kept. ")
            barcode_df.drop_duplicates(subset=["Gene"], keep="first", inplace=True)

        ic("Barcode extraction initiated...")

        for ntp_barcode in tqdm(barcode_df.itertuples()):
            gene = ntp_barcode.Gene
            barcodes = [
                f"{getattr(ntp_barcode, key)}"  # TODO: Regex grouping needed
                for key in ntp_barcode._fields
                if "Barcode" in key
            ]

            query_result = sequence_frame["Sequence"].str.contains(
                ".*".join(barcodes), regex=True  # Xarr optimization needed
            )  # TODO : This part should be optimized with n-dimensional array (Xarr)

            sequence_frame[gene] = query_result

        sequence_frame.to_parquet(
            f"{result_dir}/parquets/{chunk_number}",
            compression="snappy",
            engine="pyarrow",
            compute=True,
            write_index=True,
            write_metadata_file=True,
        )
        end_time = time.time()
        ic(f"Barcode extraction finished... in {end_time-start_time} seconds")

        return sequence_frame
    except Exception as e:
        ic(e)
        ic(traceback.format_exc())

        return -1
