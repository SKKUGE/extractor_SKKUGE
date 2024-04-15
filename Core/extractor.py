#!/usr/bin/python
# -*- coding: utf-8 -*-

# Original code is here, https://github.com/MyungjaeSong/Paired-Library.git
# This is the modified version of the program for academic uses in SKKUGE Lab

__author__ = "forestkeep21@naver.com"
__editor__ = "poowooho3@g.skku.edu"


import time
import traceback

from icecream import ic


def extractor_main(sequence_frame, ntp_barcode, logger, result_dir, chunk_number=0):
    start_time = time.time()
    try:

        gene = ntp_barcode.Gene
        barcodes = [
            f"{getattr(ntp_barcode, key)}"  # TODO: Regex grouping needed
            for key in ntp_barcode._fields
            if "Barcode" in key
        ]

        query_result = sequence_frame["Sequence"].str.contains(
            ".*".join(barcodes), regex=True  # Xarr optimization needed
        )  # TODO : This part should be optimized with n-dimensional array (Xarr)

        sequence_frame = sequence_frame.drop(columns=["Sequence"])
        sequence_frame[gene] = query_result
        sequence_frame = sequence_frame[sequence_frame[gene] == True]
        sequence_frame.to_parquet(
            f"{result_dir}/parquets/{chunk_number}",
            compression="snappy",
            engine="pyarrow",
            compute=True,
            write_index=True,
            write_metadata_file=True,
        )
        end_time = time.time()
        ic(
            f"Barcode extraction finished...{chunk_number} in {end_time-start_time} seconds"
        )

        return f"{result_dir}/parquets/{chunk_number}"

    except Exception as e:
        ic(e)
        ic(traceback.format_exc())
        logger.error(e)
        return -1
