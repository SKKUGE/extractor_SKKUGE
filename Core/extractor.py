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
        if chunk_number is None:
            raise ValueError("chunk_number is not defined")

        # ic(f"Barcode extraction initiated...{chunk_number}")

        gene = ntp_barcode.Gene
        barcodes = [
            f"{getattr(ntp_barcode, key)}"  # TODO: Regex grouping needed
            for key in ntp_barcode._fields
            if "Barcode" in key
        ]

        query_result = sequence_frame["Sequence"].str.contains(
            ".*".join(barcodes), regex=True  # Xarr optimization needed
        )
        # query_result = (
        #     [  # TODO : This part should be optimized with n-dimensional array (Xarr)
        #         sequence_frame["Sequence"].str.contains(str(barcode), regex=False)
        #         for barcode in barcodes
        #     ]
        # )

        sequence_frame = sequence_frame.drop(columns=["Sequence"])
        sequence_frame[gene] = query_result
        # sequence_frame[gene] = True
        # for q in query_result:
        #     sequence_frame[gene] &= q

        sequence_frame.to_parquet(
            f"{result_dir}/parquets/{chunk_number}",
            compression="snappy",
            engine="pyarrow",
            compute=True,
            write_index=True,
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
