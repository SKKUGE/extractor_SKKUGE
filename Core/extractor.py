#!/usr/bin/python
# -*- coding: utf-8 -*-

# Original code is here, https://github.com/MyungjaeSong/Paired-Library.git
# This is the modified version of the program for academic uses in SKKUGE Lab

__author__ = "forestkeep21@naver.com"
__editor__ = "poowooho3@g.skku.edu"


import time
import traceback

from icecream import ic


def extractor_main(
    sequence_frame, gene, barcode, logger, result_dir, sep, chunk_number=0
):
    start_time = time.time()
    try:
        if chunk_number is None:
            raise ValueError("chunk_number is not defined")

        # ic(f"Barcode extraction initiated...{chunk_number}")

        query_result = sequence_frame["Sequence"].str.contains(barcode, regex=True)
        sequence_frame[gene] = query_result

        # Reduce sparsity by dropping undetected barcode columns
        # sequence_frame["Undetected"] = False
        # sequence_frame["Undetected"] = sequence_frame["Undetected"].mask(
        #     sequence_frame.iloc[:, 2:].sum(axis=1) == 0, True
        # )

        # sequence_frame = sequence_frame.loc[sequence_frame["Undetected"] != True].drop(
        #     columns=["Undetected", "Sequence"]
        # )
        sequence_frame = sequence_frame.drop(columns=["Sequence"])

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
