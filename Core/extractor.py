#!/usr/bin/python
# -*- coding: utf-8 -*-

# Original code is here, https://github.com/MyungjaeSong/Paired-Library.git
# This is the modified version of the program for academic uses in SKKUGE Lab

__author__ = "forestkeep21@naver.com"
__editor__ = "poowooho3@g.skku.edu"


import re
import time
import traceback

from icecream import ic


def extractor_main(
    sequence_frame,
    barcode_df,
    result_dir,
    sep=",",
    logger=None,
    chunk_number=None,
    generosity=2,
):
    start_time = time.time()
    try:
        if chunk_number is None:
            raise ValueError("chunk_number is not defined")

        # pbar = tqdm(total=barcode_df.shape[0])
        # pbar.set_description(f"Barcode extraction...{chunk_number}")

        # TODO: Merge-based barcode extraction
        sequence_frame["join_key"] = 1
        barcode_df["join_key"] = 1
        sequence_frame = sequence_frame.merge(
            barcode_df, on="join_key", how="left"
        ).drop("join_key", axis=1)

        sequence_frame["match"] = sequence_frame.apply(
            lambda x: bool(re.search(x.query, x.sequence)), axis=1
        )

        # Legacy code: iteration based barcode extraction
        # for ntp in barcode_df.itertuples():
        #     gene = ntp.gene
        #     barcode = ntp.query
        #     sequence_frame[gene] = sequence_frame["sequence"].str.contains(
        #         barcode, regex=True
        #     )
        #     pbar.update(1)

        # pbar.close()

        # Drop reads that do not contain any barcode
        # sequence_frame = sequence_frame[sequence_frame["match"] > 0]
        sequence_frame["gene"] = sequence_frame["gene"].astype("category")
        sequence_frame["gene"] = sequence_frame["gene"].cat.as_known()
        pivot_sequence_frame = sequence_frame.pivot_table(
            index="id",
            columns="gene",
            values="match",
        )
        pivot_sequence_frame["detection_count"] = pivot_sequence_frame.sum(
            axis=1, numeric_only=True
        ).astype(int)
        pivot_sequence_frame = pivot_sequence_frame[
            pivot_sequence_frame.sum(axis=1) <= args.dup_threshold
        ]  # Remove reads with more than 2 barcodes
        pivot_sequence_frame = pivot_sequence_frame.drop(columns=["detection_count"])

        pivot_sequence_frame.to_parquet(
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

        del pivot_sequence_frame, sequence_frame, barcode_df
        return 0

    except Exception as e:
        ic(traceback.format_exc())
        ic(e)

        return -1
