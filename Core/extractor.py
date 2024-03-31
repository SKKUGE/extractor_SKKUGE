#!/usr/bin/python
# -*- coding: utf-8 -*-

# Original code is here, https://github.com/MyungjaeSong/Paired-Library.git
# This is the modified version of the program for academic uses in SKKUGE Lab

__author__ = "forestkeep21@naver.com"
__editor__ = "poowooho3@g.skku.edu"


import gc
import pathlib
import traceback

from icecream import ic


def extractor_main(sequence_frame, barcode_df, logger, result_dir, sep, chunk_number):
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

            if i % (2**3) == 0:  # DEBUG
                # ic(i)
                sequence_frame = sequence_frame.persist()

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
        del sequence_frame
        gc.collect()
        logger.info("Barcode extraction completed")
        logger.info("Merging parquet files...")

        return f"{result_dir}/parquets/{chunk_number}"

    except Exception as e:
        ic(e)
        ic(traceback.format_exc())
        logger.error(e)
        return -1
