# Original code is here, https://github.com/MyungjaeSong/Paired-Library.git
# This is the modified version of the program for academic uses in SKKUGE Lab

__author__ = "forestkeep21@naver.com"
__editor__ = "poowooho3@g.skku.edu"

import pathlib

import pandas as pd
import polars as pl
from pyarrow.dataset import dataset
from tqdm import tqdm


def pbar_annotation(pbar, func):
    def foo(*args, **kwargs):
        pbar.update(1)
        return func(*args, **kwargs)

    return foo


def extract_read_cnts(
    sequence_file_path: pathlib.Path,  # path-to-parquet
    barcode_file_path: pathlib.Path,  # path-to-barcode
    result_file_path: pathlib.Path,  # path-to-result
    sep=",",  # separator in barcode file
    sample_replacement_mode=False,
    logger=None,
):
    try:
        pyarrow_dataset = dataset(
            source=sequence_file_path,
            format="parquet",
        )
        seq_ldf = pl.scan_pyarrow_dataset(pyarrow_dataset)
        n_seqs = seq_ldf.select(pl.count()).collect(streaming=True)

        # Load barcode file
        try:
            barcode_table = pd.read_csv(
                barcode_file_path,
                sep=sep,
                header=None,
                dtype=str,
                engine="c",
                names=["NAME", "BARCODE"],
            )
        except Exception:
            if logger:
                logger.error("Error occured while reading barcode file")
            raise Exception("Error occured while reading barcode file")
        barcode_table.apply(
            lambda x: x.astype(str).str.upper()
        )  # Don't consider soft masking

        if not barcode_table["BARCODE"].is_unique:
            # Barcode used as a PK in the database, so duplication is not allowed
            if logger:
                logger.info("Barcode duplication detected: ")
                logger.info(
                    "Remove duplicated Barcodes... only the first one will be kept."
                )
            barcode_table.drop_duplicates(
                subset=["BARCODE"], keep="first", inplace=True
            )

        # result_df["Barcode_copy"] = result_df["Barcode"]
        # result_df = result_df.set_index("Barcode")  # TODO: tentative design
        # Load a split sequencing result using high-level I/O; validating fastq format
        # result_df["Read_counts"] = 0
        # result_df["ID"] = ""

        # seqs = skbio.io.read(
        #     sequence_file, format="fastq", verify=True, variant="illumina1.8"
        # )  # FASTQ format verification using skbio

        # seq_df = pd.DataFrame(
        #     [(seq.metadata["id"], seq._string.decode()) for seq in seqs],
        #     columns=["ID", "Sequence"],
        # )
        # seq_detection_array = np.zeros(seq_df.shape[0])

        # TODO: extract function

        # barcode_table["join"] = 1
        # barcode_table["join"] = barcode_table["join"].astype("int32")

        # seq_ldf = seq_ldf.with_columns(pl.lit(1).alias("join"))

        # result_ldf = barcode_ldf.join(seq_ldf, on="join", how="inner").drop("join")
        # result_ldf = result_ldf.with_columns(
        #     pl.col("SEQUENCE").str.contains(pl.col("BARCODE")).alias("MATCH")
        # )

        pathlib.Path.mkdir(sequence_file_path / "intermediate", exist_ok=True)
        for i, pat in enumerate(tqdm(barcode_table.BARCODE)):
            seq_ldf.with_columns(
                pl.col("SEQUENCE")
                .str.extract(
                    "(" + pat + ")",
                )
                .alias("BARCODE")
            ).collect(streaming=True).write_parquet(
                sequence_file_path / "intermediate" / f"{i}.parquet",
            )
            

        # result_ldf.groupby("BARCODE").agg(
        #     pl.sum("MATCH").alias("Read_counts"),
        #     pl.col("ID").implode().alias("ID"),
        # )

        # seq_df = seq_df.with_columns(
        #     pl.col("SEQUENCE")
        #     .str.extract(
        #         "(" + pat + ")",
        #     )
        #     .alias("BARCODE"),
        # )

        # TODO: pbar
        # for idx, row in result_df.iterrows():
        #     # query_result format
        #     #

        #     query_result = seq_df["Sequence"].str.contains(row["Barcode_copy"])
        #     # boolean indexing for fast processing
        #     result_df.at[idx, "ID"] = (
        #         seq_df.loc[query_result[query_result].index]["ID"].to_numpy().tolist()
        #     )
        #     result_df.loc[idx, "Read_counts"] = (query_result.sum(),)

        #     seq_detection_array[query_result[query_result].index] = True

        #     # TODO: Sample with replacement option
        #     # Without replacement from the sequence pool
        #     if not sample_replacement_mode:
        #         seq_df.drop(query_result[query_result].index, inplace=True, axis=0)

        # del seq_df
        # gc.collect()

        total_result_ldf = pl.scan_pyarrow_dataset(
            dataset(
                source=sequence_file_path / "intermediate",
                format="parquet",
            )
        )

        pathlib.Path.mkdir(result_file_path / sequence_file_path.name, exist_ok=True)
        total_result_ldf.collect(streaming=True).write_csv(
            result_file_path / sequence_file_path.name / "total_result.csv"
        )
        read_count_ldf = total_result_ldf.groupby("BARCODE").agg(
            pl.sum("BARCODE").alias("READ_COUNTS"),
            pl.col("ID").implode().alias("ID"),
        )

        read_count_ldf.collect(streaming=True).write_csv(
            result_file_path / sequence_file_path.name / "read_count.csv"
        )
        pass

        return 0

    except Exception as e:
        return e

    # result_df.drop("Barcode_copy", axis=1, inplace=True)
    # result_df.reset_index(inplace=True, drop=False)
    # result_df.iloc[1], result_df.iloc[-1] = result_df.iloc[-1], result_df.iloc[1]

    # def name():
    #     from datetime import datetime

    #     dt_string = datetime.now().strftime("%Y-%m-%d;%H:%M:%S")
    #     return str(f"{dt_string}")

    # result_df.to_parquet(
    #     f"{result_dir}/{name()}+{pathlib.Path(sequence_file).name}.parquet"
    # )

    # del result_df
    # gc.collect()

    # return seq_detection_array  # TODO: parquet을 읽으면 되는 것 아닌가?


def main(*args) -> pd.DataFrame:
    (sequence, barcode, result_dir, sep, sample_replacement_mode) = args[0]
    rval = extract_read_cnts(
        sequence, barcode, result_dir, sep, sample_replacement_mode
    )
    return rval
