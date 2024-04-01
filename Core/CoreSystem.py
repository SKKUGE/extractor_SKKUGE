import asyncio
import os
import pathlib
import subprocess as sp
import sys
import traceback
from types import SimpleNamespace

from psutil import cpu_count

N_PHYSICAL_CORES = cpu_count(logical=False)

os.environ["MALLOC_TRIM_THRESHOLD_"] = "65536"
from dask import config as dcf

dcf.set({"dataframe.query-planning": False})
import ctypes

import pandas as pd
import polars as pl
from dask import bag as db
from dask import dataframe as dd
from dask.distributed import Client, LocalCluster
from icecream import ic
from tqdm import tqdm


def trim_memory() -> int:
    libc = ctypes.CDLL("libc.so.6")
    return libc.malloc_trim(0)


def indel_analysis():
    ic("Indel analysis... (Not implemented yet)")
    # Use CRISPREsso API

    return 0


def binary_tree_merge(dataframes):
    if len(dataframes) == 1:
        return dataframes[0]
    mid = len(dataframes) // 2
    left = binary_tree_merge(dataframes[:mid])
    right = binary_tree_merge(dataframes[mid:])
    return left.join(right, on="ID", how="outer")


def merge_parquets(
    args,
    rvals,
):
    try:
        ic("Merging parquet files...")

        all_extraction_delayed_datasts = []
        for file in rvals:
            lazy_fragmented_parquets = pl.scan_parquet(
                f"{file}/*.parquet",
            )
            all_extraction_delayed_datasts.append(lazy_fragmented_parquets)

        lazy_combined_extraction_datasets = binary_tree_merge(
            all_extraction_delayed_datasts
        )

        # Reduce sparsity
        lazy_combined_extraction_datasets = lazy_combined_extraction_datasets.select(
            pl.exclude("^.*_right$")
        )

        return lazy_combined_extraction_datasets

    except Exception as e:
        ic(e)
        ic(traceback.format_exc())
        args.logger.error(e)
        return -1


async def get_read_count_table(lzdf):
    lzdf = lzdf.drop(columns=["Sequence"])

    ic("Remove ambiguous reads...")
    lzdf = lzdf.with_columns(ambiguity=pl.sum_horizontal(*lzdf.columns[1:]))

    lzdf = lzdf.filter(pl.col("ambiguity") <= 2)  # Remove ambiguous sequences
    lzdf = lzdf.drop("ambiguity", "ID")

    ic("Calculating read counts...")
    lzdf = lzdf.sum()

    lzdf = lzdf.collect_async(streaming=True)

    return await lzdf


async def finalize(
    lazy_combined_extraction_datasets,
    args,
    sample,
    barcode,
):
    try:
        ic("Finalizing extraction...")
        # raw_sequences = combined_extraction_datasets["Sequence"]

        lazy_combined_extraction_datasets = await get_read_count_table(
            lazy_combined_extraction_datasets
        )

        combined_extraction_datasets = lazy_combined_extraction_datasets.to_pandas().T
        combined_extraction_datasets.index.name = "Gene"
        combined_extraction_datasets.columns = ["Reads"]

        combined_extraction_datasets = combined_extraction_datasets.join(
            pd.read_csv(
                barcode,
                sep=args.sep,
                header=None,
                names=["Gene", "Barcode"],
            )
            .set_index("Gene")
            .drop(columns=["Barcode"]),
            how="outer",
        ).fillna(0)

        combined_extraction_datasets.to_csv(
            f"{args.system_structure.result_dir}/read_counts.csv",
            sep=",",
            header=True,
            index=True,
        )

        ic(f"{sample}+{barcode}: Final read count table was generated.")

        return 0
    except Exception as e:
        ic(e)
        return -1


class Helper(object):
    @staticmethod
    def mkdir_if_not(directory: pathlib.Path) -> pathlib.Path:
        """
        > If the directory doesn't exist, create it

        :param directory: The directory to create
        :type directory: str
        :return: A path object
        """
        directory.mkdir(parents=True, exist_ok=True)
        return directory

    @staticmethod
    def load_samples(directory: pathlib.Path) -> list:
        """
        It reads a file and returns a list of non-empty lines that don't start with a hash mark.

        :param directory: the directory of the samples file
        :type directory: pathlib.Path
        :return: A list of samples.
        """
        with open(directory, "r", encoding="utf-8") as file:
            lines = [line.strip("\n") for line in file.readlines()]
            sample_barcode_list = [
                line.split(",")[:2] for line in lines if line[0] != "#"
            ]

        return sample_barcode_list

    @staticmethod  # defensive
    def equal_num_samples_checker(
        proj_path: pathlib.Path, loaded_samples: list, logger
    ):
        """
        > This function checks if the number of samples in the Input folder \
            and in the User folder
        matches

        :param proj_path: pathlib.Path, loaded_samples: list, logger
        :type proj_path: pathlib.Path
        :param loaded_samples: a list of sample names
        :type loaded_samples: list
        :param logger: a logger object
        """

        if len(list(proj_path.glob("*"))) != len(loaded_samples):
            logger.warning(
                "The number of samples in the Input folder and in the User folder does not matched. Check the file list in the Input folder and the project list in the User folder."
            )

            # input_entries = [i.name for i in proj_path.glob("*")]
            # user_entries = [i for i in loaded_samples]
            logger.warning(
                f"Input folder: {len(list(proj_path.glob('*')))}, Project list samples: {len(loaded_samples)}"
            )
            # logger.warning(
            #     f"Input folder: {[i for i in input_entries if i not in user_entries]}"
            # )
            # logger.warning(
            #     f"Project list samples: {[u for u in user_entries if u not in input_entries]}"
            # )
        else:
            logger.info("The file list is correct, pass\n")


class SystemStructure(object):
    def __init__(
        self,
        user_name: str,
        project_name: str,
        base_dir: pathlib.Path = pathlib.Path.cwd() / "Data",
    ):
        from collections import defaultdict
        from typing import DefaultDict

        # https://www.notion.so/dengardengarden/s-Daily-Scrum-Reports-74d406ce961c4af78366a201c1933b66#cd5b57433eca4c6da36145d81adbbe5e
        self.user_name = user_name
        self.project_name = project_name
        self.base_dir = base_dir
        self.input_sample_organizer: DefaultDict[str, pathlib.Path] = defaultdict(
            pathlib.Path
        )
        self.input_file_organizer: DefaultDict[str, pathlib.Path] = defaultdict(
            pathlib.Path
        )
        self.output_sample_organizer: DefaultDict[str, pathlib.Path] = defaultdict(
            pathlib.Path
        )

        self.user_dir = Helper.mkdir_if_not(
            self.base_dir / f'{"User" + "/" + self.user_name}'
        )
        self.barcode_dir = Helper.mkdir_if_not(self.base_dir / "Barcodes")
        self.input_dir = Helper.mkdir_if_not(
            self.base_dir
            / f'{"Input" + "/" + self.user_name + "/" + self.project_name}'
        )
        self.project_samples_path = pathlib.Path(
            self.user_dir / f"{self.project_name}.txt"
        )
        if not self.project_samples_path.exists():
            with open(self.project_samples_path, "w", encoding="utf-8") as f:
                f.write("# Sample,Barcode\n")

        self.output_dir = Helper.mkdir_if_not(
            self.base_dir
            / f'{"Output" + "/" + self.user_name + "/" + self.project_name}'
        )  # TODO is it needed?

    def mkdir_sample(self, sample_name: str, barcode_name: str):
        # TODO
        self.input_sample_organizer[sample_name] = Helper.mkdir_if_not(
            self.input_dir / sample_name
        )
        self.output_sample_organizer[sample_name] = Helper.mkdir_if_not(
            self.output_dir / barcode_name / sample_name
        )
        self.result_dir = Helper.mkdir_if_not(self.output_sample_organizer[sample_name])
        self.parquet_dir = Helper.mkdir_if_not(self.result_dir / "parquets")
        self.full_mat_dir = Helper.mkdir_if_not(self.result_dir / "full_matrix")

        if len(os.listdir(f"{pathlib.Path.cwd() / self.parquet_dir}")) > 0:
            sp.run(
                [
                    "rm",
                    "-r",
                    f"{self.result_dir / 'parquets'}",
                ]
            )
            self.parquet_dir = Helper.mkdir_if_not(
                self.result_dir / "parquets"
            )  # Re-create the directory


# TODO: FLASh integration? https://ccb.jhu.edu/software/FLASH/
def fastp_integration():
    pass


class ExtractorRunner:
    def __init__(self, sample: str, barcode: str, args: SimpleNamespace):
        args.python = sys.executable
        # Find python executable if not specified
        args.system_structure.mkdir_sample(sample, pathlib.Path(barcode).name)
        self.sample = sample
        self.args = args

        for idx, file_path in enumerate(
            [
                p
                for p in self.args.system_structure.input_sample_organizer[
                    self.sample
                ].glob("*")
            ]
        ):
            # Load input file from input sample folder (only one file)
            if file_path.suffix in [".fastq", ".fq", ".fastq.gz", ".fq.gz"]:
                ic(
                    f"Sample name : {args.system_structure.input_sample_organizer[self.sample]}"
                )
                ic(f"File name : {file_path.stem}")
                self.args.system_structure.input_file_organizer[self.sample] = (
                    pathlib.Path.cwd() / file_path
                )
                break

            if (
                idx
                == len(
                    [
                        p
                        for p in self.args.system_structure.input_sample_organizer[
                            self.sample
                        ].glob("*")
                    ]
                )
                - 1
            ):
                raise Exception("No fastq file in the sample folder")

        # self.strInputList  => contains all splitted fastq file path; glob can be used

        self.args.system_structure.seq_split_dir = Helper.mkdir_if_not(
            self.args.system_structure.input_sample_organizer[self.sample]
            / "Split_files"
        )

        if (
            len(
                os.listdir(
                    f"{pathlib.Path.cwd() / self.args.system_structure.seq_split_dir}"
                )
            )
            > 0
        ):
            sp.run(
                [
                    "rm",
                    "-r",
                    f"{pathlib.Path.cwd() / self.args.system_structure.seq_split_dir}",
                ]
            )
            self.args.system_structure.seq_split_dir = Helper.mkdir_if_not(
                self.args.system_structure.input_sample_organizer[self.sample]
                / "Split_files"
            )  # Re-create the directory


def system_struct_checker(func):
    def wrapper(args: SimpleNamespace):

        ic("System structure check : User, Project, Input, Output")
        args.multicore = os.cpu_count() if args.multicore == 0 else args.multicore
        if os.cpu_count() < args.multicore:
            args.logger.warning(
                f"Optimal threads <= {N_PHYSICAL_CORES} : {args.multicore} is not recommended"
            )
        for key, value in sorted(vars(args).items()):
            ic(f"Argument {key}: {value}")

        ic("File num check: input folder and project list")
        Helper.equal_num_samples_checker(
            args.system_structure.input_dir, args.samples, args.logger
        )

        return func(args)

    return wrapper


@system_struct_checker
def run_pipeline(args: SimpleNamespace) -> None:
    # TODO: add parquet remove option

    from Core.extractor import extractor_main as extractor_main

    ic("Initilaizing local cluster...")

    cluster = LocalCluster(
        processes=True,
        n_workers=N_PHYSICAL_CORES,  # DEBUG
        threads_per_worker=2,
        memory_limit="4GiB",
        dashboard_address=":40928",
    )
    client = Client(cluster)
    client.amm.start()

    ic(client)
    ic(client.dashboard_link)

    loop = asyncio.get_event_loop()
    polars_results = []
    for sample, barcode in tqdm(args.samples):
        ExtractorRunner(
            sample, barcode, args
        )  # TODO: refactor its usage to avoid creating an object

        ic("Loading merged fastq file...")
        bag = db.read_text(
            args.system_structure.input_file_organizer[sample], blocksize="64MiB"
        )
        sequence_ddf = bag.to_dataframe()

        sequence_ddf = (
            sequence_ddf.to_dask_array(lengths=True)
            .reshape(-1, 4)
            .to_dask_dataframe(
                columns=["ID", "Sequence", "Separator", "Quality"],
            )
        )
        sequence_ddf = sequence_ddf.drop(columns=["Separator", "Quality"]).repartition(
            "100MB"
        )  # drop quality sequence

        ic("Save NGS reads as parquets...")
        sequence_ddf.to_parquet(
            f"{args.system_structure.seq_split_dir}",
            engine="pyarrow",
            compute=True,
        )
        ic("Parquet generation completed, load parquets")

        sequence_ddf = dd.read_parquet(
            f"{args.system_structure.seq_split_dir}",
            engine="pyarrow",
            calculate_divisions=True,
        )
        # Load barcode file
        ic("Loading barcode file...")
        barcode_row_length = sum(1 for row in open(barcode, "r"))
        chunk_size = (
            barcode_row_length // (N_PHYSICAL_CORES - 1)
            if (barcode_row_length // (N_PHYSICAL_CORES - 1)) >= 1
            else barcode_row_length
        )

        barcode_df = pd.read_csv(
            barcode,
            sep=args.sep,
            header=None,
            names=["Gene", "Barcode"],
            chunksize=chunk_size,
        )
        ic("Submitting extraction process...")
        futures = []
        for i, barcode_chunk in enumerate(barcode_df):
            futures.append(
                client.submit(
                    extractor_main,
                    sequence_ddf,
                    barcode_chunk.iloc[
                        :, [0, 1]
                    ].dropna(),  # Use only Gene and Barcode columns
                    args.logger,
                    args.system_structure.result_dir,
                    args.sep,
                    chunk_number=i,
                )
            )
        # Gather results
        ic("Gathering extraction results...")
        rvals = client.gather(futures)
        for rval in rvals:
            try:
                if rval == -1:
                    ic(traceback.format_exc())
                    raise Exception(f"extractor_main has returned with {rval}")
            except ValueError:
                ic("Parquet generation completed")
                continue

        # BUG: Memory error
        lazy_combined_extraction_datasets = merge_parquets(
            args,
            rvals,
        )

        new_length = (
            lazy_combined_extraction_datasets.select(pl.len())
            .collect(streaming=True)
            .item()
        )
        original_length = sequence_ddf.shape[0].compute()
        ic(
            new_length,
            original_length,
            new_length == original_length,
        )
        assert new_length == original_length

        del bag, sequence_ddf

        client.run(trim_memory)

        ic("Save concatenated results as parquets...")

        lazy_combined_extraction_datasets.sink_parquet(
            path=f"{args.system_structure.full_mat_dir}/full_matrix.parquet",
        )

        ic("Full result parquet generation completed, load parquets")
        lazy_combined_extraction_datasets = pl.scan_parquet(
            source=f"{args.system_structure.full_mat_dir}/*.parquet",
        )

        polars_results.append(
            loop.create_task(
                finalize(
                    lazy_combined_extraction_datasets,
                    args=args,
                    sample=sample,
                    barcode=barcode,
                )
            )
        )
        indel_analysis()  # TODO: Indel analysis using the full matrix

    loop.run_until_complete(asyncio.gather(*polars_results))
    loop.close()
    ic("All merging jobs finished...")
    ic(polars_results)
