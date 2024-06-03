import multiprocessing as mp
import os
import pathlib
import subprocess as sp
import sys
import traceback
from types import SimpleNamespace

os.environ["MALLOC_TRIM_THRESHOLD_"] = "65536"
import dask
import pandas as pd

dask.config.set({"dataframe.query-planning": False})
import ctypes

from dask import bag as db
from dask import dataframe as dd
from dask import delayed
from dask.diagnostics import ProgressBar
from dask.distributed import Client, LocalCluster, fire_and_forget, wait
from icecream import ic
from tqdm import tqdm


def trim_memory() -> int:
    libc = ctypes.CDLL("libc.so.6")
    return libc.malloc_trim(0)


def merge_parquets(
    args,
    rvals,
    sample,
    barcode,
):
    try:
        args.logger.info("Merging parquet files...")

        all_extraction_delayed_datasts = []
        for file in rvals:
            all_extraction_delayed_datasts.append(
                delayed(dd.read_parquet)(
                    file,
                    engine="pyarrow",
                    calculate_divisions=True,
                )
            )
        combined_extraction_datasets = delayed(dd.concat)(
            all_extraction_delayed_datasts, axis=1, interleave_partitions=True
        ).drop(columns=["ID"])
        # d = combined_extraction_datasets.compute() # DEBUG
        ic("Remove ambiguous reads...")
        combined_extraction_datasets = combined_extraction_datasets[
            combined_extraction_datasets.sum(axis=1, numeric_only=True) <= 2
        ]  # Remove ambiguous sequences
        combined_extraction_datasets = combined_extraction_datasets.sum(axis=0)
        combined_extraction_datasets.visualize(
            filename=f"{args.system_structure.result_dir}/read_counts.png"
        )

        ic(f"{sample}+{barcode}: Extraction future generated.")
        combined_extraction_datasets.compute().to_csv(
            f"{args.system_structure.result_dir}/read_counts.csv",
            index=True,
            single_file=True,
        )

        del all_extraction_delayed_datasts, combined_extraction_datasets

        ic(f"{sample}+{barcode}: Final read count table was generated.")
        return 0

    except Exception as e:
        ic(e)
        ic(traceback.format_exc())
        args.logger.error(e)
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
                args.logger.info(f"File name : {file_path.stem}")
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

        args.logger.info("System structure check : User, Project, Input, Output")
        args.multicore = os.cpu_count() if args.multicore == 0 else args.multicore
        if os.cpu_count() < args.multicore:
            args.logger.warning(
                f"Optimal threads <= {mp.cpu_count()} : {args.multicore} is not recommended"
            )
        for key, value in sorted(vars(args).items()):
            args.logger.info(f"Argument {key}: {value}")

        args.logger.info("File num check: input folder and project list")
        Helper.equal_num_samples_checker(
            args.system_structure.input_dir, args.samples, args.logger
        )

        return func(args)

    return wrapper


@system_struct_checker
def run_pipeline(args: SimpleNamespace) -> None:
    for sample, barcode in args.samples:
        sample = Helper.SplitSampleInfo(sample)

        extractor_runner = ExtractorRunner(sample, barcode, args)

        # Chunking
        args.logger.info("Splitting sequecnes into chunks")
        extractor_runner._split_into_chunks()

        args.logger.info("Populating command...")
        listCmd = extractor_runner._populate_command(barcode)

        args.logger.info("RunMulticore")

        # Refactor this block of code for flushing memory
        run_extractor_mp(
            listCmd,
            args.multicore,
            args.logger,
            args.verbose,
            args.system_structure.result_dir,
            sample,
        )
        sp.run(
            [
                "rm",
                "-r",
                f"{pathlib.Path.cwd() / args.system_structure.seq_split_dir}",
            ]
        )


def run_extractor_mp(
    lCmd, iCore, logger, verbose_mode: bool, result_dir: pathlib.Path, sample_name
) -> None:
    import gc
    import time

    import dask.dataframe as dd
    import numpy as np
    from tqdm import tqdm

    from extractor import main as extractor_main

    for sCmd in lCmd:
        logger.info(f"Running {sCmd} command with {iCore} cores")

    result = []
    start = time.time()
    with ProcessPoolExecutor(max_workers=iCore) as executor:
        for rval in list(tqdm(executor.map(extractor_main, lCmd), total=len(lCmd))):
            result.append(rval)
    end = time.time()
    logger.info(f"Extraction is done. {end - start}s elapsed.")

    logger.info("Generating statistics...")

    with open(f"{result_dir}/{sample_name}+read_statstics.txt", "w") as f:
        read_stat = np.concatenate([rval for rval in result], axis=0)
        detected, total_read, detection_rate = (
            read_stat.sum(),
            read_stat.shape[0],
            read_stat.sum() / read_stat.shape[0],
        )
        f.write(f"Total read: {total_read}\n")
        f.write(f"Detected read: {detected}\n")
        f.write(f"Detection rate in the sequence pool: {detection_rate}\n")

    logger.info("Generating final extraction results...")

    # TODO: asynchronous merging of parquet files
    # load multiple csv files into one dask dataframe
    # TODO : Refactor this block of code
    parquets = []
    for f in pathlib.Path(f"{result_dir}/parquets").glob("*.parquet"):
        d_parquet = dd.read_parquet(f)
        d_parquet["n_ids"] = d_parquet["ID"].apply(len, meta=("ID", "int64"))
        d_parquet = d_parquet.explode("ID")
        d_parquet["Read_counts"] = (
            d_parquet["Read_counts"] / d_parquet["n_ids"]
        )  # mutiplied and divided by the number of IDs
        parquets.append(d_parquet)
    df = dd.concat(parquets)

    # DEBUG
    df.compute().to_csv(f"{result_dir}/test.csv", index=False)

    df["RPM"] = df["Read_counts"] / df["Read_counts"].sum() * 1e6

    df.drop(["ID", "n_ids"], axis=1).groupby(
        ["Gene", "Barcode"]
    ).sum().compute().to_csv(
        f"{result_dir}/{sample_name}+extraction_result.csv", index=True
    )  # Fetch the original Read count from n_ids
    # TODO: refactor this block of code

    if verbose_mode:
        # Create NGS_ID_classification.csv

        df.drop(["Read_counts"], axis=1).dropna(subset=["ID"]).set_index(
            "ID"
        ).compute().to_csv(
            f"{result_dir}/{sample_name}+multiple_detection_test_result.csv", index=True
        )
        # Create Barcode_multiple_detection_test.csv
        df.groupby(["ID"])["Barcode"].count().compute().to_csv(
            f"{result_dir}/{sample_name}+multiple_detection_test_by_id.csv"
        )

        # Create statistics for analysis
        gc.collect()

    return
