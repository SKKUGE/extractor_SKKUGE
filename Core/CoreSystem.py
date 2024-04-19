import os
import pathlib
import re
import subprocess as sp
import sys
from datetime import datetime
from types import SimpleNamespace

from psutil import cpu_count

DEBUGGING = False
N_PHYSICAL_CORES = cpu_count(logical=False)

os.environ["MALLOC_TRIM_THRESHOLD_"] = "65536"
import pandas as pd
from dask import config as dcf

dcf.set(
    {
        "dataframe.convert-string": True,
        # "dataframe.backend": "cudf",
    }
)
import ctypes

# import dask_cudf
import dask
from dask import bag as db
from dask import dataframe as dd
from dask.distributed import Client, LocalCluster
from icecream import ic
from tqdm import tqdm

ic.configureOutput(prefix=f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} | ")
ic.configureOutput(includeContext=True)


def trim_memory() -> int:
    libc = ctypes.CDLL("libc.so.6")
    return libc.malloc_trim(0)


def match_barcode(df) -> pd.DataFrame:
    df["match"] = df.apply(lambda x: bool(re.search(x["query"], x["sequence"])), axis=1)

    return df[df["match"] == True]


def finalize(
    combined_extraction_datasets,
    sample,
    barcode,
):
    try:

        ic("Remove ambiguous reads...")
        combined_extraction_datasets = combined_extraction_datasets.set_index(
            ["ID", "Sequence"]
        )
        combined_extraction_datasets = combined_extraction_datasets[
            combined_extraction_datasets.sum(axis=1, numeric_only=True) <= 2
        ]  # Remove ambiguous sequences

        ic("Calculating read counts...")
        combined_extraction_datasets = combined_extraction_datasets.sum(axis=0)
        ic(f"{sample}+{barcode}: Summation completed.")
        return combined_extraction_datasets
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
            logger.warning(
                f"Input folder: {len(list(proj_path.glob('*')))}, Project list samples: {len(loaded_samples)}"
            )

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
        if not DEBUGGING:
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
        if not DEBUGGING:
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
    ic("Initilaizing local cluster...")

    cpu_cluster = LocalCluster(
        processes=True,
        n_workers=N_PHYSICAL_CORES,  # DEBUG
        threads_per_worker=2,
        dashboard_address=":40928",
        local_directory="./tmp",
    )
    # gpu_cluster = LocalCUDACluster(
    #     dashboard_address=":40929",
    #     local_directory="./tmp",
    # )

    cpu_client = Client(cpu_cluster)
    # gpu_client = Client(gpu_cluster)

    ic(cpu_client)
    ic(cpu_client.dashboard_link)
    # ic(gpu_client)
    # ic(gpu_client.dashboard_link)

    for sample, barcode_path in tqdm(args.samples):
        ExtractorRunner(
            sample, barcode_path, args
        )  # TODO: refactor its usage to avoid creating an object
        if not DEBUGGING:
            ic("Loading fastq...")
            bag = db.read_text(
                args.system_structure.input_file_organizer[sample], blocksize="64MiB"
            ).map(lambda x: x.strip())
            bag = cpu_client.scatter(bag).result()

            ic("Transforming fastq...")
            sequence_ddf = bag.to_dataframe()
            sequence_ddf = cpu_client.scatter(sequence_ddf).result()
            sequence_da = sequence_ddf.to_dask_array(lengths=True)
            sequence_da = sequence_da.reshape(-1, 4)
            sequence_da = cpu_client.scatter(sequence_da).result()

            ic("Repartitioning...")
            sequence_ddf = sequence_da.to_dask_dataframe(
                columns=["id", "sequence", "separator", "quality"],
            )
            sequence_ddf = sequence_ddf.drop(
                columns=["separator", "quality"]
            )  # drop quality sequence
            sequence_ddf = sequence_ddf.repartition(partition_size="4KiB")

            ic("Save raw NGS reads as parquets...")
            sequence_ddf.to_parquet(
                f"{args.system_structure.seq_split_dir}",
                engine="pyarrow",
                write_index=True,
                write_metadata_file=True,
                compute=True,
            )

            # Load barcode file
            # The barcode writing rule: 1st-column should be name, and the others are n substrings for joint matching
            ic("Loading barcode file...")
            barcode_df = pd.read_csv(
                barcode_path,
                sep=args.sep,
                header=None,
            ).fillna("")
            barcode_df.columns = ["gene"] + [
                f"barcode_{i}" for i in range(1, len(barcode_df.columns))
            ]
            # TODO: Gene uniqueness test
            if not barcode_df["gene"].is_unique:
                ic("Gene names are not unique")
                ic(barcode_df["gene"].duplicated())
                ic("Drop duplicated gene names...")
                barcode_df = barcode_df.drop_duplicates(subset="Gene")

            barcode_df.loc[
                :, ["barcode" in col for col in barcode_df.columns]
            ].transform(lambda x: x.str.upper())
            barcode_df["query"] = barcode_df.loc[
                :, ["barcode" in col for col in barcode_df.columns]
            ].apply(lambda x: ".*".join(x), axis=1)
            # barcode_df = cudf.from_pandas(barcode_df)

            ic("Parquet generation completed, load parquets")
            sequence_ddf = dd.read_parquet(
                f"{args.system_structure.seq_split_dir}",
                # engine="pyarrow",
                calculate_divisions=True,
            )

            # TODO: cartesian product of barcode and gene and calculate it
            ic("Calculating cartesian product...")
            sequence_ddf["join_key"] = 1
            barcode_df["join_key"] = 1

            delayed_sequence_ddf = dask.delayed(sequence_ddf.merge)(
                barcode_df[["gene", "query", "join_key"]], on="join_key", how="left"
            ).drop("join_key", axis=1)
            sequence_ddf = delayed_sequence_ddf.compute()
            sequence_ddf.to_parquet(
                f"{args.system_structure.result_dir}/full_matrix/",
                write_index=True,
                write_metadata_file=True,
                engine="pyarrow",
                compute=True,
            )
        # END: DEBUGGING

        sequence_ddf = dd.read_parquet(
            f"{args.system_structure.result_dir}/full_matrix/",
            engine="pyarrow",
            calculate_divisions=True,
        )
        # ic("Repartitioning...")
        # sequence_ddf = sequence_ddf.repartition(partition_size="16MiB")
        # sequence_ddf.to_parquet(
        #     f"{args.system_structure.result_dir}/full_matrix_repartitioned/",
        #     engine="pyarrow",
        #     write_index=True,
        #     write_metadata_file=True,
        #     compute=True,
        # )
        # sequence_ddf = dd.read_parquet(
        #     f"{args.system_structure.result_dir}/full_matrix_repartitioned/",
        #     engine="pyarrow",
        #     calculate_divisions=True,
        # )

        ic("Submitting extraction process...")
        sequence_ddf["match"] = False
        sequence_ddf = sequence_ddf.map_partitions(
            match_barcode,
            meta=dd.utils.make_meta(sequence_ddf),  # TODO: Remove Falses in situ?
        )

        sequence_ddf.to_parquet(
            f"{args.system_structure.result_dir}/parquets/",
            compression="snappy",
            engine="pyarrow",
            compute=True,
            write_index=True,
            write_metadata_file=True,
        )

        ic("Pivotting extraction results...")
        sequence_ddf = dd.read_parquet(
            f"{args.system_structure.result_dir}/parquets/",
            engine="pyarrow",
            calculate_divisions=True,
        )
        sequence_ddf["gene"] = sequence_ddf["gene"].astype("category")
        sequence_ddf["gene"] = sequence_ddf["gene"].cat.as_known()
        pivot_sequence_ddf = sequence_ddf.pivot_table(
            index="id",
            columns="gene",
            values="match",
        )

        pivot_sequence_ddf.columns = pivot_sequence_ddf.columns.astype(str)
        pivot_sequence_ddf = pivot_sequence_ddf[
            pivot_sequence_ddf.sum(axis=1, numeric_only=True)
            <= args.dup_threshold  # Filter out ambiguous reads
        ]

        result_series = pivot_sequence_ddf.sum(axis=0)
        result_series = result_series.compute().T
        barcode_df = barcode_df.drop(columns=["join_key"]).set_index("gene")
        barcode_df["read_count"] = result_series
        barcode_df["read_count"] = barcode_df["read_count"].fillna(0)
        barcode_df.to_csv(f"{args.system_structure.result_dir}/read_counts.csv")

    ic("All merging jobs fired. Waiting for the final result...")
