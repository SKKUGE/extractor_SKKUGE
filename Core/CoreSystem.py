import os
import pathlib
import subprocess as sp
import sys
from datetime import datetime
from types import SimpleNamespace

from psutil import cpu_count

DEBUGGING = True
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

import duckdb
from dask import bag as db
from dask.distributed import Client, LocalCluster
from icecream import ic
from tqdm import tqdm

ic.configureOutput(prefix=f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} | ")
ic.configureOutput(includeContext=True)


def trim_memory() -> int:
    libc = ctypes.CDLL("libc.so.6")
    return libc.malloc_trim(0)


def finalize(
    combined_extraction_datasets,
    sample,
    barcode,
):
    try:

        ic("Remove ambiguous reads...")
        combined_extraction_datasets = combined_extraction_datasets.set_index(["ID", "Sequence"])
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
            sample_barcode_list = [line.split(",")[:2] for line in lines if line[0] != "#"]

        return sample_barcode_list

    @staticmethod  # defensive
    def equal_num_samples_checker(proj_path: pathlib.Path, loaded_samples: list, logger):
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
        self.input_sample_organizer: DefaultDict[str, pathlib.Path] = defaultdict(pathlib.Path)
        self.input_file_organizer: DefaultDict[str, pathlib.Path] = defaultdict(pathlib.Path)
        self.output_sample_organizer: DefaultDict[str, pathlib.Path] = defaultdict(pathlib.Path)

        self.user_dir = Helper.mkdir_if_not(self.base_dir / f'{"User" + "/" + self.user_name}')
        self.barcode_dir = Helper.mkdir_if_not(self.base_dir / "Barcodes")
        self.input_dir = Helper.mkdir_if_not(
            self.base_dir / f'{"Input" + "/" + self.user_name + "/" + self.project_name}'
        )
        self.project_samples_path = pathlib.Path(self.user_dir / f"{self.project_name}.txt")
        if not self.project_samples_path.exists():
            with open(self.project_samples_path, "w", encoding="utf-8") as f:
                f.write("# Sample,Barcode\n")

        self.output_dir = Helper.mkdir_if_not(
            self.base_dir / f'{"Output" + "/" + self.user_name + "/" + self.project_name}'
        )  # TODO is it needed?

    def mkdir_sample(self, sample_name: str, barcode_name: str):
        # TODO
        self.input_sample_organizer[sample_name] = Helper.mkdir_if_not(self.input_dir / sample_name)
        self.output_sample_organizer[sample_name] = Helper.mkdir_if_not(self.output_dir / barcode_name / sample_name)
        self.result_dir = Helper.mkdir_if_not(self.output_sample_organizer[sample_name])
        self.parquet_dir = Helper.mkdir_if_not(self.result_dir / "parquets")
        # self.full_mat_dir = Helper.mkdir_if_not(self.result_dir / "full_matrix")
        # self.full_mat_dir = Helper.mkdir_if_not(self.result_dir / "full_matrix")
        if not DEBUGGING:
            if len(os.listdir(f"{pathlib.Path.cwd() / self.parquet_dir}")) > 0:
                sp.run(
                    [
                        "rm",
                        "-r",
                        f"{self.result_dir / 'parquets'}",
                    ]
                )
                self.parquet_dir = Helper.mkdir_if_not(self.result_dir / "parquets")  # Re-create the directory


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
            [p for p in self.args.system_structure.input_sample_organizer[self.sample].glob("*")]
        ):
            # Load input file from input sample folder (only one file)
            if file_path.suffix in [".fastq", ".fq", ".fastq.gz", ".fq.gz"]:
                ic(f"Sample name : {args.system_structure.input_sample_organizer[self.sample]}")
                ic(f"File name : {file_path.stem}")
                self.args.system_structure.input_file_organizer[self.sample] = pathlib.Path.cwd() / file_path
                break

            if idx == len([p for p in self.args.system_structure.input_sample_organizer[self.sample].glob("*")]) - 1:
                raise Exception("No fastq file in the sample folder")

        # self.strInputList  => contains all splitted fastq file path; glob can be used

        self.args.system_structure.seq_split_dir = Helper.mkdir_if_not(
            self.args.system_structure.input_sample_organizer[self.sample] / "Split_files"
        )
        if not DEBUGGING:
            if len(os.listdir(f"{pathlib.Path.cwd() / self.args.system_structure.seq_split_dir}")) > 0:
                sp.run(
                    [
                        "rm",
                        "-r",
                        f"{pathlib.Path.cwd() / self.args.system_structure.seq_split_dir}",
                    ]
                )
                self.args.system_structure.seq_split_dir = Helper.mkdir_if_not(
                    self.args.system_structure.input_sample_organizer[self.sample] / "Split_files"
                )  # Re-create the directory


def system_struct_checker(func):
    def wrapper(args: SimpleNamespace):

        ic("System structure check : User, Project, Input, Output")
        args.multicore = os.cpu_count() if args.multicore == 0 else args.multicore
        if os.cpu_count() < args.multicore:
            args.logger.warning(f"Optimal threads <= {N_PHYSICAL_CORES} : {args.multicore} is not recommended")
        for key, value in sorted(vars(args).items()):
            ic(f"Argument {key}: {value}")

        ic("File num check: input folder and project list")
        Helper.equal_num_samples_checker(args.system_structure.input_dir, args.samples, args.logger)

        return func(args)

    return wrapper


@system_struct_checker
def run_pipeline(args: SimpleNamespace) -> None:

    # TODO: time benchmarking
    # TODO: add parquet remove option
    ic("Initilaizing local cluster...")

    cpu_cluster = LocalCluster(
        processes=True,
        n_workers=N_PHYSICAL_CORES,  # DEBUG
        # threads_per_worker=2,
        dashboard_address=":40928",
        local_directory="./tmp",
    )

    cpu_client = Client(cpu_cluster)

    ic(cpu_client)
    ic(cpu_client.dashboard_link)

    for sample, barcode_path in tqdm(args.samples):
        start = datetime.now()
        ExtractorRunner(sample, barcode_path, args)  # TODO: refactor its usage to avoid creating an object
        # if not DEBUGGING:

        if convert_fastq_into_splitted_parquets(args, cpu_client, sample, blocksize="512KiB") < 0:
            raise Exception("FASTQ to parquet conversion failed")
        ic("Generating parquets completed")

        # END: DEBUGGING
        ic("Loading barcides...")
        barcode_df = load_barcode_file(barcode_path, args)
        # Read partitioned parquets with duckdb

        con = duckdb.connect()
        con.sql("--sql SET enable_progress_bar = true;")
        ic("Cartesian product generation...")
        cartesian_df = con.sql(
            f"SELECT * FROM read_parquet('{args.system_structure.seq_split_dir}/*.parquet') CROSS JOIN (SELECT gene, query FROM barcode_df) AS barcode_df;"
        )

        # Perform the query using DuckDB and tag truthy for each regular expression
        ic("regex search & filtering...")
        sequence_col = duckdb.ColumnExpression("sequence")
        query_col = duckdb.ColumnExpression("query")
        extraction_result = con.sql(
            """--sql
            SELECT  *, CAST(regexp_matches(sequence, query) AS INT) AS match
            FROM cartesian_df
            WHERE match > 0 ;
            """
        )

        ic("Filtering out sequences with violation...")
        filtered_ids = con.sql(
            """--sql
                SELECT id, sum(match) AS detection_count_per_sequence
                FROM  extraction_result
                GROUP BY id
                HAVING SUM(match) <= 2;
                """
        )

        # Export to parquets
        ic("Get genuine sequences for export...")
        filtered_result = con.sql(
            """--sql
            SELECT * FROM extraction_result JOIN filtered_ids ON (extraction_result.id = filtered_ids.id);
            """
        )

        ic("Export full parquet...")
        dst = f"{args.system_structure.parquet_dir}/extracted_result.parquet"
        con.sql(
            f"""--sql
            COPY filtered_result TO '{dst}' 
            (FORMAT 'parquet', COMPRESSION 'zstd', ROW_GROUP_SIZE '100_000');
            """
        )

        ic("Export barcode extraction result csv...")
        read_count_table_path = f"{args.system_structure.result_dir}/filtered_result.csv"
        con.sql(
            f"""--sql
            COPY (
            SELECT gene, SUM(match) AS read_count
            FROM read_parquet('{dst}')
            GROUP BY gene
            ) TO '{read_count_table_path}' 
            WITH (FORMAT CSV, HEADER, DELIMITER ',');
            """
        )
        end = datetime.now()
        ic(f"Elapsed time: {end - start} for {sample}+{barcode_path}")
        with open(args.system_structure.result_dir / "processing_time.txt", "w", encoding="utf-8") as f:
            f.write(f"{sample},{barcode_path}\n")
            f.write(f"# {end - start} s elapsed\n")

    ic("All merging jobs fired. Waiting for the final result...")


def convert_fastq_into_splitted_parquets(args, cpu_client, sample, blocksize="64MiB"):
    """

    Convert a FASTQ file into splitted parquet files.

    Args:
        args (object): The arguments object containing system structure information.
        cpu_client (object): The Dask client object for parallel computing.
        sample (str): The sample name.

    Returns:
        int: 0 if the conversion is successful, -1 otherwise.
    """
    try:
        ic("Loading fastq...")
        bag = db.read_text(args.system_structure.input_file_organizer[sample], blocksize=blocksize).map(
            lambda x: x.strip()
        )
        bag = cpu_client.scatter(bag).result()

        ic("Transforming fastq...")
        sequence_ddf = bag.to_dataframe()
        sequence_ddf = cpu_client.scatter(sequence_ddf).result()
        sequence_da = sequence_ddf.to_dask_array(lengths=True)
        sequence_da = sequence_da.reshape(-1, 4)
        sequence_da = cpu_client.scatter(sequence_da).result()
        sequence_ddf = sequence_da.to_dask_dataframe(
            columns=["id", "sequence", "separator", "quality"],
        )
        sequence_ddf = sequence_ddf.drop(columns=["separator"])

        ic("Save raw NGS reads as parquets...")
        sequence_ddf.to_parquet(
            f"{args.system_structure.seq_split_dir}",
            engine="pyarrow",
            # write_index=True,
            write_metadata_file=True,
            compute=True,
        )
        return 0
    except Exception as e:
        ic(e)
        return -1


def load_barcode_file(barcode_path, args):
    """
    Load barcode file and process it.

    Args:
        barcode_path (str): The path to the barcode file.
        args (object): The arguments object containing system structure information.

    Returns:
        pd.DataFrame: The processed barcode dataframe.
    """
    ic("Loading barcode file...")
    barcode_df = pd.read_csv(
        barcode_path,
        sep=args.sep,
        header=None,
    ).fillna("")
    barcode_df.columns = ["gene"] + [
        f"barcode_{i}" for i in range(1, len(barcode_df.columns))
    ]  # Support for n-th barcode

    if not barcode_df["gene"].is_unique:  # Weak assertion: gene names are unique as it is used as PK
        ic("Gene names are not unique")
        ic(barcode_df["gene"].duplicated())
        ic("Drop duplicated gene names...")
        barcode_df = barcode_df.drop_duplicates(subset="Gene")

    barcode_df.loc[:, ["barcode" in col for col in barcode_df.columns]].transform(lambda x: x.str.upper())
    barcode_df["query"] = barcode_df.loc[:, ["barcode" in col for col in barcode_df.columns]].apply(
        lambda x: ".*".join(x), axis=1
    )
    if not barcode_df["query"].is_unique:
        ic("Warning: Duplicated barcode sequences are detected. Please check the barcode file.")

    return barcode_df
