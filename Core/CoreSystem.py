import multiprocessing as mp
import os
import pathlib
import shlex
import subprocess as sp
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from types import SimpleNamespace


class Helper(object):
    @staticmethod
    def mkdir_if_not(directory: str) -> pathlib.Path:
        """
        > If the directory doesn't exist, create it

        :param directory: The directory to create
        :type directory: str
        :return: A path object
        """
        path = pathlib.Path(directory)
        path.mkdir(parents=True, exist_ok=True)
        return path

    @staticmethod
    def load_samples(dir: pathlib.Path) -> list:
        """
        It reads a file and returns a list of non-empty lines that don't start with a hash mark.

        :param dir: the directory of the samples file
        :type dir: pathlib.Path
        :return: A list of samples.
        """
        with open(dir, "r") as file:
            lines = [l.strip("\n") for l in file.readlines()]
            sample_barcode_list = [
                line.split(",")[:2] for line in lines if line[0] != "#"
            ]

        return sample_barcode_list

    @staticmethod  ## defensive
    def equal_num_samples_checker(
        proj_path: pathlib.Path, loaded_samples: list, logger
    ):
        """
        > This function checks if the number of samples in the Input folder and in the User folder
        matches

        :param proj_path: pathlib.Path, loaded_samples: list, logger
        :type proj_path: pathlib.Path
        :param loaded_samples: a list of sample names
        :type loaded_samples: list
        :param logger: a logger object
        """

        if len(list(proj_path.glob("*"))) != len(loaded_samples):
            logger.warning(
                "The number of samples in the Input folder and in the User folder does not matched."
            )

            input_entries = [i.name for i in proj_path.glob("*")]
            user_entries = [i for i in loaded_samples]
            logger.warning(
                f"Input folder: {len(list(proj_path.glob('*')))}, Project list samples: {len(loaded_samples)}"
            )
            logger.warning(
                f"Input folder: {[i for i in input_entries if i not in user_entries]}"
            )
            logger.warning(
                f"Project list samples: {[u for u in user_entries if u not in input_entries]}"
            )
        else:
            logger.info("The file list is correct, pass\n")

    @staticmethod
    def SplitSampleInfo(sample):
        # Sample\tReference\tGroup
        # logging.info("[Deprecated] Processing sample : %s" % sample)

        # sample = sample.replace(" ", r"\ ")
        return sample


# > The class creates a directory structure for a user and a project
class SystemStructure(object):
    def __init__(
        self,
        user_name: str,
        project_name: str,
    ):
        # https://www.notion.so/dengardengarden/s-Daily-Scrum-Reports-74d406ce961c4af78366a201c1933b66#cd5b57433eca4c6da36145d81adbbe5e
        self.user_name = user_name
        self.project_name = project_name
        self.project_samples_path = ""
        self.input_sample_organizer = defaultdict(pathlib.Path)
        self.input_file_organizer = defaultdict(
            pathlib.Path
        )  # .fastq.gz #TODO: FLASh integration? https://ccb.jhu.edu/software/FLASH/
        # should contain absolute path to the file
        self.output_sample_organizer = defaultdict(pathlib.Path)

        self.user_dir = Helper.mkdir_if_not("User" + "/" + self.user_name)
        self.barcode_dir = Helper.mkdir_if_not("Barcodes")
        self.input_dir = Helper.mkdir_if_not(
            "Input" + "/" + self.user_name + "/" + self.project_name
        )
        self.project_samples_path = pathlib.Path(
            "User" + "/" + self.user_name + "/" + f"{self.project_name}.txt"
        )
        if not self.project_samples_path.exists():
            with open(self.project_samples_path, "w") as f:
                f.write("")

        self.output_dir = Helper.mkdir_if_not(
            "Output" + "/" + self.user_name + "/" + self.project_name
        )

    def mkdir_sample(self, sample_name: str):
        # TODO
        self.input_sample_organizer[sample_name] = Helper.mkdir_if_not(
            self.input_dir / sample_name
        )
        self.output_sample_organizer[sample_name] = Helper.mkdir_if_not(
            self.output_dir / sample_name
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


class ExtractorRunner:
    def __init__(self, sample: str, args: SimpleNamespace):
        args.python = sys.executable
        # Find python executable if not specified
        args.system_structure.mkdir_sample(sample)
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

    def _split_into_chunks(self):
        ### Defensive : original fastq wc == split fastq wc
        # https://docs.python.org/3.9/library/subprocess.html#security-considerations
        sp.run(
            shlex.split(
                shlex.quote(
                    f'split "{self.args.system_structure.input_file_organizer[self.sample]}" -l {4 * self.args.chunk_size} -d -a 6 --additional-suffix=.fastq {self.args.system_structure.seq_split_dir}/split_'
                )
            ),
            shell=True,
            check=True,
        )

        self.args.logger.info(
            f"The number of split files:{len(list(self.args.system_structure.seq_split_dir.glob('*')))}"
        )

    def _populate_command(self, barcode):
        return [
            (
                str(pathlib.Path.cwd() / self.args.system_structure.seq_split_dir / f),
                str(
                    pathlib.Path.cwd()
                    / self.args.system_structure.barcode_dir
                    / barcode
                ),
                self.args.logger,
                f"{(pathlib.Path(self.args.system_structure.result_dir) / 'parquets').absolute()}",
                self.args.sep,
            )
            for f in sorted(os.listdir(self.args.system_structure.seq_split_dir))
            if f.endswith(".fastq")
        ]


def system_struct_checker(func):
    def wrapper(args: SimpleNamespace):
        args.multicore = os.cpu_count() if args.multicore == 0 else args.multicore
        args.logger.info("Program start")
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

        func(args)
        args.logger.info("Extraction process completed.")

    return wrapper


@system_struct_checker
def run_pipeline(args: SimpleNamespace) -> None:
    for sample, barcode in args.samples:
        sample = Helper.SplitSampleInfo(sample)

        extractor_runner = ExtractorRunner(sample, args)

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
