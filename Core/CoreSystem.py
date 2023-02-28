import logging
import subprocess as sp
import multiprocessing as mp
import shlex
import os
import pathlib
import sys
from collections import defaultdict
from types import SimpleNamespace
from concurrent.futures import ProcessPoolExecutor

import pandas as pd


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
        It reads a file and returns a list of non-empty lines that don't start with a hash

        :param dir: the directory of the samples file
        :type dir: pathlib.Path
        :return: A list of samples.
        """
        with open(dir, "r") as file:
            sample_list = [
                sample
                for sample in list(
                    filter(None, map(str.strip, file.read().split("\n")))
                )
                if sample[0] != "#"
            ]

        return sample_list

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
    def SplitSampleInfo(sample):  # Deprecated
        # Sample\tReference\tGroup
        logging.info("[Deprecated] Processing sample : %s" % sample)


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
        self.result_dir = Helper.mkdir_if_not(
            self.output_sample_organizer[sample_name]
        )


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

    def _split_into_chunks(self):

        ### Defensive : original fastq wc == split fastq wc
        # https://docs.python.org/3.9/library/subprocess.html#security-considerations
        sp.run(
            shlex.split(
                shlex.quote(
                    f"split {self.args.system_structure.input_file_organizer[self.sample]} -l {4 * self.args.chunk_size} -d -a 6 --additional-suffix=.fastq {self.args.system_structure.seq_split_dir}/split_"
                )
            ),
            shell=True,
            check=True,
        )

        self.args.logger.info(
            f"The number of split files:{len(list(self.args.system_structure.seq_split_dir.glob('*')))}"
        )

    def _populate_command(self):

        return [
            (
                str(pathlib.Path.cwd() / self.args.system_structure.seq_split_dir / f),
                str(
                    pathlib.Path.cwd()
                    / self.args.system_structure.barcode_dir
                    / self.args.barcode
                ),
                self.args.logger,
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
    for sample in args.samples:
        Helper.SplitSampleInfo(sample)

        extractor_runner = ExtractorRunner(sample, args)

        # Chunking
        args.logger.info("Splitting sequecnes into chunks")
        extractor_runner._split_into_chunks()

        args.logger.info("Populating command...")
        listCmd = extractor_runner._populate_command()

        args.logger.info("RunMulticore")
        run_extractor_mp(listCmd, args.multicore, args.logger).to_csv(
            f"{args.system_structure.result_dir}/{sample}+extraction_result.csv"
        )


def run_extractor_mp(lCmd, iCore, logger) -> pd.DataFrame:
    import time
    from extractor import main as extractor_main

    for sCmd in lCmd:
        logger.info(f"Running {sCmd} command with {iCore} cores")

    result = []
    start = time.time()
    with ProcessPoolExecutor(max_workers=iCore) as executor:
        for rval in executor.map(extractor_main, lCmd):
            result.append(rval)
    end = time.time()
    logger.info(f"Extraction is done. {end - start}s elapsed.")
    logger.info(f"All extraction subprocesses completed")
    logger.info(f"Merging extraction results...")

    df = result.pop()
    for d in result:
        df["Read_counts"] = df["Read_counts"] + d["Read_counts"]

    return df
