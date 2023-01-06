import csv
import glob
import logging
import math
import subprocess as sp
import shlex
import os
import pathlib
import pickle
import re
import sys
import pandas as pd
from collections import defaultdict
from datetime import datetime
from pdb import set_trace
from types import SimpleNamespace
from concurrent.futures import ProcessPoolExecutor

from extractor import main as extractor_main


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

    @staticmethod  ## defensive
    def CheckAllDone(strOutputProject, listSamples):
        intProjectNumInOutput = len(
            [
                i
                for i in sp.check_output("ls %s" % strOutputProject, shell=False).split(
                    "\n"
                )
                if i not in ["All_results", "Log", ""]
            ]
        )

        if intProjectNumInOutput != len(listSamples):
            logging.warning(
                "The number of samples in the output folder and in the project list does not matched."
            )
            logging.warning(
                "Output folder: %s, Project list samples: %s\n"
                % (intProjectNumInOutput, len(listSamples))
            )
        else:
            logging.info("All output folders have been created.\n")

    @staticmethod
    def CheckIntegrity(strBarcodeFile, strSeq):  ## defensive
        rec = re.compile(r"[A|C|G|T|N]")

        if ":" in strSeq:
            strSeq = strSeq.split(":")[1]

        strNucle = re.findall(rec, strSeq)
        if len(strNucle) != len(strSeq):
            logging.error(
                "This sequence is not suitable, check A,C,G,T,N are used only : %s"
                % strBarcodeFile
            )
            set_trace()
            sys.exit(1)

    @staticmethod
    def PreventFromRmMistake(strCmd):
        rec = re.compile(
            r"rm.+-rf*.+(\.$|\/$|\*$|User$|Input$|Output$)"
        )  ## This reg can prevent . / * ./User User ...
        if re.findall(rec, strCmd):
            raise Exception("%s is critical mistake! never do like this." % strCmd)


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

        self.log_dir = Helper.mkdir_if_not(self.output_dir / "Log")

        self.log = (
            self.log_dir
            / f"{str(datetime.now()).replace('-', '_').replace(':', '_').replace(' ', '_').split('.')[0]}_log.txt"
        )

        if not self.log.exists():
            with open(self.log, "w") as f:
                f.write("")

    def mkdir_sample(self, sample_name: str):
        # TODO
        self.input_sample_organizer[sample_name] = Helper.mkdir_if_not(
            self.input_dir / sample_name
        )
        self.output_sample_organizer[sample_name] = Helper.mkdir_if_not(
            self.output_dir / sample_name
        )

        self.tmp_dir = Helper.mkdir_if_not(
            self.output_sample_organizer[sample_name] / "Tmp"
        )
        self.pkl_dir = Helper.mkdir_if_not(self.tmp_dir / "Pickle")
        self.result_dir = Helper.mkdir_if_not(
            self.output_sample_organizer[sample_name] / "Result"
        )
        self.all_result_dir = Helper.mkdir_if_not(self.result_dir / "All_results")


class CoreHash(object):
    @staticmethod
    def MakeHashTable(strSeq, intBarcodeLen):
        listSeqWindow = [
            strSeq[i : i + intBarcodeLen]
            for i in range(len(strSeq))[: -intBarcodeLen - 1]
        ]
        return listSeqWindow

    @staticmethod
    def IndexHashTable(dictRef, strSeqWindow, intFirstBarcode):
        lCol_ref = dictRef[strSeqWindow]
        strBarcode = strSeqWindow
        intFirstBarcode = 1

        return (lCol_ref, strBarcode, intFirstBarcode)


class CoreGotoh(object):
    def __init__(self, strEDNAFULL="", floOg="", floOe=""):
        self.npAlnMatrix = CRISPResso2Align.read_matrix(strEDNAFULL)
        self.floOg = floOg
        self.floOe = floOe

    def GapIncentive(self, strRefSeqAfterBarcode):
        ## cripsress no incentive == gotoh
        intAmpLen = len(strRefSeqAfterBarcode)
        npGapIncentive = np.zeros(intAmpLen + 1, dtype=np.int)
        return npGapIncentive

    def RunCRISPResso2(
        self, strQuerySeqAfterBarcode, strRefSeqAfterBarcode, npGapIncentive
    ):
        listResult = CRISPResso2Align.global_align(
            strQuerySeqAfterBarcode.upper(),
            strRefSeqAfterBarcode.upper(),
            matrix=self.npAlnMatrix,
            gap_open=self.floOg,
            gap_extend=self.floOe,
            gap_incentive=npGapIncentive,
        )
        return listResult


class ExtractorRunner:
    def __init__(self, sample: str, args: SimpleNamespace):
        args.python = (
            sys.executable if args.python is None else args.python
        )  # Find python executable if not specified
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
        # intTotalLines = len(open(self.strInputFile).readlines())
        # https://docs.python.org/3.9/library/subprocess.html#security-considerations
        sp.run(
            shlex.split(
                shlex.quote(
                    f"split {self.args.system_structure.input_file_organizer[self.sample]} -l {4 * self.args.chunk_size} -d -a 3 --additional-suffix=.fastq {self.args.system_structure.seq_split_dir}/split_"
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
                self.sample,
                pathlib.Path.cwd() / self.args.system_structure.seq_split_dir / f,
                pathlib.Path.cwd()
                / self.args.system_structure.barcode_dir
                / self.args.barcode,
                
                self.args.logger,
            )
            for f in sorted(os.listdir(self.args.system_structure.seq_split_dir))
            if f.endswith(".fastq")
        ]

    def MakeOutput(self, listCmd, strSample):
        # TODO: magic keys; design changes can bring disaster
        forw = (
            listCmd[0].split("/")[-3].split(" ")[0].split(".")[0]
        )  # Barcode files are relocated to "./Barcodes" directory

        RESULT_DIR = f"./Results/temp/{forw}"
        OUTPUT_DIR = "./Results/All_results"

        if not os.path.exists(OUTPUT_DIR):
            os.makedirs(OUTPUT_DIR)

        # Get temp file names from the list of commands

        # sample_name == ./temp/sequence.fastq_1.fq_result_info.txt

        counts = glob.glob(
            os.path.join(RESULT_DIR, "*.txt")
        )  # find all temporary read count files
        findings = glob.glob(
            os.path.join(RESULT_DIR, "*.trace")
        )  # find all temporary saved sequence files

        # Initialize internal data object to add data segments
        barcode_dict_for_counts = defaultdict(int)
        list_for_findings = []

        # Sum up the counts
        for item in counts:
            # Open summary file from each dir
            with open(item, "r") as f:
                for line in f:
                    # TODO: processivity improvement with map()
                    # "Barcode" : "count"
                    line_tokens = line.split(":")
                    if len(line_tokens) > 1:
                        barcode = line_tokens[0].strip()
                        cnt = int(line_tokens[1].strip())
                        barcode_dict_for_counts[barcode] += cnt

        with open(f"{OUTPUT_DIR}/{strSample}_Summary.txt", "w") as summary:
            for barcode, count in barcode_dict_for_counts.items():
                summary.write(f"{barcode} : {count}\n")

        # Gathering the sequence findings
        if self.args.save_trace:
            for item in findings:
                with open(item, "rb") as fp:
                    pkl_obj = pickle.load(fp)
                    list_for_findings.extend(pkl_obj)

            with open(f"{OUTPUT_DIR}/{strSample}_Verbose.csv", "w") as verbose:
                csv_out = csv.writer(verbose)
                csv_out.writerow(("Barcode", "Sequence"))
                for tup in list_for_findings:
                    csv_out.writerow(tup)

        # Remove all intermediate files

        # sp.call('rm -r ./Results/temp/', shell=True)


def system_struct_checker(func):
    def wrapper(args: SimpleNamespace):

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

        # TODO: check if all files are created
        args.logger.info("Check that all folder are well created.")
        Helper.CheckAllDone(InstInitFolder.strOutputProjectDir, listSamples)  # TODO

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
        run_extractor_mp(listCmd, args.multicore, args.logger)

        # TODO
        args.logger.info("MakeOutput")
        extractor_runner.MakeOutput(listCmd, sample)


def run_extractor_mp(lCmd, iCore, logger) -> None:
    for sCmd in lCmd:
        logger.info(f"Running {sCmd} command with {iCore} cores")

    with ProcessPoolExecutor(max_workers=iCore) as executor:
        executor.map(extractor_main, lCmd)

    logger.info(f"All extraction process completed")


#LBJ editing
class ReadBarcode(object):

    def __init__(self):
        self.FilePath = ''
        self.user = ''
        self.project = ''
        self.cas_type = ''
        self.index = ''
        self.poly_t = 0

        self.BarcodeList = pd.DataFrame()

    def SelectFromExcel(self):

        df = pd.read_excel('./User/Barcode Database.xlsx', engine = 'openpyxl', sheet_name="Oligo seq") #edit to pathlib later
        out = df.loc[:, ['Nuc type', 'INDEX', 'Gene name', 'Barcode sequence']]
        self.BarcodeList = out
        return out

    def UseCSV(self):
        f_csv = open('./Input/input.csv', 'r') #edit to pathlib later
        f_read = csv.reader(f_csv)
        csv_data = [line for line in f_read]
        self.user = csv_data[1][0]              #user_name in CSV file
        self.project = csv_data[1][1]           #project_name in CSV file
        self.cas_type = csv_data[1][2]          #cas_type in CSV file
        self.index = csv_data[1][3]             #index in CSV file
        self.poly_t = csv_data[1][4]            #poly-T length in CSV file
        Barcode_temp = pd.DataFrame()
        self.BarcodeList = self.BarcodeList[self.BarcodeList['Nuc type'] == self.cas_type]  #DataFrame filter by nuclease type, ex.Cas9 or TnpB
        self.BarcodeList = self.BarcodeList[self.BarcodeList['INDEX'] == self.index]        #DataFrame filter by index

        for gene in csv_data[3:] :
            if(gene[0]==''):
                break
            Barcode_temp = pd.concat([Barcode_temp,self.BarcodeList[self.BarcodeList['Gene name'].str.contains(gene[0])]])
        
        Barcode_temp['Barcode sequence'] = 'T'*int(self.poly_t) + Barcode_temp['Barcode sequence'].astype(str)

        return Barcode_temp
