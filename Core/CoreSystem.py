import csv
import glob
import logging
import math
import multiprocessing as mp
import os
import pathlib
import pickle
import re
import subprocess as sp
import sys
from collections import defaultdict
from datetime import datetime
from pdb import set_trace
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
    def SplitSampleInfo(sample):
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

        for idx, file_path in enumerate([
            p
            for p in self.args.system_structure.input_sample_organizer[
                self.sample
            ].glob("*")
        ]):
            if file_path.suffix in [".fastq", ".fa", ".fq", ".fasta"]:
                args.logger.info(f"File name : {file_path.stem}")
                self.args.system_structure.input_file_organizer[file_path.name] = file_path
                break
            
            if idx == len([
                p
                for p in self.args.system_structure.input_sample_organizer[
                    self.sample
                ].glob("*")
            ]) - 1:
                raise Exception("No fastq file in the sample folder")
            

        # self.strInputList  => contains all splitted fastq file path; glob can be used

        self.args.system_structure.seq_split_dir = Helper.mkdir_if_not(
            self.args.system_structure.input_sample_organizer[self.sample]
            / "Split_files"
        )

    def _split_into_chunks(self):

        ### Defensive : original fastq wc == split fastq wc
        # intTotalLines = len(open(self.strInputFile).readlines())

        # TODO: split into chunks based on lines instead of sequences
        intTotalLines = int(
            sp.check_output(
                "wc -l {input_file}".format(input_file=self.strInputFile), shell=True
            ).split()[0]
        )
        intSplitNum = int(
            math.ceil(intTotalLines / float(self.args.chunk_size))
        )  ## e.g. 15.4 -> 16

        if intSplitNum == 0:
            intSplitNum = 1
        logging.info(
            "Total lines:%s, Chunk size:%s, Split number:%s"
            % (intTotalLines, self.args.chunk_size, intSplitNum)
        )

        ## Make minibatches of a single sequence file
        ## self.strInputFile : self.strInputFile : ## './Input/JaeWoo/FASTQ/JaeWoo_test_samples/Test_sample/eCas9_rep1_D4.ext.fastq'
        ## self.strInputList : self.strInputFile : ## './Input/JaeWoo/FASTQ/JaeWoo_test_samples/Test_sample/eCas9_rep1_D4.ext.txt'
        with open(self.strInputFile) as fq, open(self.strInputList, "w") as OutList:

            ## Batch processing
            for intNum in range(1, intSplitNum + 1):

                strSplitFile = (
                    self.system_structure.seq_split_dir
                    + "/{sample}_{num}.fq".format(
                        sample=os.path.basename(self.strInputFile), num=intNum
                    )
                )
                with open(strSplitFile, "w") as out:
                    OutList.write(os.path.basename(strSplitFile) + "\n")
                    intCount = 0

                    for strRow in fq:
                        intCount += 1
                        out.write(strRow)

                        if intCount == self.args.chunk_size:
                            break

        ## defensive
        # strOriginal   = sp.check_output('wc -l {input_file}'.format(input_file=self.strInputFile), shell=True)
        strSplited = sp.check_output(
            "cat {splited}/*.fq | wc -l".format(
                splited=self.system_structure.seq_split_dir
            ),
            shell=True,
        )
        # strOrigianlWc = strOriginal.split()[0]
        intSplitedWc = int(strSplited.decode(encoding="utf-8").replace("\n", ""))

        if intTotalLines != intSplitedWc:
            logging.error(
                "The number of total lines of splited file is not corresponded to origial fastq."
            )
            logging.error(
                "Original FASTQ line number : %s, Splited FASTQ line number : %s"
                % (intTotalLines, strSplited)
            )
            sys.exit(1)

    def MakeExtractorCmd(self):

        listCmd = []
        strReverse = "None"

        with open(self.strInputList) as Input:
            for strFile in Input:
                listFile = strFile.replace("\n", "").split(" ")
                strForward = self.system_structure.seq_split_dir + "/" + listFile[0]

                # if self.strPair == 'True':
                #    strReverse = self.system_structure.seq_split_dir + '/' + listFile[1]

                listCmd.append(
                    ("{python} extractor.py {forw} {barcode} {verbose}").format(
                        python=self.args.python,
                        forw=strForward,
                        barcode=self.args.system_structure.barcode,
                        verbose=self.args.save_trace,
                    )
                )
        return listCmd

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
        if mp.cpu_count() < args.multicore:
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

        # Generating a command for utilizing Python multiprocessing
        args.logger.info("MakeExtractorCmd")
        listCmd = extractor_runner.MakeExtractorCmd()

        args.logger.info("RunMulticore")
        RunMulticore(listCmd, args.multicore)  ## from CoreSystem.py

        # Need Adaptation
        args.logger.info("MakeOutput")
        extractor_runner.MakeOutput(listCmd, sample)


def AttachSeqToIndel(
    strSample, strBarcodeName, strIndelPos, strRefseq, strQueryseq, dictSub
):
    listIndelPos = strIndelPos.split("M")
    intMatch = int(listIndelPos[0])

    if "I" in strIndelPos:
        intInsertion = int(listIndelPos[1].replace("I", ""))
        strInDelSeq = strQueryseq[intMatch : intMatch + intInsertion]

    elif "D" in strIndelPos:
        intDeletion = int(listIndelPos[1].replace("D", ""))
        strInDelSeq = strRefseq[intMatch : intMatch + intDeletion]

    else:
        logging.info(
            "strIndelClass is included I or D. This variable is %s" % strIndelPos
        )
        raise Exception

    strInDelPosSeq = strIndelPos + "_" + strInDelSeq

    try:
        _ = dictSub[strSample][strBarcodeName]
    except KeyError:
        dictSub[strSample][strBarcodeName] = {}

    try:
        dictSub[strSample][strBarcodeName][strBarcodeName + ":" + strInDelPosSeq][
            "IndelCount"
        ] += 1
    except KeyError:
        dictSub[strSample][strBarcodeName][strBarcodeName + ":" + strInDelPosSeq] = {
            "IndelCount": 1
        }


def RunProgram(sCmd):
    sp.call(sCmd, shell=False)


def RunMulticore(lCmd, iCore):
    for sCmd in lCmd:
        print(sCmd)

    p = mp.Pool(iCore)
    p.map_async(RunProgram, lCmd).get()
    p.close()
