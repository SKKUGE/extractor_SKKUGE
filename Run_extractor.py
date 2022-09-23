#!/usr/bin/env python

import csv
import glob
import logging
import math
import os
import pickle
# import cPickle as pickle
import subprocess as sp
import sys
from collections import defaultdict
from optparse import OptionParser

sys.path.insert(0, os.path.dirname(os.getcwd()))
from Core.CoreSystem import InitialFolder, UserFolderAdmin, Helper, RunMulticore, CheckProcessedFiles


class clsExtractorRunner(UserFolderAdmin):
    """
    self.strOutputDir is inherited variable.

    """

    def __init__(self, strSample, options, InstInitFolder):
        UserFolderAdmin.__init__(self, strSample, options)
        self.MakeSampleFolder()

        self.strProjectFile = InstInitFolder.strProjectFile
        self.intChunkSize = options.chunk_number
        self.strPickle = options.pickle
        self.strSplit = options.split
        self.strPython = options.python
        self.barcode = options.barcode
        self.verbose = True if options.save_trace.upper() == "TRUE" else False

        ## Files needed in the FASTQ directory
        self.strFastqDir = './Input/{user}/FASTQ/{project}'.format(user=self.strUser,
                                                                   project=self.strProject)
        ## './Input/JaeWoo/FASTQ/JaeWoo_test_samples/Test_sample'
        self.strSampleDir = os.path.join(self.strFastqDir, self.strSample)

        self.strFastq_name = ''
        for strFile in os.listdir(self.strSampleDir):
            if os.path.isfile(self.strSampleDir + '/' + strFile) and strFile.split('.')[-1] == 'fastq':
                self.strFastq_name = '.'.join(strFile.split('.')[:-1])
        logging.info('File name : %s' % self.strFastq_name)

        ## './Input/JaeWoo/FASTQ/JaeWoo_test_samples/Test_sample/eCas9_rep1_D4.ext.fastq'
        self.strInputFile = os.path.join(self.strSampleDir, self.strFastq_name + '.fastq')
        ## './Input/JaeWoo/FASTQ/JaeWoo_test_samples/Test_sample/eCas9_rep1_D4.ext.txt'
        self.strInputList = os.path.join(self.strSampleDir, self.strFastq_name + '.txt')

        ## './Input/JaeWoo/FASTQ/JaeWoo_test1_samples/Test_sample/Split_files'
        self.strSplitPath = os.path.join(self.strSampleDir, 'Split_files')
        Helper.MakeFolderIfNot(self.strSplitPath)

        # TODO: no choice
        self.strPair = 'False'  # FASTQ pair: True, False

    def SplitFile(self):

        ### Defensive : original fastq wc == split fastq wc
        # intTotalLines = len(open(self.strInputFile).readlines())
        intTotalLines = int(
            sp.check_output('wc -l {input_file}'.format(input_file=self.strInputFile), shell=True).split()[0])
        intSplitNum = int(math.ceil(intTotalLines / float(self.intChunkSize)))  ## e.g. 15.4 -> 16

        if intSplitNum == 0: intSplitNum = 1
        logging.info('Total lines:%s, Chunk size:%s, Split number:%s' % (intTotalLines, self.intChunkSize, intSplitNum))

        ## Make minibatches of a single sequence file
        ## self.strInputFile : self.strInputFile : ## './Input/JaeWoo/FASTQ/JaeWoo_test_samples/Test_sample/eCas9_rep1_D4.ext.fastq'
        ## self.strInputList : self.strInputFile : ## './Input/JaeWoo/FASTQ/JaeWoo_test_samples/Test_sample/eCas9_rep1_D4.ext.txt'
        with open(self.strInputFile) as fq, \
                open(self.strInputList, 'w') as OutList:

            ## Batch processing
            for intNum in range(1, intSplitNum + 1):

                strSplitFile = self.strSplitPath + '/{sample}_{num}.fq'.format(
                    sample=os.path.basename(self.strInputFile),
                    num=intNum)
                with open(strSplitFile, 'w') as out:
                    OutList.write(os.path.basename(strSplitFile) + '\n')
                    intCount = 0

                    for strRow in fq:
                        intCount += 1
                        out.write(strRow)

                        if intCount == self.intChunkSize:
                            break

        ## defensive
        # strOriginal   = sp.check_output('wc -l {input_file}'.format(input_file=self.strInputFile), shell=True)
        strSplited = sp.check_output('cat {splited}/*.fq | wc -l'.format(splited=self.strSplitPath), shell=True)
        # strOrigianlWc = strOriginal.split()[0]
        intSplitedWc = int(strSplited.decode(encoding="utf-8").replace('\n', ''))

        if intTotalLines != intSplitedWc:
            logging.error('The number of total lines of splited file is not corresponded to origial fastq.')
            logging.error(
                'Original FASTQ line number : %s, Splited FASTQ line number : %s' % (intTotalLines, strSplited))
            sys.exit(1)

    def MakeExtractorCmd(self):

        listCmd = []
        strReverse = 'None'

        with open(self.strInputList) as Input:
            for strFile in Input:
                listFile = strFile.replace('\n', '').split(' ')
                strForward = self.strSplitPath + '/' + listFile[0]

                # if self.strPair == 'True':
                #    strReverse = self.strSplitPath + '/' + listFile[1]

                listCmd.append(('{python} extractor.py {forw} {barcode} {verbose}').format(
                    python=self.strPython,
                    forw=strForward, barcode=self.barcode, verbose=self.verbose))
        return listCmd

    def MakeOutput(self, listCmd, strSample):
        # TODO: magic keys; design changes can bring disaster
        forw = listCmd[0].split('/')[-3].split(' ')[0].split('.')[
            0]  # Barcode files are relocated to "./Barcodes" directory

        RESULT_DIR = f"./Results/temp/{forw}"
        OUTPUT_DIR = "./Results/All_results"

        if not os.path.exists(OUTPUT_DIR):
            os.makedirs(OUTPUT_DIR)

        # Get temp file names from the list of commands

        # sample_name == ./temp/sequence.fastq_1.fq_result_info.txt

        counts = glob.glob(os.path.join(RESULT_DIR, "*.txt"))  # find all temporary read count files
        findings = glob.glob(os.path.join(RESULT_DIR, "*.trace"))  # find all temporary saved sequence files

        # Initialize internal data object to add data segments
        barcode_dict_for_counts = defaultdict(int)
        list_for_findings = []

        # Sum up the counts
        for item in counts:
            # Open summary file from each dir
            with open(item, 'r') as f:
                for line in f:
                    # TODO: processivity improvement with map()
                    # "Barcode" : "count"
                    line_tokens = line.split(':')
                    if len(line_tokens) > 1:
                        barcode = line_tokens[0].strip()
                        cnt = int(line_tokens[1].strip())
                        barcode_dict_for_counts[barcode] += cnt

        with open(f'{OUTPUT_DIR}/{strSample}_Summary.txt', 'w') as summary:
            for barcode, count in barcode_dict_for_counts.items():
                summary.write(f"{barcode} : {count}\n")

        # Gathering the sequence findings
        if self.verbose:
            for item in findings:
                with open(item, 'rb') as fp:
                    pkl_obj = pickle.load(fp)
                    list_for_findings.extend(pkl_obj)

            with open(f'{OUTPUT_DIR}/{strSample}_Verbose.csv', 'w') as verbose:
                csv_out = csv.writer(verbose)
                csv_out.writerow(("Barcode", "Sequence"))
                for tup in list_for_findings:
                    csv_out.writerow(tup)

        # Remove all intermediate files

        # sp.call('rm -r ./Results/temp/', shell=True)


def Main():
    parser = OptionParser(
        'Indel search program for CRISPR CAS9 & CPF1\n<All default option> python2.7 Run_indel_searcher.py --pam_type Cas9 --pam_pos Forward')

    parser.add_option('-t', '--thread', default='1', type='int', dest='multicore',
                      help='multiprocessing number, recommendation:t<16')
    parser.add_option('-c', '--chunk_number', default='400000', type='int', dest='chunk_number',
                      help='split FASTQ, must be multiples of 4. file size < 1G recommendation:40000, size > 1G recommendation:400000')
    parser.add_option('--python', dest='python', help='The python path including the CRISPResso2')
    parser.add_option('--user', dest='user_name', help='The user name with no space')
    parser.add_option('--project', dest='project_name', help='The project name with no space')
    parser.add_option('--pickle', dest='pickle', default='False',
                      help='Dont remove the pickles in the tmp folder : True, False')
    parser.add_option('--split', dest='split', default='False',
                      help='Dont remove the split files in the input folder : True, False')
    parser.add_option('--barcode', dest='barcode', default="Barcode.txt", help='Barcode file')
    parser.add_option('--save_trace', dest='save_trace', default='True')

    options, args = parser.parse_args()

    InstInitFolder = InitialFolder(options.user_name, options.project_name, os.path.basename(__file__))
    InstInitFolder.MakeDefaultFolder()
    InstInitFolder.MakeInputFolder()

    logging.info('Program start')
    if options.multicore > 15:
        logging.warning('Optimal treads <= 15')
    logging.info(str(options))

    with open(InstInitFolder.strProjectFile) as Sample_list:

        listSamples = Helper.RemoveNullAndBadKeyword(Sample_list)
        intProjectNumInTxt = len(listSamples)

        strInputProject = f"./Input/{options.user_name}/FASTQ/{options.project_name}"

        @CheckProcessedFiles
        def RunPipeline(**kwargs):

            setGroup = set()
            for strSample in listSamples:

                tupSampleInfo = Helper.SplitSampleInfo(strSample)
                if not tupSampleInfo: continue
                strSample = tupSampleInfo

                # Locating input files in the Reference and FASTQ directory
                InstRunner = clsExtractorRunner(strSample, options, InstInitFolder)

                # Chunking
                logging.info('SplitFile')
                InstRunner.SplitFile()

                # Generating a command for utilizing Python multiprocessing
                logging.info('MakeExtractorCmd')
                listCmd = InstRunner.MakeExtractorCmd()

                logging.info('RunMulticore')
                RunMulticore(listCmd, options.multicore)  ## from CoreSystem.py

                # Need Adaptation
                logging.info('MakeOutput')
                InstRunner.MakeOutput(listCmd, strSample)

        RunPipeline(InstInitFolder=InstInitFolder,
                    strInputProject=strInputProject,
                    intProjectNumInTxt=intProjectNumInTxt,
                    listSamples=listSamples,
                    logging=logging)

    logging.info('Program end')


# END:def


if __name__ == '__main__':
    Main()
