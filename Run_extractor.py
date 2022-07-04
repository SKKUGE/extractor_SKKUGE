#!/usr/bin/env python

import logging
import math
import os
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

    def __init__(self, strSample, strRef, options, InstInitFolder):
        UserFolderAdmin.__init__(self, strSample, strRef, options, InstInitFolder.strLogPath)
        self.MakeSampleFolder()

        self.strProjectFile = InstInitFolder.strProjectFile
        self.intChunkSize = options.chunk_number
        self.strQualCutoff = options.base_quality
        self.intInsertionWin = options.insertion_window  # Insertion window 0,1,2,3,4
        self.intDeletionWin = options.deletion_window  # Deletion window 0,1,2,3,4
        self.strPamType = options.pam_type  # CRISPR type : Cpf1(2 cleavages), Cas9(1 cleavage)
        self.strPamPos = options.pam_pos  # Barcode target position : Forward (barcode + target), Reverse (target + barcode)
        self.strPickle = options.pickle
        self.strClassFASTQ = options.class_fastq
        self.strSplit = options.split
        self.strLogPath = InstInitFolder.strLogPath

        ## Files needed in the Reference directory
        self.strBarcodeFile = os.path.join(self.strRefDir, 'Barcode.txt')
        self.strReferenceSeqFile = os.path.join(self.strRefDir, 'Reference_sequence.txt')
        self.strTargetSeqFile = os.path.join(self.strRefDir, 'Target_region.txt')
        self.strRefFile = os.path.join(self.strRefDir, 'Reference.fa')

        ## The file name required for the user is 'B'arcode.txt but it may be written as 'b'arcode.txt by mistake.
        ## This part is to fix the situation as mentioned above.
        if not os.path.isfile(self.strBarcodeFile):
            if os.path.isfile(self.strRefDir + 'barcode.txt'):
                self.strBarcodeFile = self.strRefDir + 'barcode.txt'
            else:
                logging.error('Barcode path is not correct, please make sure the path correctly.')
        if not os.path.isfile(self.strReferenceSeqFile):
            if os.path.isfile(self.strRefDir + 'reference_sequence.txt'):
                self.strReferenceSeqFile = self.strRefDir + 'reference_sequence.txt'
            else:
                logging.error('Reference path is not correct, please make sure the path correctly.')
        if not os.path.isfile(self.strTargetSeqFile):
            if os.path.isfile(self.strRefDir + 'target_region.txt'):
                self.strTargetSeqFile = self.strRefDir + 'target_region.txt'
            else:
                logging.error('Target path is not correct, please make sure the path correctly.')

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

    def MakeReference(self):
        # Do always make Reference file for avoiding misuse

        with open(self.strBarcodeFile) as Barcode, \
                open(self.strTargetSeqFile) as Target, \
                open(self.strReferenceSeqFile) as Ref, \
                open(self.strRefFile, 'w') as Output:

            listBarcode = Helper.RemoveNullAndBadKeyword(Barcode)
            listTarget = Helper.RemoveNullAndBadKeyword(Target)
            listRef = Helper.RemoveNullAndBadKeyword(Ref)

            ## defensive
            assert len(listBarcode) == len(listTarget) == len(
                listRef), 'Barcode, Target and Reference must be a same row number.'

            # String pre-processing
            listName = []
            for strBar, strTar in zip(listBarcode, listTarget):
                strBar = strBar.replace('\n', '').replace('\r', '').strip().upper()
                strTar = strTar.replace('\n', '').replace('\r', '').strip().upper()

                # The strings should be composed of only A,T,C,G,N
                Helper.CheckIntegrity(self.strBarcodeFile, strBar)  ## defensive
                Helper.CheckIntegrity(self.strBarcodeFile, strTar)  ## defensive

                listName.append(strBar + ':' + strTar + '\n')

            for i, strRow in enumerate(listRef):
                strRow = strRow.replace('\r', '').strip().upper()
                Output.write('>' + listName[i] + strRow + '\n')

    def MakeExtractorCmd(self):

        listCmd = []
        strReverse = 'None'

        with open(self.strInputList) as Input:
            for strFile in Input:
                listFile = strFile.replace('\n', '').split(' ')
                strForward = self.strSplitPath + '/' + listFile[0]

                # if self.strPair == 'True':
                #    strReverse = self.strSplitPath + '/' + listFile[1]

                listCmd.append(('{python} extractor.py {forw} {reve} {ref} {pair} {GapO} {GapE}'
                                ' {Insertion_win} {Deletion_win} {PAM_type} {PAM_pos} {Qual} {outdir} {logpath}').format(
                    python=self.strPython,
                    forw=strForward, reve=strReverse, ref=self.strRefFile, pair=self.strPair,
                    GapO=self.strGapOpen, GapE=self.strGapExtend,
                    Insertion_win=self.intInsertionWin, Deletion_win=self.intDeletionWin,
                    PAM_type=self.strPamType, PAM_pos=self.strPamPos, Qual=self.strQualCutoff,
                    outdir=self.strOutSampleDir, logpath=self.strLogPath))
        return listCmd

    def RunIndelFreqCalculator(self):
        sp.call('{python} Indel_frequency_calculator.py {outdir} {sample} {logpath}'.format(python=self.strPython,
                                                                                            outdir=self.strOutSampleDir,
                                                                                            sample=self.strSample,
                                                                                            logpath=self.strLogPath),
                shell=True)
        sp.call('{python} Summary_all_trim.py {outdir} {sample} {logpath}'.format(python=self.strPython,
                                                                                  outdir=self.strOutSampleDir,
                                                                                  sample=self.strSample,
                                                                                  logpath=self.strLogPath), shell=True)
        sp.call('cp $(find ./Output/{user}/{project} -name "*.tsv") ./Output/{user}/{project}/All_results'.format(
            user=self.strUser,
            project=self.strProject), shell=True)

    def IndelNormalization(self):

        sp.call('{python} Indel_normalization.py {project_file} {user} {project}'.format(python=self.strPython,
                                                                                         project_file=self.strProjectFile,
                                                                                         user=self.strUser,
                                                                                         project=self.strProject),
                shell=True)

    def MakeOutput(self, listCmd, strSample):
        RESULT_DIR = "./Results/temp"
        OUTPUT_DIR = "./Results/All_results"

        if not os.path.exists(OUTPUT_DIR):
            os.makedirs(OUTPUT_DIR)

        # Get temp file names from the list of commands

        dirs_to_search = []
        for cmd in listCmd:
            # cmd preprocessing

            # sample_name == 0314_1.1, ...
            src_file_dir = cmd.split(' ')[2]
            split_name_token = src_file_dir.split('/')[-1].split('.')
            split_name = split_name_token[0] + '.' + split_name_token[2].split('_')[-1]
            processed_split_dir = os.path.join(RESULT_DIR, split_name)

            item = (processed_split_dir, split_name)
            dirs_to_search.append(item)

        # Initialize internal data object to add data segments
        barcode_dict = defaultdict(int)

        # Sum up the counts
        for item in dirs_to_search:
            # Open summary file from each dir
            with open(f"{item[0]}/{item[1]}_result_info.txt", 'r') as f:
                for line in f:
                    # "Barcode" : "count"
                    line_tokens = line.split(':')
                    barcode = line_tokens[0].strip()
                    cnt = int(line_tokens[1].strip())
                    barcode_dict[barcode] += cnt

        with open(f'{OUTPUT_DIR}/{strSample}_Summary.txt', 'w') as summary:
            for barcode, count in barcode_dict.items():
                summary.write(f"{barcode} : {count}\n")

        # if self.strPickle == 'False':
        #     logging.info('Delete tmp pickles')
        #     sp.call('rm {outdir}/Tmp/Pickle/*.pickle'.format(outdir=self.strOutSampleDir), shell=True)
        #
        # elif self.strSplit == 'False':
        #     logging.info('Delete splited input files')
        #     sp.call('rm {split_path}/*.fq'.format(split_path=self.strSplitPath), shell=True)


def Main():
    parser = OptionParser(
        'Indel search program for CRISPR CAS9 & CPF1\n<All default option> python2.7 Run_indel_searcher.py --pam_type Cas9 --pam_pos Forward')

    parser.add_option('-t', '--thread', default='1', type='int', dest='multicore',
                      help='multiprocessing number, recommendation:t<16')
    parser.add_option('-c', '--chunk_number', default='400000', type='int', dest='chunk_number',
                      help='split FASTQ, must be multiples of 4. file size < 1G recommendation:40000, size > 1G recommendation:400000')
    parser.add_option('-q', '--base_quality', default='20', dest='base_quality', help='NGS read base quality')
    parser.add_option('--gap_open', default='-10', type='float', dest='gap_open', help='gap open: -100~0')
    parser.add_option('--gap_extend', default='1', type='float', dest='gap_extend', help='gap extend: 1~100')
    parser.add_option('-i', '--insertion_window', default='4', type='int', dest='insertion_window',
                      help='a window size for insertions')
    parser.add_option('-d', '--deletion_window', default='4', type='int', dest='deletion_window',
                      help='a window size for deletions')
    parser.add_option('--pam_type', dest='pam_type', help='PAM type: Cas9 Cpf1')
    parser.add_option('--pam_pos', dest='pam_pos', help='PAM position: Forward Reverse')
    parser.add_option('--python', dest='python', help='The python path including the CRISPResso2')
    parser.add_option('--user', dest='user_name', help='The user name with no space')
    parser.add_option('--project', dest='project_name', help='The project name with no space')
    parser.add_option('--pickle', dest='pickle', default='False',
                      help='Dont remove the pickles in the tmp folder : True, False')
    parser.add_option('--split', dest='split', default='False',
                      help='Dont remove the split files in the input folder : True, False')
    parser.add_option('--classfied_FASTQ', dest='class_fastq', default='True',
                      help='Dont remove the ClassfiedFASTQ in the tmp folder : True, False')
    parser.add_option('--ednafull', dest='ednafull', help='The nucleotide alignment matrix')

    options, args = parser.parse_args()

    InstInitFolder = InitialFolder(options.user_name, options.project_name, os.path.basename(__file__))
    InstInitFolder.MakeDefaultFolder()
    InstInitFolder.MakeInputFolder()
    InstInitFolder.MakeOutputFolder()

    logging.basicConfig(format='%(process)d %(levelname)s %(asctime)s : %(message)s',
                        level=logging.DEBUG,
                        filename=InstInitFolder.strLogPath,
                        filemode='a')
    logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

    logging.info('Program start')
    if options.multicore > 15:
        logging.warning('Optimal treads <= 15')
    logging.info(str(options))

    with open(InstInitFolder.strProjectFile) as Sample_list:

        listSamples = Helper.RemoveNullAndBadKeyword(Sample_list)
        intProjectNumInTxt = len(listSamples)

        strInputProject = './Input/{user}/FASTQ/{project}'.format(user=options.user_name, project=options.project_name)

        @CheckProcessedFiles
        def RunPipeline(**kwargs):

            setGroup = set()
            for strSample in listSamples:

                tupSampleInfo = Helper.SplitSampleInfo(strSample)
                if not tupSampleInfo: continue
                strSample, strRef, strExpCtrl = tupSampleInfo
                setGroup.add(strExpCtrl)

                # Locating input files in the Reference and FASTQ directory
                InstRunner = clsExtractorRunner(strSample, strRef, options, InstInitFolder)

                # Chunking
                logging.info('SplitFile')
                InstRunner.SplitFile()

                # If there is no manually created "Reference.fa", make it for future processing
                # Need proper "Barcode.txt", "Reference_sequence.txt", "Target_region.txt"
                # logging.info('MakeReference')
                # InstRunner.MakeReference()

                # Generating a command for utilizing Python multiprocessing
                # This Runner program will execute Indel_searcher_crispresso_hash.py
                logging.info('MakeExtractorCmd')
                listCmd = InstRunner.MakeExtractorCmd()

                logging.info('RunMulticore')
                RunMulticore(listCmd, options.multicore)  ## from CoreSystem.py

                # Need Adaptation
                logging.info('MakeOutput')
                InstRunner.MakeOutput(listCmd, strSample)
                # logging.info('RunIndelFreqCalculator')
                # InstRunner.RunIndelFreqCalculator()

            # if setGroup == {'EXP', 'CTRL'}:
            #     InstRunner.IndelNormalization()
            # elif setGroup in [set(), set([]), set(['']), set([' '])]:
            #     pass
            # else:
            #     logging.error('The group category is not appropriate. : %s' % setGroup)
            #     logging.error('Please make sure your project file is correct.')
            #     logging.error('The group category must be Exp or Ctrl')
            #     raise Exception
            # """

        RunPipeline(InstInitFolder=InstInitFolder,
                    strInputProject=strInputProject,
                    intProjectNumInTxt=intProjectNumInTxt,
                    listSamples=listSamples,
                    logging=logging)

    logging.info('Program end')


# END:def


if __name__ == '__main__':
    Main()
