#!/usr/bin/python
# -*- coding: utf-8 -*-

# Original code is here, https://github.com/MyungjaeSong/Paired-Library.git
# This is the modified version of the program for academic uses in SKKUGE Lab

__author__ = "forestkeep21@naver.com"
__editor__ = "poowooho3@g.skku.edu"

import errno
import logging
import os
import pickle
import re
import sys
import time

BASE_DIR = os.path.dirname(sys.executable)

# for debug
BASE_DIR = os.path.dirname(os.path.realpath(__file__))


def help():
    """[Extracting nucleotide sequence from NGS result with specific barcodes]"""


def seq_validator(data):
    m = re.findall(r"^[A|a|T|t|C|c|G|g]+$", data)
    return m[0] if m else None


def count_line_in_file(file_name):
    count = 0
    for line in open(file_name, "r"):
        count += 1
    return count


def do(src_file_name, dest_file_name, sample_name, verbose=True):
    start_time = time.time()

    # 프로그램 진행율을 계산하기 위해 파일의 라인수를 센다.
    src_line_cnt = count_line_in_file(dest_file_name)
    if src_line_cnt == 0:
        print("File Not Found")
        raise
    current_cnt = 0
    extracted_line_index = []

    # 결과가 저장될 폴더 지정
    result_folder_name = os.path.join(BASE_DIR, "Results")
    # 결과가 저장될 폴더가 없다면 하나 생성
    if not os.path.exists(result_folder_name):
        os.makedirs(result_folder_name)

    # 결과가 저장될 샘플 폴더 지정 -- temp 폴더에 곧바로 저장
    sample_class = sample_name.split('.')[0]
    result_sample_dir = os.path.join(result_folder_name, f"temp/{sample_class}")

    try:
        os.makedirs(result_sample_dir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    # 총 결과 파일명 지정
    # src_file_name == '/home/dengarden/Documents/Repositories/extractor_SKKUGE/Barcode.txt'
    result_info_pkl_name = os.path.join(result_sample_dir, f"{sample_name}_result_info.trace")  # pickle path
    result_info_file_name = os.path.join(result_sample_dir, f"{sample_name}_result_info.txt")

    # file I/O -- txt
    result_info_txt = open(result_info_file_name, "w")

    # 추출할 시퀸스가 있는 파일을 읽어온다.
    data = [line.strip() for line in open(dest_file_name, "r") if seq_validator(line)]
    # 바코드가 있는 파일을 읽어온다.
    barcode_data = [line for line in open(src_file_name, "r")]

    try:
        used_data = []
        # 읽어온 바코드를 속도를 위해 모두 메모리에 올려놓고 분석을 시작한다.
        for barcode in barcode_data:
            used_lines = []

            # 바코드셋은 :를 구분자로 앞은 파일명, 뒤는 바코드로 되어있다.
            # barcode_set : ['Frag1_Pos_1_M_A', 'GCCACTGAATATAAACTT\n']
            barcode_set = barcode.split(":")
            if len(barcode_set) < 2:
                # 이름이 반드시 필요
                continue
            # 파일명에서 화이트 스페이스 삭제
            file_name = barcode_set[0].strip()
            # src_aa, dst_aa = file_name.split(sep="_")[-2], file_name.split(sep="_")[
            #     -1]  # Checking whether it is a WT barcode

            # 바코드가 valid한지 검증
            barcode = seq_validator(barcode_set[1].strip())

            # 대상이 되는 시퀸스들을 하나하나 분석한다.

            num_detected = 0

            another_mt_flag = False
            for idx, line in enumerate(data):
                # 비교를 위해 바코드, 대상 시퀸스 둘다 소문자로 변환하여 바코드가 대상 시퀸스 내에 존재하는지 검사
                if barcode.lower() in line.lower():
                    # processing buffer

                    # if src_aa == dst_aa:
                    #     # Full-sweep barcode chekcing
                    #     for barcode_complex in barcode_data:
                    #         try:
                    #             deep_barcode = barcode_complex.split(':')[1].strip()
                    #         except IndexError:
                    #             continue  # Barcode file error (example: new line character)
                    #
                    #         another_mt_flag = False
                    #         if deep_barcode == barcode:
                    #             # same barcode, pass
                    #             continue
                    #         if deep_barcode.lower() in line.lower():
                    #             # This line has another mutation, set the flag
                    #             another_mt_flag = True
                    #             break

                    if another_mt_flag:
                        # This line has another mutation, no count
                        continue

                    # 추출된 대상 시퀸스들을 pickle에 담기 위해 저장한다.
                    used_data.append((barcode, line))  # key, val
                    num_detected += 1
                    used_lines.append(idx)

            # 결과가 저장될 파일명 지정
            file_dir = os.path.join(result_sample_dir, f"{file_name}.txt")

            # 결과 파일 쓰기 시작 -- txt
            os.path.isfile(file_name)
            # writing a summary

            # barcode name 중복될 때, 찾은 시퀀스 파일이 삭제되는 문제점
            try:
                result_info_txt.write(f"{file_name}:{num_detected}\n")

            except Exception as e:
                print(e)
                print(
                    "Extraction has been done. But Making a result-info.txt is failed."
                )
                raise
            # TODO: 파일에 전부 옮겨담았다면 메모리에 올라간 전체 대상 시퀸스들에서 파일에 쓴 대상 시퀸스를 뺀다. -> bottleneck step; processing time increases about two times without it;  WT 중복이 너무 많이 발생..
            data = [i for j, i in enumerate(data) if j not in used_lines]

            # 프로그램 진행율 계산 부분
            current_cnt += 1

            if current_cnt % 100 == 0:
                progress_percentage = (float(current_cnt) / len(barcode_data)) * 100
                print("{} %".format(progress_percentage))
                print(f"{current_cnt} out of {len(barcode_data)} barcodes are processed.")

    except Exception as e:
        print(e)
        print("Extraction Failure.")
        raise

    # result_info_txt.write(f"total_reads:{src_line_cnt}")  # unnecessary
    result_info_txt.close()
    with open(result_info_pkl_name, 'wb') as fp:
        pickle.dump(used_data, fp)

    print("--- %s seconds elapsed ---" % (time.time() - start_time))


class clsParameter(object):

    def __init__(self):

        if len(sys.argv) > 1:
            self.strForwardFqPath = sys.argv[1]
            self.barcode = sys.argv[2]
            self.verbose = sys.argv[3].lower() == 'true'

        else:
            sManual = """
            Usage:

            """
            print(sManual)
            sys.exit()


if __name__ == "__main__":
    # Scarapped from Indel searcher
    InstParameter = clsParameter()

    src_file_name = os.path.join(BASE_DIR, InstParameter.barcode)

    if not os.path.isfile(src_file_name):
        print("File Not Found. Check it is in the src directory")
        raise Exception

    # print("Input sequence file name with extension: ")
    dest_file_name = InstParameter.strForwardFqPath
    # dest_file_name = os.path.join(BASE_DIR, dest_file_name)

    if not os.path.isfile(dest_file_name):
        print("File Not Found. Check it is in the src directory")
        raise Exception

    sample_name_token = InstParameter.strForwardFqPath.split('/')[-1].split('.')
    sample_name = '.'.join(sample_name_token)

    do(src_file_name, dest_file_name, sample_name, verbose=InstParameter.verbose)

    print("Extraction is completed.")

    logging.info('Program end : %s' % InstParameter.strForwardFqPath)
