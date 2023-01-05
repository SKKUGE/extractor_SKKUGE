#!/usr/bin/python
# -*- coding: utf-8 -*-

# Original code is here, https://github.com/MyungjaeSong/Paired-Library.git
# This is the modified version of the program for academic uses in SKKUGE Lab

__author__ = "forestkeep21@naver.com"
__editor__ = "poowooho3@g.skku.edu"

import os
import pickle
import re
import sys
import time
import pathlib

import pandas as pd
from tqdm import tqdm
import skbio

from Core.CoreSystem import SystemStructure

BASE_DIR = os.path.dirname(sys.executable)

# for debug
BASE_DIR = os.path.dirname(os.path.realpath(__file__))


def seq_validator(data):
    m = re.findall(r"^[A|a|T|t|C|c|G|g]+$", data)
    return m[0] if m else None


def count_line_in_file(file_name):
    count = 0
    for line in open(file_name, "r"):
        count += 1
    return count


def extract_read_cnts(
    sample_name: str,
    sequence_file: pathlib.Path,
    barcode_file: pathlib.Path,
    system_structure: SystemStructure,
) -> pd.DataFrame:
    # df index == barcode, column == read count

    # Load barcode file
    extraction_results = pd.read_csv(
        barcode_file, sep=":", header=None, names=["Gene", "Barcode"]
    ).set_index(
        "Gene"
    )  # TODO: tentative design
    extraction_results["Reads"] = 0

    # Load a splitted sequencing result using high-level I/O

    seqs = skbio.io.read(
        sequence_file, format="fastq", verify=False, variant="illumina1.8"
    )

    data = [line.strip() for line in open(sequence_file, "r") if seq_validator(line)]

    result_info_txt = open(result_info_file_name, "w")

    try:
        used_data = []
        # 읽어온 바코드를 속도를 위해 모두 메모리에 올려놓고 분석을 시작한다.

        # TODO: Possible C integration part
        for barcode in tqdm(barcode_data):
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

            # another_mt_flag = False
            for idx, line in enumerate(data):
                # 비교를 위해 바코드, 대상 시퀸스 둘다 소문자로 변환하여 바코드가 대상 시퀸스 내에 존재하는지 검사

                # TODO: debug
                # if barcode_set[0].strip().split('_')[2] == "7" and len(line) == 278:
                #     print(f">{barcode}\n{line}")

                if barcode.lower() in line.lower():
                    # processing buffer

                    # if src_aa == dst_aa:
                    #     # Full-sweep barcode checking
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

                    # if another_mt_flag:
                    # This line has another mutation, no count
                    # continue

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
            # TODO: debug
            data = [i for j, i in enumerate(data) if j not in used_lines]

    except Exception as e:
        print(e)
        print("Extraction Failure.")
        raise

    return extraction_results


def main(*args) -> pd.DataFrame:
    (sample, sequence, barcode, system_structure, logger) = args[0]

    start = time.time()
    rval = extract_read_cnts(sample, sequence, barcode, system_structure)
    end = time.time()

    logger.info(f"Extraction for {sample} is done. {end - start}s elapsed.")

    return rval
