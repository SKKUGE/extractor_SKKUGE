#!/usr/bin/env python
import logging
import os
import sys
import argparse
import pathlib
from types import SimpleNamespace
import json
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.getcwd()))

from Core.CoreSystem import (
    SystemStructure,
    Helper,
    ReadBarcode,
    run_pipeline,
)

import toolkit.FLASH_multiple_files as flash
def main():
    parser = argparse.ArgumentParser(
        prog="extractor_SKKUGE",
        description="Counting sequence reads for each barcode from NGS rawdata, tested on Python v3.9 (tentative)",
        epilog="SKKUGE_DEV, 2023-01-02 ~",
    )
    parser.add_argument(
        "-t",
        "--thread",
        default="0",
        type=int,
        dest="multicore",
        help="multiprocessing number, recommendation:t<16",
    )
    parser.add_argument(
        "-c",
        "--chunk_size",
        default="100000",
        type=int,
        dest="chunk_size",
        help="split FASTQ, indicates how many reads will be in a splitted file. file size < 1G recommendation:10000, size > 1G recommendation:100000",
    )
    parser.add_argument(
        "--barcode",
        type=pathlib.Path,
        default="Barcode.txt",
        help="Barcode file path",
    )
    parser.add_argument(
        "-u", "--user", dest="user_name", type=str, help="The user name with no space"
    )
    parser.add_argument(
        "-p",
        "--project",
        dest="project_name",
        type=str,
        help="The project name with no space",
    )
    parser.add_argument(
        "-i",
        "--ui",
        dest="ui_type",
        type=int,
        default=1,
        help="UI type, 0 : for argument mode, 1 : interactive mode(default)",
    )

    args = parser.parse_args()

    # code for interactive UI design
    if args.ui_type : 
        args.user_name = input('Enter User name : ')
        args.project_name = input('Enter Project name : ')
        args.barcode = input('Enter barcode file name : ')
    
    system_structure = SystemStructure(args.user_name, args.project_name)
    
    if args.ui_type :
        print(args.user_name, args.project_name)
        rb = ReadBarcode(args.user_name, args.project_name, args.barcode)
        barcode = rb.UseExcel()
    flash.USER = args.user_name
    flash.PROJECT = args.project_name
    flash.INPUT_FILES_PATH = system_structure.user_dir
    d_list = flash.merge()

    # Prepare logger
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    stream_handler = logging.StreamHandler()
    logger.addHandler(stream_handler)
    # TODO: follow the system directory structure

    for date in d_list:
        samples = Helper.load_samples(system_structure.input_dir / str(date) / f'{args.project_name}.txt')
        
        # Add custom arguments
        args.system_structure = system_structure
        args.samples = samples
        args.date = str(date)
        args.logger = logger

        run_pipeline(SimpleNamespace(**vars(args)))
    logger.info("Program end")


if __name__ == "__main__":
    main()
