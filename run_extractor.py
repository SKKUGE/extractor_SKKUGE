#!/usr/bin/env python
import logging
import os
import sys
import argparse
import pathlib
from types import SimpleNamespace
import json

sys.path.insert(0, os.path.dirname(os.getcwd()))

from Core.CoreSystem import (
    SystemStructure,
    Helper,
    ReadBarcode,
    run_pipeline,
)


def main():
    parser = argparse.ArgumentParser(
        prog="extractor_SKKUGE",
        description="Counting sequence reads for each barcode from NGS rawdata, tested on Python v3.9 (tentative)",
        epilog="SKKUGE_DEV, 2023-01-02 ~",
    )
    parser.add_argument(
        "-t",
        "--thread",
        default="1",
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
        "--python_path", dest="python", type=pathlib.Path, help="The python path"
    )
    parser.add_argument(
        "--barcode",
        type=pathlib.Path,
        default="Barcodes/Barcode.txt",
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
        "--save_pickle",
        action="store_true",
        help="Dont remove the pickles in the tmp folder : True, False",
    )
    parser.add_argument(
        "--save_temp",
        action="store_true",
        help="Dont remove the split files in the input folder : True, False",
    )
    parser.add_argument(
        "--save_trace", action="store_false", help="Save the trace files : True, False"
    )
    parser.add_argument("--dev", action="store_true", help="development mode")
    args = parser.parse_args()

    read_b = ReadBarcode()
    excel = read_b.SelectFromExcel()
    barcode_data = read_b.UseExcel()
    args.user_name = input('Enter User name : ')
    args.project_name = input('Enter Project name : ')
    print(args.user_name, args.project_name, barcode_data)
    system_structure = SystemStructure(args.user_name, args.project_name)

    # Prepare logger
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    stream_handler = logging.StreamHandler()
    logger.addHandler(stream_handler)
    # TODO: follow the system directory structure

    samples = Helper.load_samples(system_structure.project_samples_path)

    # Add custom arguments
    args.system_structure = system_structure
    args.samples = samples
    args.logger = logger

    run_pipeline(SimpleNamespace(**vars(args)))
    
    with open(system_structure.log_dir, "w") as f:
        json.dump(args, f, indent=4)
    logger.info("Program end")


if __name__ == "__main__":
    main()
