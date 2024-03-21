import argparse
import logging
import os
import sys
from types import SimpleNamespace
from Core.CoreSystem import (
    Helper,
    SystemStructure,
    run_pipeline,
)

sys.path.insert(0, os.path.dirname(os.getcwd()))


def main():
    """
    The main function of the extractor_SKKUGE program.

    This function parses command line arguments, sets up the system structure, logger, and samples,
    and then runs the pipeline to process the data.
    """

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
        "-u", "--user", dest="user_name", type=str, help="The user name with no space"
    )
    parser.add_argument(
        "-p",
        "--project",
        dest="project_name",
        type=str,
        help="The project name with no space\n Barcode directries are written in the project.txt in the User directory, tab separated",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="unique mutation test",
    )
    parser.add_argument(
        "--separator",
        dest="sep",
        type=str,
        help="Separator character for the barcode file. Default is ':'.",
        default=":",
    )
    # parser.add_argument(
    #     "-c",
    #     "--chunk_size",
    #     dest="chunk_size",
    #     type=int,
    #     help="barcode file chunk size",
    #     default=1,
    # )

    args = parser.parse_args()

    system_structure = SystemStructure(args.user_name, args.project_name)

    # Prepare logger
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    stream_handler = logging.StreamHandler()
    logger.addHandler(stream_handler)

    # Add custom arguments
    args.samples = Helper.load_samples(system_structure.project_samples_path)
    args.system_structure = system_structure
    args.logger = logger

    run_pipeline(SimpleNamespace(**vars(args)))
    logger.info("Program end")


if __name__ == "__main__":
    import time  #

    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))
