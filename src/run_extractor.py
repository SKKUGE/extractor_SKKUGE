import argparse
import logging
import pathlib
import uuid
from types import SimpleNamespace

from extractor.core_system import (
    run_pipeline,
)

# sys.path.insert(0, os.path.dirname(os.getcwd()))


def runner():
    session_uid = uuid.uuid4()
    parser = argparse.ArgumentParser(
        prog="extractor_SKKUGE",
        description="Counting sequence reads for each barcode from fastq files",
        epilog="SKKUGE_DEV, 2023-01-02 ~",
    )
    # parser.add_argument(
    #     "-t",
    #     "--thread",
    #     default="1",
    #     type=int,
    #     dest="multicore",
    #     help="multiprocessing number, recommendation:t<16",
    # )
    # parser.add_argument(
    #     "-c",
    #     "--chunk_size",
    #     default="100000",
    #     type=int,
    #     dest="chunk_size",
    #     help="split FASTQ, indicates how many reads will be in a splitted file. file size < 1G recommendation:10000, size > 1G recommendation:100000",
    # )
    # parser.add_argument(
    #     "-u", "--user", dest="user_name", type=str, help="The user name with no space"
    # )
    # parser.add_argument(
    #     "-p",
    #     "--project",
    #     dest="project_name",
    #     type=str,
    #     help="The project name with no space\n Barcode directries are written in the project.txt in the User directory, tab separated",
    # )
    parser.add_argument(
        "--sample_list",
        dest="sample_list",
        type=pathlib.Path,
        help="The list describing the sample name and barcode",
    )
    parser.add_argument(
        "-i",
        "--fastq_path",
        dest="fastq_path",
        type=pathlib.Path,
        help="The directory path that containing fastq files",
    )
    parser.add_argument(
        "-q",
        "--barcodes",
        dest="barcode_path",
        type=pathlib.Path,
        help="The path of barcode file",
    )
    parser.add_argument(
        "--separator",
        dest="sep",
        type=str,
        help="Separator character for the barcode file. Default is ','.",
        default=",",
    )
    parser.add_argument(
        "-o",
        "--output_path",
        dest="output_path",
        type=pathlib.Path,
        help="The directory path of output file",
    )
    parser.add_argument(
        "-t",
        "--temp",
        dest="temp_path",
        type=pathlib.Path,
        help="The directory path of temporary parquet files to be written",
        default=pathlib.Path(f"./temp/{session_uid}"),
    )

    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="unique mutation test",
    )
    parser.add_argument(
        "--sample-replacement",
        dest="sample_replacement",
        action="store_true",
        help="Sample with replacement option; Default: sample without replacement",
    )
    args, unknown = parser.parse_known_args()

    # Path check
    if (
        not args.sample_list
        and not args.fastq_path
        and not args.barcode_path
        and not args.output_path
        and not args.temp_path
    ):
        parser.print_help()
        raise FileNotFoundError("Please check the path of the input files")

    if args.sample_list and not args.sample_list.exists():
        raise FileNotFoundError(f"{args.sample_list} does not exist")
    if args.fastq_path and not args.fastq_path.exists():
        raise FileNotFoundError(f"{args.fastq_path} does not exist")
    if args.barcode_path and not args.barcode_path.exists():
        raise FileNotFoundError(f"{args.barcode_path} does not exist")
    if args.output_path and not args.output_path.exists():
        pathlib.Path.mkdir(args.output_path, parents=True)
    if args.temp_path and not args.temp_path.exists():
        pathlib.Path.mkdir(args.temp_path, parents=True)

    # system_structure = SystemStructure(args.user_name, args.project_name)

    # Prepare logger
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    stream_handler = logging.StreamHandler()
    logger.addHandler(stream_handler)
    # TODO: follow the system directory structure

    # samples_and_barcodes = Helper.load_samples(system_structure.project_samples_path)

    # Add custom arguments
    # args.system_structure = system_structure
    # args.samples = samples_and_barcodes
    args.logger = logger
    logger.info("Program start with {}".format(session_uid))
    run_pipeline(SimpleNamespace(**vars(args)))
    logger.info("Program end")


if __name__ == "__main__":
    runner()
