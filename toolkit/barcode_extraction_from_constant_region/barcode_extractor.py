# Adated to ClonTracer library

import pathlib

import dask.dataframe as dd
import dask.array as da
import dask.bag as db
import re
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from multiprocessing import cpu_count

UPSTREAM_FLANKING_SEQUENCE = "TCTTTTTACTGACTGCAGTCTGAGTCTGACAG"
DOWNSTREAM_FLANKING_SEQUENCE = "AGCAGAGCTACGCACTCTATGCTAGCTCGA"

SEARCH_PATTERN = (
    rf"({UPSTREAM_FLANKING_SEQUENCE})(.{{30}})({DOWNSTREAM_FLANKING_SEQUENCE})"
)


def pattern_matcher(pat, seq):
    result = re.search(pat, seq, flags=re.IGNORECASE)
    return result.group(2) if result else np.nan


def extract_barcodes(src_file: pathlib.Path, dest_dir: pathlib.Path) -> None:
    bags = db.read_text(str(src_file), blocksize="256MB")

    # 2*n + 1
    # id, seq, +, qual
    bags = bags.filter(lambda x: re.match(rf"^(A|T|G|C|a|t|g|c|N)", x))
    barcodes = bags.map(
        lambda x: pattern_matcher(pat=SEARCH_PATTERN, seq=x)
    ).to_dataframe(meta={0: "object"})

    with open(dest_dir / f"{src_file.stem}_stat.csv", "w") as s:
        detected, total_read = (
            barcodes.count().compute(),
            barcodes.shape[0].compute(),
        )
        detection_rate = detected / total_read

        s.write(f"File: {src_file.name}\n")
        s.write(f"Total read: {total_read}\n")
        s.write(f"Detected read: {detected}\n")
        s.write(f"Detection rate in the sequence pool: {detection_rate}\n")

    barcodes.dropna().drop_duplicates().reset_index(drop=True).compute().to_csv(
        dest_dir / f"{src_file.stem}_barcode.csv", index=True, header=False
    )


if __name__ == "__main__":
    src_files = pathlib.Path("./src/").glob("*.fastq")
    dest_dir = pathlib.Path("./dest/")

    with ProcessPoolExecutor(max_workers=cpu_count()) as executor:
        for src_file in src_files:
            # extract_barcodes(src_file, dest_dir)  # debug
            executor.submit(extract_barcodes, src_file, dest_dir)
