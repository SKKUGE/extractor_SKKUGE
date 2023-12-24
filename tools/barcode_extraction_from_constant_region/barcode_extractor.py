# Adapted to ClonTracer library

import pathlib
import re
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count

import dask.bag as db
import numpy as np

# Hyperparameters
UPSTREAM_FLANKING_SEQUENCE = "GGCCAACATGAGGATCACCCATGTCTGCAGGGC"
DOWNSTREAM_FLANKING_SEQUENCE = "GGCCAACATGAGGATCACCCATGTCTGCAGGGC"

SEARCH_PATTERN = (
    rf"({UPSTREAM_FLANKING_SEQUENCE})(.{{28}})({DOWNSTREAM_FLANKING_SEQUENCE})"
)


def pattern_matcher(pat, seq):
    result = re.search(pat, seq, flags=re.IGNORECASE)
    return result.group(2) if result else np.nan


def extract_barcodes(src_file: pathlib.Path, dest_dir: pathlib.Path) -> None:
    bags = db.read_text(str(src_file), blocksize="256MB")

    # 2*n + 1
    # id, seq, +, qual
    bags = bags.filter(lambda x: re.match(r"^(A|T|G|C|a|t|g|c|N)", x))
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
    cwd = pathlib.Path(__file__).parent.resolve()
    src_dir = cwd / "src"
    src_files = pathlib.Path(src_dir).glob("*.fastq")
    dest_dir = cwd / "dest"

    if not src_dir.exists() or not dest_dir.exists():
        src_dir.mkdir(parents=True, exist_ok=True)
        dest_dir.mkdir(parents=True, exist_ok=True)
        raise FileNotFoundError(
            'Source or destination directory has generated. Please copy-and-paste input files in "src" directory.'
        )

    with ProcessPoolExecutor(max_workers=cpu_count()) as executor:
        for src_file in src_files:
            # extract_barcodes(src_file, dest_dir)  # debug
            executor.submit(extract_barcodes, src_file, dest_dir)
