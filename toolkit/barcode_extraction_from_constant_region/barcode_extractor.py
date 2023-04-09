import pathlib

import dask.dataframe as dd
import dask.array as da
import dask.bag as db
import re
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import cpu_count

UPSTREAM_FLANKING_SEQUENCE = "TCTTTTTACTGACTGCAGTCTGAGTCTGACAG"
DOWNSTREAM_FLANKING_SEQUENCE = "AGCAGAGCTACGCACTCTATGCTAGCTCGA"

SEARCH_PATTERN = rf"({UPSTREAM_FLANKING_SEQUENCE})(.*)({DOWNSTREAM_FLANKING_SEQUENCE})"


# @dask.delayed
# def parse(block):
#     return np.genfromtxt(io.BytesIO(*block), usecols=0)


def extract_barcodes(src_file: pathlib.Path, dest_dir: pathlib.Path) -> None:
    bags = db.read_text(str(src_file), blocksize=4**10)
    arr = da.stack([bag for bag in bags], axis=0).reshape(-1, 4)

    # id, seq, +, qual
    seqs = dd.from_dask_array(arr[:, 1])

    query_result = seqs.str.extract(SEARCH_PATTERN, flags=re.IGNORECASE)
    with open(dest_dir / f"{src_file.stem}_stat.csv", "w") as s:
        detected, total_read, detection_rate, unique = (
            query_result[1].count().compute(),
            seqs.shape[0].compute(),
            query_result[1].count().compute() / seqs.shape[0].compute(),
            query_result[1].nunique().compute(),
        )
        s.write(f"File: {src_file.name}\n")
        s.write(f"Unique barcodes found: {unique}\n")
        s.write(f"Total read: {total_read}\n")
        s.write(f"Detected read: {detected}\n")
        s.write(f"Detection rate in the sequence pool: {detection_rate}\n")

    query_result.drop_duplicates(1)[1].dropna().compute().to_csv(
        dest_dir / f"{src_file.stem}_barcode.csv", index=True, header=False
    )


if __name__ == "__main__":
    src_files = pathlib.Path("./src/").glob("*.fastq")
    dest_dir = pathlib.Path("./dest/")

    with ProcessPoolExecutor(max_workers=1) as executor:
        for src_file in src_files:
            extract_barcodes(src_file=src_file, dest_dir=dest_dir)

            # executor.submit(extract_barcodes, src_file, dest_dir).result()
