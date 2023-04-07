import pathlib

UPSTREAM_FLANKING_SEQUENCE = "TCTTTTTACTGACTGCAGTCTGAGTCTGACAG"
DOWNSTREAM_FLANKING_SEQUENCE = "AGCAGAGCTACGCACTCTATGCTAGCTCGA"

SEARCH_PATTERN = rf"({UPSTREAM_FLANKING_SEQUENCE})(.*)({DOWNSTREAM_FLANKING_SEQUENCE})"


def barcode_formatter(x):
    return f"{x['i']}:{x[0]}"


def extract_barcodes(src_file: pathlib.Path, dest_dir: pathlib.Path) -> None:
    import numpy as np
    import pandas as pd

    seqs = []
    with open(src_file, "r") as f:
        # id, seq, +, qual
        lines = f.readlines()
        seqs = pd.Series(np.array(lines).reshape(-1, 4)[:, 1].tolist())

        query_result = seqs.str.contains(SEARCH_PATTERN, case=False)
        detected = seqs[query_result].drop_duplicates()
        detected = pd.DataFrame(detected)
        detected["i"] = [i for i in range(detected.shape[0])]
        detected.apply(lambda x: barcode_formatter(x), axis=1).to_csv(
            dest_dir / f"{src_file.stem}.csv", index=False, header=False
        )

        with open(dest_dir / f"{src_file.stem}_stat.csv", "w") as s:
            detected, total_read, detection_rate = (
                query_result.sum(),
                seqs.shape[0],
                query_result.sum() / seqs.shape[0],
            )
            s.write(f"Total read: {total_read}\n")
            s.write(f"Detected read: {detected}\n")
            s.write(f"Detection rate in the sequence pool: {detection_rate}\n")


if __name__ == "__main__":
    from concurrent.futures import ProcessPoolExecutor
    from multiprocessing import cpu_count

    src_files = pathlib.Path("./src/").glob("*.fastq")
    dest_dir = pathlib.Path("./dest/")

    with ProcessPoolExecutor(max_workers=1) as executor:
        for src_file in src_files:
            executor.submit(extract_barcodes, src_file, dest_dir).result()
