from Bio import SeqIO, SeqRecord, Seq
from tqdm import tqdm

import gzip
import subprocess
import shlex

# example_cmd = "/flash HKK_220314__10_1.fastq.gz HKK_220314__10_2.fastq.gz -M 400 -m 10 -O -o 0314_10"

DATE = "20221028"  # YYYYMMDD
FILENAMEPOS = 0  # relative position based on the list splitted by '_'
INPUT_FILES = """   # files separated by newline character
43_S43_L001_R2_001.fastq.gz
50_S50_L001_R2_001.fastq.gz
53_S53_L001_R2_001.fastq.gz
58_S58_L001_R2_001.fastq.gz
59_S59_L001_R2_001.fastq.gz
61_S61_L001_R2_001.fastq.gz
"""

if __name__ == "__main__":
    files = INPUT_FILES.split("\n")[1:-1]  # Additional '\n' in the head

    waiting_queue = []
    for seq in files:
        waiting_queue.append(seq)

        print(f"Target loaded successfully: {waiting_queue[-1]}")

    obj = {}
    for filename in waiting_queue:
        print(f"Processing {filename}...")
        with gzip.open(filename, "rt") as handle:
            temp = []
            for record in SeqIO.parse(handle, "fastq"):
                record.seq = record.seq.reverse_complement()
                temp.append(record)
            
            SeqIO.write(temp, f"{filename}_rev.fastq", "fastq")
            # obj[filename] = re

    # for f, records in obj.items():
