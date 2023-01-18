import subprocess
import shlex
import os
import sys
import pathlib
# example_cmd = "/flash HKK_220314__10_1.fastq.gz HKK_220314__10_2.fastq.gz -M 400 -m 10 -O -o 0314_10"

USER = ""
PROJECT = ""  # YYYYMMDD
FILENAMEPOS = 1  # relative position based on the list splitted by '_'
INPUT_FILES_PATH = pathlib.Path()
INPUT_FILES = """   # files separated by newline character
221221_HKK_1_1.fastq.gz
221221_HKK_1_2.fastq.gz
221221_HKK_2_1.fastq.gz
221221_HKK_2_2.fastq.gz
221221_HKK_3_1.fastq.gz
221221_HKK_3_2.fastq.gz
221221_HKK_4_1.fastq.gz
221221_HKK_4_2.fastq.gz
221221_HKK_5_1.fastq.gz
221221_HKK_5_2.fastq.gz
"""

# https://stackoverflow.com/questions/5389507/iterating-over-every-two-elements-in-a-list
def grouped(iterable, n):
    "s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ..."
    return zip(*[iter(iterable)] * n)

def merge():
    INPUT_FILES = [file for file in os.listdir(INPUT_FILES_PATH / PROJECT / "Input") if file.endswith(r'.fastq.gz')]
    print(INPUT_FILES_PATH)
    folders = open(INPUT_FILES_PATH / f'{PROJECT}.txt', 'w')
    waiting_queue = []
    for fwd, rev in grouped(INPUT_FILES, 2):
        waiting_queue.append((fwd, rev))

        print(f"Target loaded successfully: {waiting_queue[-1]}")

    # Generate child processes for FLASH
    for fwd, rev in waiting_queue:
        filename = fwd.split("_")[FILENAMEPOS]
        print(INPUT_FILES_PATH / PROJECT / "Input" / fwd)
        print(filename, file = folders)
        subprocess.run(
            shlex.split(
                f"mkdir -p {INPUT_FILES_PATH / PROJECT / 'Input' / filename}"
            )
        )
        cmd_input = shlex.split(
            f"./flash {INPUT_FILES_PATH / PROJECT / 'Input' / fwd} {INPUT_FILES_PATH / PROJECT / 'Input' / rev} -M 400 -m 1 -O -o {INPUT_FILES_PATH / PROJECT / 'Input' / filename / filename} -t 32"
        )

        subprocess.run(cmd_input)
    print('flash end')
    return