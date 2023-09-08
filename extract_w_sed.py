import subprocess as sp
import shlex
from collections import defaultdict
import multiprocessing as mp

import pandas as pd

BARCODE = "./Barcodes/KRAS+Barcode_220919+FOR_MISSING_VARIANTS.txt"
FASTQ = "/mnt/P41/Repositories/extractor_SKKUGE/Input/JU/FASTQ/KRAS/D3_C1/20220314_1.extendedFrags.fastq"
CORE_NUM = 32
RESULTS = defaultdict()


def execution(scmd):
    name = scmd.pop()

    with open(FASTQ, 'r') as f:
        proc1 = sp.Popen(scmd, stdin=f, stdout=sp.PIPE)
        proc2 = sp.Popen(['wc', '-l'], stdin=proc1.stdout, stdout=sp.PIPE, stderr=sp.PIPE)

        proc1.stdout.close()  # Allow proc1 to receive a SIGPIPE if proc2 exits.
        out, err = proc2.communicate()

        RESULTS[name] = int(out.strip())


def run_multicore(lcmd):
    for scmd in lcmd:
        print(name, scmd)

    ps = mp.Pool(CORE_NUM)
    async_result = ps.map_async(execution, lcmd)
    ps.close()

    return async_result.get()


if __name__ == "__main__":

    lcmd = []

    with open(BARCODE, "r") as f:
        for line in f:
            name = line.strip().split(":")[0]
            needle = line.strip().split(":")[1]

            cmd = f"grep ATG"
            scmd = shlex.split(cmd)
            scmd.append(name)
            # print(f"{scmd}")
            lcmd.append(scmd)

    # TODO Debugging parallelism
    scmd = lcmd[0]
    execution(scmd)

    run_multicore(lcmd)

    pd.DataFrame(RESULTS, index=[0]).T.to_csv("./temp/test.csv")
