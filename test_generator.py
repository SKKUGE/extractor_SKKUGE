with open("test_barcode.txt", mode='w') as barcode, open("sequence.fastq", mode='w') as seq:
    nt = ['A', 'T', 'G', 'C']

    for n in nt:
        for i in range(10000000):
            seq.write(f"{n * 300}\n")

    for n in nt:
        barcode.write(f"{n}:{n}\n")

