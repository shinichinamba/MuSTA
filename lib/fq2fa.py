#!/usr/bin/python2
"""
Usage:
    fq2fa.py fastq/a_file(.gz)
"""

import sys

def fastq2fasta(f1):
    while True:
        line = f1.readline().strip()
        if line == "":
            break
        print '>'+line[2:]
        print f1.readline().strip()
        for i in xrange(2):
            f1.readline()

def fasta2fasta(f1):
    line1 = f1.readline().strip()
    while True:
        if line1 == "":
            break
        print line1
        while True:
            line1 = f1.readline().strip()
            if line1 == "" or line1[0] == ">":
                break
            print line1

if __name__ == '__main__':
    try:
        file1 = sys.argv[1]
    except:
        print __doc__
        sys.exit(1)

    if file1[-2:] == "gz":
        import gzip
        with gzip.open(file1) as f1:
            if file1[-8:] == "fastq.gz" or file1[-5:] == "fq.gz":
                fastq2fasta(f1)
            else:
                fasta2fasta(f1)
    else:
        with open(file1) as f1:
            if file1[-5:] == "fastq" or file1[-2:] == "fq":
                fastq2fasta(f1)
            else:
                fasta2fasta(f1)
