#!/usr/bin/python2
"""
Usage:
    interleave2fa.py fasta/q_file1(.gz) fasta/q_file2(.gz)
"""

import sys

def interleave_fastq(f1, f2):
    while True:
        line = f1.readline().strip()
        if line == "":
            break
        print line
        for i in xrange(3):
            print f1.readline().strip()
        for i in xrange(4):
            print f2.readline().strip()

def interleave_fastq2fasta(f1, f2):
    while True:
        line = f1.readline().strip()
        if line == "":
            break
        print '>'+line[2:]
        print f1.readline().strip()
        for i in xrange(2):
            f1.readline()
        print '>'+f2.readline().strip()[2:]
        print f2.readline().strip()
        for i in xrange(2):
            f2.readline()

def interleave_fasta(f1, f2):
    line1 = f1.readline().strip()
    line2 = f2.readline().strip()
    while True:
        if line1 == "":
            break
        print line1
        while True:
            line1 = f1.readline().strip()
            if line1 == "" or line1[0] == ">":
                break
            print line1
        print line2
        while True:
            line2 = f2.readline().strip()
            if line2 == "" or line2[0] == ">":
                break
            print line2

if __name__ == '__main__':
    try:
        file1 = sys.argv[1]
        file2 = sys.argv[2]
    except:
        print __doc__
        sys.exit(1)

    if file1[-2:] == "gz":
        import gzip
        with gzip.open(file1) as f1:
            with gzip.open(file2) as f2:
                if file1[-8:] == "fastq.gz" or file1[-5:] == "fq.gz":
                    interleave_fastq2fasta(f1, f2)
                else:
                    interleave_fasta(f1,f2)
    else:
        with open(file1) as f1:
            with open(file2) as f2:
                if file1[-5:] == "fastq" or file1[-2:] == "fq":
                    interleave_fastq2fasta(f1, f2)
                else:
                    interleave_fasta(f1,f2)
