#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import gzip
import logging
import argparse

LOG = logging.getLogger(__name__)

__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def complement(seq):

    cdict = {"A": "T",
        "T": "A",
        "G": "C",
        "C": "G"
    }

    seq = list(seq.upper())
    nseq = ""
    for i in seq:
        nseq += cdict[i]

    return nseq


def reverse_complement(seq):

    seq = seq[::-1]

    return complement(seq)


def seq2match(seq, rc=3):

    a = ""
    b = ""
    sm = ""

    for i in seq:
        if i == a:
            b += i
        else:
            if len(b) >= rc:
                sm += "%s{%d,%d}" % (a, len(b)-1, len(b)+1)
            else:
                sm += b
            a = i
            b = i
    if len(b) >= rc:
        sm += "%s{%d, %d}" % (a, len(b)-1, len(b)+1)
    else:
        sm += b

    return sm


def read_tsv(file):

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split()


def read_barcode(file):

    r = {}

    for sample, fseq, rseq in read_tsv(file):
        r[sample] = {}
        r[sample]["R"] = [seq2match(fseq), seq2match(reverse_complement(fseq))]
        r[sample]["F"] = [seq2match(rseq), seq2match(reverse_complement(rseq))]

    return r


def read_fasta(file):

    '''Read fasta file'''
    if file.endswith(".gz"):
        fp = gzip.open(file)
    elif file.endswith(".fasta") or file.endswith(".fa"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = []
    for line in fp:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue
        if line.startswith(">"):
            line = line
            if len(seq) == 2:
                yield seq
            seq = []
            seq.append(line.split()[0])
            continue
        if len(seq) == 2:
            seq[1] += line
        else:
            seq.append(line)

    if len(seq) == 2:
        yield seq
    fp.close()


def read_fastq(file):
    '''Read fastq file'''

    LOG.info("Reading message from %r" % file)
    if file.endswith("fastq.gz") or file.endswith(".fq.gz"):
        fp = gzip.open(file)
    elif file.endswith(".fastq") or file.endswith(".fq"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    seq = []
    for line in fp:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue
        if not seq:
            seq.append(line.split()[0])
            continue

        seq.append(line)
        if len(seq) == 4:
            yield seq
            seq = []

    fp.close()


def match_barcode(sequence, barcodes):

    barcode = []

    for i in barcodes:
        barcode += re.findall(i, sequence)
        if barcode:
            break

    return barcode


def match_barcodes(sequence, barcodes, flen=500):

    r = {}
    for i in barcodes:
        fseq = match_barcode(sequence[0:flen], barcodes[i])
        rseq = match_barcode(sequence[len(sequence)-flen::], barcodes[i])
        r[i] = [fseq, rseq]

    barcode = False
    rb = ""
    fb = ""

    if (r["R"][0] and r["F"][1]):
        barcode = True
        rb = ",".join(r["R"][0])
        fb = ",".join(r["F"][1])
    elif (r["R"][1] and r["F"][0]):
        barcode = True
        rb = ",".join(r["R"][1])
        fb = ",".join(r["F"][0])
    else:
        pass

    return barcode, rb, fb


def ontlampl(file, barcode, flen=500):

    r = read_barcode(barcode)

    if (file.endswith("fastq.gz") or file.endswith(".fq.gz") or
        file.endswith(".fastq") or file.endswith(".fq")):
        fh = read_fastq(file)
        format = "fastq"
    else:
        fh = read_fastq(file)
        format = "fasta"

    fd = {}

    for line in fh:
        for i in r:
            barcode, rb, fb = match_barcodes(line[1], r[i], flen)
            if barcode:
                print("%s\t%s\t%s\t%s" % (line[0], i, rb, fb))
                name = "%s.%s" % (i, format)
                if name not in fd:
                    fd[name] = open(name, "w")
                fd[name].write("%s\n" % "\n".join(line))
                break
    for i in fd:
        fd[i].close()
    return 0


def add_hlep_args(parser):

    parser.add_argument('input', metavar='FILE', type=str,
        help='Input sequencing reads.')
    parser.add_argument('-b','--barcode', metavar='FILE', type=str, required=True,
        help='Input barcode file.')
    parser.add_argument('-fl','--flen', metavar='INT', type=int, default=500,
        help='Barcode search range, default=500.')

    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''
name:
    ontlampl.py Split ont sequencing double-ended barcodes

attention:
    ontlampl.py reads.fasta -b barcode.txt > telomere.txt
barcode.txt file format:
1	CACAAAGACACCGACAACTTTCTT	TCACACGAGTATGGAAGTCGTTCT
2	ACAGACGACTACAAACGGAATCGA	TCTATGGGTCCCAAGAGACTCGTT
version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    ontlampl(args.input, args.barcode, args.flen)


if __name__ == "__main__":

    main()
