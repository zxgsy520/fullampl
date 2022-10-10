#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import gzip
import logging
import argparse
import operator

from multiprocessing import Pool

LOG = logging.getLogger(__name__)

__version__ = "2.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_tsv(file, sep=None):

    if file.endswith(".gz"):
        fh = gzip.open(file)
    else:
        fh = open(file)

    for line in fh:
        if isinstance(line, bytes):
            line = line.decode("utf-8")
        line = line.strip()

        if not line:
            continue

        yield line.split(sep)

    fh.close()


def read_fasta(file):

    '''Read fasta file'''
    if file.endswith(".gz"):
        fp = gzip.open(file)
    elif file.endswith(".fasta") or file.endswith(".fa") or file.endswith(".faa"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    r = ""
    for line in fp:
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line = line.strip()

        if not line:
            continue
        if line.startswith(">"):
            if r:
                yield r.split("\n", 1)
            r = "%s\n" % line.strip(">")
            continue
        r += line.upper()

    if r:
        yield r.split("\n", 1)
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
            seq.append(line.strip("@").split()[0])
            continue

        seq.append(line)
        if len(seq) == 4:
            yield seq
            seq = []

    fp.close()


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


def read_bases(database):

    data = {}
    rseq = {}
    header = ""
    n = 0

    for line in read_tsv(database, "\t"):
        if line[0].startswith("#"):
            header = line
            continue
        n += 1
        data[n] = line
        seq = line[1]
        rcseq = reverse_complement(line[1])
        rseq[n] = [seq, rcseq]

    return data, rseq, header


def seq_count_bases(qseq, seqdb):

   r = {}

   for i in seqdb:
       seq, rcseq = seqdb[i]
       r[i] = qseq.count(seq)
       if seq == rcseq:
           continue
       r[i] += qseq.count(rcseq)

   return r


def run_single_match_bases(file, seqdb, thread=4, splitlen=20000):

    if file.endswith("fastq.gz") or file.endswith(".fq.gz") or file.endswith(".fastq") or file.endswith(".fq"):
        fp = read_fastq(file)
    else:
        fp = read_fasta(file)

    results = []
    n = 0
    pool = Pool(processes=thread) #并行匹配计算
    temp = ""
    tn = 0
    for line in fp:
        n += 1
        temp += "N%s" % line[1]
        tn += 1
        if tn >= splitlen:
            results.append(pool.apply_async(seq_count_bases, (temp, seqdb)))
            tn = 0
            temp = ""
    if tn >= 1:
        results.append(pool.apply_async(seq_count_bases, (temp, seqdb)))

    pool.close()
    pool.join()

    r = {}
    for temp1 in results:
        temp2 = temp1.get()
        for i in temp2:
            if i not in r:
               r[i] = 0
            r[i] += temp2[i]

    return r, n


def unimatch_bases(database, pair, single1, single2, thread=4, splitlen=20000):

    data, rseq, header = read_bases(database)

    data1 = {}
    data2 = {}
    data3 = {}
    prnum = 0.001
    srnum = 0.001

    if pair:
        data1, prnum = run_single_match_bases(pair, rseq, thread, splitlen)
    if single1 or single2:
        data2, srnum1 = run_single_match_bases(single1, rseq, thread, splitlen)
        data3, srnum2 = run_single_match_bases(single2, rseq, thread, splitlen)
        srnum = (srnum1+srnum2)/2
    header = header + ["COMBINED_COUN", "COMBINED_FREQUENCY(%)", "NOTCOMBINED_COUNT", "NOTCOMBINED_FREQUENCY(%)", "ALL_COUNT", "ALL_FREQUENCY(%)"]
    print("\t".join(header))
    for i in sorted(data.keys()):
        line = data[i]
        pn = 0
        if i in data1:
            pn = data1[i]

        sn = 0
        if i in data2:
           sn += data2[i]
        if i in data3:
           sn += data3[i]

        pnr = pn*100.0/prnum
        snr = sn*100.0/srnum
        tn = pn + sn
        tnr = pnr + snr
        temp = "{0:,}\t{1:.8f}\t{2:,}\t{3:.8f}\t{4:,}\t{5:.8f}".format(pn, pnr, sn, snr, tn, tnr)
        print("%s\t%s" % ("\t".join(line), temp))

    return 0


def add_hlep_args(parser):

    parser.add_argument("input", metavar="FILE", type=str,
        help="Input mutation library(human.lib).")
    parser.add_argument("--pair", metavar="FILE", type=str,  default="",
        help="Input the read sequence after double-end linkage(extendedFrags.fastq).")
    parser.add_argument("-s1", "--single1", metavar="FILE", type=str,  default="",
        help="Input reads with no links at both ends(notCombined_1.fastq).")
    parser.add_argument("-s2", "--single2", metavar="FILE", type=str,  default="",
        help="Input reads with no links at both endse(notCombined_2.fastq).")
    parser.add_argument("-t", "--thread", metavar="INT", type=int,  default=4,
        help="Input the number of threads to use(default=4).")
    parser.add_argument("-sl", "--splitlen", metavar="INT", type=int,  default=200000,
        help="Sets the number of sequences split per thread(default=200000).")

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
    unimatch_bases: Statistical mutation rate based on mutation database detection
attention:
    unimatch_bases human.lib --pair extendedFrags.fastq >stat_mutation.tsv
    unimatch_bases human.lib --single1 notCombined_1.fastq --single2 notCombined_2.fastq >stat_mutation.tsv
    unimatch_bases human.lib --pair extendedFrags.fastq --single1 notCombined_1.fastq --single2 notCombined_2.fastq >stat_mutation.tsv

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    unimatch_bases(args.input, args.pair, args.single1, args.single2, args.thread, args.splitlen)


if __name__ == "__main__":

    main()
