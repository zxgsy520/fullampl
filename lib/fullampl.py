#!/usr/bin/env python
#coding:utf-8

import os
import re
import sys
import logging
import argparse


LOG = logging.getLogger(__name__)

__version__ = "3.2.0"
__author__ = ("Xingguo Zhang",)
__email__ = "113178210@qq.com"
__all__ = []


def read_fasta(file):

    if file.endswith(".gz"):
        fp = gzip.open(file)
    elif file.endswith(".fasta") or file.endswith(".fa"):
        fp = open(file)
    else:
        raise Exception("%r file format error" % file)

    seqid = ''
    seq = ''

    for line in fp:
        if type(line) == type(b''):
            line = line.decode('utf-8')
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if seq!='' and seqid!='':
                yield seqid, seq
            seqid = line.strip('>')
            seq = ''
        else:
            seq += line

    if seq!='' and seqid!='':
        yield seqid, seq

    fp.close()


def produce_primer(rprimer, fprimer):

    data = {}
    primer = 'primer.fa'
    fp = open(primer, 'w')

    fp.write('>R\n%s\n>F\n%s' % (rprimer, fprimer))
    data['R'] = len(rprimer)
    data['F'] = len(fprimer)

    return primer, data


def do_blast(fasta, primer, blastn, identity):

    if not os.path.exists(blastn):
        raise Exception("blastn not found in %r" % blastn)

    cmd = '%s -task blastn-short -subject %s -query %s -perc_identity %s -max_target_seqs 2 -evalue 10 -gapopen 0 -gapextend 4 -reward 2 -penalty -2 -outfmt 6' % (blastn, primer, fasta, identity)
    result = os.popen(cmd)

    for n, line in enumerate(result):
        line = line.strip()

        if line.startswith("#") or not line:
            continue

        lines = line.split("\t")

        yield lines


def sort_seat(start, end):

    start = int(start)
    end = int(end)
    direction = "+"

    if end < start:
        start, end = end, start
        direction = "-"

    return start, end, direction


def reverse_complement(sequence):
    
    sequence = sequence[::-1].upper()
    sequence = sequence.replace('A', 't')
    sequence = sequence.replace('T', 'a')
    sequence = sequence.replace('C', 'g')
    sequence = sequence.replace('G', 'c')

    return sequence.upper()


def count_cat(plist, trim, minlen=500):

    old_ptype = ''
    old_dire = ''
    old_start = 1

    for start, end, direction, ptype in sorted(plist):
        if not old_ptype or old_ptype == ptype or old_dire == direction or (end-old_start+1) < minlen*0.75:
            old_ptype = ptype
            old_dire = direction
            old_start = start - 1
            if trim:
                old_start = end
            continue

        old_end = end
        if trim:
            old_end = start-1 

        yield old_start, old_end, old_ptype 
        old_ptype = ''


def remove_primer(fasta, rprimer, fprimer, blastn, identity, minlen, maxlen, trim):

    primer, data = produce_primer(rprimer, fprimer)
    read_dict = {}

    for line in do_blast(fasta, primer, blastn, identity):
        if int(line[3]) < data[line[1]]*0.9:
            continue

        if line[0] not in read_dict:
            read_dict[line[0]] = []
        start, end, direction  = sort_seat(line[6], line[7])
        pstart, pend, direction  = sort_seat(line[8], line[9])
        read_dict[line[0]].append([start, end, direction, line[1]])
    
    data = {}

    for seqid in read_dict:
        plist = read_dict[seqid]

        for start, end, direction in count_cat(plist, trim, minlen):
            rlen = end-start

            if rlen < minlen or rlen > maxlen:
                sys.stderr.write("%s\t%s\n" % (seqid, sorted(plist)))
                continue
            if seqid not in data:
                data[seqid] = []
            data[seqid].append([start, end, direction])

    for seqid, seq in read_fasta(fasta):
        line = seqid.split()
        seqid = line[0]
        desc = ''
        n = 0
        if len(line)>=2:
            desc = ' '.join(line[1::])
        if seqid not in data:
            continue

        for start, end, direction in data[seqid]:
            sequence = seq[start:end]
            
            if direction != "R":
                sequence = reverse_complement(sequence) 
            
            print('>%s/%s %s\n%s' % (seqid, n, desc, sequence))
            n += 1
            

def add_args(parser):

    parser.add_argument('fasta', help='')
    parser.add_argument('-R', '--rprimer', metavar='STR', type=str, default='TCCTCCGCTTATTGATATGC',
        help='Input the R-end primer sequence,default=TCCTCCGCTTATTGATATGC.')
    parser.add_argument('-F', '--fprimer', metavar='STR', type=str, default='TCCGTAGGTGAACCTGCGG',
        help='Input the R-end primer sequence,default=TCCGTAGGTGAACCTGCGG.')
    parser.add_argument('--minlen', metavar='INT', type=int, default=500,
        help='Filter the minimum read length, default=500.')
    parser.add_argument('--maxlen', metavar='INT', type=int, default=1000,
        help='Filter the maximum read length, default=1000.')
    parser.add_argument('--identity', metavar='INT', type=int, default=80,
        help='Set the identity of the comparison, default=80.')
    parser.add_argument('--no_trim', action='store_true',
        help='Input does not cut primers.')
    parser.add_argument('--blastn', metavar='FILE', type=str, default='/export/personal/software/software/blast/2.10.0/bin/blastn',
        help='Input the path of blastn.')

    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""

version: %s
contact:  %s <%s>\
    """ % (__version__, " ".join(__author__), __email__))

    parser = add_args(parser)
    args = parser.parse_args()
    if args.no_trim:
        trim = False
    else:
        trim = True

    remove_primer(args.fasta, args.rprimer, args.fprimer, args.blastn, args.identity, args.minlen, args.maxlen, trim)


if __name__ == "__main__":
    main()
