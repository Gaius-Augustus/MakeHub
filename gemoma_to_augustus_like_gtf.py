#!/usr/bin/env python3

from inspect import currentframe, getframeinfo
import argparse
import re
import os
import errno
import shutil

__author__ = "Kathairna J. Hoff"
__copyright__ = "Copyright 2020. All rights reserved."
__license__ = "Artistic Licsense"
__version__ = "1.0.5"
__email__ = "katharina.hoff@uni-greifswald.de"
__status__ = "production"

parser = argparse.ArgumentParser(
    description='Convert a GeMoMa gff3 file to AUGUSTUS-like gtf format')
parser.add_argument('-i', '--input', required=True, type=str,
                    help='File with GeMoMa predictions in gff3 format')
parser.add_argument('-o', '--output', required=True, type=str,
                    help='File with GeMoMa predictions in AUGUSTUS-like gtf format')
parser.add_argument('-d', '--outdir', required=False, type=str, default='.',
                    help="output directory to write hub to")
args = parser.parse_args()

''' Check whether args.outdir exists and is writable '''

if not os.path.isdir(args.outdir):
    frameinfo = getframeinfo(currentframe())
    print('Error in file ' + frameinfo.filename + ' at line ' +
          str(frameinfo.lineno) + ': specified argument for --outdir ' +
          "(" + args.outdir + ") is not a directory!")
    exit(1)
elif not os.access(args.outdir, os.W_OK):
    frameinfo = getframeinfo(currentframe())
    print('Error in file ' + frameinfo.filename + ' at line ' +
          str(frameinfo.lineno) + ': specified argument for --outdir ' +
          "(" + args.outdir + ") is not writable!")
    exit(1)

tmp_dir = args.outdir + "/tmp-gemoma_convert/"

try:
    os.makedirs(tmp_dir)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise

''' Function that reformats AUGUSTUS gtf format to UCSC gtf format '''

def aug2ucsc_gtf(augustus_file, ucsc_file):
    try:
        with open(augustus_file, "r") as aug_handle:
            try:
                with open(ucsc_file, "w") as ucsc_handle:
                    for line in aug_handle:
                        if re.search(r'\t\S+\tCDS\t', line) or re.search(r'\t\S+\texon\t', line) or re.search(r'\t\S+\tstart_codon\t', line) or re.search(r'\t\S+\t\d+\'-UTR\t', line):
                            if re.search(r'\S+\t\S+\t\S+\t\d+\t\d+\t\S+\t\S+\t\S+\ttranscript_id\s\"\S+\";\sgene_id\s\"\S+\";', line):
                                first_part, feature, second_part, gid, txid = re.search(
                                    r'(\S+\t\S+\t)(\S+)(\t\d+\t\d+\t\S+\t\S+\t\S+\t)(transcript_id\s\"\S+\";)\s(gene_id\s\"\S+\";)', line).groups()
                            elif re.search(r'\S+\t\S+\t\S+\t\d+\t\d+\t\S+\t\S+\t\S+\tgene_id\s\"\S+\";\stranscript_id\s\"\S+\";', line):
                                first_part, feature, second_part, txid, gid = re.search(
                                    r'(\S+\t\S+\t)(\S+)(\t\d+\t\d+\t\S+\t\S+\t\S+\t)(gene_id\s\"\S+\";)\s(transcript_id\s\"\S+\";)', line).groups()
                            if re.search(r'5\'-UTR', feature):
                                feature = "5UTR"
                            elif re.search(r'3\'-UTR', feature):
                                feature = "3UTR"
                            ucsc_handle.write(
                                first_part + feature + second_part + gid + " " +
                                txid + "\n")
            except IOError:
                frameinfo = getframeinfo(currentframe())
                print('Error in file ' + frameinfo.filename + ' at line ' +
                      str(frameinfo.lineno) + ': ' + "Failed to open file " +
                      ucsc_file + " for writing!")
                quit(1)
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " + augustus_file +
              " for reading!")
        quit(1)

def gemoma2aug_gtf(gemoma_file, aug_like_file):
    try:
        with open(aug_like_file, "w") as aug_handle:
            try:
                with open(gemoma_file, "r") as gemoma_handle:
                    for line in gemoma_handle:
                        line = line.strip()
                        if re.search(r'\tCDS\t', line) or re.search(r'\tfive_prime_UTR\t', line) or re.search(r'\tthree_prime_UTR\t', line):
                            f0, f1, f2, f3, f4, f5, f6, f7, f8 = line.split()
                            if re.match(r"ID=", f8):
                                gff3_part = f8.split(";")
                                tid_part = gff3_part[1].split("=")
                                gid = tid = tid_part[1]
                            else:
                                thismatch = re.search(r'Parent=([^;]+);?', f8)
                                gid = tid = thismatch.group(1)
                            if re.search(r'_R\d+', gid):
                                gid_groups = re.search(
                                    r'(^\S+)_R\d+', gid).groups()
                                gid = gid_groups[0]
                            if f2 == 'five_prime_UTR':
                                f2 = '5\'-UTR'
                            elif f2 == 'three_prime_UTR':
                                f2 = '3\'-UTR'
                            aug_handle.write(f0 + "\t" + f1 + "\t" + f2 + "\t" +
                                             f3 + "\t" + f4 + "\t" + f5 + "\t" +
                                             f6 + "\t" + f7 + "\t" +
                                             "gene_id \"" + gid + "\"; " +
                                             "transcript_id \"" + tid +
                                             "\";\n")
            except IOError:
                frameinfo = getframeinfo(currentframe())
                print('Error in file ' + frameinfo.filename + ' at line ' +
                      str(frameinfo.lineno) + ': ' + "Failed to open file " + gemoma_file +
                      " for reading!")
                quit(1)
    except IOError:
        frameinfo = getframeinfo(currentframe())
        print('Error in file ' + frameinfo.filename + ' at line ' +
              str(frameinfo.lineno) + ': ' + "Failed to open file " + aug_like_file +
              " for writing!")
        quit(1)


aug_file = tmp_dir + "gemoma_aug.gtf"
gemoma2aug_gtf(args.input, aug_file)
ucsc_file = args.output
aug2ucsc_gtf(aug_file, ucsc_file)

shutil.rmtree(tmp_dir)