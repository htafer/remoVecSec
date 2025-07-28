
"""
Module for detecting and processing vector contamination in genomes.
"""

import io
import sys
import subprocess
import re
from Bio import SeqIO
import logging
import argparse
import remoVecSec.removeUtils as rU


def runVecSecPipe(genomefile, dbfile):
    """
    Run vecscreen with a vector database on a genome.

    Args:
        genomefile (str): Path to genome file.
        dbfile (str): Path to vecscreen database.

    Returns:
        file-like: Output of vecscreen.
    """
    cat = subprocess.Popen(['cat', genomefile], stdout=subprocess.PIPE)
    vecsec = subprocess.Popen(['vecscreen', '-f 3', '-d ' + dbfile], stdin=cat.stdout, stdout=subprocess.PIPE)
    return io.TextIOWrapper(vecsec.stdout, encoding='utf-8')


def parseVecSec(vecsec_output, dist):
    """
    Parse the output of vecscreen.

    Args:
        vecsec_output: Output of vecscreen.
        dist (int): Max distance for merging intervals.

    Returns:
        dict: Vector intervals per contig.
    """
    vectint = {}
    for line in vecsec_output:
        decoded = line
        if re.match(r">Vector", decoded):
            contigname = re.split(r"\s+", decoded)[1]
            currentline = vecsec_output.readline()
            if re.match(r"No hits found", currentline):
                logging.info(f"skipping {decoded} since no hit\n")
                continue
            intervals = []
            currentline = vecsec_output.readline()
            while currentline:
                if re.match(r"^\s*$", currentline):
                    break
                if re.match(r"^[0-9]", currentline):
                    tempintervals = re.split(r"\s+", currentline)[:2]
                    tempintervals = [int(x) for x in tempintervals]
                    intervals.append(tempintervals)
                currentline = vecsec_output.readline()
            vectint[contigname] = rU.fuseinterval(intervals, dist)
    return vectint


def qcVector(vectint, records):
    """
    Given vector intervals and genome records, flag regions to be trimmed.

    Args:
        vectint (dict): Vector intervals.
        records (dict): Genome records.

    Returns:
        dict: Regions to be trimmed per record.
    """
    tomodify = {}
    for record in records:
        if record in vectint:
            modifications = {}
            fuseintervals = vectint[record]
            recordlength = len(records[record].seq)
            if fuseintervals and fuseintervals[0][0] < 100:
                modifications["trim5"] = fuseintervals[0][1]
            elif fuseintervals and fuseintervals[-1][1] > recordlength - 100:
                modifications["trim3"] = fuseintervals[-1][0]
            if modifications:
                tomodify[record] = modifications
    return tomodify


def main(arguments):
    parser = argparse.ArgumentParser(
        description="remoVecSec: Remove vector contamination in assembled genomes.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--genomefile', '-g', help="Genome file", type=str, required=True)
    parser.add_argument('--dbvec', '-v', help="The vecscreen database", type=str, required=False)
    parser.add_argument('--dist', '-d', help="Maximal distance for merging two intervals", type=int, required=False, default=50)
    args = parser.parse_args(arguments)
    records = SeqIO.index(args.genomefile, "fasta")
    dist = args.dist
    f = runVecSecPipe(args.genomefile, args.dbvec)
    vectint = parseVecSec(f, dist)
    rU.correctfasta(vectint, records)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
    """

    tomodify = {}

    for record in records:

        if record in vectint:
