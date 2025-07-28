
"""
Module for detecting and processing organelle (mitochondrial) contamination in genomes.
"""

import io
import sys
import subprocess
import re
from Bio import SeqIO
import argparse
import remoVecSec.removeUtils as rU


def runMitoPipe(genomefile, dbfile):
    """
    Run mitochondrial screening with a database on a genome.

    Args:
        genomefile (str): Path to genome file.
        dbfile (str): Path to organelle database.

    Returns:
        file-like: Output of screening.
    """
    blastn = subprocess.Popen([
        'blastn',
        '-query', genomefile,
        '-db', dbfile,
        '-task', 'megablast',
        '-word_size', '28',
        '-best_hit_overhang', '0.1',
        '-best_hit_score_edge', '0.1',
        '-dust', 'yes',
        '-evalue', '0.0001',
        '-perc_identity', '98.6',
        '-soft_masking', 'True',
        '-outfmt', '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
    ], stdout=subprocess.PIPE)
    awk = subprocess.Popen([
        'awk',
        '"$4>=120"'
    ], stdin=blastn.stdout, stdout=subprocess.PIPE)
    return io.TextIOWrapper(awk.stdout, encoding='utf-8')


def parseMito(contaminant_output, dist):
    """
    Parse the output of the mitochondrial screening.

    Args:
        contaminant_output: Output of blastn.
        dist (int): Max distance for merging intervals.

    Returns:
        dict: Mitochondrial intervals per contig.
    """
    contint = {}
    for line in contaminant_output:
        decoded = line
        if re.match(r"^#", decoded):
            continue
        blastData = re.split(r"\s+", decoded)
        contname = blastData[0]
        interval = [int(x) for x in blastData[6:8]]
        if contname not in contint:
            contint[contname] = []
        contint[contname].append(interval)
    for contname in contint:
        contint[contname] = rU.fuseinterval(contint[contname], dist)
    return contint


def qcMito(contint, records):
    """
    Flag records as mitochondrial if coverage >75% of sequence length.

    Args:
        contint (dict): Mitochondrial intervals.
        records (dict): Genome records.

    Returns:
        dict: Regions to be removed per record.
    """
    tomodify = {}
    for record in records:
        modifications = {}
        if record in contint:
            recordlength = len(records[record].seq)
            allcontint = contint[record]
            coverage = sum(abs(interval[0] - interval[1]) + 1 for interval in allcontint)
            if coverage / recordlength > 0.75:
                modifications["remove"] = 1
            if modifications:
                tomodify[record] = modifications
    return tomodify



def main(arguments):
    parser = argparse.ArgumentParser(
        description="remoMito: Find organelle (mitochondrial) contamination in assembled genomes.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--genomefile', '-g', help="Genome file", type=str, required=True)
    parser.add_argument('--dbmito', '-m', help="The organelle database", type=str, required=False)
    parser.add_argument('--dist', '-d', help="Maximal distance for merging two intervals", type=int, required=False, default=50)
    args = parser.parse_args(arguments)
    dist = args.dist
    records = SeqIO.index(args.genomefile, "fasta")
    f = runMitoPipe(args.genomefile, args.dbmito)
    contint = parseMito(f, dist)
    tomodify = qcMito(contint, records)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
