
"""
Module for detecting and processing contamination in assembled genomes.
"""

import io
import sys
import subprocess
import re
from Bio import SeqIO
import argparse
import remoVecSec.removeUtils as rU
import pprint

def runContaminantPipe(genomefile, dbfile):
    """
    Run contaminant screening with a contaminant database on a genome.

    Args:
        genomefile (str): Path to genome file.
        dbfile (str): Path to contaminant database.

    Returns:
        file-like: Output of contaminant screening.
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
        '-perc_identity', '90.0',
        '-outfmt', '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'
    ], stdout=subprocess.PIPE)
    awk = subprocess.Popen([
        'awk',
        '"($3>=98.0 && $4>=50)||($3>=94.0 && $4>=100)||($3>=90.0 && $4>=200)"'
    ], stdin=blastn.stdout, stdout=subprocess.PIPE)
    return io.TextIOWrapper(awk.stdout, encoding='utf-8')


def parseContaminant(contaminant_output, dist):
    """
    Parse the output of the contaminant screening.

    Args:
        contaminant_output: Output of blastn.
        dist (int): Max distance for merging intervals.

    Returns:
        dict: Contaminant intervals per contig and class.
    """
    contint = {}
    for line in contaminant_output:
        decoded = line
        if re.match(r"^#", decoded):
            continue
        blastData = re.split(r"\s+", decoded)
        contname = blastData[0]
        contclass = re.sub(r':.+', '', blastData[1])
        contclass = re.sub(r'.*VEC.*', 'VEC', contclass)
        interval = [int(x) for x in blastData[6:8]]
        if contname not in contint:
            contint[contname] = {}
        if contclass not in contint[contname]:
            contint[contname][contclass] = []
        contint[contname][contclass].append(interval)
    for contname in contint:
        allcontint = []
        for contclass in contint[contname]:
            contint[contname][contclass] = rU.fuseinterval(contint[contname][contclass], dist)
            allcontint.extend(contint[contname][contclass])
        allcontint = rU.fuseinterval(allcontint, dist)
        contint[contname]["allint"] = allcontint
    return contint


def qcContaminant(contint, records):
    """
    Given contaminant intervals and genome records, flag regions to be trimmed or removed.

    Args:
        contint (dict): Contaminant intervals.
        records (dict): Genome records.

    Returns:
        dict: Regions to be trimmed or removed per record.
    """
    tomodify = {}
    for record in records:
        if record in contint:
            recordlength = len(records[record].seq)
            modifications = {}
            allcontint = contint[record]["allint"]
            coverage = sum(abs(interval[0] - interval[1]) + 1 for interval in allcontint)
            if coverage / recordlength > 0.75:
                modifications["remove"] = 1
            contintEnd = rU.fuseinterval(allcontint, 100)
            if contintEnd and contintEnd[0][0] < 100:
                modifications["trim5"] = contintEnd[0][1]
            if contintEnd and contintEnd[-1][1] > recordlength - 100:
                modifications["trim3"] = contintEnd[-1][0]
            if modifications:
                tomodify[record] = modifications
    return tomodify



def main(arguments):
    parser = argparse.ArgumentParser(
        description="remoContaminant: Find contamination in assembled genomes.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--genomefile', '-g', help="Genome file", type=str, required=True)
    parser.add_argument('--dbcont', '-c', help="The contaminant database", type=str, required=False)
    parser.add_argument('--dist', '-d', help="Maximal distance for merging two intervals", type=int, required=False, default=50)
    args = parser.parse_args(arguments)
    dist = args.dist
    records = SeqIO.index(args.genomefile, "fasta")
    f = runContaminantPipe(args.genomefile, args.dbcont)
    contint = parseContaminant(f, dist)
    tomodify = qcContaminant(contint, records)
    pprint.pprint(tomodify)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
