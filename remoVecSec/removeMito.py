"""
This module is responsible for finding organelle in assembled genomes
"""

# import pprint
# from Bio.Seq import Seq
import io
import sys
import subprocess
import re
from Bio import SeqIO
import argparse
import removeUtils as rU
import pprint


def runMitoPipe(genomefile, dbfile):
    """
    Run contaminant screening with a contaminant database on a genome. Return a filehandle
    Args: genomfile(str) path of genomefile on which vecscreen is ran
          dbfile(str) path to the vector database used by vecscreen
    Returns a filehandle containing the output of vecscreen
    """
    # run blastn
    blastn = subprocess.Popen(['blastn',
                               '-query',  genomefile,
                               '-db',  dbfile,
                               '-task', 'megablast',
                               '-word_size', "28",
                               '-best_hit_overhang', "0.1",
                               '-best_hit_score_edge', "0.1",
                               '-dust',  "yes",
                               '-evalue', "0.0001",
                               '-perc_identity', "98.6",
                               '-soft_masking', "True",
                               '-outfmt', "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"],
                              stdout = subprocess.PIPE)
    # filter blast results
    awk = subprocess.Popen(['awk',
                             '"$4>=120"'],
                           stdin=blastn.stdout,
                           stdout=subprocess.PIPE)
    # return output of contaminant screening
    return io.TextIOWrapper(awk.stdout, encoding='utf-8')


def parseMito(contaminant_output, dist):
    """Parse the output of the mito screening
    Args: contaminant_output output of blastn

    Returns a dictionary with keys corresponding to the contigs
    containing contaminant and value corresponding to the coordinates of
    the vector inside the contig
    """
    # hash containing as key the vector, as value a list of intervals
    contint = {}
    for line in contaminant_output:
        # get line from file
        decoded = line
        # Regexp for results
        comments = re.compile("^#")
        # Regexp for finding lines starting with numbers
        # numberstart = re.compile("^[0-9]")
        # Regexp for no match
        # nohits = re.compile("No hits found")
        # Regexp for empty line
        emptyline = re.compile("^\s*$")
        # if we have a contig print / save its
        if comments.match(decoded):
            continue
        else:
            # parse results
            blastData = re.compile("\s+").split(decoded)
            # contig name
            contname = blastData[0]
            # contamination class
            interval = blastData[6:8]
            # integify parsed elements
            interval = [int(x) for x in interval]
            # classify interval boundaries
            if contname in contint:
                contint[contname].append(interval)
            else:
                contint[contname]=[]
                contint[contname].append(interval)
            
    # Fuse intervals for each contig and each contamination type
    for contname in contint:
        # allcontint contain for a given contig all merged intervals
        contint[contname] = rU.fuseinterval(
            contint[contname],
            dist)
        contint[contname]
    return contint


def qcMito(contint, records):
    """
    If the total coverage of mitochondrial matches from (3) is >75% of the
    sequence length then flag the sequence as being mitochindrial sequence
    to be excluded.
    Return: set of regions to be trimmed
    """
    # list of sequence to modify
    tomodify = {}
    # go through each record in the fasta file
    for record in records:
        modifications = {}
        if record in contint:
            # check if sum of intervals is smaller than 75% of the
            # record length
            recordlength = len(records[record].seq)
            allcontint = contint[record]
            coverage = 0
            for interval in allcontint:
                coverage = coverage + abs(interval[0]-interval[1])+1
            if coverage/recordlength > 0.75:
                modifications["remove"] = 1
            if len(modifications):
                tomodify[record] = modifications
    return tomodify


def main(arguments):
    parser = argparse.ArgumentParser(
        description="remoContaminant is a small wrapper script that allows to find\n"
        + "contamination in assembled genome. It takes as input a genome\n"
        + "file, a contaminant database and returns on stdout the location of contamination\n",
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--genomefile', '-g', help="Genome file",
                        type=str, required=True)
    parser.add_argument('--dbmito',
                        '-m',
                        help="The organelle database",
                        type=str, required=False)
    parser.add_argument('--dist',
                        '-d',
                        help="Maximal distance for merging two intervals",
                        type=int,
                        required=False,
                        default=50)
    args = parser.parse_args(arguments)
    dist = args.dist
    records = SeqIO.index(args.genomefile, "fasta")
    # Start analysis
    # Contaminant
    f = runMitoPipe(args.genomefile, args.dbmito)
    contint = parseMito(f, dist)
    tomodify = qcMito(contint, records)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
