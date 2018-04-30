"""This module is responsible for finding contamination in assembled genomes
"""

# import pprint
# from Bio.Seq import Seq
import io
import sys
import subprocess
import re
from Bio import SeqIO
import argparse
import remoVecSec.removeUtils as rU
import pprint
from functools import reduce

def runContaminantPipe(genomefile, dbfile):
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
                               '-perc_identity', "90.0",
                               '-outfmt', "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"],
                              stdout = subprocess.PIPE)
    # filter blast results
    awk = subprocess.Popen(['awk',
                             '"($3>=98.0 && $4>=50)||($3>=94.0 && $4>=100)||($3>=90.0 && $4>=200)"'],
                           stdin=blastn.stdout,
                           stdout=subprocess.PIPE)
    # return output of contaminant screening
    return io.TextIOWrapper(awk.stdout, encoding='utf-8')


def parseContaminant(contaminant_output, dist):
    """Parse the output of the contaminant screening
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
            contclass  = blastData[1]
            interval = blastData[6:8]
            # integify parsed elements
            interval = [int(x) for x in interval]
            # classify interval boundaries
            contclass = re.sub(r':.+', '', contclass)
            # Contaminant are divided into classes
            contclass = re.sub(r'.*VEC.*', 'VEC', contclass)
            # save interval in dictionary of dictionary with keys contname and contclass
            if contname in contint:
                if contclass in contint[contname]:
                    contint[contname][contclass].append(interval)
                else:
                    contint[contname][contclass] = []
                    contint[contname][contclass].append(interval)
            else:
                contint[contname] = {}
                contint[contname][contclass]=[]
                contint[contname][contclass].append(interval)
    # Fuse intervals for each contig and each contamination type
    for contname in contint:
        # allcontint contain for a given contig all merged intervals
        allcontint = []
        for contclass in contint[contname]:
            contint[contname][contclass] = rU.fuseinterval(
                contint[contname][contclass],
                dist)
            # We extend the list allcontint
            # list.
            # https://stackoverflow.com/questions/252703/difference-between-append-vs-extend-list-methods-in-python
            allcontint.extend(contint[contname][contclass])
        allcontint = rU.fuseinterval(allcontint, dist)
        contint[contname]["allint"] = allcontint
    return contint


def qcContaminant(contint, records):
    """Given a set of contaminant interval and a genome to correct flag the
    region that should be trimmed. The rules are: Contaminant matches
    from (1) are merged if they are from the same class of sequence
    (VECTOR, E.coli, IS, PHG) and they overlap or are separated by 50
    bases or less.

    If the total coverage of contaminant matches from (1) is >75% of the
    sequence length then flag the sequence as a contaminant to be
    excluded.
    If the contaminant is classed as VECTOR, E.coli, IS:./, PERM:./ or
    PHG:* and the contaminant location is within 100 bases of the the
    start or end of the sequence (or gap is the sequence is not
    contiguous), or within 100 bases of another contaminant match that is
    at an end, flag the contaminant span for trimming.
    If the contaminant is one of the above, and the match is longer than
    700 bases flag the contaminant span for trimming.
    Other matches may be false alarms. Treat them as suspect spans and
    reBLAST the hit span plus 10 Kbp of flanking sequence on each side
    against nr, HTGS, related and unrelated chromosomes (as described
    below).
    
    Return: set of regions to be trimmed
    """
    # list of sequence to modify
    tomodify = {}
    # go through each record in the fasta file
    for record in records:
        if record in contint:
            # check if sum of intervals is smaller than 75% of the
            # record length
            recordlength = len(records[record].seq)
            modifications = {}
            allcontint = contint[record]["allint"]
            coverage = 0
            for interval in allcontint:
                coverage = coverage + abs(interval[0]-interval[1])+1
            if coverage/recordlength > 0.75:
                modifications["remove"] = 1
            # Check if intervals are close to ends If the contaminant
            # is classed as VECTOR, E.coli, IS:./, PERM:./ or PHG:* and
            # the contaminant location is within 100 bases of the the
            # start or end of the sequence (or gap is the sequence is
            # not contiguous), or within 100 bases of another
            # contaminant match that is at an end, flag the contaminant
            # span for trimming.
            contintEnd = rU.fuseinterval(allcontint, 100)
            if contintEnd[0][0] < 100:
                modifications["trim5"] = contintEnd[0][1]
            if contintEnd[-1][1] > recordlength - 100:
                modifications["trim3"] = contintEnd[-1][0]
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
    parser.add_argument('--dbcont',
                        '-c',
                        help="The contaminant database",
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
    f = runContaminantPipe(args.genomefile, args.dbcont)
    contint = parseContaminant(f, dist)
    tomodify = qcContaminant(contint, records)
    pprint.pprint(tomodify)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
