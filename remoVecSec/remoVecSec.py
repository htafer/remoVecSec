"""This module is responsible for finding vector contamination in genomes
"""


# import pprint
# from Bio.Seq import Seq
import io
import sys
import subprocess
import re
from Bio import SeqIO
import logging
import argparse
import removeUtils as rU


def runVecSecPipe(genomefile, dbfile):
    """
    Run vecscreen with a vector database on a genome. Return a filehandle
    Args: genomfile(str) path of genomefile on which vecscreen is ran
          dbfile(str) path to the vector database used by vecscreen
    Returns a filehandle containing the output of vecscreen

    """
    # cat genomefile
    cat = subprocess.Popen(['cat', genomefile],
                           stdout=subprocess.PIPE)
    # pipe it to vecsec
    vecsec = subprocess.Popen(['vecscreen',
                               '-f 3',
                               '-d '+dbfile],
                              stdin=cat.stdout,
                              stdout=subprocess.PIPE)
    # return output of vecscreen
    return io.TextIOWrapper(vecsec.stdout, encoding='utf-8')


def parseVecSec(vecsec_output, dist):
    """Parse the output of vecscreen
    Args: vecsec_output output of vecsec

    Returns a dictionary with keys corresponding to the contigs
    containing vectors and value corresponding to the coordinates of
    the vector inside the contig
    """
    # hash containing as key the vector, as value a list of intervals
    vectint = {}
    # Go through file line by line
    for line in vecsec_output:
        # get line from file
        decoded = line
        # Regexp for contig
        contig = re.compile(">Vector")
        # Regexp for finding lines starting with numbers
        numberstart = re.compile("^[0-9]")
        # Regexp for no match
        nohits = re.compile("No hits found")
        # Regexp for empty line
        emptyline = re.compile("^\s*$")
        # if we have a contig print / save it
        if contig.match(decoded):
            # parse contigname, split based on regex
            contigname = re.compile("\s+").split(decoded)[1]
            # check if we have vector hit by looking at the next line
            currentline = vecsec_output.readline()
            # if we have a hit go to next line
            if nohits.match(currentline):
                logging.info("skipping ", decoded, "since no hit\n", end="")
                continue
            # if we have hits, go through them till the next emptyline
            else:
                # We have at least one interval
                # parse the contig name
                intervals=[]
                currentline = vecsec_output.readline()
                while currentline:
                    # empty line ends parsing
                    if emptyline.match(currentline):
                        break
                    # parse line
                    if numberstart.match(currentline):
                        # We got an interval
                        # get first two entries as intervals
                        tempintervals = re.compile("\s+").split(currentline)[:2]
                        # convert them to int
                        tempintervals = [int(x) for x in tempintervals]
                        # append
                        intervals.append(tempintervals)
                        # go to next line
                    currentline = vecsec_output.readline()
                # key, value of hash with intervals and contigname
                # fuseintervals
                vectint[contigname] = rU.fuseinterval(intervals, dist)
    # return hash to next function for processing
    return vectint


def qcVector(vectint, records):
    """
    Given a set of vector interval and a genome to correct remove the
    sequence of the genomes overlapping with the vector intervals, if
    the overlap occurs less than 100nts away from the contig ends
    Args: vectint list of intervals
          genome genome sequence file

    Return: corrected genome sequence. Check the warning to see if
    regions were not corrected
    """
    tomodify = {}
    for record in records:
        if record in vectint:
            modifications = {}
            # if contig in vectint
            # fuse interval
            fuseintervals = vectint[record]
            # begin or end +/- 10
            recordlength = len(records[record].seq)
            if fuseintervals[0][0] < 100:
                # correct sequence at the begining
                modifications["trim5"] = fuseintervals[0][1]
            elif fuseintervals[-1][1] > recordlength-100:
                modifications["trim3"] = fuseintervals[-1][0]
            if len(modifications):
                tomodify[record] = modifications
    return tomodify


def main(arguments):
    # Parser 
    parser = argparse.ArgumentParser(
        description="remoVecSec is a small wrapper around vecscreen that allows to remove\n"
        + "vector contamination in assembled genome. It takes as input a genome\n"
        + "file, a vecscreen database and returns on stdout the corrected genome\n"
        + "and on stderr warnings regarding vector sequences not removed.",
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--genomefile', '-g', help="Genome file",
                        type=str, required=True)
    parser.add_argument('--dbvec',
                        '-v',
                        help="The vecscreen database",
                        type=str, required=False)
    parser.add_argument('--dist',
                        '-d',
                        help="Maximal distance for merging two intervals",
                        type=int,
                        required=False,
                        default=50)
    args = parser.parse_args(arguments)
    records = SeqIO.index(args.genomefile, "fasta")
    dist = args.dist
    # Start analysis
    # Vectors
    f = runVecSecPipe(args.genomefile, args.dbvec)
    vectint = parseVecSec(f,dist)
    # Contaminant
    correctfasta(vectint, records)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))



