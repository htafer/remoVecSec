"""This module is responsible for generating the index of the genome
of interest and running the mapping process

"""
import io
import sys
import subprocess
import re
from Bio import SeqIO
import logging
import argparse


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


def parseVecSec(vecsec_output):
    """Parse the output of vecscreen
    Args: vecsec_output output of vecsec

    Returns a dictionary with keys corresponding to the contigs
    containing vectors and value corresponding to the coordinates of
    the vector inside the contig
    """
    # hash containing as key the vector, as value a list of intervals
    vectint={}
    for line in vecsec_output:
        # get line from file
        decoded = line
        # Regexp for contig
        contig = re.compile(">Vector")
        # Regexp for finding lines with starting with numbers
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
                logging.info("skipping ", decoded, "since no hit\n",end="")
                continue
            # if we have hits, go through them till the next emptyline
            else:
                # We have at least one interval
                # parse the contig name
                intervals=[]
                currentline = vecsec_output.readline()
                while currentline:
                    # empty line exnds parsing
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
                vectint[contigname] = intervals
    # return hash to next function for processing
    return vectint


def fuseinterval(intervals):
    """
    Check if two intervals are adjacent
    if they are fuse it and return only one
    interval
    Args: intervals a list of intervals
    Return a list of intervals
    """
    # if we have only one interval return it
    if len(intervals) == 1:
        return intervals
    # sort based on the coordinate of the end of the interval
    intervals = sorted(intervals,
                       key=lambda interval: interval[1])
    # go through interval, fuse the one that are adjacent
    fuseintervals=[]
    previnterval = intervals[0]
    # check adjacent intervals and merge them, return a list of intervals
    for interval in intervals[1:]:
        # adjacency: end of previous + 1 = start of current interval
        if previnterval[1]+1 == interval[0]:
            previnterval[1] = interval[1]
        else:
            fuseintervals.append(previnterval)
            previnterval = interval
    # last element has to be appended too
    fuseintervals.append(previnterval)
    # return fused intervals
    return fuseintervals


def correctfasta(vectint, genome):
    """
    Given a set of vector interval and a genome correct remove the
    sequence of the genomes overlapping with the vector intervals, if
    the overlap occurs less than 100nts away from the contig ends
    Args: vectint list of intervals
          genome genome sequence file

    Return: corrected genome sequence. Check the warning to see if
    regions were not corrected
    """
    # read fastafile
    records = SeqIO.parse(genome, "fasta")
    # open file genome"+.corrected"
    # go through each sequence in genome file
    for record in records:
        if record.id in vectint:
            # if contig in vectint
            # fuse interval
            fuseintervals = fuseinterval(vectint[record.id])
            # begin or end +/- 10
            recordlength = len(record.seq)
            for interval in fuseintervals:
                if interval[1] < 100:
                    # correct sequence at the begining
                    record.seq = record.seq[interval[1]:]
                    SeqIO.write(record,sys.stdout,"fasta")
                elif interval[0] > recordlength-100:
                    # correct sequence if at the end
                    record.seq = record.seq[:interval[0]]
                    SeqIO.write(record,sys.stdout,"fasta")
                else:
                    # in the other case, return a warning and let the
                    # user decide if the genome should be corrected or
                    # not
                    intervalstring = "{},{}".format(interval[0], interval[1])
                    logging.warning("interval "
                                    + intervalstring
                                    + " on contig "
                                    + record.id
                                    + " is too far from the contig's end and will not be processed")
        else:
            SeqIO.write(record, sys.stdout, "fasta")


def main(arguments):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--genomefile',
                        '-g',
                        help="Genome file",
                        type=str)
    parser.add_argument('--dbvec',
                        '-d',
                        help="The vecscreen database",
                        type=str)
    args = parser.parse_args(arguments)

    # f = runVecSecPipe(sys.argv[1], sys.argv[2])
    # vectint = parseVecSec(f)
    # correctfasta(vectint, sys.argv[1])
    f = runVecSecPipe(args.genomefile, args.dbvec)
    vectint = parseVecSec(f)
    correctfasta(vectint, args.genomefile)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))