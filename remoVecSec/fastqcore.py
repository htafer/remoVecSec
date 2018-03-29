"""This module is responsible for reading and processing fastq
files. Fastq pre-processing is not done here and the user should take
care of it

"""

import sys
import numpy as np


def samplefastq(filename, num_reads):
    """Read a fastq file and return a sampled fastq with the number of
    reads num_reads with the reservoir sampling algorithm. The user
    should ensure that the fastq file exists and is valid
    Args:
    num_reads(int): number of reads to sample Already checked if it
    exists
    filename(str): name of the fastq file Returns: a list of
    reads in fastq format that can be further processed
    """
    # Get the number of lines in the fastq files
    f = open(filename, 'r')
    number_of_line = sum(1 for line in f)
    f.close()
    # Check if file line number of %4
    if number_of_line % 4 != 0:
        print("number of lines in {}"
              + "is not a multiple of"
              + "4, check your FASTQ file",
              filename)
        sys.exit(0)
    # Sample the line first
    number_of_line = int(number_of_line/4)
    # Sample num_reads from the whole file
    reservoir = np.random.choice(number_of_line, num_reads, replace=False)
    # Sort reservoir
    np.ndarray.sort(reservoir)
    # Start looking for lines
    linecount = -1
    # Save sampled reads into file samplefile
    # samplefile = filename +".sample"
    samplefile = filename+".sample"
    s = open(filename+".sample", 'w')
    f = open(filename, 'r')
    for line in f:
        linecount = linecount+1
        if int(linecount/4) in reservoir:
            s.write(line)
            # samplefile.append(line)
    f.close()
    s.close()
    return samplefile




def printfastq(fastqlist):
    """Print a list containing fastqfile
    Args:
    fastqlist (list): set of reads inf fastq format stored in a list
    """
    # output on screen
    for entry in fastqlist:
        print(entry, end='')


if __name__ == '__main__':
    res = samplefastq(sys.argv[1], int(sys.argv[2]))
    printfastq(open(res,'r'))
