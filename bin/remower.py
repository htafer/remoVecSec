"""This script help removing contamination from assembled genomes
"""

import remoVecSec.removeUtils as rU
import remoVecSec.removeVec as rV
import remoVecSec.removeContaminant as rC
import remoVecSec.removeMito as rM


import argparse
import sys
from Bio import SeqIO


def main(arguments):
    # Parser 
    parser = argparse.ArgumentParser(
        description="remower is a script that allows to remove\n"
        + " contamination in assembled genome. It takes as input a genome\n"
        + "file, contamination databases returns on stdout the corrected genome\n"
        + "and on stderr warnings regarding vector sequences not removed.",
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--genomefile', '-g', help="Genome file",
                        type=str, required=True)
    parser.add_argument('--dbvec',
                        '-v',
                        help="The vecscreen database",
                        type=str, required=False)
    parser.add_argument('--dbmito',
                        '-m',
                        help="The organelle database",
                        type=str, required=False)
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
    records = SeqIO.index(args.genomefile, "fasta")
    dist = args.dist
    # Start analysis
    # Contaminant
    f = rC.runContaminantPipe(args.genomefile, args.dbcont)
    contint = rC.parseContaminant(f, dist)
    qcCont = rC.qcContaminant(contint, records)
    # Vectors
    f = rV.runVecSecPipe(args.genomefile, args.dbvec)
    vectint = rV.parseVecSec(f, dist)
    qcVect = rV.qcVector(vectint, records)
    # Mito
    f = rM.runMitoPipe(args.genomefile, args.dbmito)
    mitoint = rM.parseMito(f, dist)
    qcMit = rM.qcMito(mitoint, records)
    # Contaminant
    finalmodifs = rU.mergemodification(qcCont, qcVect, qcMit)
    # Correct Fasta
    rU.correctfasta(finalmodifs, records)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
