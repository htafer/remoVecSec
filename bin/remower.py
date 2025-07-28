
"""
remower.py: Remove contamination from assembled genomes using BLAST and VecScreen.

This script screens a genome FASTA file for common contaminants, vectors, and mitochondrial sequences,
then outputs a cleaned genome file to stdout.
"""

import remoVecSec.removeUtils as rU
import remoVecSec.removeVec as rV
import remoVecSec.removeContaminant as rC
import remoVecSec.removeMito as rM


import argparse

    # Contaminant
    parser = argparse.ArgumentParser(
        description="Remove contamination from assembled genomes using BLAST and VecScreen.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--genomefile', '-g', help="Input genome FASTA file", type=str, required=True)
    parser.add_argument('--dbvec', '-v', help="VecScreen database (adaptor/vector)", type=str)
    parser.add_argument('--dbmito', '-m', help="Organelle (mitochondrial) BLAST database", type=str)
    parser.add_argument('--dbcont', '-c', help="Contaminant BLAST database", type=str)
    parser.add_argument('--dist', '-d', help="Max distance for merging intervals", type=int, default=50)
    args = parser.parse_args(arguments)

    # Validate input
    import os
    if not os.path.isfile(args.genomefile):
        sys.stderr.write(f"Error: Genome file '{args.genomefile}' not found.\n")
        return 1
    if not (args.dbvec or args.dbmito or args.dbcont):
        sys.stderr.write("Error: At least one database (--dbvec, --dbmito, --dbcont) must be provided.\n")
        return 1

    records = SeqIO.index(args.genomefile, "fasta")
    dist = args.dist

    # Run contaminant screen
    qcCont = {}
    if args.dbcont:
        f = rC.runContaminantPipe(args.genomefile, args.dbcont)
        contint = rC.parseContaminant(f, dist)
        qcCont = rC.qcContaminant(contint, records)

    # Run vector/adaptor screen
    qcVect = {}
    if args.dbvec:
        f = rV.runVecSecPipe(args.genomefile, args.dbvec)
        vectint = rV.parseVecSec(f, dist)
        qcVect = rV.qcVector(vectint, records)

    # Run mitochondrial screen
    qcMit = {}
    if args.dbmito:
        f = rM.runMitoPipe(args.genomefile, args.dbmito)
        mitoint = rM.parseMito(f, dist)
        qcMit = rM.qcMito(mitoint, records)

    # Merge and correct
    finalmodifs = rU.mergemodification(qcCont, qcVect, qcMit)
    rU.correctfasta(finalmodifs, records)

    return 0
    finalmodifs = rU.mergemodification(qcCont, qcVect, qcMit)
    # Correct Fasta
    rU.correctfasta(finalmodifs, records)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
