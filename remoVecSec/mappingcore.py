"""This module is responsible for generating the index of the genome
of interest and running the mapping process

"""
import os
import sys
import subprocess
import re
import shlex


def generateindex(genomefile, indexing_command):
    """Generate an index for a genome based on the indexing command

    Args: genomefile (str): list of genomes sequences to
    create the indexes for.
    indexing_command (str): command to index the genomes

    Returns: none
    """
    mappers_list = ["bwa", "hisat", "STAR"]
    args = []
    for mapper in mappers_list:
        prog = re.compile(mapper)
        if prog.match(indexing_command):
            print("Indexer for mapper ", mapper, " detected")
            args = shlex.split(indexing_command)
            if mapper == "hisat":
                args.append(genomefile)  # append for input
                args.append(genomefile)  # append for output
                # print(args)
                # p = subprocess.Popen(args, check=True)   # start process
                # p.communicate()          # use it or subprocess will hang
            if mapper == "STAR":
                args.append("--genomeFastaFiles")
                args.append(genomefile)
                # print(args)
                # p = subprocess.Popen(args, check=True)
                # p.communicate()
            else:
                args.append(genomefile)
                # p= subprocess.Popen(args, check=True)
                # p.communicate()
    if args != []:
        p = subprocess.Popen(args)
        p.communicate()
        return 1

    print(indexing_command,
          "Not implemented yet, please try with hisat2, bwa or STAR")
    sys.exit(0)


def mapreads(genomefile, fastqfile, mapping_command):
    """Map the reads in the fastqfile on the genomefile
    Args: 
    genomefile (str): genome sequence
    fastqfile: (str): file containing reads
    mapping_command (str): mapping command

    Returns: none
    """
    mappers_list = ["bwa", "hisat", "STAR"]
    mapped_name = fastqfile+".sam"
    args = []
    for mapper in mappers_list:
        prog = re.compile(mapper)
        if prog.match(mapping_command):
            print("Indexer for mapper ", mapper, " detected")
            args = shlex.split(mapping_command)
            if mapper == "hisat":
                args.append("-x "+genomefile)  # append genome option
                args.append("-U "+fastqfile)   # append fastq option
                args.append("-S "+fastqfile+".sam") # append output option
                # print(args)
                # p = subprocess.Popen(args, check=True)   # start process
                # p.communicate()          # use it or subprocess will hang
                break
            if mapper == "STAR":
                # FASTQ
                args.append("--readFilesIn")
                args.append(fastqfile)
                args.append("&& mv Aligned.out.sam "+fastqfile+".sam")
                # print(args)
                # p = subprocess.Popen(args, check=True)
                # p.communicate()
                break
            else:
                args.append(genomefile)
                args.append(fastqfile)
                args.append(">"+mapped_name)
                # p= subprocess.Popen(args, check=True)
                # p.communicate()
                break

    if args != []:
        print(args)
        os.system(" ".join(args))
        return 0

    print(mapping_command, "Not implemented yet, please try with hisat2, bwa or STAR")
    sys.exit(0)


if __name__ == '__main__':
    ### Test index
    # BWA
     generateindex("../test/test.fasta", "bwa index")
    # STAR
    # generateindex("../test/test.fasta",
    #              "STAR --runMode genomeGenerate --genomeDir ../test")
    # Test mapping (sys.argv[1], sys.argv[3], sys.argv[4])
    # BWA
    # mapreads("../test/test.fasta", "../test/test.fastq", "bwa mem")
    # STAR
    # mapreads("../test/test.fasta",
    #         "../test/test.fastq",
    #         "STAR --genomeDir ../test")
    # HISAT
