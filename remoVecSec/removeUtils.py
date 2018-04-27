import pprint
import collections as co
"""This module contains utility functions
"""


def mergemodification(qcCont, qcVect, qcMit):
    """
    Some modification are overlapping.
    Merge them.
    Use defaultdicts https://docs.python.org/3/library/collections.html#defaultdict-examples
    """
    # 
    mergedmodifs = co.defaultdict(dict)
    for record in sorted(set([*qcCont,
                              *qcVect,
                              *qcMit])):
        mergedmodifs[record] = co.defaultdict(dict)
        if record in qcCont:
            for typemodif in ["trim3", "trim5", "remove"]:
                if typemodif in qcCont[record]:
                    if typemodif in mergedmodifs[record]:
                        mergedmodifs[record][typemodif].append(qcCont[record][typemodif])
                    else:
                        mergedmodifs[record][typemodif] = []
                        mergedmodifs[record][typemodif].append(qcCont[record][typemodif])

        if record in qcVect:
            for typemodif in ["trim3", "trim5", "remove"]:
                if typemodif in qcVect[record]:
                    if typemodif in mergedmodifs[record]:
                        mergedmodifs[record][typemodif].append(qcVect[record][typemodif])
                    else:
                        mergedmodifs[record][typemodif] = []
                        mergedmodifs[record][typemodif].append(qcVect[record][typemodif])

        if record in qcMit:
            if "remove" not in mergedmodifs[record]:
                mergedmodifs[record]["remove"] = []
                mergedmodifs[record]["remove"].append(qcMit[record]["remove"])

    # Reduce array of number to a single number
    # For Trim3 take min
    # For trim5 take max
    # For max dont do anything
    for record in mergedmodifs:
        if "trim5" in mergedmodifs[record]:
            mergedmodifs[record]["trim5"] = max(mergedmodifs[record]["trim5"])
        if "trim3" in mergedmodifs[record]:
            mergedmodifs[record]["trim3"] = min(mergedmodifs[record]["trim3"])
        if "remove" in mergedmodifs[record]:
            mergedmodifs[record]["remove"] = 1
    # return mergedmodifs
    return mergedmodifs


def fuseinterval(intervals, dist):
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
    fuseintervals = []
    previnterval = intervals[0]
    # check adjacent intervals and merge them, return a list of intervals
    for interval in intervals[1:]:
        # adjacency: end of previous + 1 = start of current interval
        if int(previnterval[1]+50) >= int(interval[0]):
            previnterval[1] = interval[1]
        else:
            fuseintervals.append(previnterval)
            previnterval = interval
    # last element has to be appended too
    fuseintervals.append(previnterval)
    # reverse sort the interval, so that by
    # removing the vector sequence, the coordinates are
    # still actual.
    fuseintervals = sorted(fuseintervals,
                           key=lambda interval: interval[1], reverse=True)

    # check how many intervals
    return fuseintervals

def correctfasta(vectint, records):
    """
    Given a set of vector interval and a genome to correct remove the
    sequence of the genomes overlapping with the vector intervals, if
    the overlap occurs less than 100nts away from the contig ends
    Args: vectint list of intervals
          genome genome sequence file

    Return: corrected genome sequence. Check the warning to see if
    regions were not corrected
    """


# go through each sequence in genome file
    for record in records:
        if record in vectint:
            # We have the remove keyword. Do not process sequence record
            recordseq = records[record]
            if "remove" in vectint[record]:
                continue
            if "trim3" in vectint[record]:
                # We cannot work directly on the records hash
                # duplicate the sequence, and modify it
                recordseq = recordseq[:vectint[record]["trim3"]]
            if "trim5" in vectint[record]:
                # We cannot work directly on the records hash
                # duplicate the sequence, and modify it
                recordseq = recordseq[vectint[record]["trim5"]:]
            # print modified sequence
            print(">"+record)
            print(recordseq.seq)
        else:
            # print unmodified sequence
            print(">"+record)
            print(records[record].seq)
