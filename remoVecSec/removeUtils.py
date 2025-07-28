def mergemodification(qcCont, qcVect, qcMit):

"""
Utility functions for vector and contaminant interval processing.
"""

from collections import defaultdict

def mergemodification(qcCont, qcVect, qcMit):
    """
    Merge overlapping modifications from contaminant, vector, and mitochondrial sources.

    Args:
        qcCont (dict): Contaminant modifications.
        qcVect (dict): Vector modifications.
        qcMit (dict): Mitochondrial modifications.

    Returns:
        dict: Merged modifications per record.
    """
    mergedmodifs = defaultdict(dict)
    for record in sorted(set([*qcCont, *qcVect, *qcMit])):
        mergedmodifs[record] = defaultdict(list)
        for src in (qcCont, qcVect):
            if record in src:
                for typemodif in ["trim3", "trim5", "remove"]:
                    if typemodif in src[record]:
                        mergedmodifs[record][typemodif].append(src[record][typemodif])
        if record in qcMit and "remove" in qcMit[record]:
            mergedmodifs[record]["remove"].append(qcMit[record]["remove"])

    # Reduce arrays to single values
    for record in mergedmodifs:
        if "trim5" in mergedmodifs[record]:
            mergedmodifs[record]["trim5"] = max(mergedmodifs[record]["trim5"])
        if "trim3" in mergedmodifs[record]:
            mergedmodifs[record]["trim3"] = min(mergedmodifs[record]["trim3"])
        if "remove" in mergedmodifs[record]:
            mergedmodifs[record]["remove"] = 1
    return mergedmodifs



def fuseinterval(intervals, dist):
    """
    Merge adjacent or overlapping intervals within a given distance.

    Args:
        intervals (list): List of [start, end] intervals.
        dist (int): Maximum distance to consider intervals as adjacent.

    Returns:
        list: Fused and reverse-sorted intervals.
    """
    if len(intervals) <= 1:
        return intervals
    intervals = sorted(intervals, key=lambda interval: interval[1])
    fused = []
    prev = intervals[0]
    for interval in intervals[1:]:
        if int(prev[1] + dist) >= int(interval[0]):
            prev[1] = interval[1]
        else:
            fused.append(prev)
            prev = interval
    fused.append(prev)
    return sorted(fused, key=lambda interval: interval[1], reverse=True)


def correctfasta(vectint, records):
    """
    Print corrected genome sequences after removing or trimming vector intervals.

    Args:
        vectint (dict): Intervals to remove/trim per record.
        records (dict): Genome records (e.g., from SeqIO.index).
    """
    for record in records:
        if record in vectint:
            recordseq = records[record]
            if "remove" in vectint[record]:
                continue
            if "trim3" in vectint[record]:
                recordseq = recordseq[:vectint[record]["trim3"]]
            if "trim5" in vectint[record]:
                recordseq = recordseq[vectint[record]["trim5"]:]
            if len(recordseq.seq) > 0:
                print(f">{record}")
                print(recordseq.seq)
        else:
            print(f">{record}")
            print(records[record].seq)
