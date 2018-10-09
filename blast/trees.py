from math import log, pow


def prot_dist(seq1, seq2):
    """
    Compute distance between two aligned sequences (must have the same length)
    Using Kimura's distance
    (see http://evolution.genetics.washington.edu/phylip/doc/protdist.html):
    D = - loge ( 1 - p - 0.2 p2 ).
    """
    if len(seq1) != len(seq2):
        raise ValueError(
            "Sequence lengths differ in protein distance computation %d != %d" % (
                len(seq1), len(seq2))
        )
    (diff, total) = count_diff(seq1, seq2)
    diff = float(diff) / float(total)
    dist = - log(1.0 - diff - 0.2 * pow(diff, 2))
    if dist > 0:
        return dist
    else:
        return 0.0


def nucl_dist(seq1, seq2):
    """
    Compute distance between two aligned nucleotidic sequences
    (must have the same length)
    Using Jukes Cantor method
    (see http://evolution.genetics.washington.edu/phylip/doc/dnadist.html)
    """
    if len(seq1) != len(seq2):
        raise ValueError(
            "Sequence lengths differ in protein distance computation"
        )
    (diff, total) = count_diff(seq1, seq2)
    diff = float(diff) / float(total)
    dist = -3.0 / 4.0 * log(1.0-4.0/3.0*diff)
    if dist > 0:
        return dist
    else:
        return 0.0


def count_diff(seq1, seq2):
    """
    Count number of differences between two sequences
    """
    nbdiffs = 0
    total = 0
    for i in range(0, len(seq1)):
        if seq1[i] != '-' and seq2[i] != '-':
            if seq1[i] != seq2[i]:
                nbdiffs += 1
            total += 1
    return (nbdiffs, total)
