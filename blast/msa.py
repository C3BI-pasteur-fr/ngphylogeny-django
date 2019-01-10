from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import logging

class PseudoMSA:
    """
    This class represents a pseudo Multiple Sequence alignment
    That will be build from a BLAST result and having the following
    properties:
    - columns having GAPS in the query are removed
    - subjects having several hsps are merged in the
    same sequence, filled with gaps around
    """

    query_id = ""
    query_seq = ""
    blast_prog = ""
    sequences = dict()  # Dict: key: seqname, value: sequence
    scores = dict()  # Dict: key: seqname, value: blast score
    frame = 0

    def __init__(self, query_id, query_seq, query_seq_bk, frame, blast_prog):
        """
        query : blast query sequence without gaps
        """
        self.query_id = str(query_id)
        self.query_seq = list(query_seq)
        self.query_seq_bk = query_seq_bk
        self.frame = frame
        self.blast_prog = blast_prog

    def add_hsp(self, sbjct_name, hsp):
        """
        sbjct_name: name of the subject sequence
        hsp: Bio.Blast.Record.HSP
        """
        if ((self.blast_prog != "blastx" and self.blast_prog != "tblastx")  or hsp.frame[0] == self.frame):
            seq = self.sequences.get(str(sbjct_name))
            sco = self.scores.get(str(sbjct_name))
            if seq is None:
                seq = list("-" * len(self.query_seq))
            if sco is None:
                sco = 0
            start = hsp.query_start

            if self.blast_prog=='blastx' or self.blast_prog=='tblastx':
                # Computing start on the translated query sequence
                if self.frame > 0:
                    start = (hsp.query_start - self.frame) / 3 + 1
                else:
                    start = (len(self.query_seq_bk) - hsp.query_end + 1 + self.frame) / 3 + 1

            position = start-1

            for p in range(0, len(hsp.sbjct)):
                # If no gap in the query at that position
                if hsp.query[p] != '-':
                    if hsp.sbjct[p] != '-':
                        seq[position] = hsp.sbjct[p]
                    position += 1
            self.sequences[str(sbjct_name)] = seq
            self.scores[str(sbjct_name)] = max(sco, hsp.score)
            
    def all_sequences(self):
        for id, seq in self.sequences.iteritems():
            yield (str(id), "".join(seq))

    def first_n_max_score_sequences(self, maxseqs):
        nseqs = 0
        for key, value in sorted(self.scores.iteritems(), reverse=True, key=lambda (k,v): (v,k)):
            if nseqs >= maxseqs :
                break
            else:
                seq = self.sequences.get(str(key))
                yield (str(key), "".join(seq))
                nseqs = nseqs+1
            
    def to_string(self):
        msa = ">%s\n%s\n" % (self.query_id, "".join(self.query_seq))
        for id, seq in self.sequences.iteritems():
            msa += ">%s\n%s\n" % (id, "".join(seq))
        return msa
