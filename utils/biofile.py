from Bio import SeqIO, Phylo
from Bio.Seq import Seq

from Bio.Alphabet import generic_dna
from Bio.Alphabet.IUPAC import *

import re
import magic

def detect_type(filename):
    """
    :param filename: File to read and detect the format
    :return: detected type, in [fasta, phylip, phylip-relaxed, newick, N/A]

    Tests formats using biopython SeqIO or Phylo
    """
    mimetype=magic.from_file(filename,mime=True)

    if mimetype != "text/plain" :
        return mimetype
    
    # Check Fasta Format
    try:
        nbseq = 0
        for r in SeqIO.parse(filename, "fasta"):
            nbseq += 1
        if nbseq > 0:
            return "fasta"
    except Exception:
        pass

    # Check phylip strict
    try:
        nbseq = 0
        for r in SeqIO.parse(filename, "phylip"):
            nbseq += 1
        if nbseq > 0:
            return "phylip"
    except Exception:
        pass

    # Check phylip relaxed
    try:
        nbseq = 0
        for r in SeqIO.parse(filename, "phylip-relaxed"):
            nbseq += 1
        if nbseq > 0:
            return "phylip"
    except Exception:
        pass

    # Check Newick
    try:
        nbtrees = 0
        trees = Phylo.parse(filename, 'newick')
        for t in trees:
            nbtrees += 1
        if nbtrees > 0:
            return "nhx"
    except Exception as e:
        pass

    return "txt"


def nb_sequences(filename, format):
    nbseq = 0
    length = 0
    seqaa = False
    if format == 'fasta':
        try:
            for r in SeqIO.parse(filename, "fasta"):
                tlen = len(r.seq)
                length = tlen if tlen > length else length
                nbseq += 1
                if check_aa(r.seq):
                    seqaa = True
        except Exception as e:
            print e
            pass
    elif format == 'phylip':
        try:
            for r in SeqIO.parse(filename, "phylip"):
                tlen = len(r.seq)
                length = tlen if tlen > length else length
                nbseq += 1
                if check_aa(r.seq):
                    seqaa = True
        except Exception:
            try:
                for r in SeqIO.parse(filename, "phylip-relaxed"):
                    tlen = len(r.seq)
                    length = tlen if tlen > length else length
                    nbseq += 1
                    if check_aa(r.seq):
                        seqaa = True
            except Exception:
                pass
            pass
    return (nbseq, length, seqaa)


def check_aa(sequence):
    """
    Returns True if the sequence can be considered as proteic
    """
    alphabets = [extended_protein]
    for alphabet in alphabets:
        leftover = set(str(sequence).upper()) - set(alphabet.letters)
        if not leftover:
            return True
    return False



def check_nt(sequence):
    """
    Returns True if the sequence can be considered as nucleotidic
    """
    alphabets = [ambiguous_dna, unambiguous_dna, extended_dna, ambiguous_rna, unambiguous_rna]
    for alphabet in alphabets:
        leftover = set(str(sequence).upper()) - set(alphabet.letters)
        if not leftover:
            return True
    return False

def valid_fasta(fasta_file):
    # Check uploaded file or pasted content

    mimetype=magic.from_buffer(fasta_file.read(1024),mime=True)

    print mimetype
    if mimetype != "text/plain" :
        return (0,0,False)

    fasta_file.seek(0)
    
    nbseq = 0
    length = 0
    seqaa=False
    for r in SeqIO.parse(fasta_file, "fasta"):
        tlen = len(r.seq)
        length = tlen if tlen > length else length
        nbseq += 1
        if check_aa(r.seq):
            seqaa = True
    return (nbseq, length, seqaa)



def is_fasta_one_seq(filename):
    """
    :param filename: File to read and detect the format
    :return: true if format is fasta and contains only one sequence

    Tests formats using biopython SeqIO
    """
    # Check Fasta Format
    try:
        nbseq = 0
        for r in SeqIO.parse(filename, "fasta"):
            nbseq += 1
        if nbseq == 1:
            return True
    except Exception:
        pass
    return False


def newick_clean(seqname):
    """
    Clean the sequence name to be compatible with newick format
    Try to extract species name and gene name if possible
    """
    " Removing BL_ORD_ID if any"
    seqname = re.sub(r"\s*(?i)[^\s]*\|BL_ORD_ID\|\d+\s*", "", seqname)
    species = ""
    m = re.search(r"(\[(.+?)\])", seqname)
    if m is None:
        m = re.search(r"(PREDICTED: (\w+ \w+))",seqname)
    
    if m is None:
        m = re.search(r"^[^\s]+( (\w+ \w+))",seqname)
    
    if m is not None:
        toremove= m.group(1)
        species = "_"+m.group(2)
        seqname = seqname.replace(toremove,"")
    
    species=re.sub(r"\s\(.*\)","",species)
    
    gene=""
    m = re.search(r"sp\|[^\s]*\|([\w_]+)", seqname)
    if m is not None:
        gene = m.group(1)
    else:
        m = re.findall(r"\((\w+)\)", seqname)
        if len(m) > 0:
            gene = "_"+m[0]

    out = seqname.split(" ")[0]+gene+species
    out = out.replace("[","_")
    out = out.replace("]","_")
    out = out.replace("(","_")
    out = out.replace(")","_")
    out = out.replace(",","_")
    out = out.replace(";","_")
    out = out.replace(" ","_")
    out = out.replace(":","_")
    out = re.sub(r"_+","_",out)
    out = re.sub(r"_$","",out)
    
    return out

def cleanseqname(seqname):
    """
    Cleans sequence name to be compatible with newick format
    """
    out = seqname.split(" ")[0]
    out = out.replace("[","_")
    out = out.replace("]","_")
    out = out.replace("(","_")
    out = out.replace(")","_")
    out = out.replace(",","_")
    out = out.replace(";","_")
    out = out.replace(" ","_")
    out = out.replace(":","_")
    out = re.sub(r"_+","_",out)
    out = re.sub(r"_$","",out)
    return out


def translate(sequence, frame):
    """
    It takes a sequence and translate it in the right frame.
    if frame is
     1,  2,  3: Just removes 0, 1, or 2 nt and translates
    -1, -2, -3: RevComp, then removes 0, 1, or 2 nt, and translates
    """
    my_dna = Seq(sequence, generic_dna)
    if frame < 0:
        my_dna = my_dna.reverse_complement()
    my_dna = Seq(str(my_dna)[abs(frame)-1:])
    return str(my_dna.translate())
