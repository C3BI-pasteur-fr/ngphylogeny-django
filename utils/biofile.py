from Bio import SeqIO, Phylo
from Bio.Seq import Seq

from Bio.Alphabet import generic_dna
from Bio.Alphabet.IUPAC import *

import io
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
        except Exception:
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
    # Read the whole thing upfront and normalize to text: uploaded file
    # objects (Django's InMemoryUploadedFile/TemporaryUploadedFile) yield
    # bytes, while pasted content already comes in as a text StringIO.
    # SeqIO.parse() needs the former case wrapped as text - byte lines
    # break SimpleFastaParser's EOF check (it compares readline() to ""
    # to detect the end of the file, which never matches b"").
    raw = fasta_file.read()
    if isinstance(raw, bytes):
        content = raw.decode('utf-8', errors='replace')
        mime_sample = raw[:1024]
    else:
        content = raw
        mime_sample = raw[:1024].encode('utf-8', errors='ignore')

    mimetype = magic.from_buffer(mime_sample, mime=True)

    if mimetype != "text/plain" :
        return (0,0,False)

    nbseq = 0
    length = 0
    seqaa=False
    for r in SeqIO.parse(io.StringIO(content), "fasta"):
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

    split(None)/\s (not literal " "/.split(" ")): Python's default
    whitespace handling - unlike a literal " " - is Unicode-aware, so
    this also catches e.g. a non-breaking space (U+00A0), not just the
    plain ASCII one. See sanitize_fasta_content()'s own docstring below
    for the real production bug (a NBSP silently surviving this class
    of cleanup) this same fix addresses more broadly.
    """
    out = seqname.split(None, 1)[0] if seqname.split(None, 1) else seqname
    out = out.replace("[","_")
    out = out.replace("]","_")
    out = out.replace("(","_")
    out = out.replace(")","_")
    out = out.replace(",","_")
    out = out.replace(";","_")
    out = re.sub(r"\s","_",out)
    out = out.replace(":","_")
    out = re.sub(r"_+","_",out)
    out = re.sub(r"_$","",out)
    return out


def sanitize_fasta_id(seq_id):
    """
    Normalizes a single FASTA sequence id to a plain-ASCII token every
    downstream Galaxy tool in the pipeline (MAFFT, PhyML/PhyML-SMS,
    newick_utilities' nw_display, ...) will tokenize identically.

    Real production bug this exists to prevent: a real user's uploaded
    FASTA had a non-breaking space (U+00A0) embedded inside a sequence
    id (e.g. "A0A1Q2MHV5\\xa0_1_364" - plausibly a copy-paste artifact
    from a webpage/spreadsheet). Nothing in the regular OneClick/
    Advanced/Tool submission pipeline sanitized sequence ids at all
    (cleanseqname() above is only ever called from blast/tasks.py, a
    separate code path), so it survived MAFFT/PhyML/PhyML-SMS untouched
    - those tools don't treat a NBSP as a token delimiter - into the
    output tree, where newick_utilities' own Newick parser *does* treat
    it as one, producing "ERROR: missing ')' at line 0 near '_1_364'"
    on the "Tree image" step. A more lenient parser (gotree) read the
    same tree back with no complaint - it was never a "broken tree",
    just two tools in the same pipeline silently disagreeing about
    where one token ends and the next begins.

    Same replacement set as cleanseqname() (Newick's own special
    characters, plus any Unicode whitespace - not just the literal
    ASCII space), but keeping the whole id (cleanseqname() is only ever
    used to derive a short *display* name, and deliberately truncates
    at the first whitespace instead).
    """
    out = re.sub(r"\s+", "_", seq_id)
    for ch in "()[]{}:;,":
        out = out.replace(ch, "_")
    out = re.sub(r"_+", "_", out)
    out = out.strip("_")
    return out or "seq"


def sanitize_fasta_content(content):
    """
    Rewrites every '>' header line's sequence id (the token up to the
    first whitespace) via sanitize_fasta_id() - leaving the rest of
    each header (the description) and every non-header line untouched.
    See that function's own docstring for the real bug this prevents.

    FASTA only, not PHYLIP - the far more common upload path (OneClick
    itself is documented as taking "Fasta format" - see
    templates/workflows/workflows_form.html), and the one the real
    incident behind this function was actually in; phylip's own
    differently-shaped sequence-name field was out of scope here.

    Accepts/returns either bytes or str, matching whichever the caller
    already has (an uploaded file's own chunks are bytes; pasted text
    already comes in as str) - same bytes-or-str tolerance as
    valid_fasta() itself, for the same reason.
    """
    is_bytes = isinstance(content, bytes)
    text = content.decode("utf-8", errors="replace") if is_bytes else content

    out_lines = []
    for line in text.splitlines(keepends=True):
        if line.startswith(">"):
            ending = ""
            for newline in ("\r\n", "\n", "\r"):
                if line.endswith(newline):
                    ending = newline
                    line = line[:-len(newline)]
                    break
            header = line[1:]
            # Split id from description on a *literal* ASCII space/tab
            # only, not header.split(None, ...) - the conventional
            # ">id description" FASTA boundary is exactly that, and
            # Python's default Unicode-aware split(None) would treat an
            # embedded NBSP itself as the id/description separator,
            # completely misplacing the split for the one case this
            # function exists to fix (an id with a NBSP *inside* it,
            # not between it and a real description). Any remaining
            # Unicode whitespace within the extracted id itself is what
            # sanitize_fasta_id() below then cleans up.
            match = re.match(r"([^ \t]*)([ \t].*)?$", header, re.DOTALL)
            seq_id_raw, rest = match.group(1), match.group(2) or ""
            if seq_id_raw:
                line = ">" + sanitize_fasta_id(seq_id_raw) + rest
            line += ending
        out_lines.append(line)

    result = "".join(out_lines)
    return result.encode("utf-8") if is_bytes else result


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
