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


def _guess_is_protein(seq_chars):
    """
    Best-effort nucleotide-vs-protein guess for sanitize_fasta_content()/
    sanitize_phylip_content()'s own "?" -> N/X replacement below.
    check_nt() first, not check_aa() alone or first - verified directly
    (not assumed): every letter in each nucleotide alphabet
    (unambiguous_dna.letters == "GATC", etc.) is *also* individually a
    valid amino-acid letter, since extended_protein.letters includes
    plain A/C/G/T/U among its own ambiguity codes - a pure-ACGT
    nucleotide sequence would otherwise incorrectly satisfy check_aa()
    too (protein's alphabet is by far the more permissive/inclusive of
    the two), and this function would guess "protein" for genuinely
    nucleotide input.

    "?", common alignment-gap characters, and whitespace are all
    stripped before checking (upper-cased first, same as check_aa()/
    check_nt() themselves) - both functions compare the sequence's
    *entire* character set against an alphabet's own letters, and none
    of "?"/a gap/a space is ever a member of any of them, so leaving
    them in would make every sequence that contains one fail both
    checks regardless of its real composition. Whitespace matters here
    specifically because callers pass in a PHYLIP line's own captured
    "rest" group (name/data boundary included) or several such lines
    joined together - both routinely still carry the separating
    whitespace itself, and some real PHYLIP files additionally space
    sequence data into blocks (e.g. "ACGT ACGT ACGT").

    Returns True (protein), False (nucleotide), or None if nothing's
    left to check or the remaining letters don't cleanly match either
    alphabet - callers leave "?" untouched in that case rather than
    guessing wrong.
    """
    stripped = re.sub(r"[?\-.\s]", "", seq_chars.upper())
    if not stripped:
        return None
    if check_nt(stripped):
        return False
    if check_aa(stripped):
        return True
    return None


def _unknown_base_char(is_protein):
    """None (composition couldn't be determined) -> don't replace at all."""
    if is_protein is None:
        return None
    return "X" if is_protein else "N"


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

    Also reused verbatim by sanitize_phylip_content() below - a
    sequence name's own character-cleanup rules don't actually depend
    on which format it came from, only how that format delimits the
    name from what follows it.
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
    each header (the description) untouched. See that function's own
    docstring for the real bug this prevents.

    Also replaces any "?" character found in the actual sequence data
    (never inside a header line) with "N" or "X" - whichever this
    file's own composition calls for, guessed once via
    _guess_is_protein() from every non-header character in the whole
    file (not per-sequence - a single phylogenetics input is
    overwhelmingly one homogeneous marker/alignment, never a mix of
    nucleotide and protein sequences, so one guess for the whole file
    is both simpler and no less accurate than a per-sequence one).
    "?" is a common "unknown base" placeholder in some source tools/
    formats, but isn't itself a valid symbol in any nucleotide/protein
    alphabet the tools downstream of this expect - left as "?"
    (never replaced) when the file's composition can't be confidently
    determined, rather than guessing wrong and introducing a character
    that's just as invalid as "?" was.

    FASTA only, not PHYLIP - see sanitize_phylip_content() below for
    that, and sanitize_sequence_content() for the dispatcher every real
    call site actually uses.

    Accepts/returns either bytes or str, matching whichever the caller
    already has (an uploaded file's own chunks are bytes; pasted text
    already comes in as str) - same bytes-or-str tolerance as
    valid_fasta() itself, for the same reason.
    """
    is_bytes = isinstance(content, bytes)
    text = content.decode("utf-8", errors="replace") if is_bytes else content

    is_protein = _guess_is_protein("".join(
        line for line in text.splitlines() if not line.startswith(">")))
    replacement_char = _unknown_base_char(is_protein)

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
        elif replacement_char:
            line = line.replace("?", replacement_char)
        out_lines.append(line)

    result = "".join(out_lines)
    return result.encode("utf-8") if is_bytes else result


def sanitize_phylip_content(content):
    """
    Loose/relaxed PHYLIP counterpart to sanitize_fasta_content() above -
    same real bug class (a Unicode-whitespace character, e.g. a non-
    breaking space, embedded inside a sequence name silently surviving
    into a Newick tree where it's the one tool in the pipeline that
    treats it as a token delimiter), just for PHYLIP-formatted input,
    which was never covered at all before this (every upload/paste call
    site sanitizes unconditionally regardless of detected format -
    sanitize_fasta_content() itself is a no-op on PHYLIP content, since
    none of its lines start with ">").

    Deliberately "loose": treats a sequence name as an arbitrary-length
    token up to the first run of whitespace (the same boundary rule as
    a FASTA id, and how sanitize_fasta_id() below is reused verbatim),
    not strict PHYLIP's fixed 10-character name width - real uploaded
    files vary too much in practice for that width to be worth
    enforcing here, and getting it wrong would risk corrupting sequence
    data rather than just a name.

    Only the first N lines *after* the header count as name lines - N
    (the sequence count) is the first whitespace-separated token on the
    header line (PHYLIP's own "N M" - sequence count/alignment length -
    first line, every variant's own convention, an optional trailing
    "I"/"S" interleaved/sequential flag ignored here). This is what
    correctly handles interleaved PHYLIP too, not just sequential: an
    interleaved file's later blocks are pure sequence-data continuation
    lines with *no* name field at all, and must be left completely
    untouched by sanitize_fasta_id()'s own character replacements -
    those only ever apply to the name portion of the first N lines. A
    header that can't be parsed (no leading integer) leaves the content
    untouched entirely, same as sanitize_sequence_content()'s own
    fallback below - this function only ever runs on content already
    sniffed as PHYLIP-shaped.

    Also replaces any "?" character found in the sequence-data portion
    of *every* data line (both the first N name lines' own trailing
    data, and any later interleaved continuation lines) with "N" or
    "X" - same composition guess (once for the whole file, "?"/gaps
    excluded from the check) as sanitize_fasta_content()'s own
    identical feature; see that function's own docstring for the full
    reasoning.

    Same bytes-or-str tolerance as sanitize_fasta_content().
    """
    is_bytes = isinstance(content, bytes)
    text = content.decode("utf-8", errors="replace") if is_bytes else content

    def split_ending(line):
        for newline in ("\r\n", "\n", "\r"):
            if line.endswith(newline):
                return line[:-len(newline)], newline
        return line, ""

    # Literal ASCII space/tab only for the name/sequence-data boundary,
    # not a Unicode-whitespace-aware \S/\s split - same reasoning as
    # sanitize_fasta_content()'s own id/description split: a name
    # containing an embedded NBSP (the exact case this whole function
    # exists for) would otherwise have the NBSP itself mistaken for the
    # boundary, splitting the name short *before* the character that
    # actually needs sanitizing and leaving it untouched in what this
    # treats as "rest".
    name_data_re = re.compile(r"([^ \t]*)([ \t].*)?$", re.DOTALL)

    lines = text.splitlines(keepends=True)

    header_index = None
    nseq = 0
    for i, line in enumerate(lines):
        body, _ = split_ending(line)
        if body.strip():
            header_index = i
            parts = body.split()
            try:
                nseq = int(parts[0])
            except (IndexError, ValueError):
                nseq = 0
            break

    if header_index is None or nseq <= 0:
        return content

    # First pass: gather every data line's own sequence-data portion
    # (the name-stripped "rest" for one of the first N lines, the whole
    # line for a later continuation one) to guess the file's overall
    # nucleotide/protein composition once - the same guess is then
    # applied consistently to every line in the second pass below,
    # rather than each line (mis)guessing on its own, much shorter,
    # slice of the data.
    data_chars = []
    named_lines_seen = 0
    for i in range(header_index + 1, len(lines)):
        body, _ = split_ending(lines[i])
        if not body.strip():
            continue
        if named_lines_seen < nseq:
            data_chars.append(name_data_re.match(body).group(2) or "")
            named_lines_seen += 1
        else:
            data_chars.append(body)
    replacement_char = _unknown_base_char(_guess_is_protein("".join(data_chars)))

    out_lines = list(lines)
    named_lines_seen = 0
    for i in range(header_index + 1, len(lines)):
        body, ending = split_ending(lines[i])
        if not body.strip():
            # A blank separator line (common between interleaved
            # blocks) - doesn't carry a name or data, isn't counted,
            # isn't touched.
            continue
        if named_lines_seen < nseq:
            match = name_data_re.match(body)
            name_raw, rest = match.group(1), match.group(2) or ""
            if replacement_char:
                rest = rest.replace("?", replacement_char)
            if name_raw:
                out_lines[i] = sanitize_fasta_id(name_raw) + rest + ending
            named_lines_seen += 1
        elif replacement_char:
            out_lines[i] = body.replace("?", replacement_char) + ending

    result = "".join(out_lines)
    return result.encode("utf-8") if is_bytes else result


def sanitize_sequence_content(content):
    """
    Dispatches raw uploaded/pasted content to sanitize_fasta_content()
    or sanitize_phylip_content() based on a quick peek at its own first
    non-blank line - every call site only ever has raw content in hand
    at this point, not a file on disk to run detect_type() against
    (that function needs a real path for magic/SeqIO.parse). ">" starts
    a FASTA header; two leading whitespace-separated integers starts a
    PHYLIP header (its own "N M" line). Anything else (already-invalid
    input, or a format neither sanitizer applies to) is returned
    untouched - the actual format/content validation that matters
    happens downstream (biofile.valid_fasta()/nb_sequences()), not
    here; this is purely a best-effort cleanup before that.
    """
    is_bytes = isinstance(content, bytes)
    text = content.decode("utf-8", errors="replace") if is_bytes else content

    first_line = ""
    for line in text.splitlines():
        if line.strip():
            first_line = line.strip()
            break

    if first_line.startswith(">"):
        return sanitize_fasta_content(content)

    parts = first_line.split()
    if len(parts) >= 2 and parts[0].isdigit() and parts[1].isdigit():
        return sanitize_phylip_content(content)

    return content


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
