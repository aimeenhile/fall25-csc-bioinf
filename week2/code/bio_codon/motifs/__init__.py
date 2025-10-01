# From Bio/motifs/__init__.py

"""Tools for sequence motif analysis.

Bio.motifs contains the core Motif class containing various I/O methods
as well as methods for motif comparisons and motif searching in sequences.
It also includes functionality for parsing output from the AlignACE, MEME,
and MAST programs, as well as files in the TRANSFAC format.
"""
from typing import Any, Optional, List, Dict, Tuple

from python import warnings
from python.urllib.parse import urlencode
from python.urllib.request import Request, urlopen
import numpy as np
from python.Bio.Align import Alignment
from python.Bio.Seq import Seq

from . import matrix
from . import minimal


def create(instances: List[str], alphabet: str ="ACGT") -> 'Motif':
    """Create a Motif object."""
    alignment = Alignment(instances)
    return Motif(alignment=alignment, alphabet=alphabet)


def parse(handle: Any, fmt: str, strict: bool = True) -> List['Motif']:
    """Parse an output file from a motif finding program.

    Currently supported format:
     - MINIMAL:          MINIMAL MEME output file motif

    If strict is True (default), the parser will raise a ValueError if the
    file contents does not strictly comply with the specified file format.
    """
    fmt = fmt.lower()
    if fmt == "minimal":
        from . import minimal
        return minimal.read(handle)
    else:
        raise ValueError(f"Unknown format '{fmt}'. Only 'minimal' is supported.")


def read(handle: Any, fmt: str, strict: bool = True):
    """Read a motif from a handle using the specified file-format.

    This supports the same formats as Bio.motifs.parse(), but
    only for files containing exactly one motif.  

    If the handle contains no records, or more than one record,
    an exception is raised.

    If however you want the first motif from a file containing
    multiple motifs this function would raise an exception (as
    shown in the example above). 

    Use the Bio.motifs.parse(handle, fmt) function if you want
    to read multiple records from the handle.

    If strict is True (default), the parser will raise a ValueError if the
    file contents does not strictly comply with the specified file format.
    """
    fmt = fmt.lower()
    motifs = parse(handle, fmt, strict)
    if len(motifs) == 0:
        raise ValueError("No motifs found in handle")
    if len(motifs) > 1:
        raise ValueError("More than one motif found in handle")
    motif = motifs[0]
    return motif


class Motif:
    """A class representing sequence motifs."""
    name: str
    instances: List[str]
    alphabet: str
    length: Optional[int]
    alignment: Optional[Alignment]
    counts: Any
    mask: Tuple[int, ...]
    _pseudocounts: Dict[str, float]
    _background: Dict[str, float]

    def __init__(self, alphabet: str ="ACGT", alignment: Optional[Alignment] = None, counts: Optional[Any] = None):
        """Initialize the class."""
        self.name = ""
        self.instances = []
        self._pseudocounts = {}
        self._background = {}
        self.mask = ()

        if counts is not None and alignment is not None:
            raise ValueError("Specify either counts or an alignment, don't specify both")
        elif counts is not None:
            self.alignment = None
            self.counts = matrix.FrequencyPositionMatrix(alphabet, counts)
            self.length = self.counts.length
        elif alignment is not None:
            self.length = alignment.length
            frequencies = alignment.frequencies
            for letter in alphabet:
                if letter not in frequencies:
                    frequencies[letter] = np.zeros(self.length, int)
            self.counts = matrix.FrequencyPositionMatrix(alphabet, frequencies)
            self.alignment = alignment
        else:
            self.counts = None
            self.alignment = None
            self.length = None

        self.alphabet = alphabet
        self.pseudocounts = None
        self.background = None
        self.mask = None

    def __get_mask(self):
        return self.mask

    def __set_mask(self, mask):
        if self.length is None:
            self.mask = ()
        elif mask is None:
            self.mask = tuple([1] * self.length)
        elif len(mask) != self.length:
            raise ValueError(
                "The length (%d) of the mask is inconsistent with the length (%d) of the motif"
                % (len(mask), self.length),
            )
        elif isinstance(mask, str):
            temp_mask: List[int] = []
            for char in mask:
                if char == "*":
                    temp_mask.append(1)
                elif char == " ":
                    temp_mask.append(0)
                else:
                    raise ValueError(
                        "Mask should contain only '*' or ' ' and not a '%s'" % char
                    )
            self.mask = tuple(temp_mask)
        else:
            self.mask = tuple(int(bool(c)) for c in mask)

    mask = property(__get_mask, __set_mask)
    del __get_mask
    del __set_mask

    def __get_pseudocounts(self):
        return self._pseudocounts

    def __set_pseudocounts(self, value):
        self._pseudocounts = {}
        if isinstance(value, dict):
            self._pseudocounts = {letter: float(value[letter]) for letter in self.alphabet}
        else:
            if value is None:
                value = 0.0
            self._pseudocounts = dict.fromkeys(self.alphabet, float(value))

    pseudocounts = property(__get_pseudocounts, __set_pseudocounts)
    del __get_pseudocounts
    del __set_pseudocounts

    def __get_background(self):
        return self._background

    def __set_background(self, value):
        if isinstance(value, dict):
            self._background = {letter: float(value[letter]) for letter in self.alphabet}
        elif value is None:
            self._background = dict.fromkeys(self.alphabet, 1.0)
        else:
            if not matrix._has_dna_alphabet(self.alphabet) and not matrix._has_rna_alphabet(self.alphabet):
                raise ValueError(
                    "Setting the background to a single value only works for DNA and RNA"
                    "motifs (in which case the value is interpreted as the GC content)"
                )
            T_or_U = "T" if matrix._has_dna_alphabet(self.alphabet) else "U"
            self._background["A"] = (1.0 - value) / 2.0
            self._background["C"] = value / 2.0
            self._background["G"] = value / 2.0
            self._background[T_or_U] = (1.0 - value) / 2.0
        total = sum(self._background.values())
        for letter in self.alphabet:
            self._background[letter] /= total

    background = property(__get_background, __set_background)
    del __get_background
    del __set_background

    def __getitem__(self, key: slice):
        """Return a new Motif object for the positions included in key."""

        if not isinstance(key, slice):
            raise TypeError("motif indices must be slices")
        alphabet = self.alphabet
        if self.alignment is None:
            alignment = None
            if self.counts is None:
                counts = None
            else:
                counts = {letter: self.counts[letter][key] for letter in alphabet}
        else:
            alignment = self.alignment[:, key]
            counts = None
        motif = Motif(alphabet=alphabet, alignment=alignment, counts=counts)
        motif.mask = self.mask[key]
        if alignment is None and counts is None:
            try:
                length = self.length
            except AttributeError:
                pass
            else:
                motif.length = len(range(*key.indices(length)))
        motif.pseudocounts = self.pseudocounts.copy()
        motif.background = self.background.copy()
        return motif

    @property
    def pwm(self):
        """Calculate and return the position weight matrix for this motif."""
        return self.counts.normalize(self._pseudocounts)

    @property
    def pssm(self):
        """Calculate and return the position specific scoring matrix for this motif."""
        return self.pwm.log_odds(self._background)

    def __str__(self, masked: bool =False) -> str:
        """Return string representation of a motif."""
        text: str = ""
        if self.alignment is not None:
            text += "\n".join(self.alignment)

        if masked:
            for i in range(self.length or 0):
                if self.__mask[i]:
                    text += "*"
                else:
                    text += " "
            text += "\n"
        return text

    def __len__(self) --> int:
        """Return the length of a motif.

        Please use this method (i.e. invoke len(m)) instead of referring to m.length directly.
        """
        if self.length is None:
            return 0
        else:
            return self.length

    def reverse_complement(self) -> 'Motif':
        """Return the reverse complement of the motif as a new motif."""
        alphabet = self.alphabet
        if not matrix._has_dna_alphabet(self.alphabet) and not matrix._has_rna_alphabet(self.alphabet):
            raise ValueError("Calculating reverse complement only works for DNA and RNA motifs")
        T_or_U = "T" if matrix._has_dna_alphabet(self.alphabet) else "U"

        if self.alignment is not None:
            alignment = self.alignment.reverse_complement()
            if T_or_U == "U":
                alignment.sequences = [s.replace("T", "U") for s in alignment.sequences]
            res = Motif(alphabet=alphabet, alignment=alignment)
        else:  # has counts
            new_counts = self.counts.reverse_complement()
            res = Motif(alphabet=alphabet, counts=new_counts)
            
        res.__mask = self.__mask[::-1]
        res.background = {
            "A": self.background[T_or_U],
            "C": self.background["G"],
            "G": self.background["C"],
            T_or_U: self.background["A"],
        }
        res.pseudocounts = {
            "A": self.pseudocounts[T_or_U],
            "C": self.pseudocounts["G"],
            "G": self.pseudocounts["C"],
            T_or_U: self.pseudocounts["A"],
        }
        return res

    @property
    def consensus(self) -> Seq:
        """Return the consensus sequence."""
        return self.counts.consensus

    @property
    def anticonsensus(self) -> Seq:
        """Return the least probable pattern to be generated from this motif."""
        return self.counts.anticonsensus

    @property
    def degenerate_consensus(self) -> Seq:
        """Return the degenerate consensus sequence.

        Following the rules adapted from
        D. R. Cavener: "Comparison of the consensus sequence flanking
        translational start sites in Drosophila and vertebrates."
        Nucleic Acids Research 15(4): 1353-1361. (1987).

        The same rules are used by TRANSFAC.
        """
        return self.counts.degenerate_consensus

    @property
    def relative_entropy(self) -> np.ndarray:
        """Return an array with the relative entropy for each column of the motif."""
        background = self.background
        pseudocounts = self.pseudocounts
        alphabet = self.alphabet
        counts = self.counts
        length = self.length
        values = np.zeros(length)
        if self.alignment is None:
            total = np.array(
                [
                    sum(counts[c][i] + pseudocounts[c] for c in alphabet)
                    for i in range(length)
                ]
            )
            for letter, frequencies in counts.items():
                frequencies = np.array(frequencies) + pseudocounts[letter]
                mask = frequencies > 0
                frequencies = frequencies[mask] / total[mask]
                values[mask] += frequencies * np.log2(frequencies / background[letter])
        else:
            total = np.zeros(length)
            for letter, frequencies in counts.items():
                total += np.array(frequencies) + pseudocounts[letter]
            for letter, frequencies in counts.items():
                frequencies = np.array(frequencies) + pseudocounts[letter]
                mask = frequencies > 0
                frequencies = frequencies[mask] / total[mask]
                values[mask] += frequencies * np.log2(frequencies / background[letter])
        return values

    def weblogo(self, fname, fmt="PNG", **kwds):
        """Download and save a weblogo using the Berkeley weblogo service.

        Requires an internet connection.

        The parameters from ``**kwds`` are passed directly to the weblogo server.

        Currently, this method uses WebLogo version 3.3.
        These are the arguments and their default values passed to
        WebLogo 3.3; see their website at http://weblogo.threeplusone.com
        for more information::

            'stack_width' : 'medium',
            'stacks_per_line' : '40',
            'alphabet' : 'alphabet_dna',
            'ignore_lower_case' : True,
            'unit_name' : "bits",
            'first_index' : '1',
            'logo_start' : '1',
            'logo_end': str(self.length),
            'composition' : "comp_auto",
            'percentCG' : '',
            'scale_width' : True,
            'show_errorbars' : True,
            'logo_title' : '',
            'logo_label' : '',
            'show_xaxis': True,
            'xaxis_label': '',
            'show_yaxis': True,
            'yaxis_label': '',
            'yaxis_scale': 'auto',
            'yaxis_tic_interval' : '1.0',
            'show_ends' : True,
            'show_fineprint' : True,
            'color_scheme': 'color_auto',
            'symbols0': '',
            'symbols1': '',
            'symbols2': '',
            'symbols3': '',
            'symbols4': '',
            'color0': '',
            'color1': '',
            'color2': '',
            'color3': '',
            'color4': '',

        """
        if set(self.alphabet) == set("ACDEFGHIKLMNPQRSTVWY"):
            alpha = "alphabet_protein"
        elif set(self.alphabet) == set("ACGU"):
            alpha = "alphabet_rna"
        elif set(self.alphabet) == set("ACGT"):
            alpha = "alphabet_dna"
        else:
            alpha = "auto"

        frequencies = format(self, "transfac")
        url = "https://weblogo.threeplusone.com/create.cgi"
        values = {
            "sequences": frequencies,
            "format": fmt.lower(),
            "stack_width": "medium",
            "stacks_per_line": "40",
            "alphabet": alpha,
            "ignore_lower_case": True,
            "unit_name": "bits",
            "first_index": "1",
            "logo_start": "1",
            "logo_end": str(self.length),
            "composition": "comp_auto",
            "percentCG": "",
            "scale_width": True,
            "show_errorbars": True,
            "logo_title": "",
            "logo_label": "",
            "show_xaxis": True,
            "xaxis_label": "",
            "show_yaxis": True,
            "yaxis_label": "",
            "yaxis_scale": "auto",
            "yaxis_tic_interval": "1.0",
            "show_ends": True,
            "show_fineprint": True,
            "color_scheme": "color_auto",
            "symbols0": "",
            "symbols1": "",
            "symbols2": "",
            "symbols3": "",
            "symbols4": "",
            "color0": "",
            "color1": "",
            "color2": "",
            "color3": "",
            "color4": "",
        }

        values.update({k: "" if v is False else str(v) for k, v in kwds.items()})
        data = urlencode(values).encode("utf-8")
        req = Request(url, data)
        response = urlopen(req)
        with open(fname, "wb") as f:
            im = response.read()
            f.write(im)

    def __format__(self, format_spec, **kwargs):
        """Return a string representation of the Motif in the given format.

        Currently supported formats:
         - clusterbuster: Cluster Buster position frequency matrix format
         - pfm : JASPAR single Position Frequency Matrix
         - jaspar : JASPAR multiple Position Frequency Matrix
         - transfac : TRANSFAC like files

        """
        if format_spec in ("pfm", "jaspar"):
            from Bio.motifs import jaspar

            motifs = [self]
            return jaspar.write(motifs, format_spec)
        elif format_spec == "transfac":
            from Bio.motifs import transfac

            motifs = [self]
            return transfac.write(motifs)
        elif format_spec == "clusterbuster":
            from Bio.motifs import clusterbuster

            motifs = [self]
            return clusterbuster.write(motifs, **kwargs)
        elif not format_spec:
            # Follow python convention and default to using __str__
            return str(self)
        else:
            raise ValueError("Unknown format type %s" % format_spec)

    def format(self, format_spec):
        """Return a string representation of the Motif in the given format.

        Currently supported formats:
         - clusterbuster: Cluster Buster position frequency matrix format
         - pfm : JASPAR single Position Frequency Matrix
         - jaspar : JASPAR multiple Position Frequency Matrix
         - transfac : TRANSFAC like files

        """
        return self.__format__(format_spec)


def write(motifs, fmt, **kwargs):
    """Return a string representation of motifs in the given format.

    Currently supported formats (case is ignored):
     - clusterbuster: Cluster Buster position frequency matrix format
     - pfm : JASPAR simple single Position Frequency Matrix
     - jaspar : JASPAR multiple PFM format
     - transfac : TRANSFAC like files

    """
    fmt = fmt.lower()
    if fmt in ("pfm", "jaspar"):
        from Bio.motifs import jaspar

        return jaspar.write(motifs, fmt)
    elif fmt == "transfac":
        from Bio.motifs import transfac

        return transfac.write(motifs)
    elif fmt == "clusterbuster":
        from Bio.motifs import clusterbuster

        return clusterbuster.write(motifs, **kwargs)
    else:
        raise ValueError("Unknown format type %s" % fmt)


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest(verbose=0)