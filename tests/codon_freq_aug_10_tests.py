# codon_freq_aug_10_tests.py
# total of 21 tests
# expect 4 failures and 17 passes

import doctest
from codon_freq_aug_9 import load_sequences, count_codon_frequencies, calculate_codon_percentages, calculate_codon_percentages_by_amino_acid

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def test_load_sequences():
    """
    >>> test_file = "test_sequences.fasta"

    "test_sequences.fasta" might look:
    ----------
    >seq1
    ATGACGTAG
    >seq2
    TTGATCGTT
    >seq3
    GAACTGACC
    ----------

    >>> test_sequences = load_sequences(test_file)
    >>> len(test_sequences)
    3
    >>> test_sequences[0].id
    'seq1'
    >>> test_sequences[1].seq
    Seq('ATGACGTAG') # This will fail as "TTGATCGTT" is the sequence of seq2
    """

def test_count_codon_frequencies():
    """
    >>> test_sequences = [SeqRecord(Seq("ATGACGTAG"), id="seq1"), SeqRecord(Seq("TTGATCACG"), id="seq2")]
    >>> codon_freqs = count_codon_frequencies(test_sequences)
    >>> codon_freqs["ATG"]
    1
    >>> codon_freqs["TTG"]
    1
    >>> codon_freqs["ACG"]
    1  # This will fail as the actual occurence of "ACG" is 2
    """

def test_calculate_codon_percentages():
    """
    >>> test_codon_freqs = {"ATG": 5, "TTG": 3, "GAA": 2}
    >>> codon_percentages = calculate_codon_percentages(test_codon_freqs)
    >>> codon_percentages["ATG"]
    50.0
    >>> codon_percentages["TTG"]
    30.0
    >>> codon_percentages["GAA"]
    10.0  # This will fail as the actual occurrence of "GAA" is 20%
    """

def test_calculate_codon_percentages_by_amino_acid():
    """
    >>> test_codon_freqs = {"ATG": 5, "TTT": 3, "TTC": 7, "GAA": 2, "GAG": 3}
    >>> test_codon_amino_mapping = {"Methionine": ["ATG"], "Phenylalanine": ["TTT", "TTC"], "Glutamic acid": ["GAA", "GAG"]}
    >>> codon_percentages_by_amino_acid = calculate_codon_percentages_by_amino_acid(test_codon_freqs, test_codon_amino_mapping)
    >>> codon_percentages_by_amino_acid["ATG"]
    100.0
    >>> codon_percentages_by_amino_acid["TTT"]
    30.0
    >>> codon_percentages_by_amino_acid["GAG"]
    10.0  # This will fail as the actual occurrence of "GAA" within "Glutamic acid" is 60%
    """

if __name__ == "__main__":
    doctest.testmod()

