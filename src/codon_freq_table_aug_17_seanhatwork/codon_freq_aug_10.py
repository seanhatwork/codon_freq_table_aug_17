# codon_freq_aug_10.py

from Bio import SeqIO

def load_sequences(input_file):
    """
    Loads gene sequences from a FASTA file and returns a list of sequences.

    Parameters
    ----------
    input_file : str
        The filename of a FASTA file containing the gene sequences.

    Returns
    -------
    list
        A list of gene sequences.
    """
    try:
        sequences = list(SeqIO.parse(input_file, "fasta"))
    except FileNotFoundError:
        raise FileNotFoundError("Input file '{}' not found".format(input_file))
    
    return sequences

def count_codon_frequencies(sequences):
    """
    Counts the frequency of each codon in a list of gene sequences.

    Parameters
    ----------
    sequences : list
        A list of gene sequences.

    Returns
    -------
    dict
        A dictionary containing codon frequencies.
    """
    codon_freqs = {}

    for sequence in sequences:
        gene_seq = str(sequence.seq)
        for i in range(0, len(gene_seq), 3):
            codon = gene_seq[i:i+3]
            if codon in codon_freqs:
                codon_freqs[codon] += 1
            else:
                codon_freqs[codon] = 1

    return codon_freqs

def calculate_codon_percentages(codon_freqs):
    """
    Calculates the percentage of each codon frequency out of the total.

    Parameters
    ----------
    codon_freqs : dict
        A dictionary containing codon frequencies.

    Returns
    -------
    dict
        A dictionary containing codon percentages.
    """
    total_codons = sum(codon_freqs.values())
    codon_percentages = {}

    for codon, freq in codon_freqs.items():
        codon_percentages[codon] = (freq / total_codons) * 100

    return codon_percentages

def calculate_codon_percentages_by_amino_acid(codon_freqs, codon_amino_mapping):
    """
    Calculates the percentage of each codon frequency within the same amino acid.

    Parameters
    ----------
    codon_freqs : dict
        A dictionary containing codon frequencies.
    codon_amino_mapping : dict
        A dictionary mapping codons to amino acids.

    Returns
    -------
    dict
        A dictionary containing codon percentages within each amino acid.
    """
    codon_percentages_by_amino_acid = {}

    total_freqs_by_amino_acid = {amino_acid: 0 for amino_acid in codon_amino_mapping}

    for amino_acid, codons in codon_amino_mapping.items():
        for codon_by_amino_acid in codons:
            total_freqs_by_amino_acid[amino_acid] += codon_freqs[codon_by_amino_acid]

    for amino_acid, codons in codon_amino_mapping.items():
        for codon_by_amino_acid in codons:
            codon_percentages_by_amino_acid[codon_by_amino_acid] = (codon_freqs[codon_by_amino_acid] / total_freqs_by_amino_acid[amino_acid]) * 100

    return codon_percentages_by_amino_acid

def write_csv(output_file, codon_amino_mapping, codon_percentages, codon_percentages_by_amino_acid):
    """
    Writes codon usage information to a CSV file.

    Parameters
    ----------
    output_file : str
        The filename of the CSV file to write the codon usage information to.
    codon_amino_mapping : dict
        A dictionary mapping codons to amino acids.
    codon_percentages : dict
        A dictionary containing codon percentages out of all codons.
    codon_percentages_by_amino_acid : dict
        A dictionary containing codon percentages within the same amino acid.

    Returns
    -------
    None
    """
    with open(output_file, "w") as f:
        f.write("Amino Acid, Codon, Percentage(All), Percentage(Within the same AA)\n")

        for amino_acid, codons in codon_amino_mapping.items():
            for codon_by_amino_acid in codons:
                f.write("{},{},{},{}\n".format(
                    amino_acid,
                    codon_by_amino_acid,
                    round(codon_percentages[codon_by_amino_acid], 2),
                    round(codon_percentages_by_amino_acid[codon_by_amino_acid], 2)
                ))

def main(input_file, output_file):

    codon_amino_mapping = {
        # Start codons
        "Methionine": ["ATG"],
        
        # Phenylalanine
        "Phenylalanine": ["TTT", "TTC"],
        
        # Leucine
        "Leucine": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
        
        # Isoleucine
        "Isoleucine": ["ATT", "ATC", "ATA"],
        
        # Valine
        "Valine": ["GTT", "GTC", "GTA", "GTG"],
        
        # Serine
        "Serine": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
        
        # Proline
        "Proline": ["CCT", "CCC", "CCA", "CCG"],
        
        # Threonine
        "Threonine": ["ACT", "ACC", "ACA", "ACG"],
        
        # Alanine
        "Alanine": ["GCT", "GCC", "GCA", "GCG"],
        
        # Tyrosine
        "Tyrosine": ["TAT", "TAC"],
        
        # Histidine
        "Histidine": ["CAT", "CAC"],
        
        # Glutamine
        "Glutamine": ["CAA", "CAG"],
        
        # Asparagine
        "Asparagine": ["AAT", "AAC"],
        
        # Lysine
        "Lysine": ["AAA", "AAG"],
        
        # Aspartic acid
        "Aspartic acid": ["GAT", "GAC"],
        
        # Glutamic acid
        "Glutamic acid": ["GAA", "GAG"],
        
        # Cysteine
        "Cysteine": ["TGT", "TGC"],
        
        # Tryptophan
        "Tryptophan": ["TGG"],
        
        # Arginine
        "Arginine": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
        
        # Glycine
        "Glycine": ["GGT", "GGC", "GGA", "GGG"],
        
        # Stop codons
        "Stop": ["TAA", "TGA", "TAG"],
    }

    sequences = load_sequences(input_file)
    codon_freqs = count_codon_frequencies(sequences)
    codon_percentages = calculate_codon_percentages(codon_freqs)
    codon_percentages_by_amino_acid = calculate_codon_percentages_by_amino_acid(codon_freqs, codon_amino_mapping)
    write_csv(output_file, codon_amino_mapping, codon_percentages, codon_percentages_by_amino_acid)

if __name__ == "__main__":
    main("U00096.3_gene_sequences.fasta", "U00096.3_codon_usage_amino_acid.csv")

