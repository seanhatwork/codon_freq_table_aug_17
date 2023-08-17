# Codon Frequency Calculator

The Codon Frequency Calculator is a Python package that calculates codon frequencies of given gene sequences and generates a table containing usage information. 

## Features

- Loads gene sequences from a FASTA file.
- Counts the frequency of each codon in the sequences.
- Calculates the percentage of each codon frequency out of the total.
- Calculates the percentage of each codon frequency within the same amino acid.
- Writes the codon usage information to a CSV file.

## Usage

1. Install the required dependencies by running:
   
   ```bash
   pip install biopython
   ```

2. Place your gene sequences in a FASTA file (e.g., `U00096.3_gene_sequences.fasta`).

3. Run the `codon_freq_aug_10.py` script by providing the input file (FASTA) and the output file (CSV) as arguments:

   ```bash
   python codon_freq_aug_10.py U00096.3_gene_sequences.fasta U00096.3_codon_usage_amino_acid.csv
   ```

4. The script will load the sequences, calculate codon frequencies, and generate a CSV file containing usage information for each codon, both overall and within the same amino acid.

## Input File Format

The input file should be in FASTA format.

## Output File Format

The output CSV file will have the following columns:

- Amino Acid: The corresponding amino acid for the codon.
- Codon: The codon sequence.
- Percentage(All): The percentage of the codon frequency out of all codons.
- Percentage(Within the same AA): The percentage of the codon frequency within the same amino acid.

## Dependencies

- [BioPython](https://biopython.org/): A library for working with biological sequences.

## Disclaimer

This tool is intended for research and analysis purposes and should not be considered as a substitute for professional expertise in genetics or bioinformatics.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
