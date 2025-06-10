
# --- Translation table for reverse complement ---
_complement_trans_str = str.maketrans('acgtACGT', 'TGCATGCA')

# --- Standard genetic code dictionary for DNA codon translation ---
genetic_code = {
  'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
  'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
  'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
  'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
  'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
  'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
  'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
  'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
  'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
  'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
  'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
  'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
  'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
  'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
  'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
  'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
}

# --- Function to get the reverse complement of a DNA sequence ---
def dna_rev_comp(sequence: str) -> str:
  """
  Returns the uppercase reverse complement of the input DNA sequence.

  Parameters:
  sequence (str): DNA sequence (e.g., 'atgcg')

  Returns:
  str: Reverse complement in uppercase
  """
  return sequence.translate(_complement_trans_str)[::-1]


# --- Function to translate DNA sequence to protein sequence ---
def translate_dna2aa(orf: str) -> str:
  """
  Translates a DNA sequence into a protein sequence.
  Unrecognized codons are translated as 'X'.

  Parameters:
  orf (str): Open reading frame DNA sequence (length must be divisible by 3)

  Returns:
  str: Amino acid sequence
  """
  protein = ''

  for i in range(0, (len(orf) // 3) * 3, 3):
    try:
      protein += genetic_code[orf[i:i + 3]]
    except KeyError:
      protein += 'X'  # Assign 'X' for unrecognized codons

  return protein


# --- Function to find the first occurrence of any value from a list in a string ---
def find_(string, value_list):
    """
    Finds the first occurrence of any value in a list within a string.

    Parameters:
    string (str): Input string to search (e.g. quality scores)
    value_list (list): List of characters to search for (e.g. ['#', '!', '"'])

    Returns:
    int: Index of the first match, or raises an error if no match is found
    """
    indexes = [string.find(letter) for letter in value_list]
    ind = min([index for index in indexes if index != -1])
    return ind
