import re

def calculate_gc_content(sequence):
	"""
	Calculates the GC content of a DNA sequence.
	
	Args:
		sequence (str): The DNA sequence.
	
	Returns:
		float: The GC content as a percentage.
	"""
	gc_count = sequence.upper().count('G') + sequence.upper().count('C')
	total_count = len(sequence)
	gc_content = (gc_count / total_count) * 100
	return gc_content

def find_sequence_motif(sequence, motif):
	"""
	Finds all occurrences of a motif in a DNA sequence.
	
	Args:
		sequence (str): The DNA sequence.
		motif (str): The motif to search for.
	
	Returns:
		list: A list of starting positions of the motif in the sequence.
	"""
	positions = []
	pattern = re.compile(motif)
	matches = pattern.finditer(sequence)
	for match in matches:
		positions.append(match.start())
	return positions

def transcribe_dna_to_rna(sequence):
	"""
	Transcribes a DNA sequence to RNA by replacing T with U.
	
	Args:
		sequence (str): The DNA sequence.
	
	Returns:
		str: The transcribed RNA sequence.
	"""
	rna_sequence = sequence.replace('T', 'U')
	return rna_sequence

def translate_rna_to_protein(sequence):
	"""
	Translates an RNA sequence to a protein sequence using the standard genetic code.
	
	Args:
		sequence (str): The RNA sequence.
	
	Returns:
		str: The translated protein sequence.
	"""
	codon_table = {
		'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
		'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
		'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
		'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
		'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
		'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
		'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
		'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
		'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
		'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
		'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
		'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
		'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
		'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
		'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
		'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
	}
	
	protein_sequence = ''
	codon = ''
	for base in sequence:
		codon += base
		if len(codon) == 3:
			protein_sequence += codon_table.get(codon, '')
			codon = ''
	return protein_sequence

def calculate_molecular_weight(sequence):
	"""
	Calculates the molecular weight of a DNA or protein sequence.
	
	Args:
		sequence (str): The DNA or protein sequence.
	
	Returns:
		float: The molecular weight of the sequence.
	"""
	dna_weights = {'A': 331.2, 'C': 307.2, 'G': 347.2, 'T': 322.2}
	protein_weights = {'A': 71.08, 'C': 103.14, 'D': 115.09, 'E': 129.12, 'F': 147.18, 'G': 57.05,
					   'H': 137.15, 'I': 113.17, 'K': 128.18, 'L': 113.17, 'M': 131.21, 'N': 114.11,
					   'P': 97.12, 'Q': 128.13, 'R': 156.19, 'S': 87.08, 'T': 101.11, 'V': 99.14,
					   'W': 186.21, 'Y': 163.18}
	
	sequence = sequence.upper()
	molecular_weight = 0
	
	if sequence.startswith('ATG'):
		# Assuming it's a protein sequence
		for amino_acid in sequence:
			molecular_weight += protein_weights.get(amino_acid, 0)
	else:
		# Assuming it's a DNA sequence
		for nucleotide in sequence:
			molecular_weight += dna_weights.get(nucleotide, 0)
	
	return molecular_weight


def calculate_tm(sequence):
	"""
	Calculates the melting temperature (Tm) of a DNA sequence using the nearest-neighbor method.
	
	Args:
		sequence (str): The DNA sequence.
	
	Returns:
		float: The melting temperature in degrees Celsius.
	"""
	# Melting temperature constants
	salt_correction = 16.6
	mismatch_correction = 0.41
	
	sequence = sequence.upper()
	
	# Melting temperature calculation
	tm = 4 * (sequence.count('G') + sequence.count('C')) + 2 * (sequence.count('A') + sequence.count('T'))
	tm += salt_correction * (sequence.count('K') + sequence.count('M'))
	tm -= mismatch_correction * sequence.count('S')
	
	return tm


def reverse_complement(sequence):
	"""
	Computes the reverse complement of a DNA sequence.
	
	Args:
		sequence (str): The DNA sequence.
	
	Returns:
		str: The reverse complement sequence.
	"""
	complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
	reverse_sequence = sequence[::-1]
	reverse_complement_sequence = ''.join(complement.get(base, base) for base in reverse_sequence)
	
	return reverse_complement_sequence


def calculate_hamming_distance(seq1, seq2):
	"""
	Calculates the Hamming distance between two sequences of equal length.
	
	Args:
		seq1 (str): The first sequence.
		seq2 (str): The second sequence.
	
	Returns:
		int: The Hamming distance between the sequences.
	"""
	if len(seq1) != len(seq2):
		raise ValueError("Sequences must have equal length.")
	
	hamming_distance = sum(base1 != base2 for base1, base2 in zip(seq1, seq2))
	
	return hamming_distance
