#!/usr/bin/env python
#coded by Helmut A. Carter in the laboratory of Maurice S Swanson in the month of April 2021


motif_of_interest = "YGCY"
query_sequence = "CTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG"

#tells the motif_check function which degenerate bases are equal to which canonical bases
degenerate_base_dict = {
"A":("A"),
"C":("C"),
"T":("T","U"),
"G":("G"),
"U":("T","U"),
"Y":("C","T","U"),
"R":("A","G"),
"W":("A","T","U"),
"S":("C","G"),
"K":("G","T","U"),
"M":("A","C"),
"D":("A","T","G","U"),
"V":("A","C","G"),
"H":("A","T","C","U"),
"B":("T","C","G","U"),
"N":("A","T","C","G","U")
}

#given a DNA or RNA sequence and motif sequence (assumed to be of equal length), determines whether or not the
#DNA sequence and motif sequence are equivalent and returns a boolean accordingly
def motif_check(input_str, motif_seq):
	for n in range(len(motif_seq)):
		if input_str[n] not in degenerate_base_dict[motif_seq[n]]:
			return False

	return True

#given a DNA or RNA sequence of any length, counts and returns how many times the given motif occurs in the sequence
#this is the only function you should directly call
def motif_counter(input_str, motif_seq):
	base_position = 0
	motif_count = 0
	motif_length = len(motif_seq)
	while base_position < len(input_str)-motif_length:
		if motif_check(input_str[base_position:base_position+motif_length], motif_seq) == True:
			motif_count += 1
		base_position += 1
	return motif_count


if __name__ == "__main__":
	print("This is a usage example. Modify the code or incorporate the functions into your code to use.")
	print(f"The motif {motif_of_interest} occurs in the sequence {query_sequence} {motif_counter(query_sequence, motif_of_interest)} times.")
