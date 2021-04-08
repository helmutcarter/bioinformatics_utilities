#!/usr/bin/env python3
import re
from pyfaidx import Fasta
from Bio.Seq import Seq
from os import mkdir

#settings
FDR_cutoff = 0.05  # What rMATS FDR value you want to use as a signifigance cutoff. No splicing events with values above this will be included in output.
use_count_cutoff = True  # Boolean indicating whether or not you want to exclude splicing events with less than a set number of counts.
count_cutoff = 5  # If use_count_cutoff is True, this is the cutoff for the fewest counts allowed.
SJ_flank = 100  # How many intronic bases flanking each splice junction you want the output to include.
use_JCEC = True  # Whether the program should use JCEC or JC rMATS input. JC = junction counts only, JCEC = junction counts and exon counts.
output_folder = "/path/to/your/desired/output"  # Where you want the results to go. Directory does not have to exist.
rMATS_folder = "/path/to/your/rMATS/results"  # Path to rMATS output to be processed.
path_to_fasta = "/path/to/appropriate/fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa"  # Path to appropriate primary assembly fasta for organism of interest. Must be indexed to allow fast, random access.




# This function goes through the SE rMATS output of interest, applies the cutoffs dictated in #settings, and returns a dictionary containing all the required info and sequences of interest.
def parse_rMATS(rMATS_path, genome):
	with open(rMATS_path, "r") as f:
		header_dict = {}
		rMATS_dict = {}
		found_list = []
		temp_dict = {}
		species = ""
		for line in f:
			splitline = line.strip().replace('"',"").split("\t")
			if splitline[0] == "ID":  # Determines if the current line is the header, and if it is, write each column to a dicitonary.
				for n in range(len(splitline)):
					header_dict[n] = splitline[n]
				continue
			event_id = splitline[0]

			#make temp dict and fill it with info from the current line
			temp_dict[event_id] = {}
			for index in header_dict:
				if header_dict[index] == "chr" and splitline[index][0:3] == "chr":
					temp_dict[event_id][header_dict[index]] = splitline[index][3:]
				else:
					temp_dict[event_id][header_dict[index]] = splitline[index]


			# The following code block fills the current splicing event's dictionary with the sequences for each region defined in splicing_dict: the alternative exon, and the distal and proximal ends of the upstream and downstream introns.
			temp_dict[event_id]["alt_exon"] = str(genome[temp_dict[event_id]["chr"]][int(temp_dict[event_id][splicing_dict["SE"]["alt_exon"]["5P"]]):int(temp_dict[event_id][splicing_dict["SE"]["alt_exon"]["3P"]])])
			US_exon_seq = str(genome[temp_dict[event_id]["chr"]][int(temp_dict[event_id][splicing_dict["SE"]["US_intron_5P"]["5P"]]):int(temp_dict[event_id][splicing_dict["SE"]["US_intron_5P"]["3P"]])])
			temp_dict[event_id]["US_intron_5P"] = US_exon_seq[:SJ_flank]
			temp_dict[event_id]["US_intron_3P"] = US_exon_seq[-SJ_flank:]
			DS_exon_seq = str(genome[temp_dict[event_id]["chr"]][int(temp_dict[event_id][splicing_dict["SE"]["DS_intron_5P"]["5P"]]):int(temp_dict[event_id][splicing_dict["SE"]["DS_intron_5P"]["3P"]])])
			temp_dict[event_id]["DS_intron_5P"] = DS_exon_seq[:SJ_flank]
			temp_dict[event_id]["DS_intron_3P"] = DS_exon_seq[-SJ_flank:]

			if temp_dict[event_id]["strand"] == "-":  # If the gene for the alternative exon is on the minus strand, get the referse complement of all the sequences of interest.
				for region in ("alt_exon", "US_intron_5P", "US_intron_3P", "DS_intron_5P", "DS_intron_3P"):
					temp_dict[event_id][region] = Seq(temp_dict[event_id][region]).reverse_complement()
					temp_dict[event_id][region] = str(temp_dict[event_id][region])  # Makes sure the sequence object is a string before writing it to the output dictionary.


			# Calculates the the average number of counts for the current splicing event.
			count_string = temp_dict[event_id]["IJC_SAMPLE_1"] + "," + temp_dict[event_id]["IJC_SAMPLE_2"] + "," + temp_dict[event_id]["SJC_SAMPLE_1"] + "," + temp_dict[event_id]["SJC_SAMPLE_2"]
			total = 0
			if "" in count_string.split(","):
				temp_dict[event_id]["avg_count"] = "NA"
			else:
				for count in count_string.split(","):
					total += int(count)
				temp_dict[event_id]["avg_count"] = str(total / len(count_string.split(",")))


			#adds the temp dict to the list of events for a given gene
			ensembl_gene = splitline[1]
			if ensembl_gene not in rMATS_dict:  # If the current splicing event is the first for the parent gene, adds the gene name to the dictionary so the event can belong to it.
				event_count = 0
				rMATS_dict[ensembl_gene] = {}
			identifier = ensembl_gene + "_" + str(event_count)  # In order to give each splicing event a unique identifier, successive events within the same gene are given numerically incremented names.
			rMATS_dict[ensembl_gene][identifier] = temp_dict[event_id]
			event_count += 1
	return rMATS_dict


if __name__ == "__main__":

	splicing_dict = {
					"SE":{
							"alt_exon":{"5P":"exonStart_0base","3P":"exonEnd"},
							"US_intron_5P":{"5P":"upstreamEE","3P":"exonStart_0base"},
							"US_intron_3P":{"5P":"upstreamEE","3P":"exonStart_0base"},
							"DS_intron_5P":{"5P":"exonEnd","3P":"downstreamES"},
							"DS_intron_3P":{"5P":"exonEnd","3P":"downstreamES"}}
	}  # Defines what rMATS entries correspond to regions of interest for the output.

	try:  # Creates the output directory if it does not already exist.
		mkdir(output_folder)
	except error as e:
		print(e)

	if use_JCEC == True:
		JCEC_or_JC = "JCEC"  # Decides whether to open the JCEC or JC file based on what the user supplied in #settings.
	else:
		JCEC_or_JC = "JC"

	genome = Fasta(path_to_fasta) # Instantiate the fasta access object.

	# Iterates over the 5 regions defined in splicing_dict: the alternative exon, and the distal and proximal ends of the upstream and downstream introns.
	for region in splicing_dict["SE"]:
		# Instantiates four output files for each region: one for upregulated splicing events, one for downregulated splicing events,
		# one for splicing events both up and downregulted, and one for splicing events with no significant change (to be used as background).
		with open(f"{output_folder}/{region}_up.fasta", "w+") as up_file, open(f"{output_folder}/{region}_down.fasta", "w+") as down_file, open(f"{output_folder}/{region}_both.fasta", "w+") as both_file, open(f"{output_folder}/{region}_bg.fasta", "w+") as bg_file:
			rMATS_path = f"{rMATS_folder}/SE.MATS.{JCEC_or_JC}.txt"

			rMATS_dict = parse_rMATS(rMATS_path, genome)
			for gene in rMATS_dict: # Iterates over every gene with at least one alternative splicing event in the rMATS_dict returned by parse_rMATS().
				for event in list(rMATS_dict[gene].keys()): # Iterates over every alternative splicing event in each gene in the rMATS_dict returned by parse_rMATS().
					if use_count_cutoff and float(rMATS_dict[gene][event]["avg_count"]) < count_cutoff: # Checks if a count cutoff is to be used, and applies it if appropriate based on #settings.
						if float(rMATS_dict[gene][event]["FDR"]) < FDR_cutoff: # Applies the FDR cutoff supplied in #settings.

							# The following code block sorts alternative splicing events based on delta PSI and writes it to the appropriate output file(s).
							if abs(float(rMATS_dict[gene][event]["IncLevelDifference"])) > 0.1:
								if float(rMATS_dict[gene][event]["IncLevelDifference"]) > 0:
									up_file.write(f">{event}\n{rMATS_dict[gene][event][region]}\n")
								else:
									down_file.write(f">{event}\n{rMATS_dict[gene][event][region]}\n")
								both_file.write(f">{event}\n{rMATS_dict[gene][event][region]}\n")
							else:
								bg_file.write(f">{event}\n{rMATS_dict[gene][event][region]}\n")
						else:
							bg_file.write(f">{event}\n{rMATS_dict[gene][event][region]}\n")
					else:
						bg_file.write(f">{event}\n{rMATS_dict[gene][event][region]}\n")
