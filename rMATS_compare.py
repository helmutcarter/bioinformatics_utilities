#!/usr/bin/env python3
import re
from datetime import datetime
from Bio import pairwise2
from pyfaidx import Fasta

#settings
testing = False
verbose = False
mean_counts = True
FDR_cutoff = 0.05
use_FDR_cutoff = True
count_cutoff=5
use_count_cutoff = True
absdPSI_cutoff = 0.10
use_absdPSI_cutoff = True
event_list = ["SE","MXE","A5SS","A3SS","RI"]  # List of which alternative splicing events you want to compare
intron_distances = [25]  # List of what intron distances you want to be considered for pairwise alignment
print_two = True


rMATS_folder_1 = "/orange/swanson/carter.h/ChP_Splicing_Project/Swanson_Lab_Human_ChP/seq_data/STAR_aligned_150_sjdb_multimapNmax_1/rMATS_DM1_5_12_vs_DM1_6_11_non_SS"
rMATS_folder_2 = "/orange/swanson/carter.h/ChP_Splicing_Project/Mbnl_12KO_ChP/rMATS_output_2KO_ChP_Triplicate"
path_to_fasta1 = "/orange/swanson/carter.h/genomes/human/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
path_to_fasta2 = "/orange/swanson/carter.h/genomes/mouse/ensembl/Mus_musculus.GRCm38.dna.primary_assembly.fa"
orthology_file = "/orange/swanson/carter.h/genomes/human/ensembl/human_mouse_homology_ensembl_IDs.csv"




if not verbose:
	omit_list = ["chr","strand","IncFormLen","SkipFormLen","alt_exon","DS_intron_5P","DS_intron_3P","US_intron_5P","US_intron_3P","MXE_1","MXE_2","inter_MXE_5P",
	"inter_MXE_3P","RI","short_exon","long_exon","US_exon_5P","US_exon_3P","DS_exon_5P","DS_exon_3P"]
	species2_omit_list = omit_list + ["IJC_SAMPLE_1","IJC_SAMPLE_2","SJC_SAMPLE_1","SJC_SAMPLE_2","PValue","FDR","IncLevel1","IncLevel2",
	]
else:
	omit_list = []
	species2_omit_list = []

print_two_dict = {"SE":"alt_exon",
"A5SS":"long_exon",
"A3SS":"long_exon",
"MXE":"MXE",
"RI":"RI"}


splicing_dict = {
				"SE":{
						"alt_exon":{"5P":"exonStart_0base","3P":"exonEnd","exon":True},
						"US_intron_5P":{"5P":"upstreamEE","3P":"exonStart_0base"},
						"US_intron_3P":{"5P":"upstreamEE","3P":"exonStart_0base"},
						"DS_intron_5P":{"5P":"exonEnd","3P":"downstreamES"},
						"DS_intron_3P":{"5P":"exonEnd","3P":"downstreamES"}},
				"A5SS":{
						"short_exon":{"5P":"shortES","3P":"shortEE","exon":True},
						"long_exon":{"5P":"longExonStart_0base","3P":"longExonEnd","exon":True},
						"DS_intron_5P":{"5P":"longExonEnd","3P":"flankingES"},
						"DS_intron_3P":{"5P":"longExonEnd","3P":"flankingES"}},
				"A3SS":{
						"short_exon":{"5P":"shortES","3P":"shortEE","exon":True},
						"US_intron_5P":{"5P":"flankingEE","3P":"shortES"},
						"US_intron_3P":{"5P":"flankingEE","3P":"shortES"},
						"long_exon":{"5P":"longExonStart_0base","3P":"longExonEnd","exon":True}},
				"MXE":{
						"MXE_1":{"5P":"1stExonStart_0base","3P":"1stExonEnd","exon":True},
						"MXE_2":{"5P":"2ndExonStart_0base","3P":"2ndExonEnd","exon":True},
						"inter_MXE_5P":{"5P":"1stExonEnd","3P":"2ndExonStart_0base"},
						"inter_MXE_3P":{"5P":"1stExonEnd","3P":"2ndExonStart_0base"},
						"US_intron_5P":{"5P":"upstreamEE","3P":"1stExonStart_0base"},
						"US_intron_3P":{"5P":"upstreamEE","3P":"1stExonStart_0base"},
						"DS_intron_5P":{"5P":"2ndExonEnd","3P":"downstreamES"},
						"DS_intron_3P":{"5P":"2ndExonEnd","3P":"downstreamES"}},
				"RI":{
						"RI":{"5P":"upstreamEE","3P":"downstreamES","exon":True},
						"US_exon_5P":{"5P":"upstreamES","3P":"upstreamEE"},
						"US_exon_3P":{"5P":"upstreamES","3P":"upstreamEE"},
						"DS_exon_5P":{"5P":"downstreamES","3P":"downstreamEE"},
						"DS_exon_3P":{"5P":"downstreamES","3P":"downstreamEE"}}
}

#turns rMATS file into a dictionary
def parse_rMATS(rMATS_path, genome_dict):
	AS_event_type = rMATS_path.split("/")[-1].split(".")[0]
	with open(rMATS_path, "r") as f:
		header_dict = {}
		rMATS_dict = {}
		found_list = []
		temp_dict = {}
		species = ""
		for line in f:
			splitline = line.strip().replace('"',"").split("\t")
			if splitline[0] == "ID":
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
			if use_FDR_cutoff and temp_dict[event_id]["FDR"] != "NA" and float(temp_dict[event_id]["FDR"]) > FDR_cutoff:
				continue
			for seq in splicing_dict[AS_event_type]:
				start = splicing_dict[AS_event_type][seq]["5P"]
				stop = splicing_dict[AS_event_type][seq]["3P"]
				temp_dict[event_id][seq] = genome_dict[temp_dict[event_id]["chr"]][min(int(temp_dict[event_id][start]),int(temp_dict[event_id][stop])):max(int(temp_dict[event_id][start]),int(temp_dict[event_id][stop]))]
				if "exon" not in splicing_dict[AS_event_type][seq] and len(temp_dict[event_id][seq]) > intron_distance:
					if seq[-2:] == "5P":
						temp_dict[event_id][seq] = temp_dict[event_id][seq][:intron_distance]
					elif seq[-2:] == "3P":
						temp_dict[event_id][seq] = temp_dict[event_id][seq][-intron_distance:]


				if temp_dict[event_id]["strand"] == "-":
					temp_dict[event_id][seq] = temp_dict[event_id][seq].reverse.complement


				#makes sure the sequence is a string
				temp_dict[event_id][seq] = str(temp_dict[event_id][seq])

			ensembl_gene = splitline[1]
			#get the average number of counts
			count_string = temp_dict[event_id]["IJC_SAMPLE_1"] + "," + temp_dict[event_id]["IJC_SAMPLE_2"] + "," + temp_dict[event_id]["SJC_SAMPLE_1"] + "," + temp_dict[event_id]["SJC_SAMPLE_2"]
			total = 0
			if "" in count_string.split(","):
				temp_dict[event_id]["avg_count"] = "NA"

			else:
				for count in count_string.split(","):
					total += int(count)
				temp_dict[event_id]["avg_count"] = str(total / len(count_string.split(",")))
				if use_count_cutoff and float(temp_dict[event_id]["avg_count"]) < count_cutoff:
					continue

			if temp_dict[event_id]["IncLevelDifference"] != "NA":
				temp_dict[event_id]["abs_IncLevelDifference"] = str(abs(float(temp_dict[event_id]["IncLevelDifference"])))
				if use_absdPSI_cutoff and abs(float(temp_dict[event_id]["IncLevelDifference"])) < 0.10:
					continue
			else:
				temp_dict[event_id]["abs_IncLevelDifference"] = "NA"

			#adds the temp dict to the list of events for a given gene
			if ensembl_gene not in rMATS_dict:
				event_count = 0
				rMATS_dict[ensembl_gene] = {}
				species = identify_organism(ensembl_gene)
			identifier = ensembl_gene + "_" + str(event_count)
			rMATS_dict[ensembl_gene][identifier] = temp_dict[event_id]
			event_count += 1

	print(str(len(rMATS_dict)) + " " + AS_event_type + " events passing cutoffs for " + species)
	return rMATS_dict, species

#turns an ensembl biomart file with gene orthology into a dictionary
def orthology_parse(orthology_file, species1, species2):
	orthology_dict = {}
	ensembl_species_dict = {"human":"ENSG",
							"rat":"ENSRNOG",
							"mouse":"ENSMUSG"}
	species1 = ensembl_species_dict[species1]
	species2 = ensembl_species_dict[species2]

	with open(orthology_file, "r") as f:
		orthology_matches = False
		for line in f:
			if line[0:3] != "ENS":
				continue
			splitline = line.strip().split(",")
			if species1 == splitline[0][:len(species1)] and species2 == splitline[1][:len(species2)]:
				orthology_dict[splitline[0]] = splitline[1]
				orthology_matches = True
			elif species1 == splitline[1][:len(species1)] and species2 == splitline[0][:len(species2)]:
				orthology_dict[splitline[1]] = splitline[0]
				orthology_matches = True
		if orthology_matches == False:
			print("error, orthology dict doesn't match species of rMATS files")
			print(species1)
			print(splitline[0][:len(species1)])
			print(species2)
			print(splitline[1][:len(species2)])
			print(line)

		return orthology_dict

#takes an ENSEMBL id and returns the organism
def identify_organism(gene_id):
	if gene_id[:4] == "ENSG":
		return "human"
	elif gene_id[:7] == "ENSRNOG":
		return "rat"
	elif gene_id[:7] == "ENSMUSG":
		return "mouse"

#returns the length of whichever string is longer
def longest_of_two(string1,string2):
	if len(string1) >= len(string2):
		return len(string1)
	if len(string1) < len(string2):
		return len(string2)

if __name__ == "__main__":
	orthology_dict = None
	print("reading chromosomes")
	genome_dict1 = Fasta(path_to_fasta1)
	genome_dict2 = Fasta(path_to_fasta2)
	print("done reading chromosomes")

	for intron_distance in intron_distances:
		filename_ending = f".DM1_human_ChP_outlier_contrast_vs_2KO_ChP.{intron_distance}intron.{count_cutoff}count_and_FDR_cutoff.tsv"
		for event in event_list:
			out_list = []
			rMATS_path1 = f"{rMATS_folder_1}/{event}.MATS.JCEC.txt"
			rMATS_path2 = f"{rMATS_folder_2}/{event}.MATS.JCEC.txt"
			filename = rMATS_path1 + filename_ending
			AS_event_type = rMATS_path1.split("/")[-1].split(".")[0]

			print("reading rMATS files")
			rMATS_dict1, species1 = parse_rMATS(rMATS_path1, genome_dict1)
			rMATS_dict2, species2 = parse_rMATS(rMATS_path2, genome_dict2)
			if orthology_dict == None:
				print("reading orthology file")
				orthology_dict = orthology_parse(orthology_file, species1, species2)

			print("calculating scores and writing to " + filename)

			#figure out the gene with the most splicing event in rMATS_dict2
			if verbose:
				longest_dict2_entry = 0
				for gene in rMATS_dict2:
					if len(rMATS_dict2[gene]) > longest_dict2_entry:
						longest_dict2_entry = len(rMATS_dict2[gene])

			#write the header
			header_list = []
			header_written1 = False
			header_written2 = False
			header_written3 = False
			out = open(filename, "w+")
			out.close()
			for gene1 in list(rMATS_dict1):
				for event1 in rMATS_dict1[gene1]:
					out = open(filename, "a+")

					#adds values we're interested in to the buffer to write to output
					for key in [key for key in rMATS_dict1[gene1][event1] if key not in omit_list]:
						if header_written1 == False:
							header_list.append(species1 + "_" + key)
						out_list.append(rMATS_dict1[gene1][event1][key])
					header_written1 = True
					if gene1 in orthology_dict and orthology_dict[gene1] in rMATS_dict2:
						gene2 = orthology_dict[gene1]
						for event2 in rMATS_dict2[gene2]:

							#calculates pairwise alignment scores for each set of regions
							if AS_event_type == "MXE":
								MXE_1_1_score = pairwise2.align.globalms(rMATS_dict1[gene1][event1]["MXE_1"], rMATS_dict2[gene2][event2]["MXE_1"], 2, -1, -.5, -.1, score_only=testing, one_alignment_only=True)
								MXE_2_2_score = pairwise2.align.globalms(rMATS_dict1[gene1][event1]["MXE_2"], rMATS_dict2[gene2][event2]["MXE_2"], 2, -1, -.5, -.1, score_only=testing, one_alignment_only=True)
								MXE_1_2_score = pairwise2.align.globalms(rMATS_dict1[gene1][event1]["MXE_1"], rMATS_dict2[gene2][event2]["MXE_2"], 2, -1, -.5, -.1, score_only=testing, one_alignment_only=True)
								MXE_2_1_score = pairwise2.align.globalms(rMATS_dict1[gene1][event1]["MXE_2"], rMATS_dict2[gene2][event2]["MXE_1"], 2, -1, -.5, -.1, score_only=testing, one_alignment_only=True)

								MXE_1_1_score = MXE_1_1_score[0][2] / MXE_1_1_score[0][-1]
								MXE_2_2_score = MXE_2_2_score[0][2] / MXE_2_2_score[0][-1]
								MXE_1_2_score = MXE_1_2_score[0][2] / MXE_1_2_score[0][-1]
								MXE_2_1_score = MXE_2_1_score[0][2] / MXE_2_1_score[0][-1]

								#decide which MXE matches which MXE
								if MXE_1_1_score + MXE_2_2_score > MXE_1_2_score + MXE_2_1_score:
									rMATS_dict2[gene2][event2]["5P_MXE_score"] = MXE_1_1_score
									rMATS_dict2[gene2][event2]["3P_MXE_score"] = MXE_2_2_score
									rMATS_dict2[gene2][event2]["total_score"] = MXE_1_1_score+MXE_2_2_score
								else:
									#switch the order of the MXEs for species 2
									rMATS_dict2[gene2][event2]["5P_MXE_score"] = MXE_1_2_score
									rMATS_dict2[gene2][event2]["3P_MXE_score"] = MXE_2_1_score
									rMATS_dict2[gene2][event2]["total_score"] = MXE_1_2_score+MXE_2_1_score
								rMATS_dict2[gene2][event2]["MXE_score"] = rMATS_dict2[gene2][event2]["total_score"]
							else:
								rMATS_dict2[gene2][event2]["total_score"] = 0
							for seq in splicing_dict[AS_event_type]:

								if seq == "MXE_1" or seq == "MXE_2":
									continue

								score = pairwise2.align.globalms(rMATS_dict1[gene1][event1][seq], rMATS_dict2[gene2][event2][seq], 2, -1, -.5, -.1, score_only=testing, one_alignment_only=True)
								if testing == True:
									rMATS_dict2[gene2][event2][seq + "_score"] = score / longest_of_two(rMATS_dict1[gene1][event1][seq], rMATS_dict2[gene2][event2][seq])
								else:
									if score != []:
										rMATS_dict2[gene2][event2][seq + "_score"] = score[0][2] / score[0][-1]
									else:
										rMATS_dict2[gene2][event2][seq + "_score"] = 0
								rMATS_dict2[gene2][event2]["total_score"] += rMATS_dict2[gene2][event2][seq + "_score"]


							if header_written3 == False:
								for key in [key for key in rMATS_dict2[gene2][event2] if key not in species2_omit_list]:
									header_list.append(species2 + "_" + key)

								if header_written2 == True:
									header_written3 = True
									out.write("\t".join(header_list) + "\n")
								header_written2 = True


						for event2 in sorted(rMATS_dict2[gene2].items(), key = lambda k_v: k_v[1]['total_score'], reverse = True):
							for key in [key for key in rMATS_dict2[gene2][event2[0]] if key not in species2_omit_list]:
								out_list.append(str(rMATS_dict2[gene2][event2[0]][key]))
							if not verbose:
								break
						if print_two:
							sortby = print_two_dict[AS_event_type] + "_score"

							for event2 in sorted(rMATS_dict2[gene2].items(), key = lambda k_v: k_v[1][sortby], reverse = True):
								for key in [key for key in rMATS_dict2[gene2][event2[0]] if key not in species2_omit_list]:
									out_list.append(str(rMATS_dict2[gene2][event2[0]][key]))
								if not verbose:
									break
					if header_written3:
						#have to do this a funky way otherwise the first couple lines will be indented a tab when they shouldn't be
						out_buffer = "\t".join(out_list)
						out_buffer = out_buffer.split("\n")
						for line in out_buffer:
							out.write(line.strip() + "\n")
						out_list = []
					else:
						out_list.append("\n")
					out.close()
	print("done!")
