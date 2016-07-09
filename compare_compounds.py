"""
Calculate the distance matrix for a set of compounds.
Assumes Open Babel has been installed to path.

User provides a directory of compounds in which to compare.

If a database of compounds is provided, will compare each of the input compounds
against the database. Recommend the database is in fastsearch format (fs).

It is assumed that the input compound filenames match those in the assignment file.
"""

import sys, os
import subprocess
import argparse
from collections import defaultdict

def _load_assignments(filename):
	"""
	Load the sequence - compound assignments from the input file.
	Assumes there is a header and that the first column is the compound name
	and the second is the sequence name.
	"""
	assignment = defaultdict(list)
	with open(filename, 'r') as fh:
		data = [line.strip().split() for line in fh.readlines()]
		for line in data[1:]:
			compound_name = line[0].lower()
			ancestor = line[1]
			assignment[ancestor].append(compound_name)
	return assignment

def _get_compound_filenames(directory):
	"""
	Collect the compound filenames from the input directory
	"""
	full_path = os.path.abspath(directory)
	compound_filenames = dict()
	for fn in os.listdir(directory):
		if fn.endswith("sdf"):
			compound_full_path =  full_path + "/" + fn
			name = fn.split(".")[0].lower()
			compound_filenames[name] = compound_full_path
	return compound_filenames

def _get_chemblIDs(compound_filename_dict):
	"""
	Read each compound file and get the ChemBL ID
	"""
	compound_ids = dict()
	for compound_name, compound_fn in compound_filename_dict.items():
		with open(compound_fn, 'r') as fh:
			chemID = fh.readline().strip()
			compound_ids[chemID] = compound_name.lower()
	return compound_ids

def _process_raw_distances(raw_distances_filename, chembl_id_dict, sequence_assignments_dict):
	"""
	Process the raw output file containing tanimoto distances.
	Returns a tuple - ordered list of compound names, distance matrix and result table.

	raw_distances_filename = filename of the raw tanimoto distances from babel
	chembl_id_dict = the previously constructed chembl id mapping to sequence name
	sequence_assignments_dict = the 
	"""
	chembl_id_matches_reversed = dict([(v,k) for k,v in chembl_id_dict.iteritems()])
	# Compound1 : Compound2 : Distance
	distance_dict = defaultdict(lambda : defaultdict(float))
	with open(raw_distances_filename, 'r') as fh:
		data = [line.strip().split() for line in fh.readlines()]
		for line in data:
			if not line: continue
			if line[0].startswith("NAME"):
				cur_name = line[0].split("=")[-1]
				continue
			if len(line) > 4:
				target_id = line[0][1:]
				reference_id = line[3]
				tanimoto = line[-1] # distance
				target_name = chembl_id_dict[target_id]
				reference_name = chembl_id_dict[reference_id]
				distance_dict[reference_name][target_name] = distance_dict[target_name][reference_name] = float(tanimoto)

	compound_names = distance_dict.keys()
	result_table_strs = []
	distance_matrix = []
	for reference_name in compound_names:
		row_dists = []
		for target_name in compound_names:
			distance = distance_dict[reference_name][target_name]
			row_dists.append(distance)
			for seq_name, binding_compounds in sequence_assignments_dict.iteritems():
				binds_to_both = False
				if reference_name in binding_compounds and target_name in binding_compounds:
					binds_to_both = True
				result_table_strs.append("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (chembl_id_matches_reversed[target_name], chembl_id_matches_reversed[reference_name],target_name, reference_name, distance, seq_name, binds_to_both))
		distance_matrix.append(row_dists)
	return compound_names, distance_matrix, result_table_strs


def _process_arguments(arguments):
	sequence_assignments = _load_assignments(arguments.assignment)
	compound_filenames = _get_compound_filenames(arguments.compounds)
	
	full_input_path = os.path.abspath(arguments.compounds)
	full_output_path = os.path.abspath(arguments.output)
	
	# Combine all structures into a single sdf file. This is done for simple calculation of distances.
	subprocess.call("babel %s/*.sdf %s/all_input_compounds.sdf" % 
		(full_input_path, full_output_path), shell = True)

	# Collect the chembl identifiers
	chembl_id_matches = _get_chemblIDs(compound_filenames)
	# chembl_id_matches_reversed = dict([(v,k) for k,v in chembl_id_matches.iteritems()])
	
	# Compare each compound to every other
	subprocess.call("echo > %s/raw_distances.txt" % full_output_path, shell = True)
	for compound_name, compound_fn in compound_filenames.iteritems():
		subprocess.call("echo NAME=%s >> %s/raw_distances.txt" % (compound_name, full_output_path), shell = True)
		subprocess.call("babel %s %s/all_input_compounds.sdf -ofpt -xfMACCS >>  %s/raw_distances.txt" % (compound_fn, full_output_path, full_output_path), shell = True)

	ordered_compound_names, distance_matrix, result_table_strs = _process_raw_distances("%s/raw_distances.txt" % full_output_path, chembl_id_matches, sequence_assignments)

	# Write distance matrix
	with open(full_output_path + "/distance_matrix.txt", 'w') as fh:
		print >> fh, "\t" +  "\t".join(ordered_compound_names)
		cnt = 0
		for row in distance_matrix:
			print >> fh, ordered_compound_names[cnt] + "\t" + "\t".join(map(str, row))
			cnt += 1

	# Write result table
	with open(full_output_path + "/results_table.txt", 'w') as fh:
		print >> fh, "Compound1_ID\tCompound2_ID\tCompound1_Name\tCompound2_Name\tTanimoto\tSeq\tBindsToBoth"
		for line in result_table_strs:
			print >> fh, line

	# If database provided, collect top N matching compounds for each input compound and write to disk
	if arguments.database:
		print arguments.database
		subprocess.call("mkdir %s/top_compounds" % full_output_path, shell = True)
		for compound_name, compound_fn in compound_filenames.items():
			subprocess.call('babel %s %s/top_compounds/%s_top100_compounds.sdf -s %s -at100 -xfMACCS' % (arguments.database, full_output_path, compound_name, compound_fn), shell = True)
			distances_out_str = subprocess.Popen('babel %s %s/top_compounds/%s_top100_compounds.sdf -ofpt -xfMACCS' % 
				(compound_fn, full_output_path, compound_name), shell = True, stdout=subprocess.PIPE)
			lines = [line for line in distances_out_str.stdout.read().split("\n")]
			processed_distance_strs = []
			for line in lines:
				if not line: continue
				splitted_line = line.split()
				if len(splitted_line) < 6 or not splitted_line[0].startswith(">"): continue
				target_name = splitted_line[0][1:]
				reference_name = splitted_line[3]
				tanimoto = splitted_line[-1]
				processed_distance_strs.append("\t".join([target_name, reference_name, tanimoto]))
			with open("%s/top_compounds/%s_top100_distances.txt" % (full_output_path, compound_name), 'w') as fh:
				for line in processed_distance_strs:
					print >> fh, line

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="")

	parser.add_argument('-c', '--compounds', help='Compound directory', required=True)
	parser.add_argument('-a', '--assignment', help="Sequence assignments. Tab-delimited file, one column for sequence name and another for the compound that binds to it. Assumes column header (Compound, Sequence)", required=True)
	parser.add_argument('-o', '--output', help='Output directory', required=True)
	parser.add_argument('-d', '--database', help="ChemBL database in which to search against", required=False)
	parser.add_argument('-n', '--top_n', help="Number of top scoring compound to return when scanning against database", required=False, type = int, default = 100)

	compound_dir = "/Users/julianzaugg/Documents/University/Phd/Projects/Evolutionary_Pathway/Data/ChEMBL/raine_compounds"
	output_dir = "/Users/julianzaugg/Documents/University/Phd/Projects/Evolutionary_Pathway/Data/ChEMBL/output"
	assignment_file = "/Users/julianzaugg/Documents/University/Phd/Projects/Evolutionary_Pathway/Data/ChEMBL/raine_compounds/compound_ancestor_binding_assignment.txt"
	chembl_database = "/Users/julianzaugg/Documents/University/Phd/Projects/Evolutionary_Pathway/Data/ChEMBL/chembl_21.fs"
	my_args = ["-c", compound_dir, "-a", assignment_file, "-o", output_dir, '-d', chembl_database, '-n', "100"]
	# Compound	Sequence
	args = parser.parse_args(my_args)


	_process_arguments(args)

