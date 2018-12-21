#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import os
import pathlib  # for making the folder where we store data
import platform  # for determining OS and therefore where to store data
import re
import sys
import tracemalloc
from shutil import copyfile  # for save hdf5 to work

import magic  # to detect filetype from file header
import numpy as np  # to handle arrays and NaN
import pandas as pd  # to handle dataframes and hdf5 files
import requests  # handle API requests

# to handle compressed VCFs
import lzma
import bz2
import gzip

import matplotlib  # to plot graphs
matplotlib.use("Qt5Agg")  # to make matplotlib behave nicely with PyQT5

# to read selected lines of files (reduce RAM usage for big files)
from itertools import islice

# to build graphical interface
from PyQt5 import QtCore, QtWidgets, uic
from PyQt5.QtGui import QDesktopServices
from PyQt5.QtWidgets import QApplication, QMessageBox, QDesktopWidget

# for plotting graphs
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
plt.style.use('seaborn')

import time

start_time = time.time()

# allow <Ctrl-c> to terminate the GUI
import signal

signal.signal(signal.SIGINT, signal.SIG_DFL)

# Determine where to store temporary Metallaxis Data
home_dir = os.path.expanduser('~')
if platform.system() == "Linux":
	working_directory = home_dir + "/.metallaxis/"
elif platform.system() == "Darwin":
	working_directory = home_dir + "/Library/Caches/Metallaxis/"
elif platform.system() == "Windows":
	working_directory = os.path.expandvars(r'%APPDATA%\Metallaxis\\')

# make data folder and parent folders if they don't exist
pathlib.Path(working_directory).mkdir(parents=True, exist_ok=True)

# Interface XML files
current_file_path = __file__
current_file_dir = os.path.dirname(current_file_path)
MetGUIui = os.path.join(current_file_dir, "MetallaxisGui.ui")
MetSETui = os.path.join(current_file_dir, "MetallaxisSettings.ui")
MetPROGui = os.path.join(current_file_dir, "MetallaxisProgress.ui")


# Temporary file names
global h5_output_name, annotated_h5_output_name, vcf_output_filename
h5_output_name = os.path.join(working_directory, 'input.h5')
annotated_h5_tmp_name = os.path.join(
	working_directory, 'input_annotated_tmp.h5')
annotated_h5_output_name = os.path.join(
	working_directory, 'input_annotated.h5')
vcf_output_filename = os.path.join(
	working_directory, 'vcf_output_filename.vcf')


def throw_warning_message(warning_message):
	"""
	Displays a warning dialog with a message. Accepts one string
	as argument for the message to display.
	"""
	print("Warning: " + warning_message)
	warning_dialog = QtWidgets.QMessageBox()
	warning_dialog.setIcon(QMessageBox.Warning)
	warning_dialog.setWindowTitle("Warning:")
	warning_dialog.setText(warning_message)
	warning_dialog.setStandardButtons(QMessageBox.Ok)
	# the exec means that it won't allow interaction with the GUI until the user presses "OK"
	warning_dialog.exec_()


def throw_error_message(error_message):
	"""
	Displays an error dialog with a message. Accepts one string
	as argument for the message to display.
	"""
	print("Error: " + error_message)
	error_dialog = QtWidgets.QMessageBox()
	error_dialog.setIcon(QMessageBox.Critical)
	error_dialog.setWindowTitle("Error!")
	error_dialog.setText(error_message)
	error_dialog.setStandardButtons(QMessageBox.Ok)
	error_dialog.exec_()


def decompress_vcf(type_of_compression, vcf_input_filename, headonly_bool=False, vcf_output_filename=None):
	"""
	Decompresses or not, a file in argument (accepts xz/gz/bz2), and returns
	either the head of the decompressed file, or the filename of the decompressed
	VCF depending on the provided boolean argument: "headonly_bool".
	"""
	if type_of_compression == "":
		decompressed_file_object = open(vcf_input_filename, mode="rb")
	else:
		decompressed_file_object = eval(type_of_compression).open(vcf_input_filename, mode="rb")

	if headonly_bool is True:
		decompressed_file_head = list(islice(decompressed_file_object, 100))
		decompressed_file_object.close()
		return decompressed_file_head
	else:
		with open(vcf_output_filename, "wb") as decompressed_out:
			decompressed_out.write(decompressed_file_object.read())
		decompressed_file_object.close()
		return vcf_output_filename


def set_col_to_numeric_if_isdigit(column, chunk, numeric_columns_list):
	"""
	Determines which columns of a given chunk are ints or floats, and removes
	them from a list of columns if they are neither.
	"""

	def is_number_bool(sample):
		try:
			float(sample)
		except:
			return False
		return True

	def del_col(column):
		if column in numeric_columns_list:
			numeric_columns_list.remove(column)

	for row in chunk[column]:
		row = str(row)
		if "," in row:
			del_col(column)
		if ";" in row:
			del_col(column)
		if "|" in row:
			del_col(column)
		if bool(re.match('^[0-9]+$', row)) is False:
			if is_number_bool(row) is False:
				if column in numeric_columns_list:
					numeric_columns_list.remove(column)


def load_hdf5(hdf5_filename):
	"""
	Loads a previously created .h5 file. If an analysis had already been done on a VCF
	and the analysis file was saved in hdf5 format, it can be loaded here.
	Accepts a h5 file as argument.
	"""
	if os.path.isfile(hdf5_filename):
		# global h5_input_name
		# h5_input_name = sys.argv[1]

		MetallaxisGui.loaded_vcf_lineedit.setText(os.path.abspath(hdf5_filename))
		try:
			complete_h5_file = pd.read_hdf(hdf5_filename)
		except ValueError:
			complete_h5_file = pd.read_hdf(hdf5_filename, key="df")
		return complete_h5_file

	else:
		# return error if h5 doesn't exist
		throw_error_message("Selected file does not \
		exist. You specified : " + str(hdf5_filename))


def verify_file(selected_vcf):
	"""
	Verify that given VCF is a valid file, that exists, and has a non-null filesize.
	Accepts a VCF as entry, returns a boolean: true if file exists, false if not.
	"""
	# Verify that the file exists
	if not os.path.isfile(selected_vcf):
		throw_error_message("Selected file does not \
	exist. You specified : " + str(selected_vcf))
		return False

	# Verify that the file isn't empty
	if not os.path.getsize(selected_vcf) > 0:
		throw_error_message("Selected file is empty. \
	You specified : " + str(selected_vcf))
		return False

	# Return True to continue, as VCF matched neither of our tests
	return True


def verify_vcf(decompressed_file_head):
	# verify is conform to VCFv4.1 specification:
	# The header line names the 8 fixed, mandatory columns. These columns are as follows:
	# CHROM
	# POS
	# - (Integer, Required)
	# ID
	# - No identifier should be present in more than one data record
	# REF
	# - must be one of A,C,G,T,N (case insensitive). Multiple bases are permitted
	# ALT
	# - must be one of A,C,G,T,N (case insensitive). Multiple bases are permitted
	# - or an angle-bracketed ID String (“<ID>”)
	# - or a breakend replacement string as described in the section on
	# breakends.
	# - If there are no alternative alleles, then the missing value should be used.
	# QUAL
	# - float or Integer
	# FILTER
	# INFO

	line_num = 0
	variant_num = 0
	global metadata_num  # make global as this will be used in h5_encode
	# make a list of mandatory columns in VCF we need to find for a VCF to be valid
	expected_columns = ['#CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']

	for line in decompressed_file_head:
		line_num = line_num + 1
		if line.startswith(b'#'):  # lines with # are metadata
			if line.startswith(b'#CHROM'):
				# verify header is conform to vcf 4.1 spec
				# decode byte object to utf-8 string, then split by tab to get each column
				header_line_cols = line.decode('UTF-8')
				header_line_cols = header_line_cols.rstrip()  # remove trailing/leading whitespace
				header_line_cols = header_line_cols.split("\t")
				# get index of each column name
				global chrom_col, id_col, pos_col, ref_col, alt_col, qual_col
				chrom_col = [i for i, s in enumerate(header_line_cols) if '#CHROM' in s][0]
				id_col = [i for i, s in enumerate(header_line_cols) if 'ID' in s][0]
				pos_col = [i for i, s in enumerate(header_line_cols) if 'POS' in s][0]
				ref_col = [i for i, s in enumerate(header_line_cols) if 'REF' in s][0]
				alt_col = [i for i, s in enumerate(header_line_cols) if 'ALT' in s][0]
				qual_col = [i for i, s in enumerate(header_line_cols) if 'QUAL' in s][0]

				# verify that VCF has all required columns
				if not all(column in header_line_cols for column in expected_columns):
					throw_error_message("VCF not valid: VCF doesn't not contain all required columns. Contains: " + str(
						header_line_cols))
					return False

		else:
			# if line is not part of header, split by tab to get columns
			# and verify that each column, in each line is conform to VCFv4.1
			split_line = line.decode('UTF-8').split("\t")
			for column in split_line:
				column = column.strip()

				# Verify that POS column only contains digits
				if column == split_line[pos_col]:
					if not column.isdigit():
						throw_error_message("VCF not valid: column 'POS' doesn't only contain digits: " + str(column))
						return False

				# Verify that REF column only contains ACGTN
				elif column == split_line[ref_col]:
					allowed_chars = set('ACGTN')
					if not set(column.upper()).issubset(allowed_chars):
						throw_error_message(
							"VCF not valid: column 'REF' doesn't only contain A,C,G,T,N: " + str(column))
						return False

				# Verify that ALT column only contains ACGTN or <ID>
				elif column == split_line[alt_col]:
					if not column.startswith("<") and column.endswith(">"):
						allowed_chars = set('ACGTN')
						if not set(column).issubset(allowed_chars):
							throw_error_message(
								"VCF not valid: column 'ALT' doesn't only contain A,C,G,T,N or <ID>: " + str(column))
							return False

				# Verify that QUAL column only contains number (float or int) or "."
				elif column == split_line[qual_col]:
					if column.isdigit():
						break
					elif column == ".":
						break
					else:
						allowed_chars = set('123456789.')
						if not set(column).issubset(allowed_chars):
							try:
								float(column)

							except ValueError:
								throw_error_message(
									"VCF not valid: column 'QUAL' doesn't only contain digits: " + str(column))
								return False
			variant_num += 1

		# set global variable to be total number of lines of metadata (will be used to extract just the data in h5_encode()
		metadata_num = int(line_num - variant_num)

	# return Error or warning depending on variant numbers
	if variant_num == 0:
		throw_error_message("VCF is empty, there are no variants at all in this vcf, please use a different vcf")
		return False
	elif variant_num < 5:
		throw_error_message("VCF contains too few variants to analyse, please use a different vcf")
		return False
	elif 5 < variant_num < 30:
		throw_warning_message("VCF contains very few variants, only rudimentary statistics can be performed")
		return True
	else:
		# if more than 35 variants then VCF is fine, return without any alert
		return True


def parse_vcf(vcf_input_filename):
	"""
	Takes a VCF in input, runs both file and VCF verifications, decompresses it, and extracts metadata  and statistics
	Accepts a string of a VCF filename in input.
	Returns a tuple of a metadata dictionary, a statistics dictionary, and the filename of the decompressed VCF.
	"""
	# Verify that the selected VCF is a valid file (ie. that it exists, and has a non-null size)
	MetallaxisGui.progress_bar(3, "Verifying VCF: verifying that file is valid")
	file_is_valid = verify_file(vcf_input_filename)
	if not file_is_valid:
		return

	vcf_filetype = magic.from_file(vcf_input_filename)
	# Decompress selected vcf
	if "XZ" in vcf_filetype:
		MetallaxisGui.detected_filetype_label.setText("xz compressed VCF")
		decompressed_file_head = decompress_vcf("lzma", vcf_input_filename, headonly_bool=True)

	elif "bzip2" in vcf_filetype:
		MetallaxisGui.detected_filetype_label.setText("bz2 compressed VCF")
		decompressed_file_head = decompress_vcf("bz2", vcf_input_filename, headonly_bool=True)

	elif "gzip" in vcf_filetype:
		MetallaxisGui.detected_filetype_label.setText("gz compressed VCF")
		decompressed_file_head = decompress_vcf("gzip", vcf_input_filename, headonly_bool=True)

	elif "Variant Call Format" in vcf_filetype:
		MetallaxisGui.detected_filetype_label.setText("uncompressed VCF")
		decompressed_file_head = decompress_vcf("", vcf_input_filename, headonly_bool=True)

	else:
		throw_error_message("Selected file must be a VCF file")
		return

	# now we have a returned decompressed file object verify if
	# contents are valid vcf
	MetallaxisGui.progress_bar(8, "Verifying VCF: verifying that file is a valid VCF")
	vcf_is_valid = verify_vcf(decompressed_file_head)
	if not vcf_is_valid:
		return

	MetallaxisGui.progress_bar(9, "Decompressing VCF")
	if "XZ" in vcf_filetype:
		decompressed_file = decompress_vcf("lzma", vcf_input_filename, vcf_output_filename=vcf_output_filename)

	elif "bzip2" in vcf_filetype:
		decompressed_file = decompress_vcf("bz2", vcf_input_filename, vcf_output_filename=vcf_output_filename)

	elif "gzip" in vcf_filetype:
		decompressed_file = decompress_vcf("gzip", vcf_input_filename, vcf_output_filename=vcf_output_filename)

	elif "Variant Call Format" in vcf_filetype:
		decompressed_file = decompress_vcf("", vcf_input_filename, vcf_output_filename=vcf_output_filename)

	MetallaxisGui.loaded_vcf_lineedit.setText(os.path.abspath(vcf_input_filename))

	# Calculate counts of different variant types
	def add_to_dict_iterator(dictionary, key, iterator_value):
		"""
		iterates a key in a dictionary. for a given dictionary name and key name it will add iterator_value
		to the value of the key, unless it doesn't exist in which case it will create it.
		"""
		if key not in dictionary:
			# add empty list to dict as key
			dictionary[key] = 0
			dictionary[key] = dictionary[key] + iterator_value
		else:
			# append value to list
			dictionary[key] = dictionary[key] + iterator_value

	variant_stats = {}
	length_of_all_indels = 0
	variant_stats["Total_SNP_Count"] = 0
	variant_stats["Total_Indel_Count"] = 0
	global list_chromosomes, ALT_Types
	list_chromosomes = set()  # use set instead of list so it won't store duplicate values
	ALT_Types = set()  # use set instead of list so it won't store duplicate values

	# If VCF only has SNP then collect data on distribution of different nucleotides
	alt_types_only_snp = True
	with open(vcf_output_filename) as decompressed_out:
		for line in decompressed_out:
			if not line.startswith('#'):
				line = line.split("\t")
				alt = str(line[alt_col])
				if set(alt).issubset(set('ACTG')) and len(alt) != 1:
					alt_types_only_snp = False

	with open(vcf_output_filename) as decompressed_out:
		for line in decompressed_out:
			if not line.startswith('#'):
				line = line.split("\t")
				list_chromosomes.add(line[chrom_col])
				alt = str(line[alt_col])
				# if our alt types are only SNP or transposable elements
				# then count them
				if alt_types_only_snp is True:
					ALT_Types.add(alt)
					add_to_dict_iterator(
						variant_stats, line[alt_col] + "_Alt_Count", 1)

				if len(line[ref_col]) == len(line[alt_col]):
					variant_stats["Total_SNP_Count"] += 1
					add_to_dict_iterator(
						variant_stats, line[chrom_col] + "_Chrom_SNP_Count", 1)
					add_to_dict_iterator(
						variant_stats, line[chrom_col] + "_Chrom_Variant_Count", 1)
				else:
					variant_stats["Total_Indel_Count"] += 1
					add_to_dict_iterator(
						variant_stats, line[chrom_col] + "_Chrom_Indel_Count", 1)
					add_to_dict_iterator(
						variant_stats, line[chrom_col] + "_Chrom_Variant_Count", 1)
					length_of_all_indels += len(line[alt_col])

	total_chrom_snp_count, total_chrom_indel_count = 0, 0
	for key, value in variant_stats.items():
		if "_Chrom_SNP_Count" in key:
			total_chrom_snp_count += value
		if "_Chrom_Indel_Count" in key:
			total_chrom_indel_count += value

	if length_of_all_indels > 0 and variant_stats["Total_Indel_Count"] > 0:
		variant_stats["Avg_Indel_Length"] = float(
			length_of_all_indels / variant_stats["Total_Indel_Count"])
		variant_stats["Avg_Indel_Length"] = round(variant_stats["Avg_Indel_Length"], 3)
	variant_stats["Avg_SNP_per_Chrom"] = int(
		total_chrom_snp_count / len(list_chromosomes))
	variant_stats["Avg_Indel_per_Chrom"] = int(
		total_chrom_indel_count / len(list_chromosomes))
	variant_stats["Avg_Variant_per_Chrom"] = int(
		(total_chrom_snp_count + total_chrom_indel_count) / len(list_chromosomes))

	variant_stats["List_Chromosomes"] = list_chromosomes

	# if alt_types isn't empty then add to statistics dictionary
	if ALT_Types != set():
		variant_stats["ALT_Types"] = ALT_Types

	# Extract Metadata from VCF
	metadata_dict = {}
	# Match groups either side of an "=", after a "##". e.g. filename=xyz, source=tangram, etc.
	regex_metadata = re.compile('(?<=##)(.*?)=(.*$)')

	MetallaxisGui.progress_bar(10, "Extracting VCF metadata")
	with open(vcf_output_filename) as decompressed_out:
		vcf_line_nb, metadata_line_nb = 0, 0

		for line in decompressed_out:
			if line.startswith('##'):
				metadata_tag = str(regex_metadata.search(line).group(1))
				metadata_result = str(regex_metadata.search(line).group(2))
				# classify uppercase metadata (e.g. "INFO", "FILTER") differently
				if metadata_tag.isupper():
					metadata_type = metadata_tag
				else:
					metadata_type = "basic"
					# truncate long metadata to avoid database errors and distorting the GUI
					if len(metadata_tag) > 20:
						metadata_tag = metadata_tag[:20] + "..."
					if len(metadata_result) > 95:
						metadata_result = metadata_result[:95] + "...<truncated due to length>"

				metadata_dict_entry = [metadata_type, metadata_tag, metadata_result]
				if not metadata_dict_entry in metadata_dict.values():
					metadata_dict[metadata_line_nb] = metadata_dict_entry
				metadata_line_nb += 1

	return metadata_dict, variant_stats, decompressed_file


def h5_encode(selected_vcf, decompressed_file, variant_stats, metadata_dict):
	"""
	accepts as input a raw vcf file, and optionally, dictionaries for variant counts (variant_stats) and metadata
	(metadata_dict). Then encodes them as tables into a HDF5 file. Returns filename of created HDF5.
	"""
	h5_file = pd.HDFStore(h5_output_name, mode='w')

	# write each entry from metadata_dict to a new table in a HDF5 file
	for metadata_line_nb in metadata_dict:
		metadata_tag = str(metadata_dict[metadata_line_nb][1])
		metadata_result = str(metadata_dict[metadata_line_nb][2])
		if not metadata_tag.isupper():
			metadata_line = {'Tag': metadata_tag, 'Result': metadata_result}
			metadata_line = pd.DataFrame(
				metadata_line, index=[metadata_line_nb])
			h5_file.append("metadata", metadata_line, data_columns=True,
			               min_itemsize=150, complib=complib, complevel=int(complevel))

	# write each entry from  to a new table in a HDF5 file
	h5_stat_table_index = 0
	for key, value in variant_stats.items():
		key = str(key)
		value = str(value)
		if len(key) > 40:
			key = key[:40] + "..."
		if len(value) > 200:
			value = value[:200] + "..."
		variant_stats_line = {'Tag': key, 'Result': value}
		variant_stats_line = pd.DataFrame(
			variant_stats_line, index=[h5_stat_table_index])

		h5_file.append("stats", variant_stats_line, data_columns=True,
		               min_itemsize=250, complib=complib, complevel=int(complevel))
		h5_stat_table_index += 1

	chunked_vcf_len = sum(1 for row in open(decompressed_file, 'r'))

	annotate_nb = 0
	annotate_percent = 35

	# PARSE INFO COLUMNS
	# the info column of a vcf is long and hard to read if
	# displayed as is, but it is composed of multiple key:value tags
	# that we can parse as new columns, making them filterable.
	global info_cols_to_add
	info_cols_to_add = set()
	chunked_vcf = pd.read_csv(selected_vcf,
	                          delim_whitespace=True,
	                          # skip rows that started with "#"
	                          skiprows=range(0, metadata_num - 1),
	                          # use chunk size that was set in settings
	                          chunksize=int(h5chunksize),
	                          low_memory=False,
	                          # make default data type an object (ie. string)
	                          dtype=object)

	# Run through every chunk of the whole VCF and get the key out of the
	# key:value pair and add to a set (so no duplicates will be added) that
	# will later become new column names
	for chunk in chunked_vcf:
		for line in chunk["INFO"]:
			if ";" in line:
				line_split = line.split(";")
				for col in line_split:
					col_to_add = col.split("=")[0]
					info_cols_to_add.add(col_to_add)

	global table_column_names
	chunked_vcf = pd.read_csv(selected_vcf,
	                          delim_whitespace=True,
	                          skiprows=range(0, metadata_num - 1),
	                          chunksize=int(h5chunksize),
	                          low_memory=False,
	                          # make default data type an object (ie. string)
	                          dtype=object)

	for chunk in chunked_vcf:
		# set the new info column names to be empty by default
		for col in info_cols_to_add:
			chunk[col] = "."

		line_nb = 0
		for line in chunk["INFO"]:
			# split the INFO column by ; which separates the different
			# key:value pairs
			if ";" in line:
				line_split = line.split(";")
				# get both sides of the = to get the key and the value
				for col in line_split:
					col_split = col.split("=")
					key_to_add = col_split[0]
					# col_split will be greater than 1 if there is an = sign
					# a = means there is a key:value pair to be extracted
					if len(col_split) > 1:
						data_to_add = col_split[1]
						chunk[key_to_add].values[line_nb] = data_to_add
					# if there is no = sign then there is no key:value pair just
					# a tag so set it tag as a boolean column
					else:
						chunk[key_to_add].values[line_nb] = "True"
			line_nb += 1

		# Rename column so we get 'CHROM' not '#CHROM' from chunk.keys()
		chunk.rename(columns={'#CHROM': 'CHROM'}, inplace=True)

		# set variable table_column_names to have both normal VCF
		# columns and the parsed columns from INFO column
		table_column_names = list(chunk.keys()) + list(info_cols_to_add)

		# SET COLS WITH NUMBERS TO BE NUMERIC TYPE
		# set columns that only contain numbers to be numeric dtype
		# otherwise they are just strings, and can't be used with a
		# dash separated filter

		# make a list from all column names, we will later remove the
		# non-numeric columns from the list
		global numeric_columns
		numeric_columns = table_column_names

		chunk = chunk.replace(".", np.NaN)

		# run the function for every column to remove non-numeric columns
		for column in chunk.keys():
			set_col_to_numeric_if_isdigit(column, chunk, numeric_columns)

		# convert the remaining columns in "numeric_columns" list to numeric datatype
		for column in numeric_columns:
			chunk[column] = pd.to_numeric(chunk[column])

		# only update progress bar every 20 lines to avoid performance hit
		# of refreshing the whole GUI at every line
		if annotate_nb % 20 == 0:
			annotate_progress = annotate_percent + \
			                    (annotate_nb / chunked_vcf_len) * (43 - annotate_percent)
			MetallaxisGui.progress_bar(
				float(annotate_progress), "Encoding H5 database")

		# append the modified chunk to the h5 database without indexing
		h5_file.append("df", chunk, index=False, data_columns=True,
		               min_itemsize=80, complib=complib, complevel=int(complevel))
		annotate_nb += 1

	# index columns explicitly, now that we have finished adding data to it
	MetallaxisGui.progress_bar(46, "Indexing H5 database")
	h5_file.create_table_index(
		"df", columns=table_column_names, optlevel=9, kind='full')
	h5_file.close()
	return h5_output_name


def annotate_h5(h5_output_name):
	"""
	Uses EBI's VEP API to annotate lines where ID is valid. This is still highly experimental,
	it only works about 50% of the time and on small VCF files.

	Accepts as input, filename of HDF5 store created by h5_encode.
	Returns filename of annotated HDF5 store, if completes without errors.
	"""
	ebi_rest_api = "https://rest.ensembl.org"
	ext = "/vep/human/id/"
	headers = {"Content-Type": "application/json", "Accept": "application/json"}

	all_hdf5_ids = []
	hdf5_input_read = pd.read_hdf(h5_output_name, key="df")
	hdf5_input_read = hdf5_input_read[hdf5_input_read.ID.notnull()]

	if len(hdf5_input_read) == 0:
		throw_warning_message("No ids in VCF, aborting annotation.")
		raise IndexError  # raise error to quit annotation and use non-annotated h5

	for line in hdf5_input_read.itertuples():
		# turn tuple into a list, but exclude the first item of list
		# because its the h5 index and not part of our original data
		line = list(line)[1:]
		# get the id from id_column and add it to all_hdf5_ids list
		line_id = line[id_col]
		# ids are strings so if we ave a valid non-null string add it as an ID to send to EBI
		if type(line_id) == str:
			if line_id != "." and line_id != "nan":
				all_hdf5_ids.append(line_id)  # list of all ids in hdf5 input

	if len(all_hdf5_ids) == 0:
		# return an error if we have no IDs so that metallaxis can use the non-annotated h5 instead
		throw_error_message("No valid ids to annotate in vcf")
		raise IndexError

	annotated_h5_tmp = pd.HDFStore(annotated_h5_tmp_name, mode='w')
	annotated_row_ids = []

	annotate_nb = 0
	# divide our all_hdf5_ids list, into multiple sublists, and call API with sublists to avoid
	# reaching the API's limit
	api_call_number = 1
	maximum_apis_to_send = 50
	for i in range(0, len(all_hdf5_ids), maximum_apis_to_send):
		# update progress bar dynamically based off current API call number out of total number of API calls to make
		MetallaxisGui.progress_bar(51, "Annotate H5: Making API Call to EBI (" + str(api_call_number) +
		                           "/" + str(int(round((len(all_hdf5_ids) / maximum_apis_to_send), 0) + 1)) + ")")
		api_ids_sublist = all_hdf5_ids[i:i + maximum_apis_to_send]

		formatted_api_ids_sublist = ','.join('"{0}"'.format(ind_id) for ind_id in api_ids_sublist)
		formatted_api_ids_sublist = '{ "ids" : [' + formatted_api_ids_sublist + '] }'
		api_call = requests.post(ebi_rest_api + ext, headers=headers, data=formatted_api_ids_sublist)

		if api_call.ok:
			data = api_call.text
			# save sublist's API response as a JSON
			json_output = os.path.join(working_directory, 'VEP_API_' + str(api_call_number) + '.json')
			with open(json_output, 'w') as json_file:
				json_file.write(data)

			api_call_number += 1

		else:
			if api_call_number < 2:  # if not even one of the API calls worked
				throw_error_message("API call failed")
				api_call.raise_for_status()  # raises error to quit this try loop so Metallaxis can use non-annotated h5
			else:
				# at least one the API calls worked, so continue using that annotation
				throw_warning_message("The API call failed, only doing partial annotation")

	# specify the columns we want to look for in the VEP JSONs
	columns_to_annotate = ['impact', 'consequence_terms', 'gene_id', 'gene_symbol',
	                       'biotype', 'distance', 'gene_symbol_source', 'transcript_id',
	                       'cdna_start', 'cdna_end', 'most_severe_consequence']
	# calculate what annotation_columns will be made so they can all be set to default
	# values even if one of the json makes no mention of such a column. we defined
	# columns_to_annotate but if there are multiple values for the defined columns
	# then new columns will be generated which will have to be filled with default values
	# this avoids there being a "invalid combinate" error on appending the data to the h5

	# Loop through all the JSONs, adding columns that match our columns_to_annotate list to our annotation_columns set. We do this because sometimes there are multiple values for each specified column and we want to know about all of them before we extract that information
	MetallaxisGui.progress_bar(53, "Annotate H5: Parsing downloaded JSONs for annotations")
	annotation_columns = set()
	# loop through all the saved JSONs
	for json_number in range(1, api_call_number):
		json_file = os.path.join(working_directory, 'VEP_API_' + str(json_number) + '.json')

		# open the json file as a JSON object
		with open(json_file) as json_input:
			json_raw_data = json.load(json_input)

		# loop through all the dictionaries embedded in the JSON object
		for data in json_raw_data:
			# if returned JSON has a transcript_consequence extract its
			# data and add to h5 file
			if 'transcript_consequences' in data:
				for subdata in data['transcript_consequences']:
					for column in columns_to_annotate:
						if column in subdata.keys():
							if isinstance(subdata[column], list) and len(subdata[column]) > 1:
								for each_col in range(1, len(subdata[column])):
									annotation_columns.add((column + "_" + str(each_col)).upper())
							else:
								annotation_columns.add(column.upper())
						else:
							annotation_columns.add(column.upper())

	MetallaxisGui.progress_bar(54, "Annotate H5: writing annotated chunk to temporary h5")
	annotate_percent = 55
	# for each json file generated from the API call, parse the previously found annotation
	# columns from the JSON to our dataframe
	for json_number in range(1, api_call_number):
		json_file = os.path.join(working_directory, 'VEP_API_' + str(json_number) + '.json')

		# open the json file as a JSON object
		with open(json_file) as json_input:
			json_raw_data = json.load(json_input)

		# loop through all the dictionaries embedded in the JSON object
		for data in json_raw_data:
			# if returned JSON has a transcript_consequence extract its
			# data and add to h5 file
			if 'transcript_consequences' in data:
				for subdata in data['transcript_consequences']:
					# only update progress bar every 20 lines to avoid performance hit
					# from doing it every line
					if annotate_nb % 20 == 0:
						annotate_len = (len(all_hdf5_ids) * len(data))
						annotate_progress = annotate_percent + (annotate_nb / annotate_len) * (70 - annotate_percent)
						MetallaxisGui.progress_bar(float(annotate_progress),
						                           "Annotate H5: writing annotated chunk to temporary h5")

					# get h5 row from id and add JSON data to it and add it to a different h5
					annotation_id = data['id']
					annotation_alt = subdata['variant_allele']

					# get unannotated row with id that match annotated row and have same ALT
					selected_h5_row = pd.read_hdf(h5_output_name, key="df", where=(
							"ID=='" + annotation_id + "' and ALT=='" + annotation_alt + "'"))

					if not selected_h5_row.empty:
						# get annotation information if it exists else insert a "."
						def get_annotation_for_columns(columns):
							for column in columns:
								if column in subdata.keys():
									# for a given column (e.g. biotype, impact) there can be
									# a list of values instead of just one value, to find out
									# weather to add multiple columns or just set one we test
									# if its a list and has more than 1 value
									if isinstance(subdata[column], list) and len(subdata[column]) > 1:
										for each_col in range(1, len(subdata[column])):
											selected_h5_row[str(column + "_" + str(each_col)).upper()] = \
												subdata[column][each_col]
											annotation_columns.add((column + "_" + str(each_col)).upper())
									else:
										selected_h5_row[str(column).upper()] = subdata[column]
										annotation_columns.add(column.upper())

						# initialise all the found columns
						for column in annotation_columns:
							selected_h5_row[column] = "."

						# run the function to replace the default "." with the real value if it exists
						get_annotation_for_columns(columns_to_annotate)

						# Convert all columns to strings
						for column in selected_h5_row:
							selected_h5_row[column] = selected_h5_row[column].apply(str)

						# Append the row to h5 database
						annotated_h5_tmp.append("df", selected_h5_row, index=False, data_columns=True, min_itemsize=100,
						                        complib=complib, complevel=int(complevel))

						annotated_row_ids.append(annotation_id)
						annotate_nb += 1

	# read through vcf again, and load rows with no annotation-friendly IDs to a new h5 database
	MetallaxisGui.progress_bar(70, "Annotate H5: writing non-annotated chunk to temporary h5")

	# get all rows from non-annoted h5, except which have IDs in annotated_row_ids
	non_annotated_h5_row = pd.read_hdf(h5_output_name, key="df")
	non_annotated_h5_row = non_annotated_h5_row[~non_annotated_h5_row.ID.isin(annotated_row_ids)]

	if not non_annotated_h5_row.empty:
		for col in annotation_columns:
			non_annotated_h5_row[col] = "."

		for column in non_annotated_h5_row:
			non_annotated_h5_row[column] = non_annotated_h5_row[column].apply(str)

		annotated_h5_tmp.append("df", non_annotated_h5_row, index=False, data_columns=True,
		                        min_itemsize=100, complib=complib, complevel=int(complevel))

	annotated_h5_tmp.close()

	# Read the just written h5 file replacing all "." with NaN which is read as NULL instead of a string.
	# convert columns containing only numbers to numeric datatype, and remove columns with no data for any row
	annotated_numeric_columns = list(numeric_columns) + list(annotation_columns)
	print(annotated_numeric_columns)

	annotated_tmp_read_obj = pd.read_hdf(annotated_h5_tmp_name, key="df")
	annotated_tmp_read_obj = annotated_tmp_read_obj.replace(".", np.NaN)

	for column in annotated_tmp_read_obj.keys():
		set_col_to_numeric_if_isdigit(column, annotated_tmp_read_obj, annotated_numeric_columns)

	for column in annotated_numeric_columns:
		annotated_tmp_read_obj[column] = pd.to_numeric(annotated_tmp_read_obj[column])

	# remove columns that only contained null values
	annotated_tmp_read_obj = annotated_tmp_read_obj.dropna(axis=1, how='all')

	# save appended dataframe to new h5 file and remove the temporary one
	annotated_tmp_read_obj.to_hdf(annotated_h5_output_name, key='df', mode='w', format='table', data_columns=True,
	                              complib=complib, complevel=int(complevel))

	if os.path.exists(annotated_h5_tmp_name):
		os.remove(annotated_h5_tmp_name)

	# embed metadata and statistics in new H5
	h5_file = pd.HDFStore(annotated_h5_output_name)
	for metadata_line_nb in metadata_dict:
		metadata_tag = str(metadata_dict[metadata_line_nb][1])
		metadata_result = str(metadata_dict[metadata_line_nb][2])
		if not metadata_tag.isupper():
			metadata_line = {'Tag': metadata_tag, 'Result': metadata_result}
			metadata_line = pd.DataFrame(
				metadata_line, index=[metadata_line_nb])
			h5_file.append("metadata", metadata_line, data_columns=True,
			               min_itemsize=150, complib=complib, complevel=int(complevel))

	h5_stat_table_index = 0
	for var_counts_key, var_counts_value in var_counts.items():
		var_counts_key = str(var_counts_key)
		var_counts_value = str(var_counts_value)
		if len(var_counts_key) > 40:
			var_counts_key = var_counts_key[:40] + "..."
		if len(var_counts_value) > 60:
			var_counts_value = var_counts_value[:60] + "..."
		var_counts_line = {'Tag': var_counts_key, 'Result': var_counts_value}
		var_counts_line = pd.DataFrame(
			var_counts_line, index=[h5_stat_table_index])

		h5_file.append("stats", var_counts_line, data_columns=True,
		               min_itemsize=100, complib=complib, complevel=int(complevel))
		h5_stat_table_index += 1

	h5_file.create_table_index("df", columns=table_column_names, optlevel=9, kind='full')
	h5_file.close()

	return annotated_h5_output_name


# Build graphical interface constructed in XML
gui_window_object, gui_base_object = uic.loadUiType(MetGUIui)


class MetallaxisGuiClass(gui_base_object, gui_window_object):
	"""
	Class that constructs the PyQt graphical interface
	"""

	def __init__(self):
		super(gui_base_object, self).__init__()
		self.setupUi(self)
		self.setWindowTitle("Metallaxis")
		# initialise progress bar
		self.MetallaxisProgress = MetallaxisProgress()
		self.MetallaxisProgress.show()

		# Center GUI on screen
		qt_rectangle = self.frameGeometry()
		center_point = QDesktopWidget().availableGeometry().center()
		qt_rectangle.moveCenter(center_point)

		# start measuring memory
		tracemalloc.start()
		global normal_mem
		normal_mem = tracemalloc.take_snapshot()

		self.progress_bar(1, "setting up gui")
		# buttons on interface
		self.open_vcf_button.clicked.connect(self.select_and_parse)
		# menus on interface
		self.actionOpen_VCF.triggered.connect(self.select_and_parse)
		self.actionSave_as_HDF5.triggered.connect(self.save_hdf5)
		self.actionQuit.triggered.connect(self.close)

		# Link "Github Page" button on menu to its URL
		def open_github(url):
			url = "https://github.com/SL-LAIDLAW/Metallaxis"
			QDesktopServices.openUrl(QtCore.QUrl(url))

		self.actionGithub_Page.triggered.connect(open_github)

		# Link "Documentation" button on menu to its URL
		def open_docs(url):
			url = "https://metallaxis.readthedocs.io"
			QDesktopServices.openUrl(QtCore.QUrl(url))

		self.actionMetallaxis_Documentation.triggered.connect(open_docs)

		def open_about_tab():
			self.tabWidget.setCurrentIndex(3)

		# set first tab as default
		self.tabWidget.setCurrentIndex(0)
		# on click "About", open about tab
		self.actionAbout.triggered.connect(open_about_tab)

		self.MetallaxisSettings = MetallaxisSettings()

		def show_settings_window():
			self.MetallaxisSettings.show()

		self.actionSettings.triggered.connect(show_settings_window)
		self.MetallaxisSettings.annotate_species_comboBox.addItems(['Other', 'Human'])

		# convert gui settings to global variables
		global complevel, complib, h5chunksize
		complevel = self.MetallaxisSettings.compression_level_spinBox.text()
		complib = self.MetallaxisSettings.compression_comboBox.currentText()
		h5chunksize = self.MetallaxisSettings.vcf_chunk_size.text()

		# on changing Species combobox in Settings, run the changed_species_combobox
		# function that'll enable or disable the "annotation" checkbox
		self.MetallaxisSettings.annotate_species_comboBox.currentTextChanged.connect(self.changed_species_combobox)

	def progress_bar(self, percent, message):
		"""
		Function that updates the on screen progress bar. Accepts two arguments, one
		for % progression, the other for the current status to display.
		"""
		self.MetallaxisProgress.progressbar_message.setText(message)
		percent = round(percent, 2)
		self.MetallaxisProgress.progressbar_progress.setValue(percent)

		snapshot2 = tracemalloc.take_snapshot()
		mem_usage = sum(stat.size for stat in snapshot2.statistics('lineno'))
		mem_usage_mb = round(mem_usage / 1000000, 2)
		mem_usage_mb = ("%s" % (mem_usage_mb))

		time_secs = str(round((time.time() - start_time), 2))

		self.MetallaxisProgress.progressbar_ram_usage.setText(mem_usage_mb)
		self.MetallaxisProgress.progressbar_time.setText(time_secs)
		print(str(percent) + "% : " + message + " | Time: " + time_secs + " \
		| RAM (MB): " + mem_usage_mb)
		MetallaxisApp.processEvents()  # refresh the GUI

	def empty_qt_layout(self, qt_layout_name):
		"""
		Empties the Qt layout provided in argument. This allows us to load new files without
		restarting the application. Accepts name of layout as argument.
		"""
		while 1:
			layout_widget = qt_layout_name.takeAt(0)
			if not layout_widget:
				break
			layout_widget.widget().deleteLater()

	def save_hdf5(self):
		"""
		Open a dialog where the user can chose where to place hdf5 savefile
		"""
		save_dialog = QtWidgets.QFileDialog()
		save_dialog.setAcceptMode(save_dialog.AcceptSave)
		save_folder = save_dialog.getSaveFileName(self, 'Save Analayis as HDF5', filter="*.h5")[0]
		if self.MetallaxisSettings.annotation_checkbox.isChecked():
			copyfile(annotated_h5_output_name, save_folder)
		else:
			copyfile(h5_output_name, save_folder)

	def changed_species_combobox(self):
		"""
		Deactivates checkbox for VCF annotation if a species other than 'Human'
		is selected, as only human VCFs can be annotated without changing EBI
		VEP API settings.
		"""
		vcf_species = self.MetallaxisSettings.annotate_species_comboBox.currentText()
		if vcf_species != 'Human':
			self.MetallaxisSettings.annotation_checkbox.setChecked(False)
			self.MetallaxisSettings.annotation_checkbox.setEnabled(False)
			self.MetallaxisSettings.annotate_vcf_label.setEnabled(False)
		else:
			self.MetallaxisSettings.annotation_checkbox.setEnabled(True)
			self.MetallaxisSettings.annotate_vcf_label.setEnabled(True)

	def select_file(self):
		"""
		Opens a file dialog where the user can chose an input file.
		"""
		select_dialog = QtWidgets.QFileDialog()
		select_dialog.setAcceptMode(select_dialog.AcceptSave)
		selected_vcf = select_dialog.getOpenFileName(self, filter="VCF Files (*.vcf \
			*.vcf.xz *.vcf.gz *.vcf.bz2) ;;Metallaxis HDF5 Files(*.h5) ;;All Files(*.*)")
		selected_vcf = selected_vcf[0]
		# if the user cancels the select_file() dialog, then run select again
		while selected_vcf == "":
			throw_error_message("No selected file")
			return False
		return selected_vcf

	def select_and_parse(self):
		"""
		Runs select_file() to get an input and based off extension, either runs the parse_VCF() or displays the contents
		of the h5 file directly. Accepts no arguments, and returns nothing. This function exists solely to call other
		functions, as menu items in PyQt can only call one function.
		"""
		selected_file = self.select_file()
		if selected_file is False:
			return  # If the user cancels the file picker
		h5_only = False
		if selected_file.endswith(".h5"):
			complete_h5_file = load_hdf5(selected_file)
			h5_only = True
			self.write_h5_data_to_interface(complete_h5_file, selected_file)
		else:
			selected_vcf = selected_file

		# reopen progress bar for loading new file
		self.MetallaxisProgress = MetallaxisProgress()
		self.MetallaxisProgress.show()

		if not h5_only:
			# get metadata and variant counts from vcf
			metadata_dict, var_counts, decompressed_file = parse_vcf(selected_vcf)

			# convert vcf into a hdf5 object
			h5_file = h5_encode(selected_vcf, decompressed_file, var_counts, metadata_dict)

			# Read H5 for actual table populating
			complete_h5_file = pd.read_hdf(h5_file, key="df")

			# populate interface with information from hdf5
			self.write_h5_data_to_interface(complete_h5_file, h5_file)

			# get annotation data
			if self.MetallaxisSettings.annotation_checkbox.isChecked():
				try:
					complete_h5_file = annotate_h5(h5_file)
					complete_h5_file = pd.read_hdf(complete_h5_file, key="df")
				except:
					throw_warning_message("Annotation did not succeed, proceeding with non-annotated h5")
					self.MetallaxisSettings.annotation_checkbox.setChecked(False)

		# populate table
		self.populate_table(complete_h5_file)

	def filter_table(self):
		"""
		Filters table based on chosen filters. Comma separated, dash separated, and single filters exist. This function
		reads the chosen filter from the interface, requests rows matching those filters from the HDF5 store and
		populates table with that data.
		Accepts no arguments (retrieves all information from interface), and returns no value.
		"""
		selected_filter = self.filter_box.currentText()
		filter_text = self.filter_lineedit.text()
		# remove leading / trailing whitespace from request
		filter_text = re.sub(r"\s+", "", filter_text)
		filter_text = filter_text.upper()

		if "-" in filter_text and "," in filter_text:
			throw_error_message("Please only use either comma separated values or a dash separated range")
			return

		elif "-" in filter_text:
			if selected_filter not in numeric_columns:
				throw_error_message(
					"Can only filter a dash separated range on numeric columns: " + str(numeric_columns))
				return

			split_filter_text = filter_text.split("-")
			# Filter out Null values to avoid corrupting item count
			split_filter_text = filter(None, split_filter_text)
			split_filter_text = list(split_filter_text)
			if len(split_filter_text) == 2:
				self.filter_text.setText(
					"Filtering to show " + selected_filter + " from " + str(split_filter_text[0]) + " to " + str(
						split_filter_text[1]))

				if split_filter_text[0] > split_filter_text[1]:
					filter_condition = selected_filter + ">=" + split_filter_text[1] + " & " + selected_filter + "<=" + \
					                   split_filter_text[0]
				elif split_filter_text[0] < split_filter_text[1]:
					filter_condition = selected_filter + ">=" + split_filter_text[0] + " & " + selected_filter + "<=" + \
					                   split_filter_text[1]
				else:
					filter_condition = selected_filter + "==" + split_filter_text[0]

				if self.MetallaxisSettings.annotation_checkbox.isChecked():
					filtered_h5_table = pd.read_hdf(annotated_h5_output_name, key="df", where=filter_condition)
				else:
					filtered_h5_table = pd.read_hdf(h5_output_name, key="df", where=filter_condition)
			else:
				throw_error_message("Please only enter 2 values separated by a dash")
				return

		elif "," in filter_text:
			split_filter_text = filter_text.split(",")
			# Filter out Null values to avoid corrupting item count
			split_filter_text = filter(None, split_filter_text)
			split_filter_text = list(split_filter_text)
			nb_filters = len(split_filter_text)
			if nb_filters >= 1:
				self.filter_text.setText("Filtering to show " + selected_filter + ": " + str(split_filter_text))
				filter_condition = selected_filter + " in " + str(split_filter_text)

				if self.MetallaxisSettings.annotation_checkbox.isChecked():
					filtered_h5_table = pd.read_hdf(annotated_h5_output_name, key="df", where=filter_condition)
				else:
					filtered_h5_table = pd.read_hdf(h5_output_name, key="df", where=filter_condition)

			else:
				self.filter_text.setText(" ")
				throw_error_message("Please enter 2 or more values separated by a comma")
				return

		elif filter_text == "":
			self.filter_text.setText("No Filter Selected")
			if self.MetallaxisSettings.annotation_checkbox.isChecked():
				filtered_h5_table = pd.read_hdf(annotated_h5_output_name, key="df")
			else:
				filtered_h5_table = pd.read_hdf(h5_output_name, key="df")
		else:
			self.filter_text.setText("Filtering to show " + selected_filter + ": " + str(filter_text))
			filter_condition = selected_filter + "==\"" + filter_text + "\""
			if self.MetallaxisSettings.annotation_checkbox.isChecked():
				filtered_h5_table = pd.read_hdf(annotated_h5_output_name, key="df", where=filter_condition)
			else:
				filtered_h5_table = pd.read_hdf(h5_output_name, key="df", where=filter_condition)

		self.populate_table(filtered_h5_table)

	def write_h5_data_to_interface(self, complete_h5_file, h5_input_name):
		"""
		function that clears the interface if it already has data,
		then runs the annotate_h5() populate_table() functions. Requires h5
		read object, and h5 filename.
		"""
		# activate widgets that are disabled before VCF is chosen
		self.loaded_vcf_lineedit.setEnabled(True)
		self.loaded_vcf_label.setEnabled(True)
		self.meta_detected_filetype_label.setEnabled(True)
		self.metadata_area_label.setEnabled(True)
		self.viewer_tab_table_widget.setEnabled(True)
		self.filter_table_btn.setEnabled(True)
		self.filter_label.setEnabled(True)
		self.filter_lineedit.setEnabled(True)
		self.filter_box.setEnabled(True)

		# get column numbers for ID, POS, etc.
		self.progress_bar(47, "Extracting column data")

		column_names = list(complete_h5_file.keys())
		global chrom_col, id_col, pos_col, ref_col, alt_col, qual_col
		chrom_col = [i for i, s in enumerate(column_names) if 'CHROM' in s][0]
		id_col = [i for i, s in enumerate(column_names) if 'ID' in s][0]
		pos_col = [i for i, s in enumerate(column_names) if 'POS' in s][0]
		ref_col = [i for i, s in enumerate(column_names) if 'REF' in s][0]
		alt_col = [i for i, s in enumerate(column_names) if 'ALT' in s][0]
		qual_col = [i for i, s in enumerate(column_names) if 'QUAL' in s][0]

		# Clear metadata layouts
		self.empty_qt_layout(self.dynamic_metadata_label_results)
		self.empty_qt_layout(self.dynamic_metadata_label_tags)
		# clear statistics & plot layouts
		self.empty_qt_layout(self.dynamic_stats_value_label)
		self.empty_qt_layout(self.dynamic_stats_key_label)
		self.empty_qt_layout(self.stat_plot_layout)

		self.filter_table_btn.clicked.connect(self.filter_table)

		metadata_dict = {}
		for line in pd.read_hdf(h5_input_name, key="metadata").itertuples():
			# turn tuple into a list, but exclude the first item of list
			# because its the h5 index and not part of our original data
			line = list(line)[1:]
			metadata_tag = str(line[0])
			metadata_result = str(line[1])
			metadata_dict[metadata_tag] = metadata_result
			self.dynamic_metadata_label_tags.addWidget(
				QtWidgets.QLabel(metadata_tag, self))
			self.dynamic_metadata_label_results.addWidget(
				QtWidgets.QLabel(metadata_result, self))

		var_counts = {}
		for line in pd.read_hdf(h5_input_name, key="stats").itertuples():
			line = list(line)[1:]
			var_counts_key = str(line[0])
			var_counts_value = str(line[1])
			var_counts[var_counts_key] = var_counts_value

			# increase readability by showing a space instead of _ on interface (e.g. "1 Chrom Indel Count"
			# insetad of 1_Chrom_Indel_Count
			var_counts_key = var_counts_key.replace("_", " ")

			# add statistics dictionary data to interface, as labels
			new_label = QtWidgets.QLabel(str(var_counts_key), self)
			new_label.setWordWrap(True)
			self.dynamic_stats_key_label.addWidget(new_label)

			new_label = QtWidgets.QLabel(str(var_counts_value), self)
			new_label.setWordWrap(True)
			self.dynamic_stats_value_label.addWidget(new_label)

		if "ALT_Types" in var_counts:
			ALT_Types = eval(var_counts["ALT_Types"])

		self.progress_bar(49, "Plotting Statistics")

		# plot piechart of proportions of SNP/Indel
		if 'Total_SNP_Count' in var_counts:
			total_figure = plt.figure()
			graph = total_figure.add_subplot(111)
			graph.pie([var_counts['Total_SNP_Count'], var_counts['Total_Indel_Count']], labels=['SNP', 'Indels'],
			          autopct='%1.1f%%')
			# set x and y axes to be equal to get a perfect circle as a piechart
			graph.axis('equal')
			plt.title('Proportion of SNP/Indels')
			total_figure.tight_layout()
			graph.legend()
			self.stat_plot_layout.addWidget(FigureCanvas(total_figure))

		if "ALT_Types" in var_counts:
			# plot piechart of proportions of types of ALT
			# get the value for each ALT_Types key in order, per type of Alt so it can be graphed
			alt_values_to_plot = []
			for alt in ALT_Types:
				dict_key = alt + "_Alt_Count"
				alt_values_to_plot.append(var_counts[dict_key])

			total_figure = plt.figure()
			graph = total_figure.add_subplot(111)
			graph.pie(alt_values_to_plot, labels=ALT_Types, autopct='%1.1f%%')
			# set x and y axes to be equal to get a perfect circle as a piechart
			graph.axis('equal')
			plt.title('Proportion of different mutations')
			total_figure.tight_layout()
			graph.legend()
			self.stat_plot_layout.addWidget(FigureCanvas(total_figure))

		# plot piechart of proportions of types of ALT
		# get the nb of mutations for each chromosome
		if "List_Chromosomes" in var_counts:
			values_to_plot = []
			global list_chromosomes  # we're editing a global so it needs to be declared global again
			list_chromosomes = eval(var_counts['List_Chromosomes'])
			list_chromosomes = list(list_chromosomes)
			for chrom in list_chromosomes:
				dict_key = chrom + "_Chrom_Variant_Count"
				if dict_key in var_counts:
					values_to_plot.append(var_counts[dict_key])

			if values_to_plot != []:
				total_figure = plt.figure()
				graph = total_figure.add_subplot(111)
				# convert items in values_to_plot to int, so that matplotlib orders them correctly
				values_to_plot = [int(x) for x in values_to_plot]
				graph.bar(list(list_chromosomes), values_to_plot)
				plt.title('Distribution of Mutations by Chromosome')
				plt.xlabel('Chromosome')
				plt.ylabel('Number of Variants')
				total_figure.tight_layout()
				self.stat_plot_layout.addWidget(FigureCanvas(total_figure))

		self.chrom_selection_stat_comboBox.addItems(list_chromosomes)
		# setup variants by position graph for first chromosome in list
		if list_chromosomes[0]:
			self.changed_chrom_stat_combobox(list_chromosomes[0])
		# if selected chromosome changes, create new graph for variants by position for chosen chromosome
		self.chrom_selection_stat_comboBox.currentTextChanged.connect(self.changed_chrom_stat_combobox)

	def changed_chrom_stat_combobox(self, chrom=None):
		"""
		This function is run on changing the combobox on the statistics pane. On changing chromosome it
		will generate a new matplotlib graph for the chosen chromosome (which can be given as an
		argument to override combobox choice).
		"""
		# if no optional argument is provided then read chrom selection, from combobox
		if chrom == None:
			chrom = self.chrom_selection_stat_comboBox.currentText()
		chrom_data_subset_variants = []
		chrom_data_subset_ranges = []
		# empty layout from previous selection
		self.empty_qt_layout(self.chrom_stat_plot_layout)

		# avoid type errors by making chrom an int if its in numeric_columns
		if 'CHROM' in numeric_columns:
			chrom = int(chrom)

		# filter complete_h5_file to only results from chosen chromosome
		chrom_data = complete_h5_file[(complete_h5_file['CHROM'] == chrom)]
		min_pos = chrom_data['POS'].min()
		max_pos = chrom_data['POS'].max()
		# calculate the size of the chromosome based on smallest and largest POS values
		chrom_size = max_pos - min_pos
		# divide that size up so we can see variants by each part of the chromosome
		chrom_size_10 = int(chrom_size / 12)
		for i in range(min_pos, max_pos, chrom_size_10):
			if i != min_pos:  # don't count first one as there will be no data
				chrom_data_subset_range = str(i - chrom_size_10) + "-" + str(i)
				chrom_data_subset_filter = (chrom_data['POS'] >= (i - chrom_size_10)) & (chrom_data['POS'] <= i)
				chrom_data_subset_vars = len(chrom_data.loc[chrom_data_subset_filter])

				chrom_data_subset_ranges.append(chrom_data_subset_range)
				chrom_data_subset_variants.append(chrom_data_subset_vars)

		# plot data into bar plots
		total_figure = plt.figure()
		graph = total_figure.add_subplot(111)
		graph.bar(chrom_data_subset_ranges, chrom_data_subset_variants)
		plt.title('Distribution of Variants by Position in Chr ' + str(chrom))
		plt.xticks(rotation=70)
		plt.xlabel('Position in Chr ' + chrom)
		plt.ylabel('Number of Variants')
		total_figure.tight_layout()
		self.chrom_stat_plot_layout.addWidget(FigureCanvas(total_figure))

	def populate_table(self, selected_h5_data):
		if selected_h5_data is None:
			throw_error_message("Can't Populate Table: was passed a None object. Verify input H5 or VCF is not currupt")
			return

		# clear current table
		self.viewer_tab_table_widget.setRowCount(0)
		self.viewer_tab_table_widget.setColumnCount(0)

		# create empty table with correct length and width
		table_length = selected_h5_data.shape[0]
		table_width = selected_h5_data.shape[1]

		self.viewer_tab_table_widget.setRowCount(table_length)
		self.viewer_tab_table_widget.setColumnCount(table_width)

		column_names = list(selected_h5_data.keys())

		# clear chrom_filter_box
		self.filter_box.clear()
		self.filter_text.setText(" ")

		# set filter_box to list column_names
		self.filter_box.addItems(column_names)
		# read header labels from the dataframe keys
		self.viewer_tab_table_widget.setHorizontalHeaderLabels(column_names)

		# populate table
		vcf_line_nb = 0
		annotate_percent = 80
		for line in selected_h5_data.itertuples():
			if table_length >= 1000:
				# only update progress bar every 300 lines to avoid performance hit
				# from doing it every line
				if vcf_line_nb % 300 == 0:
					annotate_progress = annotate_percent + (vcf_line_nb / table_length) * (100 - annotate_percent)
					self.progress_bar(float(annotate_progress), "Populating Table from H5")

			# turn tuple into a list, but exclude the first item of list
			# because its the h5 index and not part of our original data
			line = list(line)[1:]

			vcf_field_nb = 0
			for vcf_field in line:
				# we replaced "." with np.NaN earlier so replace it back now
				if vcf_field is np.NaN:
					vcf_field = "."
				elif str(vcf_field) == "nan":
					vcf_field = "."
				self.viewer_tab_table_widget.setItem(
					vcf_line_nb, vcf_field_nb, QtWidgets.QTableWidgetItem(str(vcf_field)))
				vcf_field_nb += 1
			vcf_line_nb += 1
		self.progress_bar(100, "Populating Table...done")
		# close progress bar when file is completely loaded
		self.MetallaxisProgress.close()


progress_window_object, progress_base_object = uic.loadUiType(MetPROGui)


class MetallaxisProgress(progress_base_object, progress_window_object):
	"""
	Status bar that shows progress when opening files with Metallaxis
	"""

	def __init__(self):
		super(progress_base_object, self).__init__()
		self.setupUi(self)
		self.setWindowTitle("Metallaxis")


settings_window_object, settings_base_object = uic.loadUiType(MetSETui)


class MetallaxisSettings(settings_base_object, settings_window_object):
	"""
	Settings window that allows changing of default options.
	"""

	def __init__(self):
		super(settings_base_object, self).__init__()
		self.setupUi(self)
		self.setWindowTitle("Metallaxis Settings")
		self.compression_comboBox.addItems(['zlib', 'blosc', 'bzip2', 'lzo'])
		self.vcf_chunk_size.setText("5000")
		# Center settings pannel on screen
		qt_rectangle = self.frameGeometry()
		center_point = QDesktopWidget().availableGeometry().center()
		qt_rectangle.moveCenter(center_point)
		# set minimum size to size it takes up
		self.setMinimumSize(self.sizeHint())


if __name__ == '__main__':
	MetallaxisApp = QApplication(sys.argv)
	MetallaxisGui = MetallaxisGuiClass()

	# get annotation by default
	MetallaxisGui.MetallaxisSettings.annotation_checkbox.setChecked(True)

	# get input file
	MetallaxisGui.progress_bar(2, "Parsing arguments")
	global h5_only
	h5_only = False

	# if we load a h5 file from a previous analysis then we skip the usual analysis and verifications
	if len(sys.argv) == 2 and sys.argv[1].endswith(".h5"):
		h5_only = True
		complete_h5_file = load_hdf5(sys.argv[1])
		MetallaxisGui.write_h5_data_to_interface(complete_h5_file, sys.argv[1])

	elif len(sys.argv) == 2:  # if we give it vcf or compressed vcf
		selected_vcf = os.path.abspath(sys.argv[1])

	elif len(sys.argv) == 1:  # if we don't give any args then open file picker
		selected_file = MetallaxisGui.select_file()
		if selected_file is False:
			throw_error_message("No file selected, qutting Metallaxis")
			sys.exit(1)

		elif selected_file.endswith(".h5"):
			complete_h5_file = load_hdf5(selected_file)
			h5_only = True
			MetallaxisGui.write_h5_data_to_interface(complete_h5_file, selected_file)
		else:
			selected_vcf = selected_file

	else:  # if we give more than 1 arg
		print("Error: Metallaxis can only take one argument, a vcf file")
		exit(1)

	# if we loaded a VCF file, do the verification and analysis
	if not h5_only:

		# get metadata and variant counts from vcf
		metadata_dict, var_counts, decompressed_file = parse_vcf(selected_vcf)

		# convert vcf into a hdf5 object
		h5_file = h5_encode(selected_vcf, decompressed_file, var_counts, metadata_dict)

		# Read H5 for actual table populating
		complete_h5_file = pd.read_hdf(h5_file, key="df")

		# populate interface with information from hdf5
		MetallaxisGui.write_h5_data_to_interface(complete_h5_file, h5_file)

		# get annotation data
		if MetallaxisGui.MetallaxisSettings.annotation_checkbox.isChecked():
			if h5_only is not True:
				try:
					complete_h5_file = annotate_h5(h5_file)
					complete_h5_file = pd.read_hdf(complete_h5_file, key="df")
				except:
					throw_warning_message("Annotation did not succeed, proceeding with non-annotated h5")
					MetallaxisGui.MetallaxisSettings.annotation_checkbox.setChecked(False)

	# actions that are to be done once we have our h5 file (if we loaded a h5 file then it starts here)
	# populate table
	MetallaxisGui.populate_table(complete_h5_file)

	# show GUI
	MetallaxisGui.show()

	# exit program on quitting the GUI
	sys.exit(MetallaxisApp.exec_())
