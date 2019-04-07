#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import requests
import os
import yaml  # for reading settings file
import pathlib  # for making the folder where we store data
import platform  # for determining OS and therefore where to store data
import re
import sys
from shutil import copyfile  # for save analysis

import magic  # to detect filetype from file header
import numpy as np  # to handle arrays and NaN
import pandas as pd  # to handle dataframes
import sqlite3  # handle sqlite db
import wget  # to download from FTP

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
from PyQt5.QtSvg import QSvgWidget

# Import SVG Drawing Classes
from metallaxis import SVGClasses

# for plotting graphs
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt

plt.style.use('seaborn')

# allow <Ctrl-c> to terminate the GUI
import signal

signal.signal(signal.SIGINT, signal.SIG_DFL)

# Determine where to store temporary Metallaxis Data
home_dir = os.path.expanduser('~')
if platform.system() == "Linux":
	config_directory = home_dir + "/.metallaxis/"
elif platform.system() == "Darwin":
	config_directory = home_dir + "/Library/Caches/Metallaxis/"
elif platform.system() == "Windows":
	config_directory = os.path.expandvars(r'%APPDATA%\Metallaxis\\')

# make data folder and parent folders if they don't exist
pathlib.Path(config_directory).mkdir(parents=True, exist_ok=True)

config_file = os.path.join(config_directory, "metallaxis_config.yaml")

# Interface XML files
current_file_path = __file__
current_file_dir = os.path.dirname(current_file_path)
MetGUIui = os.path.join(current_file_dir, "MetallaxisGui.ui")
MetSETui = os.path.join(current_file_dir, "MetallaxisSettings.ui")
MetPROGui = os.path.join(current_file_dir, "MetallaxisProgress.ui")

# Annotation executables
snpsift_jar = os.path.join(current_file_dir, "annotation/SnpSift.jar")
snpeff_jar = os.path.join(current_file_dir, "annotation/SnpEff.jar")


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


def read_config(config_file):
	with open(config_file, 'r') as configf:
		try:
			config = yaml.safe_load(configf)
			return config
		except yaml.YAMLError as exc:
			throw_error_message(exc)
			return False


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


def is_number_bool(sample):
	try:
		float(sample)
	except:
		return False
	return True


def set_col_to_numeric_if_isdigit(column, chunk, numeric_columns_list):
	"""
	Determines which columns of a given chunk are ints or floats, and removes
	them from a list of columns if they are neither.
	"""

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
		if "True" in row:
			del_col(column)
		if "False" in row:
			del_col(column)
		if bool(re.match('^[0-9]+$', row)) is False:
			if is_number_bool(row) is False:
				if column in numeric_columns_list:
					numeric_columns_list.remove(column)


def load_sqlite(sqlite_filename):
	"""
	Loads a previously created .sqlite file. If an analysis had already been done on a VCF
	and the analysis file was saved in sqlite format, it can be loaded here.
	Accepts a sqlite3 file as argument.
	"""
	if os.path.isfile(sqlite_filename):
		# Verify file exists and return its read object
		MetallaxisGui.loaded_vcf_lineedit.setText(os.path.abspath(sqlite_filename))
		loaded_db_connection = sqlite3.connect(sqlite_filename)
		complete_sqlite_file = pd.read_sql("SELECT * FROM df", loaded_db_connection)
		return complete_sqlite_file

	else:
		throw_error_message("Selected file does not \
		exist. You specified : " + str(sqlite_filename))


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
	global metadata_num  # make global as this will be used in database_encode
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

		# set global variable to be total number of lines of metadata (will be used to extract just the data in database_encode()
		metadata_num = int(line_num - variant_num)

	# return Error or warning depending on variant numbers
	if variant_num == 0:
		throw_error_message("VCF is empty, there are no variants at all in this vcf, please use a different vcf")
		return False
	elif variant_num < 5:
		# TODO: uncomment the lines below, and remove "return TRUE"
		return True
	# throw_error_message("VCF contains too few variants to analyse, please use a different vcf")
	# return False
	elif 5 < variant_num < 30:
		# TODO: Uncomment the below as is just to speed up development that I removed GUI
		# throw_warning_message("VCF contains very few variants, only rudimentary statistics can be performed")
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
	with open(decompressed_file) as decompressed_out:
		for line in decompressed_out:
			if not line.startswith('#'):
				line = line.split("\t")
				alt = str(line[alt_col])
				if set(alt).issubset(set('ACTG')) and len(alt) != 1:
					alt_types_only_snp = False

	with open(decompressed_file) as decompressed_out:
		for line in decompressed_out:
			if not line.startswith('#'):
				line = line.split("\t")

				# Convert CHR 1 to be CHR 01 for correct sorting
				my_chrom = line[chrom_col]
				if is_number_bool(my_chrom):
					if int(my_chrom) < 10:
						my_chrom = '0' + str(my_chrom)
				else:
					my_chrom = str(my_chrom)
				line[chrom_col] = my_chrom
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
					add_to_dict_iterator(variant_stats, line[chrom_col] + "_Chrom_SNP_Count", 1)
					add_to_dict_iterator(variant_stats, line[chrom_col] + "_Chrom_Variant_Count", 1)
				else:
					variant_stats["Total_Indel_Count"] += 1
					add_to_dict_iterator(variant_stats, line[chrom_col] + "_Chrom_Indel_Count", 1)
					add_to_dict_iterator(variant_stats, line[chrom_col] + "_Chrom_Variant_Count", 1)
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
	with open(decompressed_file) as decompressed_out:
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


def already_annotated(vcf_file):
	with open(vcf_file, 'r') as vcf_read_obj:
		for line in vcf_read_obj:
			if line.startswith("##INFO=<ID=ANN"):
				return True


def annotate_vcf(vcf_file):
	clinvar_url = "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
	clinvar_index_url = "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi"

	dbSNP_url = "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/00-All.vcf.gz"
	dbSNP_index_url = "ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh38p7/VCF/00-All.vcf.gz.tbi"

	dbSNP_index_path = os.path.join(config['working_dir'], 'dbsnp.vcf.gz.tbi')
	if not os.path.isfile(dbSNP_index_path):
		wget.download(dbSNP_index_url, out=dbSNP_index_path)

	dbSNP_path = os.path.join(config['working_dir'], 'dbsnp.vcf.gz')
	if not os.path.isfile(dbSNP_path):
		MetallaxisGui.progress_bar(35, "Downloading Annotation databases (this will take a while)")
		wget.download(dbSNP_url, out=dbSNP_path)

	clinvar_index_path = os.path.join(config['working_dir'], 'clinvar.vcf.gz.tbi')
	if not os.path.isfile(clinvar_index_path):
		wget.download(clinvar_index_url, out=clinvar_index_path)

	clinvar_path = os.path.join(config['working_dir'], 'clinvar.vcf.gz')
	if not os.path.isfile(clinvar_path):
		MetallaxisGui.progress_bar(37, "Downloading Annotation databases (this will take a while)")
		wget.download(clinvar_url, out=clinvar_path)

	MetallaxisGui.progress_bar(5, "Running annotation on VCF (this will take some time)")

	annotate_cmd = "java -Xmx" + config['max_memory'] + "G -jar " + snpeff_jar + " " + config[
		'genome_version'] + " " + vcf_file + " > " + annotated_vcf_output_filename
	os.system(annotate_cmd)

	annotate_cmd = "java -Xmx" + config[
		'max_memory'] + "G -jar " + snpsift_jar + " annotate " + dbSNP_path + " " + annotated_vcf_output_filename + " > " + annotated_vcf_output_filename + "1"
	os.system(annotate_cmd)

	annotate_cmd = "java -Xmx" + config[
		'max_memory'] + "G -jar " + snpsift_jar + " annotate " + clinvar_path + " " + annotated_vcf_output_filename + "1 > " + annotated_vcf_output_filename + "2"
	os.system(annotate_cmd)

	os.remove(annotated_vcf_output_filename)
	os.remove(annotated_vcf_output_filename + "1")
	os.rename(annotated_vcf_output_filename + "2", annotated_vcf_output_filename)
	return annotated_vcf_output_filename


def database_encode(decompressed_file, variant_stats, metadata_dict):
	"""
	accepts as input a raw vcf file, and dictionaries for variant counts (variant_stats) and metadata
	(metadata_dict). Then encodes them as tables into a database.
	Returns filename of created sqlite database.
	"""
	sqlite_output = sqlite3.connect(sqlite_output_name)
	cursor = sqlite_connection.cursor()
	cursor.execute("DROP TABLE IF EXISTS df;")
	cursor.execute("DROP TABLE IF EXISTS stats;")
	cursor.execute("DROP TABLE IF EXISTS metadata;")
	cursor.execute("DROP TABLE IF EXISTS chrom_genes;")
	sqlite_output.commit()

	# write each entry from metadata_dict to a new "metadata" table in database
	for metadata_line_nb in metadata_dict:
		metadata_tag = str(metadata_dict[metadata_line_nb][1])
		metadata_result = str(metadata_dict[metadata_line_nb][2])
		if not metadata_tag.isupper():
			metadata_line = {'Tag': metadata_tag, 'Result': metadata_result}
			metadata_line = pd.DataFrame(
				metadata_line, index=[metadata_line_nb])

			metadata_line.to_sql('metadata', sqlite_output, if_exists='append', index=False)

	# write each entry to a new "stats" table in database
	stat_table_index = 0
	for key, value in variant_stats.items():
		key = str(key)
		value = str(value)
		if len(key) > 40:
			key = key[:40] + "..."
		if len(value) > 200:
			value = value[:200] + "..."
		variant_stats_line = {'Tag': key, 'Result': value}
		variant_stats_line = pd.DataFrame(
			variant_stats_line, index=[stat_table_index])

		variant_stats_line.to_sql('stats', sqlite_output, if_exists='append', index=False)
		stat_table_index += 1

	# chunked_vcf_len = sum(1 for row in open(decompressed_file, 'r'))

	# PARSE INFO COLUMNS
	# the info column of a vcf is long and hard to read if
	# displayed as is, but it is composed of multiple key:value tags
	# that we can parse as new columns, making them filterable.
	global info_cols_to_add
	info_cols_to_add = set()
	chunked_vcf = pd.read_csv(decompressed_file,
	                          sep="\t",
	                          # skip rows that started with "#"
	                          skiprows=range(0, metadata_num - 1),
	                          # use chunk size that was set in settings
	                          chunksize=int(config['vcf_chunk_size']),
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

	anno_info_cols_to_add = []
	with open(decompressed_file, 'r') as vcf_read_obj:
		for line in vcf_read_obj:
			if not line.startswith('#'):
				break

			if line.startswith('##INFO=<ID=ANN'):
				line = line.split("Functional annotations: '")[1]
				line = line.split("|")
				for col in line:
					col = col.replace('"', "")
					col = col.replace('>', "")
					col = col.replace("'", '')
					col = re.sub('^\s+', '', col)
					col = re.sub('\s+$', '', col)
					anno_info_cols_to_add.append(col)

	chunked_vcf = pd.read_csv(decompressed_file,
	                          sep="\t",
	                          skiprows=range(0, metadata_num - 1),
	                          chunksize=int(config['vcf_chunk_size']),
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

		# Get rid of INFO column now that it exists as multiple columns
		chunk = chunk.drop(columns=['INFO'])

		if len(anno_info_cols_to_add) > 0:
			for col in anno_info_cols_to_add:
				chunk[col] = "."

		line_nb = 0
		if "ANN" in chunk:
			for line in chunk["ANN"]:
				line_split = line.split("|")
				for col_num in range(0, len(anno_info_cols_to_add)):
					current_col = anno_info_cols_to_add[col_num]
					if len(line_split) > col_num:
						chunk[current_col].values[line_nb] = line_split[col_num]
				line_nb += 1

			chunk = chunk.drop(columns=['ANN'])

		# Rename column so we get 'CHROM' not '#CHROM' from chunk.keys()
		chunk.rename(columns={'#CHROM': 'CHROM'}, inplace=True)

		# To fix sorting problems
		chunk['CHROM'] = chunk['CHROM'].replace('1', '01')
		chunk['CHROM'] = chunk['CHROM'].replace('2', '02')
		chunk['CHROM'] = chunk['CHROM'].replace('3', '03')
		chunk['CHROM'] = chunk['CHROM'].replace('4', '04')
		chunk['CHROM'] = chunk['CHROM'].replace('5', '05')
		chunk['CHROM'] = chunk['CHROM'].replace('6', '06')
		chunk['CHROM'] = chunk['CHROM'].replace('7', '07')
		chunk['CHROM'] = chunk['CHROM'].replace('8', '08')
		chunk['CHROM'] = chunk['CHROM'].replace('9', '09')

		# SET COLS WITH NUMBERS TO BE NUMERIC TYPE
		# set columns that only contain numbers to be numeric dtype
		# otherwise they are just strings, and can't be used with a
		# dash separated filter

		# make a list from all column names, we will later remove the
		# non-numeric columns from the list
		global numeric_columns
		numeric_columns = list(chunk.keys())
		#
		chunk = chunk.replace(".", np.NaN)

		# run the function for every column to remove non-numeric columns
		for column in chunk.keys():
			set_col_to_numeric_if_isdigit(column, chunk, numeric_columns)

		# convert the remaining columns in "numeric_columns" list to numeric datatype
		for column in numeric_columns:
			chunk[column] = pd.to_numeric(chunk[column])

		chunk.to_sql('df', sqlite_output, if_exists='append', index=False)

	# for col in anno_info_cols_to_add:
	# 	sql_request = "DELETE FROM df WHERE %s IS NULL OR trim(%s) = '';" % (col, col)
	# 	conn.execute(sql_request)
	# 	conn.commit()
	# 	conn.close()
	return sqlite_output


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

		# Setup inital GUI
		self.graphicsView.setMaximumHeight(0)

		# Center GUI on screen
		qt_rectangle = self.frameGeometry()
		center_point = QDesktopWidget().availableGeometry().center()
		qt_rectangle.moveCenter(center_point)

		# buttons on interface
		self.open_vcf_button.clicked.connect(self.select_and_parse)
		# menus on interface
		self.actionOpen_VCF.triggered.connect(self.select_and_parse)
		self.actionSave_Analysis.triggered.connect(self.save_analysis)
		self.actionQuit.triggered.connect(self.close)

		# Link "Github Page" button on menu to its URL
		def open_github(url):
			url = "https://github.com/SL-LAIDLAW/Metallaxis"
			QDesktopServices.openUrl(QtCore.QUrl(url))

		self.actionGithub_Page.triggered.connect(open_github)

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

	def progress_bar(self, percent, message):
		"""
		Function that updates the on screen progress bar. Accepts two arguments, one
		for % progression, the other for the current status to display.
		"""
		self.MetallaxisProgress.progressbar_message.setText(message)
		percent = round(percent, 2)
		self.MetallaxisProgress.progressbar_progress.setValue(percent)
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

	def save_svg(self):
		"""
		Open a dialog where the user can chose where to place to save svg of variant view
		"""
		save_dialog = QtWidgets.QFileDialog()
		save_dialog.setAcceptMode(save_dialog.AcceptSave)
		save_folder = save_dialog.getSaveFileName(self, 'Save variant view as vector image', filter="*.svg")[0]
		copyfile(svg_output_name, save_folder)

	def save_analysis(self):
		"""
		Open a dialog where the user can chose where to place to save analysis as sqlite
		"""
		save_dialog = QtWidgets.QFileDialog()
		save_dialog.setAcceptMode(save_dialog.AcceptSave)
		save_folder = save_dialog.getSaveFileName(self, 'Save Analayis as database', filter="*.sqlite")[0]
		copyfile(sqlite_output_name, save_folder)

	def select_file(self):
		"""
		Opens a file dialog where the user can chose an input file.
		"""
		select_dialog = QtWidgets.QFileDialog()
		select_dialog.setAcceptMode(select_dialog.AcceptSave)
		selected_vcf = select_dialog.getOpenFileName(self, filter="VCF Files (*.vcf \
			*.vcf.xz *.vcf.gz *.vcf.bz2) ;;Metallaxis Database Files(*.sqlite) ;;All Files(*.*)")
		selected_vcf = selected_vcf[0]
		# if the user cancels the select_file() dialog, then run select again
		while selected_vcf == "":
			throw_error_message("No selected file")
			return False
		return selected_vcf

	def select_and_parse(self, cli_arg=False):
		"""
		Runs select_file() to get an input and based off extension, either runs the parse_VCF() or displays the contents
		of the database file directly. Accepts no arguments, and returns nothing. This function exists solely to call other
		functions, as menu items in PyQt can only call one function.
		"""
		self.MetallaxisProgress.show()

		if not cli_arg:
			selected_file = self.select_file()
			if selected_file is False:
				return  # If the user cancels the file picker
		else:
			selected_file = cli_arg

		load_session = False
		if selected_file.endswith(".sqlite"):
			load_session = True
			loaded_database = load_sqlite(selected_file)
			self.write_database_to_interface(loaded_database)
		else:
			selected_vcf = selected_file

		# reopen progress bar for loading new file
		self.MetallaxisProgress = MetallaxisProgress()
		self.MetallaxisProgress.show()

		if not load_session:

			# get metadata and variant counts from vcf
			if self.MetallaxisSettings.annotation_checkbox.isChecked():
				if not already_annotated(selected_vcf):
					selected_vcf = annotate_vcf(selected_vcf)

			# parse vcf, convert to a database, and write database data to interface
			metadata_dict, var_counts, decompressed_file = parse_vcf(selected_vcf)
			global db_connection
			db_connection = database_encode(decompressed_file, var_counts, metadata_dict)
			loaded_database = pd.read_sql("SELECT * FROM df", db_connection)
			self.write_database_to_interface(loaded_database)

		# populate table
		self.populate_table(loaded_database)

	def hide_graphics_view(self):
		self.graphicsView.setMaximumHeight(0)

	def show_graphics_view(self):
		self.graphicsView.setMaximumHeight(16777215)

	def toggle_graphics_view(self):
		current_height = self.graphicsView.maximumHeight()
		if current_height == 0:
			self.show_graphics_view()
		else:
			self.hide_graphics_view()

	def reload_generate_variant_graphic(self):
		self.generate_variant_graphic(True)

	def generate_variant_graphic(self, read_pos_input=False):
		self.show_graphics_view()
		current_row = self.viewer_tab_table_widget.currentRow()
		# if no row selected then stop function
		if current_row == -1:
			return

		current_chr = self.viewer_tab_table_widget.item(current_row, 0).text()
		current_pos = int(self.viewer_tab_table_widget.item(current_row, 1).text())
		current_alt = self.viewer_tab_table_widget.item(current_row, 4).text()

		# SVG Setup
		varScene = SVGClasses.Scene('variant_scene')
		line_length = 650
		varScene.add(SVGClasses.Line((50, 100), (line_length + 50, 100)))

		my_query = "SELECT * FROM df where CHROM =='" + str(current_chr) + "'"
		chrom_data = pd.read_sql(my_query, db_connection)

		def get_default_min_max():
			# Find default positions
			if (current_pos / 2) < 2500000:
				min_pos = (current_pos / 2)
				max_pos = current_pos + (current_pos / 2)
			elif (current_pos - 2500000) > chrom_data['POS'].min():
				min_pos = (current_pos - (2500000 - 1))
				max_pos = (current_pos + (2500000 - 1))
			else:
				min_pos = chrom_data['POS'].min()
				max_pos = min_pos + (5000000 - 1)

			min_pos = int(float(min_pos))
			max_pos = int(float(max_pos))
			self.graphics_max_pos_textin.setText(str(max_pos))
			self.graphics_min_pos_textin.setText(str(min_pos))
			return min_pos, max_pos

		if not read_pos_input:
			min_pos, max_pos = get_default_min_max()
		else:
			max_pos = int(self.graphics_max_pos_textin.text())
			min_pos = int(self.graphics_min_pos_textin.text())
			if min_pos == 0 or max_pos == 0:
				min_pos, max_pos = get_default_min_max()
			# inverse, if user put the two wrong way around
			if max_pos < min_pos:
				self.graphics_max_pos_textin.setText(str(min_pos))
				self.graphics_min_pos_textin.setText(str(max_pos))
				min_pos = int(self.graphics_max_pos_textin.text())
				max_pos = int(self.graphics_min_pos_textin.text())

		# reduce length of request if it will fail due to size
		if (max_pos - min_pos) > (5000000 - 1):
			max_pos = min_pos + (5000000 - 1)

		# if ENSEMBL gene location not already found, then fetch it
		chrom_genes_empty_bool = pd.read_sql(
			"SELECT name FROM sqlite_master WHERE type='table' AND name='chrom_genes';", db_connection).empty
		if chrom_genes_empty_bool:
			min_db_pos, max_db_pos = 0, 0  # make sure annotation runs
		else:
			min_db_pos = pd.read_sql(
				'SELECT DISTINCT MIN(start) FROM chrom_genes WHERE chrom = ' + current_chr + ' AND gene_id <> "";',
				db_connection)
			max_db_pos = pd.read_sql(
				'SELECT DISTINCT MAX(end) FROM chrom_genes WHERE chrom = ' + current_chr + ' AND gene_id <> "";',
				db_connection)
			min_db_pos = list(min_db_pos.iterrows())[0][0]
			max_db_pos = list(max_db_pos.iterrows())[0][0]

		if not min_db_pos <= current_pos <= max_db_pos:

			# if chrom '01' then flatten to '1' so works with the API
			if is_number_bool(current_chr):
				current_chr = str(int(current_chr))

			# TODO: ask for organsim in latin in settings
			api_request_str = "https://rest.ensembl.org/overlap/region/%s/%s:%d-%d?feature=gene;feature=transcript;\
			feature=cds;feature=exon" % (config['organism'], current_chr, min_pos, max_pos)

			def request_ensembl_gene_pos():
				api_request = requests.get(api_request_str, headers={"Content-Type": "application/json"})
				if not api_request.ok:
					throw_error_message(
						"ENSEMBLE rejected API Call, please try a smaller nucleotide range, or verify internet access")
					return False
				else:
					return api_request

			ensembl_gene_pos_req = request_ensembl_gene_pos()
			if ensembl_gene_pos_req:
				sqlite_output = sqlite3.connect(sqlite_output_name)

				# if we got any data back from ENSEMBL
				if len(ensembl_gene_pos_req.json()) > 0:
					ensembl_gene_pos_df = pd.DataFrame.from_dict(ensembl_gene_pos_req.json())
					ensembl_gene_pos_df['chrom'] = current_chr
					ensembl_gene_pos_df.to_sql('chrom_genes', sqlite_output, if_exists='append', index=True)

		# get percentage of line length that variant is to show up
		def rel_position_on_line(position):
			rel_pos = position - min_pos
			rel_pos = float(float(rel_pos) / float(max_pos - min_pos))
			rel_pos = rel_pos * line_length
			rel_pos = rel_pos + 50
			return rel_pos

		alleles_to_draw = []
		my_query = "SELECT external_name,start,end FROM (SELECT DISTINCT gene_id, external_name,start,end FROM chrom_genes WHERE gene_id <> '' and start >= %d and end <= %d);" % (
		min_pos, max_pos)
		alleles_list = pd.read_sql(my_query, db_connection)
		if not alleles_list.empty:
			allele_nb = 0
			for index, line in alleles_list.iterrows():

				start_pos_on_line = rel_position_on_line(line['start'])
				if start_pos_on_line < 50:
					start_pos_on_line = 50

				end_pos_on_line = rel_position_on_line(line['end'])
				if end_pos_on_line > (line_length + 50):
					end_pos_on_line = line_length + 50

				new_allele = SVGClasses.Allele(rel_position_on_line(line['start']), end_pos_on_line,
				                               line['external_name'], allele_nb)

				# add allele to list thatll be added to SVG
				alleles_to_draw.append(new_allele)

				# this means alleles get all the colors
				allele_nb += 1
				if allele_nb > 4:
					allele_nb = 0

		te_pos = rel_position_on_line(current_pos)

		self.empty_qt_layout(self.graphicsView_layout)
		self.graphics_chr_label.setText(str(current_chr))

		if current_alt.startswith('<INS'):
			varScene.add(SVGClasses.TE(te_pos, "ins"))
		if current_alt.startswith('<DEL'):
			varScene.add(SVGClasses.TE(te_pos, "del"))
		else:
			varScene.add(SVGClasses.TE(te_pos))

		if len(alleles_to_draw) > 0:
			for allele in alleles_to_draw:
				if allele.getWidth() < 7:
					alleles_to_draw.remove(allele)

			for allele in alleles_to_draw:
				if len(alleles_to_draw) > 6:
					allele = allele.removeName()
				varScene.add(allele)

		varScene.write_svg(svg_output_name)
		self.graphicsView_layout.addWidget(QSvgWidget(svg_output_name))

	def filter_table(self):
		"""
		Filters table based on chosen filters. Comma separated, dash separated, and single filters exist. This function
		reads the chosen filter from the interface, requests rows matching those filters from the database and
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
			split_filter_text = filter_text.split("-")
			# Filter out Null values to avoid corrupting item count
			split_filter_text = filter(None, split_filter_text)
			split_filter_text = list(split_filter_text)
			split_filter_text = ["'" + str(split_filter_text[i]) + "'" for i in range(0, len(split_filter_text))]
			if len(split_filter_text) == 2:
				self.filter_text.setText(
					"Filtering to show " + selected_filter + " from " + str(split_filter_text[0]) + " to " + str(
						split_filter_text[1]))

				if split_filter_text[0] > split_filter_text[1]:
					filter_condition = selected_filter + ">=" + split_filter_text[
						1] + " and " + selected_filter + "<=" + \
					                   split_filter_text[0]
				elif split_filter_text[0] < split_filter_text[1]:
					filter_condition = selected_filter + ">=" + split_filter_text[
						0] + " and " + selected_filter + "<=" + \
					                   split_filter_text[1]
				else:
					filter_condition = selected_filter + "=" + split_filter_text[0]


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
				split_filter_text = str(split_filter_text).replace('[', '')
				split_filter_text = str(split_filter_text).replace(']', '')
				split_filter_text = str(split_filter_text).replace(',', ' OR ' + selected_filter + "=")
				filter_condition = selected_filter + "=" + split_filter_text

			else:
				self.filter_text.setText(" ")
				throw_error_message("Please enter 2 or more values separated by a comma")
				return

		elif filter_text == "":
			self.filter_text.setText("No Filter Selected")
			filtered_table = pd.read_sql_query("SELECT * from df", sqlite_connection)
			self.populate_table(filtered_table)
			return

		else:
			self.filter_text.setText("Filtering to show " + selected_filter + ": " + str(filter_text))
			filter_condition = "'" + selected_filter + "'=='" + filter_text + "'"

		filtered_table = pd.read_sql_query("SELECT * from df where " + filter_condition, sqlite_connection)
		self.populate_table(filtered_table)

	def write_database_to_interface(self, loaded_database):
		"""
		function that clears the interface if it already has data,
		then runs the populate_table() functions.
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
		self.view_variant_btn.setEnabled(True)
		self.chrom_selection_stat_comboBox.setEnabled(True)
		self.chrom_selection_label.setEnabled(True)

		# get column numbers for ID, POS, etc.
		self.progress_bar(47, "Extracting column data")

		column_names = list(loaded_database.keys())
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
		# self.empty_qt_layout(self.dynamic_stats_value_label)
		# self.empty_qt_layout(self.dynamic_stats_key_label)
		self.empty_qt_layout(self.stat_plot_layout)

		self.filter_table_btn.clicked.connect(self.filter_table)
		self.view_variant_btn.clicked.connect(self.generate_variant_graphic)
		self.graphics_hide_view_btn.clicked.connect(self.hide_graphics_view)
		self.graphics_reload_btn.clicked.connect(self.reload_generate_variant_graphic)
		self.export_svg_toolbtn.clicked.connect(self.save_svg)

		metadata_sql_result = pd.read_sql_query("SELECT DISTINCT Tag,Result FROM metadata", sqlite_connection)
		for i in range(0, len(metadata_sql_result)):
			self.dynamic_metadata_label_tags.addWidget(QtWidgets.QLabel(metadata_sql_result['Tag'][i], self))
			self.dynamic_metadata_label_results.addWidget(QtWidgets.QLabel(metadata_sql_result['Result'][i], self))

		var_counts = {}
		stats_sql_result = pd.read_sql_query("SELECT * FROM stats", sqlite_connection)
		for i in range(0, len(stats_sql_result)):
			var_counts_key = stats_sql_result['Tag'][i]
			var_counts_value = stats_sql_result['Result'][i]
			var_counts[var_counts_key] = var_counts_value

		self.progress_bar(49, "Plotting Statistics")

		if "ALT_Types" in var_counts:
			ALT_Types = eval(var_counts["ALT_Types"])
			# plot piechart of proportions of types of ALT
			# get the value for each ALT_Types key in order, per type of Alt so it can be graphed
			alt_values_to_plot = []
			for alt in ALT_Types:
				dict_key = alt + "_Alt_Count"
				alt_values_to_plot.append(eval(var_counts[dict_key]))

			alt_types_clean_label = [str(alt_type).replace('<', '').replace('>', '') for alt_type in ALT_Types]
			total_figure = plt.figure()
			graph = total_figure.add_subplot(111)
			graph.pie(alt_values_to_plot, labels=alt_types_clean_label, autopct='%1.1f%%')
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
				graph_df = pd.Series(values_to_plot, index=list_chromosomes)
				graph_df = graph_df.sort_index()
				# graph_df graph_df.sort_values(column=""
				graph.bar(graph_df.index, graph_df.values, tick_label=graph_df.index)
				plt.title('Distribution of Mutations by Chromosome')
				plt.xlabel('Chromosome')
				plt.ylabel('Number of Variants')
				total_figure.tight_layout()
				self.stat_plot_layout.addWidget(FigureCanvas(total_figure))
				self.chrom_selection_stat_comboBox.addItems(graph_df.index)

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

		# filter loaded_database to only results from chosen chromosome
		my_query = "SELECT * FROM df where CHROM =='" + str(chrom) + "'"
		chrom_data = pd.read_sql(my_query, db_connection)
		# chrom_data = loaded_database[(loaded_database['CHROM'] == chrom)]
		min_pos = chrom_data['POS'].min()
		max_pos = chrom_data['POS'].max()
		# calculate the size of the chromosome based on smallest and largest POS values
		chrom_size = max_pos - min_pos

		if chrom_size is np.NaN:
			chrom_size = 0

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
		plt.xlabel('Position in Chr ' + str(chrom))
		plt.ylabel('Number of Variants')
		total_figure.tight_layout()
		self.chrom_stat_plot_layout.addWidget(FigureCanvas(total_figure))

	def populate_table(self, selected_data):
		if selected_data is None:
			throw_error_message(
				"Can't Populate Table: was passed a None object. Verify that input databbase or VCF is not corrupt")
			return

		# clear current table
		self.viewer_tab_table_widget.setRowCount(0)
		self.viewer_tab_table_widget.setColumnCount(0)

		# create empty table with correct length and width
		table_length = selected_data.shape[0]
		table_width = selected_data.shape[1]

		self.viewer_tab_table_widget.setRowCount(table_length)
		self.viewer_tab_table_widget.setColumnCount(table_width)

		column_names = list(selected_data.keys())

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
		for line in selected_data.itertuples():
			if table_length >= 1000:
				# only update progress bar every 300 lines to avoid performance hit
				# from doing it every line
				if vcf_line_nb % 300 == 0:
					annotate_progress = annotate_percent + (vcf_line_nb / table_length) * (100 - annotate_percent)
					self.progress_bar(float(annotate_progress), "Populating Table")

			# turn tuple into a list, but exclude the first item of list
			# because its the pandas index and not part of our original data
			line = list(line)[1:]

			vcf_field_nb = 0
			for vcf_field in line:
				# we replaced "." with np.NaN earlier so replace it back now
				if vcf_field is np.NaN:
					vcf_field = "."
				elif str(vcf_field) == "nan":
					vcf_field = "."
				elif vcf_field is "":
					vcf_field = "."
				elif vcf_field is None:
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

		# Initiate Buttons
		self.change_wd_btn.clicked.connect(self.set_working_dir)
		self.save_settings_btn.clicked.connect(self.save_settings)

		# Center settings pannel on screen
		qt_rectangle = self.frameGeometry()
		center_point = QDesktopWidget().availableGeometry().center()
		qt_rectangle.moveCenter(center_point)
		# set minimum size to size it takes up
		self.setMinimumSize(self.sizeHint())

	def set_working_dir(self):
		select_dialog = QtWidgets.QFileDialog()
		select_dialog.setAcceptMode(select_dialog.AcceptSave)
		working_directory = select_dialog.getExistingDirectory(directory=config_directory,
		                                                       caption='Select folder to be working directory')
		self.working_directory_lineedit.setText(working_directory)

	def save_settings(self):
		config = {}
		config['working_dir'] = self.working_directory_lineedit.text()
		config['vcf_chunk_size'] = self.vcf_chunk_size.text()
		# config['auto_annotate'] = self.vcf_chunk_size
		config['max_memory'] = self.max_memory_lineedit.text()
		config['genome_version'] = self.genome_version_lineEdit.text()
		config['organism'] = self.organism_lineedit.text().replace(' ', '_')

		# Write settings to YAML file
		with open(config_file, 'w') as configf:
			configf.write(yaml.safe_dump(config, default_flow_style=False))
		self.close()


if __name__ == '__main__':
	MetallaxisApp = QApplication(sys.argv)
	MetallaxisGui = MetallaxisGuiClass()

	# get annotation by default
	MetallaxisGui.MetallaxisSettings.annotation_checkbox.setChecked(True)

	# If config file doesn't exist setup one
	if not os.path.isfile(config_file):
		throw_warning_message("No Metallaxis settings were found, asking user to set them")

		# Open settings window for settings to be changed
		MetallaxisGui.MetallaxisSettings.working_directory_lineedit.setText(config_directory)
		MetallaxisGui.MetallaxisSettings.exec_()

	# Read settings into config dictionary
	config = read_config(config_file)

	# Temporary file names
	global sqlite_output_name, vcf_output_filename
	svg_output_name = os.path.join(config['working_dir'], 'variant_pic.svg')
	sqlite_output_name = os.path.join(config['working_dir'], 'database.sqlite')
	sqlite_connection = sqlite3.connect(sqlite_output_name, isolation_level=None)
	vcf_output_filename = os.path.join(config['working_dir'], 'vcf_output_filename.vcf')
	annotated_vcf_output_filename = os.path.join(config['working_dir'], 'vcf_annot_filename.vcf')

	if len(sys.argv) == 2:
		MetallaxisGui.select_and_parse(sys.argv[1])

	# show GUI
	MetallaxisGui.show()

	# exit program on quitting the GUI
	sys.exit(MetallaxisApp.exec_())
