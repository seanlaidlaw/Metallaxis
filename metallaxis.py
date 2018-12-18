#!/usr/bin/env python3
# -*- coding: utf-8 -*-

## Importer des modules de Python Standard
import sys
import re
import os
from shutil import copyfile  # for save hdf5 to work
import tracemalloc
import json
import requests
# pour gerer les VCF compressé
import lzma
import bz2
import gzip

import pathlib  # for making the folder where we store data
import platform  # for determining OS and therefore where to store data

# Importer les modules de tierce-partie
import magic  # pour detecter type de fichier
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Qt5Agg")

# pour lire uniquement certains lignes des fichiers (reduit conso RAM)
from itertools import islice

# Pour l'interface graphique
from PyQt5 import QtCore, QtGui, QtWidgets, uic
from PyQt5.QtCore import QUrl
from PyQt5.QtGui import QDesktopServices
from PyQt5.QtWidgets import QApplication, QMainWindow, QDialog, QMessageBox, QProgressBar, QTableWidget, QLabel, QDesktopWidget

# for plotting graphs
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt

import time
start_time = time.time()

# allow <Ctrl-c> to terminate the GUI
import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)




# Determine where to store Metallaxis Data
home_dir = os.path.expanduser('~')
if platform.system() == "Linux":
	working_directory = home_dir + "/.metallaxis/"
elif platform.system() == "Darwin":
	working_directory = home_dir + "/Library/Caches/Metallaxis/"
elif platform.system() == "Windows":
	working_directory = os.path.expandvars(r'%APPDATA%\Metallaxis\\')

# make data folder and parent folders if they don't exist
pathlib.Path(working_directory).mkdir(parents=True, exist_ok=True)

# Tempory file names
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
	Generer une dialogue d'avertissement avec l'argument comme message
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
	Generer une dialogue d'alerte avec l'argument comme message d'alerte
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
	Décompresse le fichier d'entrée (si compressé avec with xz/gz/bz2), et
	retourne un objet exploitable: soit le head du fichier, soit le fichier
	entier selon l'option headonly_bool.
	"""
	# la commande open() ouvre juste le fichier, il ne le charge pas en
	# entier ou même de tout en RAM. Ça c'est fait quand on la lis en
	# list plus bas, où on ne lira que les 100 premiers lignes pour eviter
	# de charger de très gros fichiers en RAM
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
	if os.path.isfile(hdf5_filename):
		# global h5_input_name
		# h5_input_name = sys.argv[1]

		MetallaxisGui_object.loaded_vcf_lineedit.setText(os.path.abspath(hdf5_filename))
		try:
			complete_h5_file = pd.read_hdf(hdf5_filename)
		except ValueError:
			complete_h5_file = pd.read_hdf(hdf5_filename, key="df")
		return complete_h5_file

	else:
		# return error if h5 doesn't exist
		throw_error_message("ERROR: Selected file does not \
		exist. You specified : " + str(hdf5_filename))


# Charge l'interface graphique construit en XML
gui_window_object, gui_base_object = uic.loadUiType("MetallaxisGui.ui")

class MetallaxisGui(gui_base_object, gui_window_object):
	"""
	Classe qui construit l'interface graphique Qt sur lequel repose Metallaxis
	"""
	def __init__(self):
		super(gui_base_object, self).__init__()
		self.setupUi(self)
		self.setWindowTitle("Metallaxis")
		# initalise progress bar
		self.MetallaxisProgress = MetallaxisProgress()
		self.MetallaxisProgress.show()

		# Center GUI on screen
		qtRectangle = self.frameGeometry()
		centerPoint = QDesktopWidget().availableGeometry().center()
		qtRectangle.moveCenter(centerPoint)

		# start mesuring memory
		tracemalloc.start()
		global normal_mem
		normal_mem = tracemalloc.take_snapshot()


		self.progress_bar(1,"setting up gui")
		# boutons sur interface
		self.open_vcf_button.clicked.connect(self.select_and_process)
		# menus sur interface
		self.actionOpen_VCF.triggered.connect(self.select_and_process)
		self.actionSave_as_HDF5.triggered.connect(self.save_hdf5)
		self.actionQuit.triggered.connect(self.close)

		# Relier bouton "Github Page" du menu avec l'URL
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
		"""docstring for progress_bar"""
		self.MetallaxisProgress.progressbar_message.setText(message)
		percent = round(percent, 2)
		self.MetallaxisProgress.progressbar_progress.setValue(percent)

		snapshot2 = tracemalloc.take_snapshot()
		mem_usage = sum(stat.size for stat in snapshot2.statistics('lineno'))
		mem_usage_mb = round(mem_usage/1000000, 2)
		mem_usage_mb =("%s" % (mem_usage_mb))

		time_secs = str(round((time.time() - start_time),2))

		self.MetallaxisProgress.progressbar_ram_usage.setText(mem_usage_mb)
		self.MetallaxisProgress.progressbar_time.setText(time_secs)
		print(str(percent) + "% : " + message + " | Time: "+ time_secs + " \
		| RAM (MB): " + mem_usage_mb)
		MetallaxisApp.processEvents()  # refresh the gui


	def empty_qt_layout(self, qt_layout_name):
		"""
		Vide le Qt layout afin qu'on puisse changer de fichier VCF sans
		garder l'information du dernier sur l'interface.

		Accepte le nom d'un Qt layout comme argument.
		"""
		while 1:
			layout_widget = qt_layout_name.takeAt(0)
			if not layout_widget:
				break
			layout_widget.widget().deleteLater()

	def save_hdf5(self):
		"""
		Ouvre une dialogue pour que l'utilisateur puisse choisir un endroit
		ou enregister le fichier hdf5
		"""
		save_dialog = QtWidgets.QFileDialog()
		save_dialog.setAcceptMode(save_dialog.AcceptSave)
		save_folder = save_dialog.getSaveFileName(self, 'Save Analayis as HDF5',filter="*.h5")[0]
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
		Ouvre une dialogue pour que l'utilisateur puisse choisir un fichier
		VCF, et ensuite il le passe à la fonction decompress_vcf()
		"""
		select_dialog = QtWidgets.QFileDialog()
		select_dialog.setAcceptMode(select_dialog.AcceptSave)
		selected_vcf = select_dialog.getOpenFileName(self, filter="VCF Files (*.vcf \
			*.vcf.xz *.vcf.gz *.vcf.bz2) ;;Metallaxis HDF5 Files(*.h5) ;;All Files(*.*)")
		selected_vcf = selected_vcf[0]
		# if the user cancels the select_file() dialog, then run select again
		while selected_vcf == "":
			throw_error_message("ERROR: No selected file")
			return False
		return selected_vcf

	def select_and_process(self):
		selected_file = self.select_file()
		if selected_file is False:
			return
		h5_only = False
		if selected_file.endswith(".h5"):
			complete_h5_file = load_hdf5(selected_file)
			h5_only = True
			self.write_h5_to_interface(complete_h5_file, h5_only, selected_file)
		else:
			selected_vcf = selected_file

		# reopen progress bar for loading new file
		self.MetallaxisProgress = MetallaxisProgress()
		self.MetallaxisProgress.show()

		if not h5_only:
			# get metadata and variant counts from vcf
			metadata_dict, var_counts, decompressed_file = self.process_vcf(selected_vcf)

			# convert vcf into a hdf5 object
			h5_file = self.h5_encode(selected_vcf, decompressed_file, var_counts=var_counts, metadata_dict=metadata_dict)

			# Read H5 for actual table populating
			complete_h5_file = pd.read_hdf(h5_file, key="df")

			# populate interface with information from hdf5
			self.write_h5_to_interface(complete_h5_file, h5_only, h5_file, var_counts=var_counts, metadata_dict=metadata_dict, selected_vcf=selected_vcf)

			# get annotation data
			if self.MetallaxisSettings.annotation_checkbox.isChecked():
				if h5_only is not True:
					try:
						complete_h5_file = self.annotate_h5(complete_h5_file, selected_vcf)
						complete_h5_file = pd.read_hdf(complete_h5_file, key="df")
					except:
						throw_warning_message("Annotation did not succeed, proceeding with non-annotated h5")
						self.MetallaxisSettings.annotation_checkbox.setChecked(False)

		# populate table
		self.populate_table(complete_h5_file)


	def filter_table(self):
		"""
		Selectionne les données correspondant au filtre choisi, et les envoient
		au fonction populate_table pour remplacer les lignes du tableau avec
		uniquement les données qui correspondent à notre filtre.
		"""
		selected_filter = self.filter_box.currentText()
		filter_text = self.filter_lineedit.text()
		# enleve espace blanc de la requete
		filter_text = re.sub(r"\s+", "", filter_text)
		filter_text = filter_text.upper()

		if "-" in filter_text and "," in filter_text:
			throw_error_message("Please only use either comma separated values or a dash separated range")
			return

		elif "-" in filter_text:
			if selected_filter not in numeric_columns:
				throw_error_message("Can only filter a dash separated range on numeric columns: " + str(numeric_columns))
				return

			split_filter_text = filter_text.split("-")
			# filtrer les valeurs nulles ou strings vides pour pas gener le comptage des item
			split_filter_text = filter(None, split_filter_text)
			split_filter_text = list(split_filter_text )
			if len(split_filter_text) == 2:
				self.filter_text.setText("Filtering to show " + selected_filter + " from "+ str(split_filter_text[0]) + " to " + str(split_filter_text[1]))

				if split_filter_text[0] > split_filter_text[1]:
					filter_condition = selected_filter + ">=" + split_filter_text[1] + " & " + selected_filter + "<=" + split_filter_text[0]
				elif split_filter_text[0] < split_filter_text[1]:
					filter_condition = selected_filter + ">=" + split_filter_text[0] + " & " + selected_filter + "<=" + split_filter_text[1]
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
			# filtrer les valeurs nulles ou strings vides pour pas gener le comptage des item
			split_filter_text  = filter(None, split_filter_text)
			split_filter_text  = list(split_filter_text)
			nb_filters = len(split_filter_text)
			if nb_filters >= 1:
				self.filter_text.setText("Filtering to show "  + selected_filter + ": " + str(split_filter_text))
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
			filter_condition = selected_filter+"==\""+filter_text+"\""
			if self.MetallaxisSettings.annotation_checkbox.isChecked():
				filtered_h5_table = pd.read_hdf(annotated_h5_output_name, key="df", where=filter_condition)
			else:
				filtered_h5_table = pd.read_hdf(h5_output_name, key="df", where=filter_condition)

		self.populate_table(filtered_h5_table)


	def verify_file(self,selected_vcf):
		"""
		Verifer si le fichier donnée existe, et n'est pas vide.
		Retourne True si selected_vcf est un fichier valide.
		"""
		self.progress_bar(3,"Verifying VCF: verifying that file is valid")

		# verifier que le fichier existe
		if not os.path.isfile(selected_vcf):
			throw_error_message("ERROR: Selected file does not \
		exist. You specified : " + str(selected_vcf))
			return False

		# verifier que le fichier n'est pas vide
		if not os.path.getsize(selected_vcf) > 0:
			throw_error_message("ERROR: Selected file is empty. \
		You specified : " + str(selected_vcf))
			return False

		# retourne True pour continuer
		return True


	def verify_vcf(self,decompressed_file_head):
		self.progress_bar(8,"Verifying VCF: verifying that VCF is valid")
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

		# also verify no lines start with a # after end of header
		line_num = 0
		variant_num = 0
		global metadata_num, qual_is_numeric
		qual_is_numeric = True
		for line in decompressed_file_head:
			line_num = line_num + 1
			if line.startswith(b'#'):
				# isoler le ligne du header avec les colonnes
				if line.startswith(b'#CHROM'):
					# verify header is conform to vcf 4.1 spec
					# decode byte object to utf-8 string and split by tab
					header_line_cols = line.decode('UTF-8').split("\t")
					# get index of each column
					global chrom_col, id_col, pos_col, ref_col, alt_col, qual_col
					chrom_col = [i for i, s in enumerate(header_line_cols) if '#CHROM' in s][0]
					id_col = [i for i, s in enumerate(header_line_cols) if 'ID' in s][0]
					pos_col = [i for i, s in enumerate(header_line_cols) if 'POS' in s][0]
					ref_col = [i for i, s in enumerate(header_line_cols) if 'REF' in s][0]
					alt_col = [i for i, s in enumerate(header_line_cols) if 'ALT' in s][0]
					qual_col = [i for i, s in enumerate(header_line_cols) if 'QUAL' in s][0]
					# verifier que l'entete contient tous les colonnes obligatoires du VCF 4.1
					if not all(x in header_line_cols  for x in ['#CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']):
						throw_error_message("ERROR: VCF not valid: VCF doesn not contain all required columns")
						return False

			else:
				# if line is not part of header, split by tab to get columns
				# and verify each column for each line is VCFv4.1 conforming
				split_line = line.decode('UTF-8').split("\t")
				for col in split_line:
					col = col.strip()
					# verifier que le pos ne conteint que des chiffres
					if col == split_line[pos_col]:
						if not col.isdigit():
							throw_error_message("ERROR: VCF not valid: column 'POS' doesn't only contain digits: " + str(col))
							return False

					elif col == split_line[ref_col]:
						allowed_chars = set('ACGTN')
						# verifier que REF ne contient pas autre chose que ACGTN
						# (les seules characters authorisés selon VCF 4.1)
						if not set(col.upper()).issubset(allowed_chars):
							throw_error_message("ERROR: VCF not valid: column 'REF' doesn't only contain A,C,G,T,N: " + str(col))
							return False

					elif col == split_line[alt_col]:
						# verifier que ALT ne contient pas autre chose que ACGTN
						# ou un identifiant entre <> (les seules characters
						# authorisés selon VCF 4.1)
						if not col.startswith("<") and col.endswith(">"):
							allowed_chars = set('ACGTN')
							if not set(col).issubset(allowed_chars):
								throw_error_message("ERROR: VCF not valid: column 'ALT' doesn't only contain A,C,G,T,N or <ID>: " + str(col))
								return False

					elif col == split_line[qual_col]:
						# verifier que le QUAL ne contient que un entier
						# ou un float ou un "."
						if col.isdigit():
							break
						elif col == ".":
							qual_is_numeric = False
							break
						else:
							allowed_chars = set('123456789.')
							if not set(col).issubset(allowed_chars):
								try:
									float(col)

								except ValueError:
									throw_error_message("ERROR: VCF not valid: column 'QUAL' doesn't only contain digits: " + str(col))
									return False
				variant_num += 1

			metadata_num = int(line_num - variant_num)

		if variant_num == 0:
			throw_error_message("ERROR: VCF is empty, there are no variants at all in this vcf, please use a different vcf")
			return False
		elif variant_num < 5:
			throw_error_message("ERROR: VCF contains too few variants to analyse, please use a different vcf")
			return False
		elif variant_num > 5 and variant_num < 30:
			throw_warning_message("Warning: VCF contains very few variants, only rudementary statistics can be perfomred")
			return True
		else:
			# si plus de 35 variants retouner sans alerte
			return True



	def process_vcf(self, vcf_input_filename):
		"""
		Effectue les verifications et analyses sur le VCF choisi
		"""
		# verifier que le fichier est valide, on verra s'il est un
		# vcf valide apres décompression
		file_is_valid = self.verify_file(vcf_input_filename)
		if not file_is_valid:
			return


		# TODO: replace arg_file_type with something more intuitive like file_type
		arg_file_type = magic.from_file(vcf_input_filename)
		# Decompresse fichiers selectionées
		if "XZ" in arg_file_type:
			self.detected_filetype_label.setText("xz compressed VCF")
			decompressed_file_head = decompress_vcf("lzma", vcf_input_filename, headonly_bool=True)

		elif "bzip2" in arg_file_type:
			self.detected_filetype_label.setText("bz2 compressed VCF")
			decompressed_file_head  = decompress_vcf("bz2", vcf_input_filename, headonly_bool=True)

		elif "gzip" in arg_file_type:
			self.detected_filetype_label.setText("gz compressed VCF")
			decompressed_file_head  = decompress_vcf("gzip", vcf_input_filename, headonly_bool=True)

		elif "Variant Call Format" in arg_file_type:
			self.detected_filetype_label.setText("uncompressed VCF")
			decompressed_file_head  = decompress_vcf("", vcf_input_filename, headonly_bool=True)

		else:
			throw_error_message("Error: Selected file must be a VCF file")
			return

		# now we have a returned decompressed file object verify if
		# contents are valid vcf
		vcf_is_valid = self.verify_vcf(decompressed_file_head)
		if not vcf_is_valid:
			return


		self.progress_bar(9,"Decompressing VCF")
		# si le vcf est valide alors decompressons tout le fichier
		if "XZ" in arg_file_type:
			decompressed_file = decompress_vcf("lzma", vcf_input_filename, vcf_output_filename=vcf_output_filename)

		elif "bzip2" in arg_file_type:
			decompressed_file = decompress_vcf("bz2", vcf_input_filename, vcf_output_filename=vcf_output_filename)

		elif "gzip" in arg_file_type:
			decompressed_file = decompress_vcf("gzip", vcf_input_filename, vcf_output_filename=vcf_output_filename)

		elif "Variant Call Format" in arg_file_type:
			decompressed_file = decompress_vcf("", vcf_input_filename, vcf_output_filename=vcf_output_filename)

		self.loaded_vcf_lineedit.setText(os.path.abspath(vcf_input_filename))

		# Calculate counts of different variant types
		def add_to_dict_iterator(dictionary, key, iterator_value):
			"""
			appends value to list at dictionary[key], but if key doesn't
			already exist then it adds it
			"""
			if key not in dictionary:
				# add empty list to dict as key
				dictionary[key] = 0
				dictionary[key] = dictionary[key] + iterator_value
			else:
				# append value to list
				dictionary[key] = dictionary[key] + iterator_value


		var_counts = {}
		length_of_all_indels = 0
		var_counts["Total_SNP_Count"] = 0
		var_counts["Total_Indel_Count"] = 0
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
					# if our alt types are only SNP or transposable élements
					# then count them
					if alt_types_only_snp is True:
						ALT_Types.add(alt)
						add_to_dict_iterator(
							var_counts, line[alt_col]+"_Alt_Count", 1)

					if len(line[ref_col]) == len(line[alt_col]):
						var_counts["Total_SNP_Count"] += 1
						add_to_dict_iterator(
							var_counts, line[chrom_col]+"_Chrom_SNP_Count", 1)
						add_to_dict_iterator(
							var_counts, line[chrom_col]+"_Chrom_Variant_Count", 1)
					else:
						var_counts["Total_Indel_Count"] += 1
						add_to_dict_iterator(
							var_counts, line[chrom_col]+"_Chrom_Indel_Count", 1)
						add_to_dict_iterator(
							var_counts, line[chrom_col]+"_Chrom_Variant_Count", 1)
						length_of_all_indels += len(line[alt_col])


		total_chrom_snp_count, total_chrom_indel_count = 0, 0
		for key, value in var_counts.items():
			if "_Chrom_SNP_Count" in key:
				total_chrom_snp_count += value
			if "_Chrom_Indel_Count" in key:
				total_chrom_indel_count += value


		if length_of_all_indels > 0 and var_counts["Total_Indel_Count"] > 0:
			var_counts["Avg_Indel_Length"] = float(
				length_of_all_indels / var_counts["Total_Indel_Count"])
			var_counts["Avg_Indel_Length"] = round(var_counts["Avg_Indel_Length"], 3)
		var_counts["Avg_SNP_per_Chrom"] = int(
			total_chrom_snp_count / len(list_chromosomes))
		var_counts["Avg_Indel_per_Chrom"] = int(
			total_chrom_indel_count / len(list_chromosomes))
		var_counts["Avg_Variant_per_Chrom"] = int(
			(total_chrom_snp_count + total_chrom_indel_count) / len(list_chromosomes))


		if ALT_Types != set():
			var_counts["ALT_Types"] = str(ALT_Types).strip("{").strip("}")


		# Extract Metadata from VCF
		metadata_dict = {}
		# Match des groupes des deux cotés du "=", apres un "##"
		regex_metadata = re.compile('(?<=##)(.*?)=(.*$)')

		self.progress_bar(10,"Extracting VCF metadata")
		# parser vcf decompressé dans plusieurs dictionnaires
		with open(vcf_output_filename) as decompressed_out:
			vcf_line_nb, metadata_line_nb = 0, 0
			for line in decompressed_out:
				if line.startswith('##'):
					metadata_tag = str(regex_metadata.search(line).group(1))
					metadata_result = str(regex_metadata.search(line).group(2))
					# dans un premier temps on va pas s'interesser pas aux
					# metadonéées 'en majiscules' ("INFO" "FILTER" "ALT")
					# peuvent être regroupés ensemble dans un tableau
					if metadata_tag.isupper():
						metadata_type = metadata_tag
					else:
						metadata_type = "basic"
						# too long metadata distorts the GUI window and causes
						# database errors so truncate if its too long
						if len(metadata_tag) > 20:
							metadata_tag = metadata_tag[:20] + "..."
						if len(metadata_result) > 95:
							metadata_result = metadata_result[:95] + "...<truncated due to length>"

					metadata_dict_entry = [metadata_type, metadata_tag, metadata_result]
					if not metadata_dict_entry in metadata_dict.values():
						metadata_dict[metadata_line_nb] = metadata_dict_entry
					metadata_line_nb += 1

		return (metadata_dict, var_counts, decompressed_file)


	def h5_encode(self, selected_vcf, decompressed_file=None, var_counts=None, metadata_dict=None):
		"""
		Lis en entier le vcf selectionné, bloc par bloc et le met dans un
		fichier HDF5.
		"""
		h5_file = pd.HDFStore(h5_output_name ,mode='w')

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

		chunked_vcf_len = sum(1 for row in open(decompressed_file, 'r'))
		chunked_vcf = pd.read_csv(selected_vcf,
		                          delim_whitespace=True,
		                          skiprows=range(0, metadata_num - 1),
		                          chunksize=int(h5chunksize),
		                          low_memory=False,
		                          # make default data type an object (ie. string)
		                          dtype=object)
		annotate_nb = 0
		annotate_percent = 35

		# PARSE INFO COLUMNS
		# the info column of a vcf is long and hard to read if
		# displayed as is, but it is composed of multiple key:value tags
		# that we can parse as new columns, making more visual sense
		global info_cols_to_add
		info_cols_to_add = set()
		chunked_vcf = pd.read_csv(selected_vcf,
		                          delim_whitespace=True,
		                          skiprows=range(0, metadata_num - 1),
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

			# set global variable table_column_names to have both normal table
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
				self.progress_bar(
					float(annotate_progress), "Encoding H5 database")

			# append the modified chunk to the h5 database without indexing
			h5_file.append("df", chunk, index=False, data_columns=True,
			               min_itemsize=80, complib=complib, complevel=int(complevel))
			annotate_nb += 1

		# index columns explicitly, now that we have finished adding data to it
		self.progress_bar(46, "Indexing H5 database")
		h5_file.create_table_index(
			"df", columns=table_column_names, optlevel=9, kind='full')
		h5_file.close()
		return h5_output_name


	def annotate_h5(self, selected_h5_data, selected_vcf):
		"""
		Uses EBI's VEP API to annotate lines where ID is valid
		"""
		ebi_rest_api = "https://rest.ensembl.org"
		ext = "/vep/human/id/"
		headers = {"Content-Type": "application/json", "Accept": "application/json"}

		api_ids = []
		selected_h5_data = pd.read_hdf(h5_output_name, key="df", where="ID!='.'")

		for line in selected_h5_data.itertuples():
			# turn tuple into a list, but exclude the first item of list
			# because its the h5 index and not part of our original data
			line = list(line)[1:]
			# get the id from id_column and add it to api_ids list
			line_id = line[id_col]
			# ids are strings so if we ave a valid non-null string add it as an ID to send to EBI
			if type(line_id) == str:
				if line_id != "." and line_id != "nan":
					api_ids.append(line_id)

		if len(api_ids) == 0:
			# return an error if we have no IDs so that metallaxis can use the non-annotated h5 instead
			throw_error_message("No valid ids to annotate in vcf")
			raise IndexError

		annotated_h5_tmp = pd.HDFStore(annotated_h5_tmp_name, mode='w')
		annotated_row_ids = []

		annotate_nb = 0
		# divide api_ids list into sublists, and make API calls based on that to avoid reaching the API's max length
		api_call_number = 1
		api_max_at_once = 50
		for i in range(0, len(api_ids), api_max_at_once):
			self.progress_bar(51, "Annotate H5: Making API Call to EBI (" + str(api_call_number) +
			                  "/" + str(int(round((len(api_ids)/api_max_at_once),0))) + ")")
			api_ids_sublist = api_ids[i:i + api_max_at_once]

			formatted_api_ids_sublist = ','.join('"{0}"'.format(ind_id) for ind_id in api_ids_sublist)
			formatted_api_ids_sublist = '{ "ids" : [' + formatted_api_ids_sublist + '] }'
			api_call = requests.post(ebi_rest_api + ext, headers=headers, data=formatted_api_ids_sublist)

			if not api_call.ok:
				throw_error_message("ERROR: API call failed")
				api_call.raise_for_status()
			# vep_json = json.loads(api_call.text)

			data = api_call.text
			# save sublist's API response as a JSON
			json_output = os.path.join(working_directory, 'VEP_API_' + str(api_call_number) + '.json')
			with open(json_output, 'w') as json_file:
				json_file.write(data)

			api_call_number += 1

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
		self.progress_bar(53, "Annotate H5: Parsing downloaded JSONs for annotations")
		annotation_columns = set()
		# loop through all the saved JSONs
		for json_number in range(1, api_call_number):
			json_file = os.path.join(working_directory, 'VEP_API_' + str(json_number) + '.json')

			# open the json file as a JSON object
			with open(json_file) as json_input:
				json_raw_data = json.load(json_input)

			# loop through all the dictionaries embeded in the JSON object
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

		self.progress_bar(54, "Annotate H5: writing annotated chunk to temporary h5")
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
							annotate_len = (len(api_ids) * len(data))
							annotate_progress = annotate_percent + (annotate_nb / annotate_len) *  (70 - annotate_percent)
							self.progress_bar(float(annotate_progress), "Annotate H5: writing annotated chunk to temporary h5")

						# get h5 row from id and add JSON data to it and add it to a different h5
						api_id = data['id']
						# see if using input instead of id makes a diff
						api_alt = subdata['variant_allele']

						# get unannotated row with id that match annotated row and have same ALT
						selected_h5_row = pd.read_hdf(h5_output_name, key="df", where=(
								"ID=='" + api_id + "' and ALT=='" + api_alt + "'"))

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

							# initalise all the found columns
							for column in annotation_columns:
								selected_h5_row[column] = "."

							# run the function to replace the default "." with the real value if it exists
							get_annotation_for_columns(columns_to_annotate)

							# Convert all columns to strings
							for column in selected_h5_row:
								selected_h5_row[column] = selected_h5_row[column].apply(str)

							# Append the row to h5 database
							annotated_h5_tmp.append("df", selected_h5_row, index=False, data_columns=True,
							                        min_itemsize=100, complib=complib, complevel=int(complevel))
							annotated_row_ids.append(api_id)
							annotate_nb += 1


		# read through vcf again, and load rows with no annotation-friendly IDs to a new h5 database
		annotate_nb = 0
		annotate_percent = 70
		chunked_vcf_len = sum(1 for row in open(selected_vcf, 'r'))
		self.progress_bar(70, "Annotate H5: writing non-annotated chunk to temporary h5")
		for id in api_ids:
			if id not in annotated_row_ids:
				non_annotated_h5_row = pd.read_hdf(h5_output_name, key="df", where=(
						"ID=='" + id + "'"))
				if not non_annotated_h5_row.empty:
					for col in annotation_columns:
						non_annotated_h5_row[col] = "."

					for column in non_annotated_h5_row:
						non_annotated_h5_row[column] = non_annotated_h5_row[column].apply(str)

					print("getting row")
					annotated_h5_tmp.append("df", non_annotated_h5_row, index=False, data_columns=True,
					                        min_itemsize=80, complib=complib, complevel=int(complevel))
					print(non_annotated_h5_row)
					print("got row")

					# only update progress bar every 20 lines to avoid
					# performance hit from doing it every line
					if annotate_nb % 20 == 0:
						annotate_progress = annotate_percent + (annotate_nb / chunked_vcf_len) * (80 - annotate_percent)
						self.progress_bar(float(annotate_progress), "Annotate H5: writing non-annotated chunk to temporary h5")
					annotate_nb += 1

		annotated_h5_tmp.close()

		print("tmp file closed")

		# Read the just written h5 file replacing all "." with NaN which is read as NULL instead of a string, convert
		# columns containing only numbers to numeric datatype, and remove columns with no data for any row
		annotated_numeric_columns = list(numeric_columns) + list(annotation_columns)
		print(annotated_numeric_columns)

		print("start read")
		annotated_tmp_read_obj = pd.read_hdf(annotated_h5_tmp_name, key="df")
		annotated_tmp_read_obj = annotated_tmp_read_obj.replace(".", np.NaN)
		print("replaced empties")

		for column in annotated_tmp_read_obj.keys():
			set_col_to_numeric_if_isdigit(column, annotated_tmp_read_obj, annotated_numeric_columns)

		for column in annotated_numeric_columns:
			annotated_tmp_read_obj[column] = pd.to_numeric(annotated_tmp_read_obj[column])

		# remove columns that only contained null values
		annotated_tmp_read_obj = annotated_tmp_read_obj.dropna(axis=1, how='all')

		# save appended dataframe to new h5 file and remove the temporary one
		annotated_tmp_read_obj.to_hdf(annotated_h5_output_name, key='df', mode='w', format='table', data_columns=True, complib=complib, complevel=int(complevel))


		if os.path.exists(annotated_h5_tmp_name):
			os.remove(annotated_h5_tmp_name)

		# embed metadat and statistics in new H5
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



	def write_h5_to_interface(self, complete_h5_file, h5_only, h5_input_name, var_counts=None, metadata_dict=None, selected_vcf=None):
		"""
		function that clears the interface if it already has data,
		then runs the annotate_h5() populate_table() functions. Requires h5
		read object, h5 only boolean and optionally takes selected_vcf string
		if annotation is chosen.
		"""
		# active les widgets qui sont desactivés tant qu'on a pas de VCF selectioné
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
		self.progress_bar(47,"Extracting column data")
		print( complete_h5_file )
		print(complete_h5_file.keys())

		column_names = list(complete_h5_file.keys())
		global chrom_col, id_col, pos_col, ref_col, alt_col, qual_col
		chrom_col = [i for i, s in enumerate(column_names) if 'CHROM' in s][0]
		id_col = [i for i, s in enumerate(column_names) if 'ID' in s][0]
		pos_col = [i for i, s in enumerate(column_names) if 'POS' in s][0]
		ref_col = [i for i, s in enumerate(column_names) if 'REF' in s][0]
		alt_col = [i for i, s in enumerate(column_names) if 'ALT' in s][0]
		qual_col = [i for i, s in enumerate(column_names) if 'QUAL' in s][0]

		# effacer espace metadonées (utile si on charge un fichier apres un autre)
		self.empty_qt_layout(self.dynamic_metadata_label_results)
		self.empty_qt_layout(self.dynamic_metadata_label_tags)
		# effacer aussi espace statistiques
		self.empty_qt_layout(self.dynamic_stats_value_label)
		self.empty_qt_layout(self.dynamic_stats_key_label)

		self.filter_table_btn.clicked.connect(self.filter_table)

		if not h5_only:
			for metadata_line_nb in metadata_dict:
				metadata_tag = metadata_dict[metadata_line_nb][1]
				metadata_result = metadata_dict[metadata_line_nb][2]
				if not metadata_tag.isupper():
					# Generer dynamiquement du texte pour le titre et resultat pour
					# chaque type de metadonnée non-majiscule
					self.dynamic_metadata_label_tags.addWidget(
						QtWidgets.QLabel(metadata_tag, self))
					self.dynamic_metadata_label_results.addWidget(
						QtWidgets.QLabel(metadata_result, self))

			for key, value in var_counts.items():
				key=key.replace("_"," ")
				key=key+":"
				self.dynamic_stats_key_label.addWidget(
					QtWidgets.QLabel(str(key), self))
				self.dynamic_stats_value_label.addWidget(
					QtWidgets.QLabel(str(value), self))

		else:
			for line in pd.read_hdf(h5_input_name , key="metadata").itertuples():
				# turn tuple into a list, but exclude the first item of list
				# because its the h5 index and not part of our original data
				line = list(line)[1:]
				metadata_tag = line[0]
				metadata_result = line[1]
				self.dynamic_metadata_label_tags.addWidget(
					QtWidgets.QLabel(metadata_tag, self))
				self.dynamic_metadata_label_results.addWidget(
					QtWidgets.QLabel(metadata_result, self))

			for line in pd.read_hdf(h5_input_name , key="stats").itertuples():
				var_counts_key = line[0]
				var_counts_value = line[1]
				self.dynamic_stats_key_label.addWidget(
					QtWidgets.QLabel(str(var_counts_key), self))
				self.dynamic_stats_value_label.addWidget(
					QtWidgets.QLabel(str(var_counts_value), self))


		self.progress_bar(49, "Plotting Statistics")

		# clear plot layout to avoid dupication on file reload
		self.empty_qt_layout(self.stat_plot_layout)

		# plot piechart of proportions of SNP/Indel
		total_figure = plt.figure()
		graph = total_figure.add_subplot(111)
		graph.pie([var_counts['Total_SNP_Count'], var_counts['Total_Indel_Count']], labels=['SNP', 'Indels'], autopct='%1.1f%%')
		# installer les axes de x et y sont egals pour assurer le rond
		graph.axis('equal')
		plt.title('Proportion of SNP/Indels')
		total_figure.tight_layout()
		graph.legend()
		self.stat_plot_layout.addWidget(FigureCanvas(total_figure))

		# plot piechart of proportions of types of ALT
		# get the value for each ALT_Types key in order, per type of Alt so it can be graphed
		alt_values_to_plot = []
		for alt in ALT_Types:
			dict_key = alt + "_Alt_Count"
			alt_values_to_plot.append(var_counts[dict_key])

		total_figure = plt.figure()
		graph = total_figure.add_subplot(111)
		graph.pie(alt_values_to_plot, labels=ALT_Types, autopct='%1.1f%%')
		# installer les axes de x et y sont egals pour assurer le rond
		graph.axis('equal')
		plt.title('Proportion of different mutations')
		total_figure.tight_layout()
		graph.legend()
		self.stat_plot_layout.addWidget(FigureCanvas(total_figure))

		# plot piechart of proportions of types of ALT
		# get the nb of mutations for each chromosome
		values_to_plot = []
		global list_chromosomes  # we're editing a global so it needs to be declared global again
		list_chromosomes = sorted(list_chromosomes, key=str)
		for chrom in list_chromosomes:
			dict_key = chrom + "_Chrom_Variant_Count"
			values_to_plot.append(var_counts[dict_key])

		total_figure = plt.figure()
		graph = total_figure.add_subplot(111)
		graph.bar(list(list_chromosomes), values_to_plot)
		plt.title('Distribution of mutations by Chromosome')
		plt.xlabel('Chromosome')
		plt.ylabel('Number of Variants')
		total_figure.tight_layout()
		graph.legend()
		self.stat_plot_layout.addWidget(FigureCanvas(total_figure))




	def populate_table(self, selected_h5_data):
		# effacer tableau actuelle
		self.viewer_tab_table_widget.setRowCount(0)
		self.viewer_tab_table_widget.setColumnCount(0)

		# Creer Tableau vide avec bonne nombre de colonnes et lignes
		table_length = selected_h5_data.shape[0]
		table_width = selected_h5_data.shape[1]

		self.viewer_tab_table_widget.setRowCount(table_length)
		self.viewer_tab_table_widget.setColumnCount(table_width)

		global column_names
		column_names = list(selected_h5_data.keys())

		# effacer chrom_filter_box
		self.filter_box.clear()
		self.filter_text.setText(" ")

		# set filter_box to list column_names
		self.filter_box.addItems(column_names)
		# read header labels from the dataframe keys
		self.viewer_tab_table_widget.setHorizontalHeaderLabels(column_names)

		# Remplir Tableau
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
		self.progress_bar(100,"Populating Table...done")
		# close progress bar when file is completely loaded
		self.MetallaxisProgress.close()




progress_window_object, progress_base_object = uic.loadUiType("MetallaxisProgress.ui")
class MetallaxisProgress(progress_base_object, progress_window_object):
	"""
	Status bar that shows progress when opening files with Metallaxis
	"""
	def __init__(self):
		super(progress_base_object, self).__init__()
		self.setupUi(self)
		self.setWindowTitle("Metallaxis")


settings_window_object, settings_base_object = uic.loadUiType("MetallaxisSettings.ui")
class MetallaxisSettings(settings_base_object, settings_window_object):
	"""
	La fenetre de paramètres qui permet de modifier les paramètres de defaut
	"""
	def __init__(self):
		super(settings_base_object, self).__init__()
		self.setupUi(self)
		self.setWindowTitle("Metallaxis Settings")
		self.compression_comboBox.addItems(['zlib', 'blosc', 'bzip2', 'lzo'])
		self.vcf_chunk_size.setText("5000")
		# Center settings pannel on screen
		qtRectangle = self.frameGeometry()
		centerPoint = QDesktopWidget().availableGeometry().center()
		qtRectangle.moveCenter(centerPoint)


if __name__ == '__main__':
	MetallaxisApp = QApplication(sys.argv)
	MetallaxisGui_object = MetallaxisGui()

	# obtenir vcf d'entrée
	MetallaxisGui_object.progress_bar(2, "Parsing arguments")
	global h5_only
	h5_only = False

	# if we load a h5 file from a previous analysis then we skip the usual analysis and verifications
	if len(sys.argv) == 2 and sys.argv[1].endswith(".h5"):
		h5_only = True
		complete_h5_file = load_hdf5(sys.argv[1])
		MetallaxisGui_object.write_h5_to_interface(complete_h5_file, h5_only, sys.argv[1])

	elif len(sys.argv) == 2:  # if we give it vcf or compressed vcf
		selected_vcf = os.path.abspath(sys.argv[1])

	elif len(sys.argv) == 1:  # if we don't give any args then open file picker
		selected_file = MetallaxisGui_object.select_file()
		if selected_file is False:
			throw_error_message("ERROR: No file selected, qutting Metallaxis")
			sys.exit(1)

		elif selected_file.endswith(".h5"):
			complete_h5_file = load_hdf5(selected_file)
			h5_only = True
			MetallaxisGui_object.write_h5_to_interface(complete_h5_file, h5_only, selected_file)
		else:
			selected_vcf = selected_file

	else:  # if we give more than 1 arg
		print("Error: Metallaxis can only take one argument, a vcf file")
		exit(1)

	# if we loaded a VCF file, do the verification and analysis
	if not h5_only:

		# get metadata and variant counts from vcf
		metadata_dict, var_counts, decompressed_file = MetallaxisGui_object.process_vcf(selected_vcf)

		# convert vcf into a hdf5 object
		h5_file = MetallaxisGui_object.h5_encode(selected_vcf, decompressed_file,  var_counts=var_counts, metadata_dict=metadata_dict)

		# Read H5 for actual table populating
		complete_h5_file = pd.read_hdf(h5_file, key="df")

		# populate interface with information from hdf5
		MetallaxisGui_object.write_h5_to_interface(complete_h5_file, h5_only, h5_file, var_counts=var_counts, metadata_dict=metadata_dict, selected_vcf=selected_vcf)


		# get annotation data
		if MetallaxisGui_object.MetallaxisSettings.annotation_checkbox.isChecked():
			if h5_only is not True:
				try:
					complete_h5_file = MetallaxisGui_object.annotate_h5(complete_h5_file, selected_vcf)
					complete_h5_file = pd.read_hdf(complete_h5_file, key="df")
				except:
					throw_warning_message("Annotation did not succeed, proceeding with non-annotated h5")
					MetallaxisGui_object.MetallaxisSettings.annotation_checkbox.setChecked(False)


	# actions that are to be done once we have our h5 file (if we loaded a h5 file then it starts here)
	# populate table
	MetallaxisGui_object.populate_table(complete_h5_file)

	# show GUI
	MetallaxisGui_object.show()

	# exit program on quitting the GUI
	sys.exit(MetallaxisApp.exec_())
