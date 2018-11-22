#!/usr/bin/env python3
# -*- coding: utf-8 -*-

## Importer des modules de Python Standard
import sys
import re
import os
# pour gerer les VCF compressé
import lzma
import bz2
import gzip

# Importer les modules de tierce-partie
import magic  # pour detecter type de fichier
import pandas as pd
# pour lire uniquement certains lignes des fichiers (reduit conso RAM)
from itertools import islice
from PyQt5 import QtCore, QtGui, QtWidgets, uic
from PyQt5.QtWidgets import QApplication, QMainWindow, QDialog, QMessageBox, QTableWidget, QLabel
from PyQt5.QtGui import QDesktopServices
from PyQt5.QtCore import QUrl

def decompress_vcf(compression_type, selected_vcf, headonly):
	"""
	Décompresse le fichier d'entrée (si compressé avec with xz/gz/bz2), et
	retourne un objet exploitable: soit le head du fichier, soit le fichier
	entier selon l'option headonly.
	"""
	# la commande open() ouvre juste le fichier, il ne le charge pas en
	# entier ou même de tout en RAM. Ça c'est fait quand on la lis en
	# list plus bas, où on ne lira que les 100 premiers lignes pour eviter
	# de charger de très gros fichiers en RAM
	if compression_type == "":
		decompressed_file_object = open(selected_vcf, mode="rb")
	else:
		decompressed_file_object = eval(compression_type).open(selected_vcf, mode="rb")

	if headonly is True:
		decompressed_file_head = list(islice(decompressed_file_object, 100))
		decompressed_file_object.close()
		return decompressed_file_head
	else:
		with open("decompressed_vcf_output.vcf", "wb") as decompressed_out:
			decompressed_out.write(decompressed_file_object.read())
		decompressed_file_object.close()
		return decompressed_file_object



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
		# boutons sur interface
		self.open_vcf_button.clicked.connect(self.select_and_process)
		# menus sur interface
		self.actionOpen_VCF.triggered.connect(self.select_and_process)
		self.actionQuit.triggered.connect(self.close)

		# Relier bouton "Github Page" du menu avec l'URL
		def open_github(url):
			url = "https://github.com/SL-LAIDLAW/Metallaxis"
			QDesktopServices.openUrl(QtCore.QUrl(url))

		self.actionGithub_Page.triggered.connect(open_github)

		def open_about_tab():
			self.tabWidget.setCurrentIndex(4)
		# set first tab as default
		self.tabWidget.setCurrentIndex(0)
		# on click "About", open about tab
		self.actionAbout.triggered.connect(open_about_tab)

		self.MetallaxisSettings = MetallaxisSettings()
		def show_settings_window():
			self.MetallaxisSettings.show()

		self.actionSettings.triggered.connect(show_settings_window)

		# selectionner vcf d'entrée si c'est pas fourni
		if len(sys.argv) == 1:
			selected_vcf = self.select_and_process()
		elif len(sys.argv) == 2:
			# obtenir le chemin absolue afin d'être dans les memes conditions
			# que si on le selectionnait
			selected_vcf = os.path.abspath(sys.argv[1])
			self.process_vcf(selected_vcf)
		else:
			print("Error: Metallaxis can only take one argument, a vcf file")
			exit(1)


	def throw_warning_message(self, warning_message):
		"""
		Generer une dialogue d'avertissement avec l'argument comme message
		"""
		print("Warning: " + warning_message)
		warning_dialog = QtWidgets.QMessageBox()
		warning_dialog.setIcon(QMessageBox.Warning)
		warning_dialog.setWindowTitle("Warning:")
		warning_dialog.setText(warning_message)
		warning_dialog.setStandardButtons(QMessageBox.Ok)
		warning_dialog.exec_()

	def throw_error_message(self, error_message):
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
			self.throw_error_message("Please only use either comma separated values or a dash separated range")
			return

		# TODO: make this functionnable for int columns (like pos)
		elif "-" in filter_text:
			split_filter_text = filter_text.split("-")
			# filtrer les valeurs nulles ou strings vides pour pas gener le comptage des item
			split_filter_text   = filter(None, split_filter_text)
			split_filter_text   = list(split_filter_text )
			if len(split_filter_text) == 2:
				self.filter_text.setText("Filtering to show " + selected_filter + " from "+ str(split_filter_text[0]) + " to " + str(split_filter_text[1]))
				self.throw_error_message("Sorry, dash separated filters aren't currently working")
				return
			else:
				self.filter_text.setText(" ")
				self.throw_error_message("Please only enter 2 values separated by a dash")
				return

		elif "," in filter_text:
			split_filter_text = filter_text.split(",")
			# filtrer les valeurs nulles ou strings vides pour pas gener le comptage des item
			split_filter_text  = filter(None, split_filter_text)
			split_filter_text  = list(split_filter_text)
			nb_filters = len(split_filter_text)
			if nb_filters >= 2:
				self.filter_text.setText("Filtering to show "  + selected_filter + ": " + str(split_filter_text))
				filter_condition = ""
				filter_iterator = 1
				for each_filter in split_filter_text:
					comma_filter = selected_filter + "=='" + each_filter + "'"
					if filter_iterator < nb_filters:
						comma_filter = comma_filter  + " or "
					filter_condition = filter_condition + comma_filter
					filter_iterator += 1

				filtered_h5_table = pd.read_hdf('input.h5', where=filter_condition)
			else:
				self.filter_text.setText(" ")
				self.throw_error_message("Please enter 2 or more values separated by a comma")
				return

		elif filter_text == "":
			self.filter_text.setText("No Filter Selected")
			filtered_h5_table = pd.read_hdf('input.h5')
		else:
			self.filter_text.setText("Filtering to show " + selected_filter + ": " + str(filter_text))
			filter_condition = selected_filter+"==\""+filter_text+"\""
			filtered_h5_table = pd.read_hdf('input.h5', where=filter_condition)

		self.populate_table(filtered_h5_table)


	def verify_file(self,selected_vcf):
		"""
		Verifer si le fichier donnée existe, et n'est pas vide.
		Retourne True si selected_vcf est un fichier valide.
		"""

		# verifier que le fichier existe
		if not os.path.isfile(selected_vcf):
			self.throw_error_message("ERROR: Selected file does not \
		exist. You specified : " + str(selected_vcf))
			return False

		# verifier que le fichier n'est pas vide
		if not os.path.getsize(selected_vcf) > 0:
			self.throw_error_message("ERROR: Selected file is empty. \
		You specified : " + str(selected_vcf))
			return False

		# retourne True pour continuer
		return True


	def verify_vcf(self,decompressed_file_head):
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
		global metadata_num
		for line in decompressed_file_head:
			line_num = line_num + 1
			if line.startswith(b'#'):
				# isoler le ligne du header avec les colonnes
				if line.startswith(b'#CHROM'):
					# verify header is conform to vcf 4.1 spec
					# decode byte object to utf-8 string and split by tab
					header_line_cols = line.decode('UTF-8').split("\t")
					# get index of each column
					pos_col = [i for i, s in enumerate(header_line_cols) if 'POS' in s][0]
					ref_col = [i for i, s in enumerate(header_line_cols) if 'REF' in s][0]
					alt_col = [i for i, s in enumerate(header_line_cols) if 'ALT' in s][0]
					qual_col = [i for i, s in enumerate(header_line_cols) if 'QUAL' in s][0]
					# verifier que l'entete contient tous les colonnes obligatoires du VCF 4.1
					if not all(x in header_line_cols  for x in ['#CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']):
						self.throw_error_message("ERROR: VCF not valid: VCF doesn not contain all required columns")
						return False

			else:
				split_line = line.decode('UTF-8').split("\t")
				for col in split_line:
					col = col.strip()
					# verifier que le pos ne conteint que des chiffres
					if col == split_line[pos_col]:
						# allowed_chars = set('123456789')
						# if not set(col).issubset(allowed_chars):
						if not col.isdigit():
							self.throw_error_message("ERROR: VCF not valid: column 'POS' doesn't only contain digits: " + str(col))
							return False

					elif col == split_line[ref_col]:
						allowed_chars = set('ACGTN')
						# verifier que REF ne contient pas autre chose que ACGTN
						# (les seules characters authorisés selon VCF 4.1)
						if not set(col.upper()).issubset(allowed_chars):
							self.throw_error_message("ERROR: VCF not valid: column 'REF' doesn't only contain A,C,G,T,N: " + str(col))
							return False

					elif col == split_line[alt_col]:
						# verifier que ALT ne contient pas autre chose que ACGTN
						# ou un identifiant entre <> (les seules characters
						# authorisés selon VCF 4.1)
						if not col.startswith("<") and col.endswith(">"):
							allowed_chars = set('ACGTN')
							if not set(col).issubset(allowed_chars):
								self.throw_error_message("ERROR: VCF not valid: column 'ALT' doesn't only contain A,C,G,T,N or <ID>: " + str(col))
								return False

					elif col == split_line[qual_col]:
						# verifier que le QUAL ne contient que un entier
						# ou un float ou un "."
						if col.isdigit():
							return True
						elif col == ".":
							return True
						else:
							allowed_chars = set('123456789.')
							if not set(col).issubset(allowed_chars):
								try:
									float(col)
									return True
								except ValueError:
									self.throw_error_message("ERROR: VCF not valid: column 'QUAL' doesn't only contain digits: " + str(col))
									return False
								return False
				variant_num += 1

			metadata_num = int(line_num - variant_num)

		if variant_num == 0:
			self.throw_error_message("ERROR: VCF is empty, there are no variants at all in this vcf, please use a different vcf")
			return False
		elif variant_num < 5:
			self.throw_error_message("ERROR: VCF contains too few variants to analyse, please use a different vcf")
			return False
		elif variant_num > 5 and variant_num < 30:
			self.throw_warning_message("Warning: VCF contains very few variants, only rudementary statistics can be perfomred")
			return True
		else:
			# si plus de 35 variants retouner sans alerte
			return True


	def select_vcf(self):
		"""
		Ouvre une dialogue pour que l'utilisateur puisse choisir un fichier
		VCF, et ensuite il le passe à la fonction decompress_vcf()
		"""
		select_dialog = QtWidgets.QFileDialog()
		select_dialog.setAcceptMode(select_dialog.AcceptSave)
		selected_vcf = select_dialog.getOpenFileName(self, filter="VCF Files (*.vcf \
			*.vcf.xz *.vcf.gz *.vcf.bz2) ;;All Files(*.*)")
		selected_vcf = selected_vcf[0]
		return selected_vcf


	def process_vcf(self, selected_vcf):
		"""
		Effectue les verifications et analyses sur le VCF choisi
		"""
		# verifier que le fichier est valide, on verra s'il est un
		# vcf valide apres décompression
		file_is_valid = self.verify_file(selected_vcf)
		if not file_is_valid:
			return


		# TODO: replace arg_file_type with something more intuitive like file_type
		arg_file_type = magic.from_file(selected_vcf)
		# Decompresse fichiers selectionées
		if "XZ" in arg_file_type:
			self.detected_filetype_label.setText("xz compressed VCF")
			decompressed_file_head = decompress_vcf("lzma", selected_vcf, headonly=True)

		elif "bzip2" in arg_file_type:
			self.detected_filetype_label.setText("bz2 compressed VCF")
			decompressed_file_head  = decompress_vcf("bz2", selected_vcf, headonly=True)

		elif "gzip" in arg_file_type:
			self.detected_filetype_label.setText("gz compressed VCF")
			decompressed_file_head  = decompress_vcf("gzip", selected_vcf, headonly=True)

		elif "Variant Call Format" in arg_file_type:
			self.detected_filetype_label.setText("uncompressed VCF")
			decompressed_file_head  = decompress_vcf("", selected_vcf, headonly=True)
		else:
			self.throw_error_message("Error: Selected file must be a VCF file")
			return

		# now we have a returned decompressed file object verify if
		# contents are valid vcf
		vcf_is_valid = self.verify_vcf(decompressed_file_head)
		if not vcf_is_valid:
			return


		# si le vcf est valide alors decompressons tout le fichier
		if "XZ" in arg_file_type:
			decompressed_file= decompress_vcf("lzma", selected_vcf, headonly=False)

		elif "bzip2" in arg_file_type:
			decompressed_file = decompress_vcf("bz2", selected_vcf, headonly=False)

		elif "gzip" in arg_file_type:
			decompressed_file = decompress_vcf("gzip", selected_vcf, headonly=False)

		elif "Variant Call Format" in arg_file_type:
			decompressed_file = decompress_vcf("", selected_vcf, headonly=False)

		h5_file = self.h5_encode(selected_vcf)


		# active les widgets qui sont desactivés tant qu'on a pas de VCF selectioné
		self.loaded_vcf_lineedit.setText(selected_vcf)
		self.loaded_vcf_lineedit.setEnabled(True)
		self.loaded_vcf_label.setEnabled(True)
		self.meta_detected_filetype_label.setEnabled(True)
		self.metadata_area_label.setEnabled(True)
		self.viewer_tab_table_widget.setEnabled(True)
		self.filter_table_btn.setEnabled(True)
		self.filter_label.setEnabled(True)
		self.filter_lineedit.setEnabled(True)
		self.filter_box.setEnabled(True)

		# effacer espace metadonées (utile si on charge un fichier apres un autre)
		self.empty_qt_layout(self.dynamic_metadata_label_results)
		self.empty_qt_layout(self.dynamic_metadata_label_tags)
		self.viewer_tab_table_widget.setRowCount(0)  # supprime tout les lignes
		# effacer chrom_filter_box
		self.filter_box.clear()
		self.filter_text.setText(" ")


		metadata_dict = {}
		# Match des groupes des deux cotés du "=", apres un "##"
		regex_metadata = re.compile('(?<=##)(.*?)=(.*$)')

		# parser vcf decompressé dans plusieurs dictionnaires
		with open("decompressed_vcf_output.vcf") as decompressed_out:
			vcf_line_nb, metadata_line_nb = 0, 0
			# # TODO: replace metadata functions with calling metadata from h5 once i store it there
			for line in decompressed_out:
				if line.startswith('##'):
					metadata_tag = regex_metadata.search(line).group(1)
					metadata_result = regex_metadata.search(line).group(2)
					# dans un premier temps on va pas s'interesser pas aux
					# metadonéées 'en majiscules' ("INFO" "FILTER" "ALT")
					# peuvent être regroupés ensemble dans un tableau
					if metadata_tag.isupper():
						metadata_type = metadata_tag
					else:
						metadata_type = "basic"
					metadata_dict_entry = [metadata_type, metadata_tag, metadata_result]
					if not metadata_dict_entry in metadata_dict.values():
						metadata_dict[metadata_line_nb] = metadata_dict_entry
					metadata_line_nb += 1

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

		# Read H5 for actual table populating
		complete_h5_file = pd.read_hdf('input.h5')
		self.populate_table(complete_h5_file)



	def populate_table(self, selected_h5_data):
		# effacer tableau actuelle
		self.viewer_tab_table_widget.setRowCount(0)
		self.viewer_tab_table_widget.setColumnCount(0)

		# Creer Tableau vide avec bonne nombre de colonnes et lignes
		table_length = selected_h5_data.shape[0]
		table_width = selected_h5_data.shape[1]

		self.viewer_tab_table_widget.setRowCount(table_length)
		self.viewer_tab_table_widget.setColumnCount(table_width)

		self.viewer_tab_table_widget.setHorizontalHeaderLabels(column_names)
		# set filter_box to list column_names
		self.filter_box.addItems(column_names)
		self.filter_table_btn.clicked.connect(self.filter_table)

		# Remplir Tableau
		vcf_line_nb = 0
		for line in selected_h5_data.itertuples():
			# turn tuple into a list, but exclude the first item of list
			# because its the h5 index and not part of our original data
			line = list(line)[1:]

			vcf_field_nb = 0
			for vcf_field in line:
				self.viewer_tab_table_widget.setItem(
					vcf_line_nb, vcf_field_nb, QtWidgets.QTableWidgetItem(vcf_field))
				vcf_field_nb += 1
			vcf_line_nb += 1


	def select_and_process(self):
		selected_vcf = self.select_vcf()
		self.process_vcf(selected_vcf)


	def h5_encode(self, selected_vcf):
		"""
		Lis en entier le vcf selectionné, bloc par bloc et le met dans un
		fichier HDF5.
		"""
		complevel = self.MetallaxisSettings.compression_level_spinBox.text()
		complib = self.MetallaxisSettings.compression_comboBox.currentText()
		h5chunksize = self.MetallaxisSettings.vcf_chunk_size.text()
		h5_file = pd.HDFStore('input.h5',mode='w')
		chunked_vcf = pd.read_csv(selected_vcf,
							delim_whitespace=True,
							skiprows= range(0,metadata_num-1),
							chunksize=int(h5chunksize),
							low_memory=False,
							dtype=object, # make default data type an object (ie. string)
							)

		for chunk in chunked_vcf:
			# Rename Columns
			chunk.rename(columns = {'#CHROM':'CHROM'}, inplace = True)
			# TODO: remove index from here and add after all the appending
			h5_file.append('df',chunk, index=True, data_columns=True, min_itemsize=80,complib=complib,complevel=int(complevel))
			global column_names
			column_names = list(chunk.keys())

		h5_file.close()


settings_window_object, settings_base_object = uic.loadUiType("MetallaxisSettings.ui")
class MetallaxisSettings(settings_base_object, settings_window_object):
	"""
	La fenetre de pparamètres qui permet de modifier les paramètres de defaut
	"""
	def __init__(self):
		super(settings_base_object, self).__init__()
		self.setupUi(self)
		self.setWindowTitle("Metallaxis Settings")
		self.compression_comboBox.addItems(['zlib', 'blosc', 'bzip2', 'lzo'])
		self.vcf_chunk_size.setText("5000")


# Si le script est executé directement, lance l'interface graphique
if __name__ == '__main__':
	MetallaxisApp = QApplication(sys.argv)
	MetallaxisGui_object = MetallaxisGui()
	MetallaxisGui_object.show()
	sys.exit(MetallaxisApp.exec_())
