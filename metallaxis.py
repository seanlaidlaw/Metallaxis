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
import allel  # pour convertir vcf en h5
import h5py  # pour lire les fichiers h5
# pour lire uniquement certains lignes des fichiers (reduit conso RAM)
from itertools import islice
from PyQt5 import QtCore, QtGui, QtWidgets, uic
from PyQt5.QtWidgets import QApplication, QDialog, QMessageBox, QTableWidget, QLabel
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



	# # transforme les fichiers decompressés en fichiers HDF5 afin qu'on puisse
	# # analyser les VCF très gros et qui sont normalement trop gros pour pouvoir
	# # stocker en RAM
	# allel.vcf_to_hdf5("decompressed_vcf_output.vcf", "input_file.h5", overwrite=True)
	# h5_input = h5py.File("input_file.h5", mode="r")
	# return h5_input

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
		Parses user instructions to filter table and displays what it ill filter
		"""
		#TODO actually filter the table not just tell user what it will do one day
		selected_filter = self.filter_box.currentText()
		filter_text = self.filter_lineedit.text()
		filter_text = re.sub(r"\s+", "", filter_text)
		filter_text = filter_text.upper()

		if "-" in filter_text and "," in filter_text:
			self.throw_error_message("Please only use either comma separated values or a dash separated range")

		elif "-" in filter_text:
			dash_split_filter_text = filter_text.split("-")
			# filtrer les valeurs nulles ou strings vides pour pas gener le comptage des item
			dash_split_filter_text  = filter(None, dash_split_filter_text)
			dash_split_filter_text  = list(dash_split_filter_text)
			if len(dash_split_filter_text) == 2:
				self.filter_text.setText("Filtering to show " + selected_filter + " from "+ str(dash_split_filter_text[0]) + " to " + str(dash_split_filter_text[1]))
			else:
				self.filter_text.setText(" ")
				self.throw_error_message("Please only enter 2 values separated by a dash")

		elif "," in filter_text:
			comma_split_filter_text = filter_text.split(",")
			# filtrer les valeurs nulles ou strings vides pour pas gener le comptage des item
			comma_split_filter_text  = filter(None, comma_split_filter_text)
			comma_split_filter_text  = list(comma_split_filter_text)
			if len(comma_split_filter_text) > 0:
                            self.filter_text.setText("Filtering to show "  + selected_filter + ": " + str(comma_split_filter_text))
			else:
				self.filter_text.setText(" ")
				self.throw_error_message("Please only enter 2 values separated by a comma")

		else:
			self.filter_text.setText("Filtering to show " + selected_filter + ": " + str(filter_text))


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


# Obtenir Metadonnées à partir du header du fichier VCF:
##source=Tangram
##ALT=<ID=INS:ME:AL,Description="Insertion of ALU element">
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
		metadata_dict = {}
		# Match des groupes des deux cotés du "=", apres un "##"
		regex_metadata = re.compile('(?<=##)(.*?)=(.*$)')
		with open("decompressed_vcf_output.vcf") as decompressed_out:
			# determiner taille du table
			table_width, table_length = 0, 0
			for line in decompressed_out:
				if not line.startswith('#'):
					if table_width < len(line.split("\t")):
						table_width = len(line.split("\t"))
					table_length += 1

		self.viewer_tab_table_widget.setRowCount(table_length)
		self.viewer_tab_table_widget.setColumnCount(table_width)

		# parser vcf decompressé dans plusieurs dictionnaires
		with open("decompressed_vcf_output.vcf") as decompressed_out:
			vcf_line_nb, metadata_line_nb = 0, 0
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
				elif line.startswith('#'):
					line = line.strip()
					line = line.strip("#")
					column_names = line.split("\t")
					self.viewer_tab_table_widget.setHorizontalHeaderLabels(column_names)
				else:
					line = line.strip()
					vcf_field_nb = 0
					for vcf_field in line.split("\t"):
						vcf_field = vcf_field.strip()
						self.viewer_tab_table_widget.setItem(
							vcf_line_nb, vcf_field_nb, QtWidgets.QTableWidgetItem(vcf_field))
						vcf_field_nb += 1
					vcf_line_nb += 1

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

		# set filter_box to list column_names
		self.filter_box.addItems(column_names)
		self.filter_table_btn.clicked.connect(self.filter_table)


	def select_and_process(self):
		selected_vcf = self.select_vcf()
		self.process_vcf(selected_vcf)



# Si le script est executé directement, lance l'interface graphique
if __name__ == '__main__':
	MetallaxisApp = QApplication(sys.argv)
	MetallaxisGui_object = MetallaxisGui()
	MetallaxisGui_object.show()
	sys.exit(MetallaxisApp.exec_())
