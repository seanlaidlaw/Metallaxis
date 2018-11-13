#!/usr/bin/env python3
# -*- coding: utf-8 -*-

## Importer des modules de Python Standard
import sys
import re
# pour gerer les VCF compressé
import lzma
import bz2
import gzip

# Importer les modules de tierce-partie
import magic  # pour detecter type de fichier
import allel  # pour convertir vcf en h5
import h5py  # pour lire les fichiers h5
from PyQt5 import QtCore, QtGui, QtWidgets, uic
from PyQt5.QtWidgets import QApplication, QDialog, QMessageBox, QTableWidget, QLabel



def decompress_vcf(compression_type, selected_vcf):
	"""
	Décompresse le fichier d'entrée (si compressé avec with xz/gz/bz2), et
	retourne un fichier H5
	"""
	if compression_type == "":
		decompressed_arg_file = open(selected_vcf, mode="rb")
	else:
		decompressed_arg_file = eval(compression_type).open(selected_vcf, mode="rb")
	with open("decompressed_vcf_output.vcf", "wb") as decompressed_out:
		decompressed_out.write(decompressed_arg_file.read())
	decompressed_arg_file.close()

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
		self.open_vcf_button.clicked.connect(self.select_vcf)
		# menus sur interface
		self.actionOpen_VCF.triggered.connect(self.select_vcf)


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

		# Detecte type de fichier et decompresse si compressé
		try:
			# utilise "magic" pour determiner si compressé ou non car .vcf
			# dans le nom ne veut pas forcement dire le bon format
			arg_file_type = magic.from_file(selected_vcf)
		except FileNotFoundError:
			# catch erreur si fichier n'existe pas
			self.throw_error_message("ERROR: Selected file does not \
				exist. You specified : " + str(selected_vcf))
			return

		# Decompresse fichiers selectionées en fichiers h5
		if "XZ" in arg_file_type:
			self.detected_filetype_label.setText("xz compressed VCF")
			decompress_vcf("lzma", selected_vcf)

		elif "bzip2" in arg_file_type:
			self.detected_filetype_label.setText("bz2 compressed VCF")
			decompress_vcf("bz2", selected_vcf)

		elif "gzip" in arg_file_type:
			self.detected_filetype_label.setText("gz compressed VCF")
			decompress_vcf("gzip", selected_vcf)

		elif "Variant Call Format" in arg_file_type:
			self.detected_filetype_label.setText("uncompressed VCF")
			decompress_vcf("", selected_vcf)
		else:
			self.throw_error_message("Error: Selected file must be a VCF file")
			return 1

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




# Si le script est executé directement, lance l'interface graphique
if __name__ == '__main__':
	MetallaxisApp = QApplication(sys.argv)
	MetallaxisGui_object = MetallaxisGui()
	MetallaxisGui_object.show()
	sys.exit(MetallaxisApp.exec_())
