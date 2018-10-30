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
from PyQt5.QtWidgets import QApplication, QDialog, QMessageBox, QLabel



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
		super(gui_base_object ,self).__init__()
		self.setupUi(self)
		self.setWindowTitle("Metallaxis")
		# boutons sur interface
		self.open_vcf_button.clicked.connect(self.select_vcf)
		# menus sur interface
		self.actionOpen_VCF.triggered.connect(self.select_vcf)


	def throw_error_message(self,error_message):
		"""
		Generer une dialogue d'alerte avec l'argument comme message d'alerte
		"""
		print("Error: " + error_message)
		error_dialog  = QtWidgets.QMessageBox()
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


	def select_vcf(self):
		"""
		Ouvre une dialogue pour que l'utilisateur puisse choisir un fichier
		VCF, et ensuite il le passe à la fonction decompress_vcf()
		"""
		select_dialog = QtWidgets.QFileDialog()
		select_dialog.setAcceptMode(select_dialog.AcceptSave)
		selected_vcf = select_dialog.getOpenFileName(self,filter=
                    "VCF Files (*.vcf *.vcf.xz *.vcf.gz *.vcf.bz2)\
                    ;;All Files(*.*)")
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

		# Decompresse fichiers selectionées en fichiers h5
		if "XZ" in arg_file_type:
				print("Detected Filetype: xz compressed VCF")
				h5_input = decompress_vcf("lzma",selected_vcf)

		elif "bzip2" in arg_file_type:
				print("Detected Filetype: gz compressed VCF")
				h5_input = decompress_vcf("bz2",selected_vcf)

		elif "gzip" in arg_file_type:
				print("Detected Filetype: gz compressed VCF")
				h5_input = decompress_vcf("gzip",selected_vcf)

		elif "Variant Call Format" in arg_file_type:
				print("Detected filetype: uncompressed VCF")
				h5_input = decompress_vcf("",selected_vcf)
		else:
				self.throw_error_message("Error: Selected file must be a VCF file")
				return 1
		# active les widgets qui sont desactivés tant qu'on a pas de VCF selectioné
		self.loaded_vcf_lineedit.setText(selected_vcf)
		self.loaded_vcf_lineedit.setEnabled(True)
		self.loaded_vcf_label.setEnabled(True)
		self.metadata_area_label.setEnabled(True)

		# effacer espace metadonées (utile si on charge un fichier apres un autre)
		self.empty_qt_layout(self.dynamic_metadata_label_results)
		self.empty_qt_layout(self.dynamic_metadata_label_tags)


		# Obtenir Metadonnées à partir du header du fichier VCF:
		##source=Tangram
		##ALT=<ID=INS:ME:AL,Description="Insertion of ALU element">
		##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">

		metadata_dict = {}
		# Match des groupes des deux cotés du "=", apres un "##"
		regex_metadata = re.compile('(?<=##)(.*?)=(.*$)')
		with open("decompressed_vcf_output.vcf") as decompressed_out:
			metadata_line_nb = 0
			for line in decompressed_out:
				if line.startswith('##'):
					metadata_tag = regex_metadata.search(line).group(1)
					metadata_result = regex_metadata.search(line).group(2)

					metadata_dict_entry = [metadata_tag, metadata_result]
					metadata_dict[metadata_line_nb] = metadata_dict_entry
					metadata_line_nb = metadata_line_nb + 1

		for metadata_line_nb in metadata_dict:
			metadata_tag = metadata_dict[metadata_line_nb][0]
			metadata_result = metadata_dict[metadata_line_nb][1]
			if not metadata_tag.isupper():
				# Generer dynamiquement du texte pour le titre et resultat pour
				# chaque type de metadonnée non-majiscule
				self.dynamic_metadata_label_tags.addWidget(QtWidgets.QLabel(metadata_tag, self))
				self.dynamic_metadata_label_results.addWidget(QtWidgets.QLabel(metadata_result, self))


# Si le script est executé directement, lance l'interface graphique
if __name__ == '__main__':
	MetallaxisApp = QApplication(sys.argv)
	MetallaxisGui_object = MetallaxisGui()
	MetallaxisGui_object.show()
	sys.exit(MetallaxisApp.exec_())
