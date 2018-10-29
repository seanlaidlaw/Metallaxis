#!/usr/bin/env python3
# -*- coding: utf-8 -*-

## Importer des modules de Python Standard
import sys
# pour gerer les VCF compressé
import lzma
import bz2
import gzip

# Importer les modules de tierce-partie
import magic  # pour detecter type de fichier
import allel  # pour convertir vcf en h5
import h5py  # pour lire les fichiers h5
from PyQt5 import QtCore, QtGui, QtWidgets, uic
from PyQt5.QtWidgets import QApplication, QDialog



def decompress_vcf(compression_type, selected_vcf):
	"""
	Décompresse le fichier d'entrée (si compressé avec with xz/gz/bz2), et
	retourne un fichier H5
	"""
	if compression_type == "":
		decompressed_arg_file = open(selected_vcf, mode="rb")
	else:
		decompressed_arg_file = eval(compression_type).open(selected_vcf, mode="rb")
	decompressed_out = open("decompressed_vcf_output.vcf", "wb")
	decompressed_out.write(decompressed_arg_file.read())
	decompressed_arg_file.close()
	decompressed_out.close()
	# transforme les fichiers decompressés en fichiers HDF5 afin qu'on puisse
	# analyser les VCF très gros et qui sont normalement trop gros pour pouvoir
	# stocker en RAM
	allel.vcf_to_hdf5("decompressed_vcf_output.vcf", "input_file.h5", overwrite=True)
	h5_input = h5py.File("input_file.h5", mode="r")
	return h5_input

# Charge l'interface graphique construit en XML
gui_window_object, gui_base_object = uic.loadUiType("MetallaxisGui.ui")

class MetallaxisGui(gui_base_object, gui_window_object):
	"""
	Classe qui construit l'interface graphique Qt sur lequel repose Metallaxis
	"""
	def __init__(self):
		super(gui_base_object ,self).__init__()
		self.setupUi(self)
		self.open_vcf_button.clicked.connect(self.select_vcf)

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
		self.loaded_vcf_lineedit.setText(selected_vcf)

		# Detecte type de fichier et decompresse si compressé
		try:
				# utilise "magic" pour determiner si compressé ou non car .vcf
				# dans le nom ne veut pas forcement dire le bon format
				arg_file_type = magic.from_file(selected_vcf)
		except FileNotFoundError:
				# catch erreur si fichier n'existe pas
				print("ERROR: File given in argument does not exist. You specified : " + str(selected_vcf))
				exit(1)

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
				print("Error: Argument must be a VCF file")
				exit(1)  # quitte le logiciel avec code 1 (erreur)



# Si le script est executé directement, lance l'interface graphique
if __name__ == '__main__':
    MetallaxisApp = QApplication(sys.argv)
    MetallaxisGui_object = MetallaxisGui()
    MetallaxisGui_object.show()
    sys.exit(MetallaxisApp.exec_())
