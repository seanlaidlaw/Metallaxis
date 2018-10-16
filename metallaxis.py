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


def decompress_vcf(compression_type, selected_vcf):
	"""
	Décompresse le fichier d'entrée (si compressé avec with xz/gz/bz2), et
        retourne un fichier VCF décompressé
	"""
	if compression_type == "":
		decompressed_arg_file = open(selected_vcf, mode="rb")
	else:
		decompressed_arg_file = eval(compression_type).open(selected_vcf, mode="rb")
	decompressed_out = open("decompressed_vcf_output.vcf", "wb")
	decompressed_out.write(decompressed_arg_file.read())
	decompressed_arg_file.close()
	decompressed_out.close()
	return "decompressed_vcf_output.vcf"


# make sure script is only given one arg and that its a vcf
if len(sys.argv) != 2:
	print("Error: Can only take 1 argument : an optionally compressed vcf file")
	exit(1)  # exit with error

# Detecte type de fichier et decompresse si compressé
try:
	# utilise "magic" pour determiner si compressé ou non car .vcf
	# dans le nom ne veut pas forcement dire le bon format
	selected_vcf = sys.argv[1]
	arg_file_type = magic.from_file(selected_vcf)
except FileNotFoundError:
	# catch erreur si fichier n'existe pas
	print("ERROR: File given in argument does not exist. You specified : " + str(selected_vcf))
	exit(1)

# Decompress compressed files to normal vcf text file
if "XZ" in arg_file_type:
	print("Detected Filetype: xz compressed VCF")
	input_vcf = decompress_vcf("lzma",selected_vcf)

elif "bzip2" in arg_file_type:
	print("Detected Filetype: gz compressed VCF")
	input_vcf = decompress_vcf("bz2",selected_vcf)

elif "gzip" in arg_file_type:
	print("Detected Filetype: gz compressed VCF")
	input_vcf = decompress_vcf("gzip",selected_vcf)

elif "Variant Call Format" in arg_file_type:
	print("Detected filetype: uncompressed VCF")
	input_vcf = decompress_vcf("",selected_vcf)
else:
	print("Error: Argument must be a VCF file")
	exit(1)  # quitte le logiciel avec code 1 (erreur)
