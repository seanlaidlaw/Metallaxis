#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
	name="Metallaxis",
	version="0.9",
	license='GPLv3',
	author="Sean Laidlaw",
	author_email="seanlaidlaw95@gmail.com",
	url='https://github.com/SL-LAIDLAW/Metallaxis',
	description="A graphical python-based VCF viewer with optional annotation",
	packages=['metallaxis'],
	package_data={'': ['*.ui', 'annotation/*']},
	include_package_data=True,
	install_requires=[
		'python-magic',
		'pandas',
		'numpy',
		'PyQt5',
		'requests',
		'wget',
		'matplotlib'
	],
	classifiers=[
		# Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
		'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
		'Development Status :: 4 - Beta',
		'Operating System :: OS Independent',
		'Programming Language :: Python :: 3.6',
		'Topic :: Scientific/Engineering :: Bio-Informatics'
	]
)
