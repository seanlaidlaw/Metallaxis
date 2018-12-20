[![Known Vulnerabilities](https://snyk.io/test/github/SL-LAIDLAW/Metallaxis/badge.svg?targetFile=requirements.txt)](https://snyk.io/test/github/SL-LAIDLAW/Metallaxis?targetFile=requirements.txt)
[![Maintainability](https://api.codeclimate.com/v1/badges/636053f63e1587622300/maintainability)](https://codeclimate.com/github/SL-LAIDLAW/Metallaxis/maintainability)
[![Documentation Status](https://readthedocs.org/projects/metallaxis/badge/?version=latest)](https://metallaxis.readthedocs.io/en/latest/?badge=latest)

[![Python 3.6](https://img.shields.io/badge/python-3.6-blue.svg)](https://www.python.org/downloads/release/python-360/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# Metallaxis

Metallaxis is a Python graphical interface for viewing and annotating VCF
files. On loading a VCF file or compressed variant (vcf.gz,vcf.xz) it will make
an API call to EMSEMBL's [VEP](https://www.ensembl.org/vep) with the IDs of the
VCF.

On the interface basic statistics and graphs are displayed in the
statistics pane, and a Table pane is also avaliable on which variants can be
filtered by data from VCF, from annotation, or parsed form the INFO column.

## Features
- INFO column splitting into columns that can be sorted
- Filtering of all VCF columns
- Automatic annotation from VEP (provided VCF is human)
- Automatically generated statistics and graphs
- Savable analysis as a HDF5 data store

## Authors
Sean Laidlaw & Qiqi He

## Requirements:
Python:
- Python 3.6

Libraries
- python-magic : 0.4.15
- pandas : 0.23.4
- numpy : 1.15.4
- tables : 3.4.4
- PyQt5 : 5.11.2
- requests : 2.20.1
- matplotlib : 3.0.2


## Installation

### Pip
Use the package manager [pip](https://pip.pypa.io/en/stable/) to install
Metallaxis.

```bash
pip install Metallaxis
```

### From Source
Metallaxis can also be installed from source, by running a git clone and then
running the make file which will install dependancies and install as a python
module.

```bash
git clone https://github.com/SL-LAIDLAW/Metallaxis
cd Metallaxis && make
```


## Usage
To open the GUI and get started, run Metallaxis as a module:
```bash
python3 -m metallaxis
```


## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License
[GPLv3](https://choosealicense.com/licenses/gpl-3.0/)
