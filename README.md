[![CI](https://github.com/delphes/delphes/actions/workflows/ci.yml/badge.svg)](https://github.com/delphes/delphes/actions/workflows/ci.yml)
[![DOI](https://zenodo.org/badge/21390046.svg)](https://zenodo.org/badge/latestdoi/21390046)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/delphes.svg)](https://anaconda.org/conda-forge/delphes)

# Delphes-BlindCalorimeter

Delphes is a C++ framework, performing a fast multipurpose detector response simulation.

More details can be found on the Delphes website https://delphes.github.io .

This is a fork of Delphes that introduces a blind calorimeter feature, allowing users to define insensitive bins in the calorimeter. This extension is useful for studying detector performance and carrying out customized analyses within Delphes.

# Overview of New Feature

The following feature is added to the **`SimpleCalorimeter`** module: 

- A new parameter **`InsensitiveEtaPhiBins`** can be defined.  
- Bins specified in this parameter will **not register any particles, hits, or tracks**.  

An example usage of this module can be found in the **`delphes_card_CMS_blind.tcl`** card.  


# Quick start with Delphes

Commands to get the code:

```
wget http://cp3.irmp.ucl.ac.be/downloads/Delphes-3.5.0.tar.gz

tar -zxf Delphes-3.5.0.tar.gz
```

Commands to compile the code:

```
cd Delphes-3.5.0

make
```

Finally, we can run Delphes:

```
./DelphesHepMC3
```

Command line parameters:

```
./DelphesHepMC3 config_file output_file [input_file(s)]
  config_file - configuration file in Tcl format
  output_file - output file in ROOT format,
  input_file(s) - input file(s) in HepMC format,
  with no input_file, or when input_file is -, read standard input.
```

For more detailed documentation, please visit https://delphes.github.io/workbook

