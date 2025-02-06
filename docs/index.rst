.. MSCI documentation master file, created by
   sphinx-quickstart.

============================================
MSCI: Mass Spectrometry Content Information Python Library
============================================

MSCI is a Python package designed for the assessment of peptide fragmentation spectra information content. It helps researchers identify indistinguishable peptides in a given proteome by analyzing spectral similarity scores. 

.. image:: INTRODUCTION.png
   :alt: MSCI workflow illustration

Contents
--------

.. toctree::
   :maxdepth: 2

   installation
   usage
   API
   examples
   contributing

Installation
============

Prerequisites:

- Python 3.8 - 3.11
- Anaconda
- Matchms
- Pyteomics
- OpenMS
- NumPy, Pandas, Biopython

To install MSCI, run:

.. code-block:: bash

   pip install msci

Alternatively, use Anaconda:

.. code-block:: bash

   conda create -n msci_env python=3.9
   conda activate msci_env
   pip install msci

API
===

The MSCI package offers functionalities for:

- **Data Import:** Load proteomes and spectral libraries.
- **Spectra Prediction & Processing:** Predict peptide spectra and filter fragments.
- **Spectra Grouping:** Group peptides based on m/z and iRT values.
- **Similarity Measurement:** Compute spectral similarity using different scoring functions.
- **Output & Visualization:** Export similarity results and generate fragmentation plots.

For full API documentation, see :doc:`API`.

Usage
=====

Example workflow:

For a full tutorial, visit our Colab notebook:
https://colab.research.google.com/drive/1ny97RNgvnpD7ZrHW8TTRXWCAQvIcavkk

Graphical User Interface (GUI)
==============================

MSCI includes a web-based GUI for non-programmers, accessible at:
https://msci--proteomicsunit.streamlit.app/

You can use the GUI for:

- **Peptide Twin Analysis**: Identify indistinguishable peptide pairs.
- **Peptide Twin Checker**: Search for similar peptides in precomputed databases.

Contributing
============

We welcome contributions! Visit our GitHub repository for guidelines:
https://github.com/zahrael97/MSCI

License
=======

MSCI is released under the MIT License.
