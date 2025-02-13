MSCI
=======================

Introduction
------------

Peptide identification by mass spectrometry relies on the interpretation of fragmentation spectra based on:
- The **m/z pattern** (mass-to-charge ratio),
- **Relative intensities** of detected fragments,
- **Retention time (RT)**.

Given a proteome, we explored how many peptides generate highly similar fragmentation spectra with current MS methods.
MSCI is a Python package built to assess the **information content of peptide fragmentation spectra**. 



Main features
-----
The MSCI package offers functionalities for:

- **Data Import**: Load proteomes and spectral libraries.
- **Spectra Prediction & Processing**: Predict peptide spectra and filter fragments.
- **Spectra Grouping**: Group peptides based on m/z and iRT values.
- **Similarity Measurement**: Compute spectral similarity using different scoring functions.
- **Output & Visualization**: Export similarity results and generate fragmentation plots.

For full API documentation, see `MSCI Documentation <https://msci.readthedocs.io>`_.

Usage
-----

Example workflow:

For a full tutorial, visit our Colab notebook:  
`MSCI Colab Notebook <https://colab.research.google.com/drive/1ny97RNgvnpD7ZrHW8TTRXWCAQvIcavkk>`_.

License
-------

MSCI is released under the MIT License.



