MSCI Documentation
=================

.. contents:: Table of Contents
   :depth: 2
   :local:

Preprocessing Data
==================

This module provides functionality to generate variable-length peptides from protein sequences
and extract peptides from a FASTA file.

Functions
---------

generate_variable_length_peptides(protein_sequence, min_length=8, max_length=11)
-------------------------------------------------------------------------------

Generates all possible peptides of varying lengths from a given protein sequence.

**Parameters:**
    - **protein_sequence** (*str*) -- The protein sequence from which peptides are generated.
    - **min_length** (*int, optional*) -- The minimum length of peptides to generate (default: 8).
    - **max_length** (*int, optional*) -- The maximum length of peptides to generate (default: 11).

**Returns:**
    - **list** -- A list of generated peptides of varying lengths.

**Example:**
    .. code-block:: python

        peptides = generate_variable_length_peptides("ABCDEFG", min_length=3, max_length=5)
        print(peptides)
        # Output: ['ABC', 'BCD', 'CDE', 'DEF', 'EFG', 'ABCD', 'BCDE', 'CDEF', 'DEFG', 'ABCDE', 'BCDEF', 'CDEFG']


extract_peptides_from_fasta(fasta_path, min_length=8, max_length=11)
--------------------------------------------------------------------

Reads a FASTA file and extracts peptides from each protein sequence.

**Parameters:**
    - **fasta_path** (*str*) -- The path to the FASTA file.
    - **min_length** (*int, optional*) -- The minimum length of peptides to extract (default: 8).
    - **max_length** (*int, optional*) -- The maximum length of peptides to extract (default: 11).

**Returns:**
    - **list** -- A list of peptides extracted from the protein sequences in the FASTA file.

**Example:**
    .. code-block:: python

        peptides = extract_peptides_from_fasta("example.fasta", min_length=3, max_length=5)
        print(peptides)
        # Output: ['ABC', 'BCD', 'CDE', 'DEF', 'EFG', ...]

keep_top_n_peaks(spectrum, n)
------------------------------

Filters a spectrum to retain only the top `n` most intense peaks.

**Parameters:**
    - **spectrum** (*object*) -- A spectrum object containing mass-to-charge ratio (m/z) peaks.
    - **n** (*int*) -- The number of top peaks to retain.

**Returns:**
    - **object** -- A spectrum object with only the top `n` peaks.

**Example:**
    .. code-block:: python

        filtered_spectrum = keep_top_n_peaks(spectrum, n=5)
        print(filtered_spectrum)

filter_spectra_by_top_peaks(input_file_path, output_file_path, n_peaks)
------------------------------------------------------------------------

Reads a pickled list of spectra, processes each spectrum to keep only the top `n` peaks, and saves the results.

**Parameters:**
    - **input_file_path** (*str*) -- Path to the input pickle file containing spectra.
    - **output_file_path** (*str*) -- Path to save the processed spectra as a pickle file.
    - **n_peaks** (*int*) -- The number of top peaks to retain in each spectrum.

**Returns:**
    - **list** -- A list of processed spectra with only the top `n` peaks.


reading MS spectra
---------------------

This module provides functionality to read and process mass spectrometry files, including MSP, MGF, and MZML formats.

Functions
~~~~~~~~

read_msp_file(filename)
--------------------

Reads an MSP file and returns a DataFrame containing the spectra information.

:Parameters:
    - **filename** (*str*) -- The path to the MSP file

:Returns:
    **pandas.DataFrame** with columns:
        - **Name** -- The name of the spectrum
        - **MW** -- Mass/charge of the spectrum
        - **iRT** -- Indexed retention time

read_mgf_file(filename)
--------------------

Reads an MGF file and returns a list of spectra data.

:Parameters:
    - **filename** (*str*) -- The path to the MGF file

:Returns:
    **list** of dictionaries containing:
        - mz_values
        - intensities
        - MW
        - RT

read_mzml_file(filename)
--------------------

Reads an MZML file and returns a list of processed spectrum data.

:Parameters:
    - **filename** (*str*) -- The path to the MZML file

:Returns:
    **list** of dictionaries containing processed spectrum data

read_ms_file(filename)
--------------------

Determines the file format and calls the appropriate function to read the mass spectrometry file.

:Parameters:
    - **filename** (*str*) -- The path to the mass spectrometry file

:Returns:
    **pandas.DataFrame** or **list** depending on the file format

Grouping MS1 Module
--------------------

This module provides functions for grouping MS1 peptides based on mass-to-charge ratio (m/z) and indexed retention time (iRT) using k-d tree data structures and tolerance calculations.

Functions
~~~~~~~~

make_data_compatible(index_df)
--------------------

Converts a DataFrame into a list of tuples compatible with further processing.

:Parameters:
    - **index_df** (*pandas.DataFrame*) -- DataFrame containing mass spectrometry data with columns ``MW`` and ``iRT``

:Returns:
    **list** of tuples in format ``(index, MW, iRT)``

within_ppm(pair, ppm_tolerance1, ppm_tolerance2)
--------------------

Checks if two peptide pairs are within specified tolerances.

:Parameters:
    - **pair** (*tuple*) -- Two peptide tuples ``((index1, MW1, iRT1), (index2, MW2, iRT2))``
    - **ppm_tolerance1** (*float*) -- PPM tolerance for m/z values
    - **ppm_tolerance2** (*float*) -- Absolute tolerance for iRT values

:Returns:
    **bool** -- True if within tolerances, False otherwise

within_tolerance(pair, tolerance1, tolerance2)
--------------------

Checks if peptide pairs are within absolute tolerances.

:Parameters:
    - **pair** (*tuple*) -- Two peptide tuples ``((index1, MW1, iRT1), (index2, MW2, iRT2))``
    - **tolerance1** (*float*) -- Absolute tolerance for m/z values
    - **tolerance2** (*float*) -- Absolute tolerance for iRT values

:Returns:
    **bool** -- True if within tolerances, False otherwise

find_combinations_kdtree(data, tolerance1, tolerance2, use_ppm=True)
--------------------

Uses k-d tree for efficient querying of valid peptide combinations.

:Parameters:
    - **data** (*list*) -- Peptide data tuples ``(index, MW, iRT)``
    - **tolerance1** (*float*) -- Tolerance for m/z values
    - **tolerance2** (*float*) -- Tolerance for iRT values
    - **use_ppm** (*bool*) -- Use PPM tolerance if True, absolute if False

:Returns:
    **list** of valid peptide pairs

Similarity Module
---------------

This module calculates similarity between mass spectrometry spectra using various methods.

Functions and Classes
~~~~~~~~~~~~~~~~~~~

ndotproduct(x, y, m=0, n=0.5, na_rm=True)
--------------------

Calculates normalized dot product between spectra.

:Parameters:
    - **x** (*pandas.DataFrame*) -- First spectrum (m/z and intensities)
    - **y** (*pandas.DataFrame*) -- Second spectrum (m/z and intensities)
    - **m** (*float*) -- M/z values exponent (default: 0)
    - **n** (*float*) -- Intensity values exponent (default: 0.5)
    - **na_rm** (*bool*) -- Remove missing values (default: True)

:Returns:
    **float** -- Normalized dot product

nspectraangle(x, y, m=0, n=0.5, na_rm=True)
--------------------

Calculates normalized spectral angle between spectra.

:Parameters:
    - **x** (*pandas.DataFrame*) -- First spectrum (m/z and intensities)
    - **y** (*pandas.DataFrame*) -- Second spectrum (m/z and intensities)
    - **m** (*float*) -- M/z values exponent (default: 0)
    - **n** (*float*) -- Intensity values exponent (default: 0.5)
    - **na_rm** (*bool*) -- Remove missing values (default: True)

:Returns:
    **float** -- Normalized spectral angle

joinPeaks(tolerance=0, ppm=0)
--------------------

Class that joins peaks from two spectra based on m/z and intensity values.

:Parameters:
    - **tolerance** (*float*) -- Absolute tolerance for m/z matching
    - **ppm** (*float*) -- PPM tolerance for m/z matching

Methods:
    - **match(x, y)** -- Matches peaks from two spectra
        - Parameters: Two DataFrames with m/z and intensities
        - Returns: Tuple of matched peaks DataFrames

Mutation Module
-------------

Tools for processing proteins by simulating peptide digestion and introducing mutations.

Classes and Functions
~~~~~~~~~~~~~~~~~~~

ProteinMutator
--------------------

.. class:: ProteinMutator(proteome_file, mutations_file, output_dir, digestion_method)

    Handles protein mutations and peptide generation.

    :Parameters:
        - **proteome_file** (*str*) -- Path to FASTA proteome file
        - **mutations_file** (*str*) -- Path to TSV mutations file
        - **output_dir** (*str*) -- Output directory path
        - **digestion_method** (*callable*) -- Function returning peptide list

    Methods:
        - **load_proteome()** -- Loads proteome sequences
        - **load_mutations()** -- Loads mutation data
        - **process_protein(target_protein_accession)** -- Processes single protein
        - **process_all_proteins()** -- Processes all proteins

tryptic_digest(sequence)
--------------------

Simulates tryptic digestion of protein sequence.

:Parameters:
    - **sequence** (*str*) -- Protein sequence

:Returns:
    **list** -- Resulting peptides

Example Usage
-----------

.. code-block:: python

    from MSCI.Preprocessing.Koina import PeptideProcessor
    from MSCI.Grouping_MS1.Grouping_mw_irt import process_peptide_combinations
    from MSCI.Preprocessing.read_msp_file import read_msp_file
    from MSCI.Similarity.spectral_angle_similarity import process_spectra_pairs
    from MSCI.data.digest import parse_fasta_and_digest, tryptic_digest, peptides_to_csv
    from matchms.importing import load_from_msp
    import random
    import numpy as np
    import pandas as pd

    # Parse FASTA and perform digestion
    result = parse_fasta_and_digest(
        "https://github.com/proteomicsunitcrg/MSCI/blob/main/tutorial/sp_human_2023_04.fasta",
        digest_type="trypsin"
    )
    peptides_to_csv(result, "random_tryptic_peptides.txt")

    # Initialize and process peptides
    processor = PeptideProcessor(
        input_file="random_tryptic_peptides.txt",
        collision_energy=30,
        charge=2,
        model_intensity="Prosit_2020_intensity_HCD",
        model_irt="Prosit_2019_irt"
    )
    processor.process('random_tryptic_peptides.msp')

    # Load and process spectra
    File = 'random_tryptic_peptides.msp'
    spectra = list(load_from_msp(File))
    mz_tolerance = 1
    irt_tolerance = 5

    mz_irt_df = read_msp_file(File)
    Groups_df = process_peptide_combinations(mz_irt_df, mz_tolerance, irt_tolerance, use_ppm=False)

    # Process similarity pairs
    Groups_df.columns = Groups_df.columns.str.strip()
    index_array = Groups_df[['index1','index2']].values.astype(int)
    result = process_spectra_pairs(index_array, spectra, mz_irt_df, tolerance=0, ppm=10)
    result.to_csv("output.csv", index=False)

    # Plot spectra comparison
    import matplotlib.pyplot as plt
    print(mz_irt_df.iloc[19])
    print(mz_irt_df.iloc[36])
    spectra[19].plot_against(spectra[36])
    plt.savefig('spectra_comparison.png')
