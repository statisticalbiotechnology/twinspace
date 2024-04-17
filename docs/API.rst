API
==========



read_msp_file
---------------------
The function parses .msp files, capturing key parameters - the spectrum name, molecular weight (MW), and index of relative retention time (iRT). 

The acquired data is structured into a pandas DataFrame.

.. code::

   >>>  from MSCI.Preprocessing import read_msp_file
   >>>  mz_irt_df = read_msp_file(filename)

.. 
Arguments:

filename: A string representing the path to the MSP file to be read.

read_mzml_file
---------------------
The function parses .mzml files, capturing key parameters - the spectrum name, molecular weight (MW), and retention time (RT). 

The acquired data is structured into a pandas DataFrame.

.. code::

   >>>  from MSCI.Preprocessing import read_mzml_file
   >>>  mz_irt_df = read_mzml_file(filename)

.. 
Arguments:

filename: A string representing the path to the mzml file to be read.

Grouping_mw_rt 
-----



.. code::

   >>> 


joinPeaks
-----
The joinPeaks class is designed for matching peaks between two spectra based on a specified tolerance and parts per million (ppm).
The match method merges two spectra, sorts them by m/z values, and assigns indices to peaks. 
Peaks within the specified tolerance are matched and assigned the same index. The method then filters out duplicated indices with NaN intensities and updates m/z values where intensity is missing. 
The final matched spectra are returned as two separate data frames.
.. code::

   >>> joiner = joinPeaks(tolerance=value1, ppm=value2)
   >>> matched_x, matched_y = joiner.match(x, y)

process_combin 
-----



.. code::

   >>> 




create_mutated_sequence
-----



.. code::

   >>> 


Generate_hla
---------------------



.. code::

   >>>


