Installation Guide
==================

This guide will help you install the MSCI package and its dependencies.

Clone the Repository
--------------------
If you prefer to work with the latest code directly from the repository, you can clone it using the following command:

.. code-block:: bash

    git clone https://github.com/proteomicsunitcrg/MSCI.git

Installing MSCI via pip
-----------------------
You can install MSCI directly using pip, which will also handle the necessary dependencies.

.. code-block:: bash

    pip install MSCI==0.2.0

If you're working in a Jupyter notebook or a Colab environment, you may need to adjust the Python path to include the MSCI module:

.. code-block:: python

    import sys
    sys.path.append('/path/to/MSCI')


**Note**: If you already have `matchms` installed and your environment prompts you to restart the session, select 'Cancel' to avoid any disruptions.

Working with the Repository
---------------------------
After cloning the repository, navigate to the MSCI directory:

.. code-block:: bash

    cd MSCI

This step ensures that you are in the correct directory before running any scripts or modules from the MSCI package.

Installing in Google Colab
--------------------------
If you're using Google Colab, you can execute the following code to install the MSCI package and its dependencies:

.. code-block:: python

    pip install MSCI==0.2.0

Remember to cancel any restart prompts from Colab to avoid disrupting your environment.

Conclusion
----------
Once installed, you can start using the MSCI package to process mass spectrometry data, perform peptide grouping, and calculate similarity scores.
