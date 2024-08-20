#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

# Define your package's dependencies
requirements = ['Click>=7.0', ]

# Define your test dependencies
test_requirements = ["matchms", " scipy"]

setup(
    author="Zahra ELHAMRAOUI",
    author_email='zahra.elhamraoui@crg.eu',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="MSCI assesses peptide fragmentation spectra information content.",
    entry_points={
        'console_scripts': [
            # Define any command-line scripts here, if applicable
        ],
    },
    install_requires=requirements,
    license="MIT license",
    include_package_data=True,
    keywords='MSCI',
    name='MSCI',
    packages=find_packages(include=['MSCI', 'MSCI.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/proteomicsunitcrg/MSCI',
    version='0.1.0',
    zip_safe=False,
)
