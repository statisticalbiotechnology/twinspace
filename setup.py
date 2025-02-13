#!/usr/bin/env python

from setuptools import setup, find_packages

with open("README.rst", "r", encoding="utf-8") as fh:
    long_description = fh.read()

requirements = [
    'Click>=7.0',
    'streamlit',
    'matchms',
    'scipy'
]

setup(
    author="Zahra ELHAMRAOUI",
    author_email='zahra.elhamraoui@crg.eu',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.8',
    ],
    description="MSCI assesses peptide fragmentation spectra information content.",
    long_description=long_description,  # ✅ Adding long description
    long_description_content_type="text/markdown",  # ✅ Specifies markdown format
    install_requires=requirements,
    license="MIT license",
    include_package_data=True,
    keywords='MSCI',
    name='MSCI',
    packages=find_packages(include=['MSCI', 'MSCI.*']),
    test_suite='tests',
    tests_require=requirements,
    url='https://github.com/proteomicsunitcrg/MSCI',
    version='0.2.0',
    zip_safe=False,
)
