#!/usr/bin/env python

from distutils.core import setup

setup(
    name='gms-gnomad-annot',
    version='0.1',
    description='Tool to add gnomAD values to a GMS annotation file',
    author='Kilannin Krysiak',
    author_email='kkrysiak@wustl.edu',
    url='https://github.com/griffithlab/gnomad_annotation',
    entry_points={
        "console_scripts":[
            "annotate-gnomad = annotator",
        ]
    },
    install_requires=[
        'marisa_trie',
    ],
    classifiers=[
        'Development Status :: 4 - Beta',

        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',

        "Programming Language :: Python :: 3.5"
    ],
)
