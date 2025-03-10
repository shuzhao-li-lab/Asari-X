Asari-X
=======
Repository-Scale eXposome data mining, part of the Asari Suite.

Introduction
============
Asari-X is a tool for detecting signatures of exposomic compounds in large metabolomics datasets. This includes:

- environmental pollutants
- drugs
- dietary ingredients
- microbial products

These compounds are often not present at high frequency or high concentration in human samples, therefore requiring separate computational approaches for detection and monitoring. These are different from the Asari data processing routines used for LC and GC-MS pre-processing.

Asari-X is a screening tool with an emphasis on high-performance signature detection. The goal is to have timely interpretation and analysis even on increasingly larger datasets.

Installation
============

Asari-X can be installed using pip:

<ADD INSTRUCTIONS>

Or from source. Installation on Anaconda is not supported. 

Workflows and Algorithms
========================

The Asari-X workflow consists of two steps: 

1. exposome signature generation (optional)
2. signature search and scoring

In step 1, compound signatures can optionally be generated 
via prediction from parent compounds using reaction rules 
or alternatively you can move directly to step 2 by using our 
empirical signature database or your own. 

In step 2, compound signatures are extracted from either: a) centroided
mzML files or b) asari processed feature tables. Option A offers maximal sensitivity but B will provide information on a feature level. 

Signature Generation
====================

# To Be Added


Signature Search and Scoring
============================

# mzML Search

# Asari Feature Table Search

Data Input Types (GC or LC)
===========================

# Need to think about this more


Use 
===
Help information:

`python3 -m xasari.main`

To process all mzML files under directory mydir/projectx_dir:

`python3 -m xasari.main --ref referenceDatabase.json --input mydir/projectx_dir --para mydir/projectx_dir/parameters.json`

