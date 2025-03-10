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

`pip3 install asarix`

Alternatively, you can install from source using `pip3 install <asarix_package_dir>`.

Anaconda is not explicitly supported.

Installation will add two scripts on your local path to interact with Asari-X. First is the CLI which is accessibly by running `asarix`, the second is a basic GUI which can be reached at `asarix-gui`. GUI is provided for convenience; however, is currently not officially supported yet. 

Usage and Quick Start
=====================

The Asari-X workflow consists of two steps: 

1. exposome signature generation (optional)
2. signature search and scoring

In step 1, compound signatures can optionally be generated 
via prediction from parent compounds using reaction rules 
or alternatively you can move directly to step 2 by using our 
empirical signature database or your own. 

In step 2, compound signatures are extracted from either: a) centroided
mzML files or b) asari processed feature tables. Option A offers maximal sensitivity but B will provide information on a feature level. 

See below for instructions on each step.

What is a Signature
===================

A signature is the searchable unit in Asari-X. A signature minimally has an m/z value to search for and can be related to a parent compound or metabolite via metadata fields. Below are examples of signatures:

Minimal Example
-
<EXAMPLE 1>

Generated Example
-
This will be the resulting format from using the signature generator

<EXAMPLE 2>

Complex Example
-
An example of a qualified signature that can be passed manually.

<EXAMPLE 3>


Signature Generation
====================

Signatures can be manually generated and provided to the search step (see below) or first signatures can be generated de novo from input compounds and reactions. Currently, the reaction logic is very simple and formula based only i.e., it does not consider the feasibility of the reactions based on structure. The number of reaction steps is bounded by an integer, the reaction depth, which specifies the maximum number of reactions that can be chained together. 

Note that many reaction chains yield degenerate or redundant combinations.

Signatures are generated using the 'build_signatures' command, this requires also passing a list of reactions, compounds, and a target to save the resulting signatures to. 

Here is an example:

`python3 ./asarix/main.py build_signatures -i . -r "<PATH_TO_RXN.json>" -c "<PATH_TO_CPD.json>" -s "<SIGNATURE_PATH.json>"`

This will react each cpd in PATH_TO_CPD.json with each permutation of reactions with a magnitude of reaction_depth or less in PATH_TO_RNX.json and save them at the specified SIGNATURE_PATH.json file for use in signature search and scoring. 

Signature Search
================
