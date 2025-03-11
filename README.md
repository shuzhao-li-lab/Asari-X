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
- etc.

These compounds are often not present at high frequency or high concentration in human samples, therefore requiring separate computational approaches for detection and monitoring. These are different from the Asari data processing routines used for LC and GC-MS pre-processing.

Asari-X is a screening tool with an emphasis on high-performance signature detection. The goal is to have timely interpretation and analysis even on increasingly larger datasets.

Installation
============

Asari-X can be installed using pip (once released on PyPi):

`pip3 install asarix`

Alternatively, you can install from source using `pip3 install <asarix_package_dir>`.

Anaconda is not explicitly supported. 

This alpha version of Asari-X requires a hacked version of JMS to work properly. A future update of JMS will eliminate the need for this hack. The workaround is to find the installed location of JMS on your system and replace it with the JMS_hack directory provided in source. THIS WILL BE UNNECESSARY IN A FUTURE VERSION!

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
```
{
    "data": [
        {
            "neutral_formula_mass": 499.937494,
            "neutral_formula": "C8HF17O3S",
            "reactions": "",
            "uuid": "65afb228-569d-4154-ae2a-6e6f1f4a98dc"
        },
        {
            "neutral_formula_mass": 498.953478,
            "neutral_formula": "C8H2F17NO2S",
            "reactions": "",
            "uuid": "0fc5ebc1-205a-4692-8018-e5ccfca81dd5"
        },
    ],
}
```
Generated Example
-
This will be the resulting format from using the signature generator

```
{
    "data": [
        {
            "neutral_formula_mass": 499.937494,
            "neutral_formula": "C8HF17O3S",
            "reactions": "",
            "uuid": "65afb228-569d-4154-ae2a-6e6f1f4a98dc"
        },
        {
            "neutral_formula_mass": 498.953478,
            "neutral_formula": "C8H2F17NO2S",
            "reactions": "",
            "uuid": "0fc5ebc1-205a-4692-8018-e5ccfca81dd5"
        },
    ],
    "metadata": "generated_automatically"
}
```

Complex Example
-
An example of a qualified signature that can be passed manually.

(Not Tested yet for Asari-X)

```
{
"data": [
    {
        "interim_id": "Omeprazole_C17H19N3O3S",
        "Database_referred": [
            "bowen2023"
        ],
        "identity": [
            "Omeprazole"
        ],
        "MS1_pseudo_Spectra": [
            {
            "Compound ID": "Omeprazole",
            "Compound Name": "Omeprazole",
            "Molecular Formula": "C17H19N3O3S",
            "Phase I/II": "N/A",
            "m/z": "346.1221",
            "RT (secs)": "33.9",
            "Ion Form": "[M+H]+",
            "MSI Confidence Level": "2",
            "Schymanski Confidence Level": "2a"
            },
            ...
        ]
        }
    ]
"metadata": "example empCpd format for signatures (untested)"
}
```


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

Signature search can be performed on mzML files in a scan-level search without feature detection or on Asari feature tables (to come). Scan-level search has the potential to be more sensitive but with higher false positive rates whereas feature table mining is predicated on the feature detection algorithm correctly peak picking the xenobiotic features. 

Scan Level Searching
-

Scan-level searching does not perform traditional feature detection, rather, in scan-level mode, the presence of the same mz value in consecutive scans is considered evidence of an analyte with that m/z value even if the signal would not be considered a feature in Asari or other tools. During searching, scans with the signature m/z values will be found, and recorded for future scoring (see next step). Separating search from scoring allows more complex worfklows in the future where searches can be scored using different criteria.

To run a scan-level search, use the following command:

`python3 ./asarix/main.py mzml_search -i <mzml_directory> -s <signatures_for_search.json>`

The default search assumes a mass accuracy of 10 ppm, to change this, either pass `--mz_tolerance_ppm=<ppm_tol>` or `-z=<ppm_tol>`.

Feature Level Searching
-

<TO IMPLEMENT>

Signature Scoring
=================

Scoring is the final step in an Asari-X analysis. The result are scores per signature that are then mapped back to their parent compounds. The resulting scores can then be used by downstream tools. 

Scan Level Scoring
-

Each observed m/z value in the mass spectrometry data represents an ion or a noise signal detected during that scan. For a single m/z value we cannot determine if that m/z value and its accompanying intensity value is due to a true ion signal or due to transient background noise. However, we can distinguish these at the experiment level by looking for consecutive scan sets and scoring them. 

The consecutive scan set (CSS) is the scorable unit in Asari-X for scan-level scoring. The idea is that true ion signals should occur in multiple, consecutive scans, whereas noise signal will be random and not occur in consecutive scans. A CSS is defined as set of scans containing an m/z value of interest where the max scan delta between any two scans is capped by a gap_tolerance which by default is 0 (i.e., no gaps allowed in CSSs).

Step-by-step, this is how scoring is done.

1. For a given signature, collect all scans with that m/z value (from search)
2. For each m/z value, calculate the frequency of that m/z value as a percentage of scans. 
    * Note that a very common m/z value will have a frequency close to 1.0, while a very rare m/z value will have a frequency closer to 0.0. Thus, consecutive scans with a rare m/z value is a more unlikely event than consecutive scans with a common m/z value and this likelihood can be calculated statistically. (see next steps in scoring).
3. For the set of scans with a given m/z value, using a provided 'gap tolerance' find all sets of consecutive scans with the m/z value for that signature. 
    * The gap tolerance allows us to have missing scans but still be considered consecutive. A potential cause of this is when signals are close to the detection limit. Larger gap tolerance values will find more consecutive scan sets but they will be of lower quality.
4. For the CSS calculate the following:
    * The left and right base of the CSS in retention time.
    * The apex of the CSS based on intensity and reported in retention time.
    * The SNR of the region based on the ratio of the apex intensity to the minimum intensity of the CSS, larger values require there to be a more definitive peak.
    * The intensity integral of the CSS.
5. Now calculate the probability the scan sets:
    * For a CSS of magnitude N, with an m/z value observed K times over M total scans in the sample calculate: np.product([[K - 1 - i]/[M - 1 - i]] for i in range(N)). In essence, we are allowing the first scan in the CSS to be free, as the definition of a CSS implies that the first m/z value is not selected at random, then simply assuming that the presence of that m/z value in the following scans is completely random and totally independent. This is not the case for actual signals, but a very low probability implies that the m/z values are not randomly distributed and that it is 
    highly unlikely to be noise. 
    * Report the 1-probability from above as the score of the CSS. 
    * Many CSSs will receive scores close to 0 on account of the fact that most m/z values are relatively rare and do co-occur, indicating that most signal is MS is actually from non-random sources.
6. Now we can add an additional scoring layer to take advantage of isotopologue patterns expected for actual signals as follows:
    * Assume that the monoisotopic version of the signature is the most abundant form
    * Sort remaining isotopologues by decreasing NAP, JMS_hack adds this functionality.
    * For each isotopologue, including the monoisotopic version, generate its intensity vector for the set of scans it shares with the monoisotopic version. 
    * Note that each isotopologue had its own CSS that was scored previously with its own probability. We can, for each isotopologue, weight that score by its correlation with the monoisotopic version and sum all such scores. This yields a total score for the signature that captures information from all isotopologues. The minimum required correlation can be set by a user-provided parameter but defaults to 0.5, isotopologues with correlation below this are discarded. 

To score an mzML search result, run the following command:

`python3 ./asarix/main.py mzml_search_score -i <input_directory>`

Feature Level Scoring
- 

TO IMPLEMENT 

Downstream Analyses using Asari-X Results
=========================================

TO ADD

Future Directions
=================

Asari-X is released as an alpha meaning there will be sharp edges and bugs; however, the core algorithms are completed. Any feedback or comments is appreciated and may be addressed by opening a GitHub issue or Q&A. 

Next steps include:
1. Implementation of feature-level search and scoring
2. Remove JMS hack dependency, will require updated JMS
3. Better scoring - not sure what form this will take but surely there is room for improvement
4. More examples for use cases
5. Writing the manuscript