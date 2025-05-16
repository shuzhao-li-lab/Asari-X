"""
Asari-X Exposure Prediction Module

This module provides functionality to generate mass spectrometry signatures, search mzML files using these signatures, score scan results, and predict exposure categories based on these scores.

Dependencies:
- asarix
- numpy
- sklearn
- argparse
- json
"""

import json
import os
import numpy as np
import argparse
from asarix.signature_generator import SignatureGenerator
from asarix.scan_search import mzML_Searcher
from asarix.scan_score import mzML_Search_Scorer
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score

def generate_signatures(compounds, reactions, reaction_depth=3):
    """Generate mass spectrometry signatures from given compounds and reactions.

    Args:
        compounds (list): List of compound definitions.
        reactions (list): List of reaction definitions.
        reaction_depth (int, optional): Depth of reaction exploration. Defaults to 3.

    Returns:
        dict: Generated signatures.
    """
    SG = SignatureGenerator.from_compounds_reactions(compounds, reactions)
    SG.generate_signatures(reaction_depth=reaction_depth)
    return SG.signatures

def search_signatures(mzml_input_dir, signatures, mz_tolerance_ppm=5):
    """Search mzML files for the provided signatures.

    Args:
        mzml_input_dir (str): Directory containing mzML files.
        signatures (dict): Signatures to search for.
        mz_tolerance_ppm (int, optional): Mass-to-charge ratio tolerance in ppm. Defaults to 5.
    """
    XS = mzML_Searcher(signatures, mzml_input_dir, mz_tolerance_ppm)
    XS.search()

def score_scans(mzml_input_dir, snr_cutoff=3, scan_cutoff=5):
    """Score scan results based on signal-to-noise ratio and scan count cutoffs.

    Args:
        mzml_input_dir (str): Directory containing mzML files.
        snr_cutoff (int, optional): Signal-to-noise ratio cutoff. Defaults to 3.
        scan_cutoff (int, optional): Scan count cutoff. Defaults to 5.
    """
    scan_files = mzML_Search_Scorer.filter_inputs(mzml_input_dir, extension_filter=".scans_ASARIX.json")
    SS = mzML_Search_Scorer(snr_cutoff, scan_cutoff, scan_files)
    SS.score()

def prepare_data(score_files, target_formula):
    """Prepare scored data for predictive modeling by extracting scores for the target formula.

    Args:
        score_files (list): List of file paths containing scores.
        target_formula (str): Target formula to extract scores.

    Returns:
        numpy.ndarray: Array of scores formatted for modeling.
    """
    scores = []
    for file in score_files:
        with open(file) as fh:
            data = json.load(fh)
            formula_scores = {sig['neutral_formula']: sig['score'] for sig in data['signature_map']}
            scores.append(formula_scores.get(target_formula, 0))
    return np.array(scores).reshape(-1, 1)

def predict_exposure(score_data, labels):
    """Predict exposure categories using a Random Forest classifier.

    Args:
        score_data (numpy.ndarray): Prepared data array for prediction.
        labels (list): List of exposure labels for training.

    Returns:
        tuple: Predictions and model accuracy.
    """
    clf = RandomForestClassifier(n_estimators=100)
    clf.fit(score_data, labels)
    predictions = clf.predict(score_data)
    accuracy = accuracy_score(labels, predictions)
    return predictions, accuracy

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Asari-X Exposure Prediction")
    parser.add_argument('--compounds', required=True, help='Path to compounds.json file')
    parser.add_argument('--reactions', required=True, help='Path to reactions.json file')
    parser.add_argument('--mzml_dir', required=True, help='Directory containing mzML files')
    parser.add_argument('--labels', required=True, help='Path to labels.json file containing exposure labels')
    parser.add_argument('--target_formula', required=True, help='Target formula for analysis')

    args = parser.parse_args()

    compounds = json.load(open(args.compounds))['data']
    reactions = json.load(open(args.reactions))['data']
    labels = json.load(open(args.labels))['labels']

    signatures = generate_signatures(compounds, reactions)
    search_signatures(args.mzml_dir, signatures)
    score_scans(args.mzml_dir)

    score_files = mzML_Search_Scorer.filter_inputs(args.mzml_dir, extension_filter=".scores.json")
    score_data = prepare_data(score_files, args.target_formula)

    predictions, accuracy = predict_exposure(score_data, labels)
    print(f"Predictions: {predictions}")
    print(f"Accuracy: {accuracy}")