"""
A collection of helper functions for throughout Asari-X.

Most importantly the consecutive scan set generation is defined here.
"""

from collections import defaultdict

def logo():
    """
    By making this a function, we can reuse it elsewhere.
    """
    return  """ 
###############################
#                             #
#           ASARI-X           #
#                             #
###############################
"""

def sum_formula_dicts(formula_dicts):
    """
    Given an interable full of dicts, reduce them by summing the values of shared keys
    and storing in a new dictionary.

    Each dict is assumed to be a key value pair where the value is an integer.

    Args:
        formula_dicts (list): list of dicts with values that are integers

    Returns:
        dict: summed values for all keys across all dicts
    """
    _d = defaultdict(int)
    for wd in formula_dicts:
        for key in wd:
            assert isinstance(wd[key], int)
            _d[key] += wd[key]
    return _d

# todo - maybe this should not be here?
def consecutive_scans(scans, max_gap=2, min_group_size=2):
    """
    Given a list of integers representing scan numbers, find all such subsets of the
    list of scans such that there is a continuous set of scans separated by no more than 
    max_gap between any two scans and a minimum number of min_group_size consecutive scans. 

    Args:
        scans (list): list of scans, scans are integers
        max_gap (int, optional): minimum gap between consecutive scans. Defaults to 2.
        min_group_size (int, optional): minimum size of CSS to be returned. Defaults to 2.

    Returns:
        list: list of lists representing consecutive scan sets
    """
    if scans.shape[0]:
        groups = [[scans[0]]]
        for scan in scans[1:]:
            if scan == groups[-1][-1] + 1:
                groups[-1].append(scan)
            elif scan == groups[-1][-1]:
                pass
            else:
                current_gap = 0
                gap_filled = False
                while current_gap < max_gap:
                    current_gap += 1
                    if scan == groups[-1][-1] + 1 + current_gap:
                        groups[-1].append(scan)
                        gap_filled = True
                        break
                if not gap_filled:
                    groups.append([scan])
        return [g for g in groups if len(g) >= min_group_size]
    return []