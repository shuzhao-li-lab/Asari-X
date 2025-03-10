from collections import defaultdict
import numpy as np

def sum_formula_dicts(formula_dicts):
    _d = defaultdict(int)
    for wd in formula_dicts:
        for key in wd:
            _d[key] += wd[key]
    return _d

def consecutive_scans(scans, max_gap=2, min_group_size=2, cutoff=10000):
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