import json
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from utils import consecutive_scans
from scipy.stats import spearmanr
import tqdm

DEBUG = False
CUTOFF = 0
SNR_CUTOFF = 2.5

class ScanScorer():
    def __init__(self):
        self.frequencies = {}
        self.max_scans = {}

    @staticmethod
    def digest_signatures(sig_dict, cutoff=CUTOFF):
        return {
            s: {
                "scans": np.array([x[0] for x in d if x[1] > cutoff]),
                "intensities": [x[1] for x in d if x[1] > cutoff],
                "masses": [x[2] for x in d if x[1] > cutoff],
                "times": [x[3] for x in d if x[1] > cutoff],
            } for s, d in sig_dict.items()
        }
    
    @staticmethod
    def topo_sort_signatures(sig_dict):
        _s = {} 
        for sig in sig_dict.keys():
            parent, adduct_iso = sig.split("$")
            adduct_iso, order = adduct_iso.split(";")
            if "," in adduct_iso:
                adduct = ','.join(adduct_iso.split(",")[:-1])
            else:
                adduct = adduct_iso
                _ = ''
            key = parent + "_" + adduct
            if key not in _s:
                _s[key] = {}
            _s[key][int(order)] = sig
        _t = {}
        for key, order_dict in _s.items():
            _t[key] = []
            i = 0
            while True:
                if i in order_dict:
                    _t[key].append(order_dict[i])
                else:
                    break
                i += 1
            if not _t[key]:
                del _t[key]
        return _t

    @staticmethod
    def score_signatures(file):
        scores = {"scores": {}}
        sig_dict = json.load(open(file))
        digested = ScanScorer.digest_signatures(sig_dict['hits'])
        topo_sig = ScanScorer.topo_sort_signatures(sig_dict['hits'])
        for signature in topo_sig.keys():
            S = ScanScorer.score_signature(signature, topo_sig, digested, sig_dict["max_scan"])
            for k, v in S.items():
                if v[0] > 0:
                    if signature not in scores["scores"]:
                        scores["scores"][signature] = []
                    scores["scores"][signature].append({
                        "left_base": k[0],
                        "apex": k[1],
                        "right_base": k[2],
                        "score": v[0],
                        "scans": v[1],
                        "freq": v[2],
                        "integral": v[3],
                        "mz": v[4]
                    })
        with open(file.replace(".scans.json", ".scores.json"), 'w+') as out_fh:
            for k, v in sig_dict.items():
                if k != 'hits':
                    scores[k] = v
            scores = ScanScorer.consolidate_sig_scores(scores)
            json.dump(scores, out_fh, indent=4)

    @staticmethod
    def score_signature(signature, topo_sig, digested, max_scans):
        scores = {}
        maps = []
        
        for iso in topo_sig[signature]:
            map = {
                "I": dict(zip(digested[iso]["scans"], digested[iso]["intensities"])),
                "M": dict(zip(digested[iso]["scans"], digested[iso]["masses"])),
                "T": dict(zip(digested[iso]["scans"], digested[iso]["times"]))
            }
            maps.append(map)


        scan_sets = [consecutive_scans(digested[iso]["scans"]) for iso in topo_sig[signature]]
        filter_empty = []
        for ss in scan_sets:
            if ss:
                filter_empty.append(ss)
            else:
                break
        scan_sets = filter_empty

        ion_counts = [len(digested[iso]["scans"]) for iso in topo_sig[signature]]
        if scan_sets and scan_sets[0]:
            for index in np.ndindex(tuple([len(s) for s in scan_sets])):
                working_scan_sets = [scan_sets[i][j] for i,j in enumerate(index)]
                for i in range(len(working_scan_sets)):
                    if i > 0:
                        if len(working_scan_sets[i]) > len(working_scan_sets[i-1]):
                            continue
                #setup
                left_base = maps[0]["T"][min(working_scan_sets[0])]
                right_base = maps[0]["T"][max(working_scan_sets[0])]
                apex = (None, 0)
                for k, v in maps[0]["I"].items():
                    if k in working_scan_sets[0]:
                        if v > apex[1]:
                            apex = (k, v)
                apex = maps[0]["T"][apex[0]]
                if (left_base, apex, right_base) not in scores:
                    scores[(left_base, apex, right_base)] = (0, 0)

                working_intensities = [[maps[i]["I"][s] for s in working_scan_sets[i]] for i in range(len(working_scan_sets))]

                #subscores
                working_probs = [ScanScorer.cluster_prob(wss, count, max_scans) for wss, count in zip(working_scan_sets, ion_counts)]
                snr_filters = [ScanScorer.cluster_snr(wi, SNR_CUTOFF) for wi in working_intensities]

                correlations = []
                for i, _ in enumerate(working_scan_sets):
                    for_i = dict(zip(working_scan_sets[i], working_intensities[i]))
                    for_corr_calc = [for_i.get(wss_0_sno, 0) for wss_0_sno in working_scan_sets[0]]
                    corr = spearmanr(for_corr_calc, working_intensities[0]).statistic
                    if np.isnan(corr):
                        corr = 0
                    correlations.append(corr)
                score = []
                integral = 0
                if snr_filters[0]:
                    for i, (_, prob, corr) in enumerate(zip(snr_filters, working_probs, correlations)):
                        corr_prob = prob
                        #corr_prob = min((prob * np.product(num_tests)), 1)
                        corr = corr * (corr > .5)
                        score.append((1-corr_prob) * corr)
                        if corr: 
                            integral += np.sum(working_intensities[i])
                total_score = np.sum(score)
                if total_score > scores[(left_base, apex, right_base)][0]:
                    scores[(left_base, apex, right_base)] = (
                        total_score, 
                        len(working_scan_sets[0]), 
                        ion_counts[0] / max_scans, 
                        int(integral), float(np.mean(list(maps[0]["M"].values()))))    
        return scores

    @staticmethod
    def consolidate_sig_scores(scores):
        sigscores = {}
        sigmap = scores["sigmap"]
        extra_sigs = {}
        sigmap.update(extra_sigs)
        for sig, psuedo_feature_list in scores["scores"].items():
            sig = sig.replace("_M", "$M") + ";0"
            best_score = 0
            for pseudo_feature in psuedo_feature_list:
                for uuid in sigmap[sig]:
                    if uuid not in sigscores:
                        sigscores[uuid] = 0
                    sigscores[uuid] += pseudo_feature["score"]
        new_signatures = []
        for _d in scores["signature_map"]:
            if _d['uuid'] in sigscores:
                _d["score"] = sigscores[_d["uuid"]]
                new_signatures.append(_d)
            else:
                _d["score"] = 0
        scores['signature_map'] = new_signatures
        del scores['sigmap']

        return scores

            
    @staticmethod
    def cluster_prob(cluster, total_instances, max_scans):
        probs = []
        for k, _ in enumerate(cluster[1:]):
            numerator = total_instances - 1 - k
            denominator = max_scans - 1 - k
            probs.append(numerator / denominator)
        P = np.product(probs)
        return P
    
    @staticmethod
    def cluster_snr(intensities, snr_cutoff, apex_mode="max", baseline_mode="min"):
        modes = {
            "max": np.max,
            "median": np.median,
            "mean": np.mean,
            "min": np.min
        }
        return modes[apex_mode](intensities) > modes[baseline_mode](intensities) * snr_cutoff


if __name__ == '__main__':
    import sys
    import multiprocessing as mp
    mp.freeze_support()
    S = ScanScorer()
    files = [x for x in os.listdir(sys.argv[1]) if x.endswith(".scans.json")]
    with mp.Pool(8) as workers:
        r = list(tqdm.tqdm(workers.imap(S.score_signatures, files), total=len(files)))
