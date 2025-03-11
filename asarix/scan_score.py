"""
This module implements the methods needed to perform scan scoring using the 
consecutive scan set method described in README.md.
"""

import os
import logging
import json
import multiprocessing as mp
import numpy as np
import tqdm

from scipy.stats import spearmanr
from asarix.utils import consecutive_scans


class mzML_Search_Scorer():
    """
    This object implements the mzML search scoring
    """
    def __init__(self, snr_cutoff=None, scan_cutoff=None, scan_files=None):
        self.frequencies = {}
        self.max_scans = {}
        self.snr_cutoff = snr_cutoff
        self.scan_cutoff = scan_cutoff
        self.scan_files = scan_files[:1]
        assert self.snr_cutoff > 0, "snr_cutoff must be positive"
        assert self.scan_cutoff >= 0, "scan_cutoff must be non-negative"

    @staticmethod
    def filter_inputs(input, extension_filter="ASARIX.json"):
        """
        This function takes a provide input file or directory, infer if it is a 
        directory or file and generate an appropriate iterable containing all 
        absolute paths to all files with the expected extension specified by 
        extension_filter. 

        Args:
            input (str): path to mzML or directory with mzML files within it.
            extension_filter (str, optional): return files matching this extension. Defaults to "mzML".

        Raises:
            Warning: raise if no mzML files were found

        Returns:
            list: list of absolute paths to matching files
        """
        logging.info(f"Processing {input} for {extension_filter} files")
        if os.path.isdir(input):
            logging.info(f"{input} appears to be directory")
            pass_filter = []
            for f in os.listdir(input):
                f = f.rstrip()
                if f.endswith(extension_filter):
                    print("\t", f)
                    pass_filter.append(os.path.join(os.path.abspath(input), f))
            logging.info(f"{len(pass_filter)} were found in directory")
            return pass_filter
        elif os.path.isfile(input):
            logging.info(f"{input} appears to be a file")
            if input.endswith(extension_filter):
                return [os.path.abspath(input)]
        else:
            logging.warn(f"{input} could not be processed")
            raise Warning("Could not infer mzML file locations")


    @staticmethod
    def from_params(params):
        """
        Instantiate from Asari-X params

        Args:
            params (dict): Asari-X 

        Returns:
            mzML_Search_Scorer: _description_
        """

        scan_files = mzML_Search_Scorer.filter_inputs(params['input'], extension_filter=".scans_ASARIX.json")
        return mzML_Search_Scorer(params['snr_cutoff'], params['scan_cutoff'], scan_files)
    
    def score(self):
        jobs = [(x, self.snr_cutoff, self.scan_cutoff) for x in self.scan_files]
        with mp.Pool(mp.cpu_count()) as workers:
            r = list(tqdm.tqdm(workers.imap(self.score_signatures_wrapped, jobs), total=len(self.scan_files)))

    def score_signatures_wrapped(self, job):
        return self.score_signatures(job[0], job[1], job[2])

    @staticmethod
    def digest_signatures(sig_dict, scan_cutoff):
        return {
            s: {
                "scans": np.array([x[0] for x in d if x[1] > scan_cutoff]),
                "intensities": [x[1] for x in d if x[1] > scan_cutoff],
                "masses": [x[2] for x in d if x[1] > scan_cutoff],
                "times": [x[3] for x in d if x[1] > scan_cutoff],
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
    def score_signatures(file, snr_cutoff, scan_cutoff):
        scores = {"scores": {}}
        sig_dict = json.load(open(file))
        digested = mzML_Search_Scorer.digest_signatures(sig_dict['hits'], scan_cutoff)
        topo_sig = mzML_Search_Scorer.topo_sort_signatures(sig_dict['hits'])
        for signature in topo_sig.keys():
            S = mzML_Search_Scorer.score_signature(signature, topo_sig, digested, sig_dict["max_scan"], snr_cutoff)
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
        with open(file.replace(".scans_ASARIX.json", ".scores.json"), 'w+') as out_fh:
            for k, v in sig_dict.items():
                if k != 'hits':
                    scores[k] = v
            scores = mzML_Search_Scorer.consolidate_sig_scores(scores)
            json.dump(scores, out_fh, indent=4)

    @staticmethod
    def score_signature(signature, topo_sig, digested, max_scans, snr_cutoff):
        print(signature)
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
                working_probs = [mzML_Search_Scorer.cluster_prob(wss, count, max_scans) for wss, count in zip(working_scan_sets, ion_counts)]
                snr_filters = [mzML_Search_Scorer.cluster_snr(wi, snr_cutoff) for wi in working_intensities]

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