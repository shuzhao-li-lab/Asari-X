import json
import pymzml
import os
from collections import defaultdict

class ScanLevelSearcher():
    def __init__(self, signatures, blanks, samples, KCD):
        self.signatures = []
        self.load_signatures(signatures)
        self.blanks = blanks
        self.samples = samples
        self.KCD = KCD

    def load_signatures(self, signatures):
        if isinstance(signatures, str) and signatures.endswith(".json"):
            with open(signatures) as sig_fh:
                signatures = json.load(sig_fh)
        if isinstance(signatures, list) and isinstance(signatures[0], dict):
            self.signatures = signatures

    def mine_file(self, 
                  infile):
        all_hits = []
        signature_map = defaultdict(set)
        scan_no = 0
        modes = set()
        search_mz_single = self.KCD.search_mz_single
        try:
            experiment = pymzml.run.Reader(infile)
            specs = [(None, 'pos' if spec['positive scan'] else 'neg', spec.scan_time_in_minutes() * 60, spec.mz, spec.i) for spec in experiment if spec.ms_level==1]
            for scan_no, (_, spec_mode, scan_time, spec_mzs, spec_is) in enumerate(specs):
                modes.add(spec_mode)
                for m, i in zip(spec_mzs, spec_is):
                    result = search_mz_single(m, mode=spec_mode, mz_tolerance_ppm=10)
                    for t in result:
                        all_hits.extend([(t['interim_id'] + '$' + t['ion_relation'], scan_no, int(i), m, scan_time)])
                        for cpd in t['compounds']:
                            signature_map[t['interim_id'] + '$' + t['ion_relation']].add(cpd['uuid'])
        except:
            pass
        feature_dict = {
            "sigmap": {k: list(v) for k,v in signature_map.items()},
            "sample": infile,
            "max_scan": scan_no,
            "hits": self.hits_to_feature_dict(all_hits),
            "signature_map": self.signatures
        }
        if list(modes):
            feature_dict["mode"] = list(modes)[0] if len(modes) == 1 else "multiple"
        else:
            feature_dict["mode"] = None
        return feature_dict
    
    @staticmethod
    def export_feature_dict(feature_dict):
        s = 'formula_mass@ion\tscan_numbers\tintensity\n'
        for k,v in feature_dict.items():
            s += '\t'.join([k, ','.join([str(x) for x in v[0]]), ','.join([str(x) for x in v[1]])] ) + '\n'
        out_path = os.path.join(".", feature_dict["sample"].replace('.mzML', '.xasari'))
        with open(out_path, 'w+') as out_fh:
            out_fh.write(s)

    @staticmethod
    def hits_to_feature_dict(hits):
        feature_dict = {}
        for hit in hits:
            if hit[0] not in feature_dict:
                feature_dict[hit[0]] = []
            feature_dict[hit[0]].append(hit[1:])
        return feature_dict
    
    @staticmethod
    def save_scan_data(feature_dict):
        out_path = os.path.join(".", os.path.basename(feature_dict["sample"]).replace('.mzML', '.scans.json'))
        with open(out_path, 'w+') as out_fh:
            json.dump(feature_dict, out_fh, indent=4)
        

    def search(self):
        for file in self.blanks + self.samples:
            feature_dict = self.mine_file(file)
            self.save_scan_data(feature_dict)

if __name__ == '__main__':
    pass