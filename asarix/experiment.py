import json
import os
import pandas as pd
import tqdm
from scan_search import ScanLevelSearcher
import multiprocessing as mp
from jms.dbStructures import knownCompoundDatabase

class Experiment():
    def __init__(self):
        self.signatures = []
        self.all_samples = []
        self.blanks = []
        self.study_samples = []
        self.feature_table = None

    def prepare_database(self, signatures):
        if signatures.endswith(".json"):
            with open(signatures) as sig_fh:
                signatures = json.load(sig_fh)
        if isinstance(signatures, list) and isinstance(signatures[0], dict):
            self.signatures = signatures
        self.KCD = knownCompoundDatabase()
        self.KCD.mass_index_list_compounds(self.signatures)
        self.KCD.build_emp_cpds_index(primary_only=True, include_C13=True)
        self.export_empCpds()

    def export_empCpds(self):
        outfile = os.path.join(".", 'empCpds.json')
        with open(outfile, 'w', encoding='utf-8') as f:
            json.dump(self.KCD.mass_indexed_compounds, f, ensure_ascii=False, indent=2)

    def register_samples(self, directory):
        abs_path_dir = os.path.abspath(directory)
        for file in os.listdir(abs_path_dir):
            abs_path_file = os.path.join(abs_path_dir, file)
            if abs_path_file.endswith(".mzML"):
                self.all_samples.append(abs_path_file)

    def register_feature_table(self, feature_table):
        self.feature_table = pd.read_csv(feature_table, sep="\t")

    def partition_samples(self, blank_substring="blank"):
        if self.all_samples:
            for sample in self.all_samples:
                if blank_substring in os.path.basename(sample).lower():
                    self.blanks.append(sample)
                else:
                    self.study_samples.append(sample)
        elif self.feature_table:
            for sample in self.feature_table.columns[11:]:
                if blank_substring in sample.lower():
                    self.blanks.append(sample)
                else:
                    self.study_samples.append(sample)

    def mine_scans(self):
        jobs = [[sample, self.signatures, self.KCD] for sample in self.all_samples]
        with mp.Pool(4) as workers:
            r = list(tqdm.tqdm(workers.imap(self.mine_scans_single, jobs), total=len(self.all_samples)))


    @staticmethod
    def mine_scans_single(job):
        sample, signatures, KCD = job
        S = ScanLevelSearcher(signatures, [], [sample], KCD)
        S.search()

    def mine_features(self):
        pass


if __name__ == '__main__':
    import sys
    E = Experiment()
    E.prepare_database(sys.argv[1])
    E.export_empCpds()
    E.register_samples(sys.argv[2])
    E.partition_samples("blank")
    E.mine_scans()