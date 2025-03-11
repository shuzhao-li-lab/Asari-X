import json
import pymzml
import os
from collections import defaultdict
from jms.dbStructures import knownCompoundDatabase
import tqdm

import logging
logging.getLogger(__name__)
#todo - seems that the logging does not always work correctly, see build_KCD
#todo - the logs are being redirected to khipu.log...


class mzML_Searcher():
    def __init__(self, signatures, mzml_files, ppm):
        self.signatures = signatures
        self.mzml_files = mzml_files[:30]
        self.KCD = self.build_KCD()
        self.ppm = ppm

    def build_KCD(self):
        """
        For the set of provided signatures, build the the knownCompoundDatabase
        (KCD) to allow for the search to occur. 

        Returns:
            knownCompoundDatabase: KCD for the signatures
        """
        logging.info(f"building KCD from signatures")
        KCD = knownCompoundDatabase()
        KCD.mass_index_list_compounds(self.signatures)
        KCD.build_emp_cpds_index(primary_only=True, include_C13=True)
        return KCD
    
    def search(self):
        """
        This method will execute the search_file function on each mzML_file
        in parallel (eventually). This is a wrapper essentially. 

        All Asari-X searchers implement a search method, allowing a future 
        abstract base class implementation, and a common interface.
        """
        for file in tqdm.tqdm(self.mzml_files, desc="searching mzML"):
            logging.info(f"searching {file}")
            feature_dict = self.search_file(file)
            self.save_scan_data(feature_dict)

    def search_file(self, file):
        """
        For a given mzML file, search for all scans with mz values in KCD. 

        As the KCD was configured from signatures.

        Args:
            file (_type_): _description_

        Returns:
            _type_: _description_
        """
        infile = file
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
                    result = search_mz_single(m, mode=spec_mode, mz_tolerance_ppm=self.ppm)
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
    def hits_to_feature_dict(hits):
        feature_dict = {}
        for hit in hits:
            if hit[0] not in feature_dict:
                feature_dict[hit[0]] = []
            feature_dict[hit[0]].append(hit[1:])
        return feature_dict
    
    @staticmethod
    def export_feature_dict(feature_dict):
        s = 'formula_mass@ion\tscan_numbers\tintensity\n'
        for k,v in feature_dict.items():
            s += '\t'.join([k, ','.join([str(x) for x in v[0]]), ','.join([str(x) for x in v[1]])] ) + '\n'
        out_path = os.path.join(".", feature_dict["sample"].replace('.mzML', '.asarix'))
        with open(out_path, 'w+') as out_fh:
            out_fh.write(s)

    def save_scan_data(self, feature_dict):
        out_path = os.path.join(".", os.path.abspath(feature_dict["sample"]).replace('.mzML', '.scans_ASARIX.json'))
        with open(out_path, 'w+') as out_fh:
            logging.info(f"saving scan data to {out_path}")
            json.dump(feature_dict, out_fh, indent=4)
        

    @staticmethod
    def filter_inputs(input, extension_filter="mzML"):
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
            for f in [x for x in os.listdir(input) if x.endswith(extension_filter)]:
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
        This method instantiates an instance of an mzML searcher from Asari-X params. 

        Args:
            params (dict): Asari-X params

        Returns:
            mzML_Searcher: a configured mzML_Searcher instance
        """
        logging.info(f"creating mzML_searcher from params, input: {params['input']}")
        mzml_files = mzML_Searcher.filter_inputs(params['input'])
        return mzML_Searcher(params['signatures'], mzml_files, params['mz_tolerance_ppm'])