import csv
import os
import io
import json
import time
from mass2chem.formula import calculate_formula_mass

with io.open(os.path.join(os.path.dirname(__file__), "other/FDA_drug_library.csv"),'r',encoding='utf-8',errors='ignore') as infile, \
     io.open('/tmp/FDA_temp.txt','w',encoding='ascii',errors='ignore') as outfile:
    for line in infile:
        print(*line.split(), file=outfile)

compounds = []
for entry in csv.DictReader(open('/tmp/FDA_temp.txt')):
    try:
        compounds.append(
            {
                "primary_id": entry['CatalogNumber'],
                "primary_db": 'FDA_drug_library',
                "name": entry["Item Name"],
                "neutral_formula": entry['Formula'],
                "neutral_formula_mass": calculate_formula_mass(entry['Formula']),
                "charged_formula": None,
                "charge": 0,
                "charged_formula_mass": None
            }
        )
    except:
        pass
json.dump({
    "compounds": compounds,
    "metadata": {
        "source": 'L1021_Discovery_Probe_FDA_Drug_Database',
        "date": time.time(),
        "notes": "missing some compounds with elements not covered in mass2chem"
    }
}, open(os.path.join(os.path.abspath(os.path.dirname(__file__)), "cpds/FDA_drugs.json"), 'w+'), indent=4)
