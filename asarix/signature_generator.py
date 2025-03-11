import uuid
import tqdm
import os
import json
from utils import sum_formula_dicts
from mass2chem.formula import parse_chemformula_dict, dict_to_hill_formula, calculate_formula_mass, calculate_mass
from collections import defaultdict
from itertools import product

import logging
logger = logging.getLogger(__name__)

class SignatureGenerator():
    def __init__(self, compounds, reactions, signatures):
        self.compounds, self.uuid_map = self.__initialize_compounds(compounds)
        self.reactions  = self.__initialize_reactions(reactions)
        self.signatures = signatures
        #self.__at_creation = [compounds, reactions, signatures]
        assert isinstance(self.compounds, (list, type(None)))
        assert isinstance(self.reactions, (list, type(None)))
        assert isinstance(self.signatures, (list, type(None)))

    @staticmethod
    def __initialize_compounds(compounds):
        if compounds:
            logging.info(f"total of {len(compounds)} compounds found in input.")
            new_cpds, uuid_map = [], {}
            for compound in compounds:
                if 'neutral_formula' in compound and compound['neutral_formula']:
                    cpd = {}
                    cpd["parent"] = None
                    cpd["reactions"] = []
                    cpd["uuid"] = str(uuid.uuid4())
                    for k, v in compound.items():
                        cpd[k] = v
                    new_cpds.append(cpd)
                    uuid_map[cpd['uuid']] = cpd
                else:
                    print(f"error parsing {compound}")
            logging.info(f"total of {len(new_cpds)} compounds parsed correctly.")
            return new_cpds, uuid_map
        else:
            return None, None
    
    @staticmethod
    def __initialize_reactions(reactions):
        # placeholder function
        if reactions:
            logging.info(f"total of {len(reactions)} reactions found in input.")
        return reactions
    
    @staticmethod
    def __initialize_signatures(signatures):
        # placeholder function
        if signatures:
            logging.info(f"total of {len(signatures)} reactions found in input.")
        return signatures

    @staticmethod
    def from_signatures(signatures, compounds=None, reactions=None):
        return SignatureGenerator(compounds, reactions, signatures)

    @staticmethod
    def from_compounds_reactions(compounds, reactions):
        #assert compounds is not None and reactions is not None
        logging.info("creating signature generator from compounds and reactions")
        return SignatureGenerator(compounds, reactions, None)

    def cartesian_product_reactions(self, reaction_depth):
        logging.info(f"generating product of reactions to depth {reaction_depth}")
        all_reactions = [{"reaction_name": "", "formula_dict": {}}]
        for i in tqdm.tqdm(range(reaction_depth), position=0, leave=False, desc="Iterating Depth"):
            for combo in tqdm.tqdm(product(self.reactions, repeat=i), position=1, leave=False, desc="Iterating Reaction Permutations"):
                combined_formula = sum_formula_dicts([r["formula_dict"] for r in combo])
                combined_rxn_name = "+".join(r["reaction_name"] for r in combo)
                if combined_rxn_name:
                    all_reactions.append({
                        "reaction_name": combined_rxn_name,
                        "formula_dict": combined_formula
                    })
        logging.info(f"generated {len(all_reactions)} reaction permutations")
        return all_reactions

    def generate_signatures(self, reaction_depth=3):
        logging.info(f"generating signatures to depth {reaction_depth}")
        reaction_depth += 1
        products = []
        all_rxns = self.cartesian_product_reactions(reaction_depth)
        for cpd in tqdm.tqdm(self.compounds, position=0, leave=False, desc="Permuting Compounds:"):
            for rxn in tqdm.tqdm(all_rxns, position=1, leave=False, desc="Permuting Reactions:"):
                p_formula = sum_formula_dicts([parse_chemformula_dict(cpd['neutral_formula']), rxn['formula_dict']])
                if min(p_formula.values()) > 0:
                    products.append({
                        "neutral_formula_mass": calculate_mass(p_formula),
                        "neutral_formula": dict_to_hill_formula(p_formula),
                        "reactions": rxn['reaction_name'],
                        "uuid": str(uuid.uuid4())
                    })
        logging.info(f"generation produced {len(products)} signatures")
        self.signatures = products
        self.uuid_map.update({p['uuid']: p for p in products})
        
    def save_signatures(self, signature_path):
        logging.info(f"saving signatures to {signature_path}")
        assert self.signatures is not None, "cannot save null signatures"
        #assert not os.path.exists(signature_path), "will not overwrite existing signatures!"
        with open(signature_path, 'w+') as signature_path_fh:
            json.dump({"data": self.signatures, "metadata": "generated_automatically"}, signature_path_fh, indent=4)
