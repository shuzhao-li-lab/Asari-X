import uuid
import tqdm
import json
import numpy as np
from utils import sum_formula_dicts
from mass2chem.formula import parse_chemformula_dict, dict_to_hill_formula, calculate_formula_mass, calculate_mass
from collections import defaultdict
from itertools import product

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
            return new_cpds, uuid_map
        else:
            return None, None
    
    @staticmethod
    def __initialize_reactions(reactions):
        # placeholder function
        return reactions
    
    @staticmethod
    def __initialize_signatures(signatures):
        # placeholder function
        return signatures

    @staticmethod
    def from_signatures(signatures, compounds=None, reactions=None):
        return SignatureGenerator(compounds, reactions, signatures)

    @staticmethod
    def from_compounds_reactions(compounds, reactions):
        #assert compounds is not None and reactions is not None
        return SignatureGenerator(compounds, reactions, None)

    def cartesian_product_reactions(self, reaction_depth):
        all_reactions = [{
            "reaction_name": '',
            "formula_dict": {}
        }]
        num_reactions = len(self.reactions)
        for i in tqdm.tqdm(range(reaction_depth), position=0, leave=False, desc="Reaction Depth:"):
            for index in tqdm.tqdm(np.ndindex(tuple([num_reactions for _ in range(i)])), position=1, leave=False, desc="Enumerating Rxn Combos"):
                to_combine = [self.reactions[j] for j in index]
                combined_formula = sum_formula_dicts([r['formula_dict'] for r in to_combine])
                combined_rxn_name = "+".join([r['reaction_name'] for r in to_combine])
                if combined_rxn_name:
                    all_reactions.append({
                        "reaction_name": combined_rxn_name,
                        "formula_dict": combined_formula
                    })
        return all_reactions

    def generate_signatures(self, reaction_depth=3):
        assert self.signatures is None, "signatures already generated or provided"    
        self.signatures = []
        reactions = self.cartesian_product_reactions(reaction_depth)
        

    def save_signatures(self, signature_path):
        assert self.signatures is not None, "cannot save null signatures"


class Reactor():
    def __init__(self):
        self.uuid_map = {}
        self.compounds = []
        self.reactions = []
        self.products = []

#    def set_compounds(self, compounds):
#        for compound in compounds:
##            if 'neutral_formula' in compound and compound['neutral_formula']:
#                cpd = {}
#                cpd["parent"] = None
#                cpd["reactions"] = []
#                cpd["uuid"] = str(uuid.uuid4())
#                for k, v in compound.items():
#                    cpd[k] = v
#                self.compounds.append(cpd)
#                self.uuid_map[cpd['uuid']] = cpd
    
#    def set_reactions(self, reactions):
#        self.reactions = [r for r in reactions if r['frequency'] > 20]

    def generate_reactions(self, reaction_depth):
        all_reactions = [{
            "reaction_name": '',
            "formula_dict": {}
        }]
        num_reactions = len(self.reactions)
        for i in tqdm.tqdm(range(reaction_depth), position=0, leave=False, desc="Reaction Depth:"):
            for index in tqdm.tqdm(np.ndindex(tuple([num_reactions for _ in range(i)])), position=1, leave=False, desc="Enumerating Rxn Combos"):
                to_combine = [self.reactions[j] for j in index]
                combined_formula = sum_formula_dicts([r['formula_dict'] for r in to_combine])
                combined_rxn_name = "+".join([r['reaction_name'] for r in to_combine])
                if combined_rxn_name:
                    all_reactions.append({
                        "reaction_name": combined_rxn_name,
                        "formula_dict": combined_formula
                    })
        return all_reactions


    def react_all(self, reaction_depth):
        reaction_depth += 1
        products = []
        all_rxns = self.generate_reactions(reaction_depth)
        for cpd in tqdm.tqdm(self.compounds, position=0, leave=False, desc="Compounds:"):
            for rxn in tqdm.tqdm(all_rxns, position=1, leave=False, desc="Reactions:"):
                p_formula = sum_formula_dicts([parse_chemformula_dict(cpd['neutral_formula']), rxn['formula_dict']])
                if min(p_formula.values()) > 0:
                    products.append({
                        "neutral_formula_mass": calculate_mass(p_formula),
                        "neutral_formula": dict_to_hill_formula(p_formula),
                        "reactions": rxn['reaction_name'],
                        "uuid": str(uuid.uuid4())
                    })
        self.products = products
        self.uuid_map.update({p['uuid']: p for p in products})


    def save_signatures(self, out_path):
        with open(out_path, 'w+') as out_fh:
            json.dump(self.products, out_fh, indent=True)



if __name__ == '__main__':
    import json
    import sys
    R = Reactor()
    R.set_compounds(json.load(open(sys.argv[1]))["compounds"])
    R.set_reactions([x for x in json.load(open(sys.argv[2]))["reactions"] if x['frequency'] > 10])
    R.react_all(int(sys.argv[4]))
    R.save_signatures(sys.argv[3] + ".json" if not sys.argv[3].endswith(".json") else sys.argv[3])
    

