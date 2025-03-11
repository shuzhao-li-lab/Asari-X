"""
The signature generator class takes a set of compounds and reactions and generates signatures for search. 

This is used during mzML scan searching currently and will be required for feature level searching in the future. 

"""
import json
import uuid
import logging
from itertools import product

import tqdm

from asarix.utils import sum_formula_dicts
from mass2chem.formula import parse_chemformula_dict, dict_to_hill_formula, calculate_formula_mass, calculate_mass

logger = logging.getLogger(__name__)

class SignatureGenerator():
    """
    SignatureGenerator creates signatures for search. 

    You can create a signature generator eitehr from existing signatures or from
    a set of input compounds and reactions. 
    """
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
        """
        Compounds may require additional data wrangling to be useful. This method cleans
        up the input compounds for use. 

        This also assigns the compounds a unique UUID, allowing the mapping of signatures to 
        compounds and reactions easily.

        Args:
            compounds (list): list of input compounds

        Returns:
            tuple: two element tuple with the filtered cpds and the uuid map of them or None, None if parsing failed
        """
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
        """
        This is placeholder to process reactions in the future

        Args:
            signatures (list): list of reactions

        Returns:
            list: the input list of reactions
        """
        if reactions:
            logging.info(f"total of {len(reactions)} reactions found in input.")
        return reactions
    
    @staticmethod
    def __initialize_signatures(signatures):
        """
        This is placeholder to process signatures in the future

        Args:
            signatures (list): list of signatures

        Returns:
            list: the input list of signatures
        """
        # placeholder function
        if signatures:
            logging.info(f"total of {len(signatures)} signatures found in input.")
        return signatures

    @staticmethod
    def from_signatures(signatures):
        """
        Create a SignatureGenerator from signatures. This prevents de novo signature
        generation as signatures are provided in this case from the user.

        Args:
            signatures (list): list of signatures

        Returns:
            SignatureGenerator: SignatureGenerator object with signatures already generated
        """
        return SignatureGenerator(None, None, signatures)

    @staticmethod
    def from_compounds_reactions(compounds, reactions):
        """
        Create a SignatureGenerator from compounds and reactions. This is needed for de novo
        signature generation as signatures can be produced from compounds and reactions.

        Args:
            compounds (list): list of compounds
            reactions (list): list of reactions

        Returns:
            SignatureGenerator: SignatureGenerator object ready for de novo signature generation
        """
        #assert compounds is not None and reactions is not None
        logging.info("creating signature generator from compounds and reactions")
        return SignatureGenerator(compounds, reactions, None)

    def cartesian_product_reactions(self, reaction_depth):
        """
        When we need to reaction M compounds with N reactions up to K times, it is easier to simply 
        calculate the cartesian product of the reactions by calculating the cartesian product of N 
        with repeat of i for i between 0 and K, thus generating the powerset of all reaction combinations. 

        We can then react each reaction combination with each compound later to generate signatures. 

        Args:
            reaction_depth (int): the maximum depth to which reactions should be permuted.

        Returns:
            list: powerset of reaction permutations out to depth
        """
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
        """
        With a configured SignatureGenerator loaded with compounds and reactions, generate signatures
        and save them into the SignatureGenerator. 

        Once generated, the signatures are assigned to a unique UUID and stored in the UUID map of the 
        object. 

        Args:
            reaction_depth (int, optional): the maximum depth to which reactions should be permuted.. Defaults to 3.
        """
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
        """
        Save the signatures in the object to specified signature path.

        Args:
            signature_path (str): path to which we should save the signature path
        """
        logging.info(f"saving signatures to {signature_path}")
        assert self.signatures is not None, "cannot save null signatures"
        #assert not os.path.exists(signature_path), "will not overwrite existing signatures!"
        with open(signature_path, 'w+') as signature_path_fh:
            json.dump({"data": self.signatures, "metadata": "generated_automatically"}, signature_path_fh, indent=4)
