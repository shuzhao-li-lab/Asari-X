import argparse
import json
import logging
import os

from default_parameters import PARAMETERS
from utils import logo
from signature_generator import SignatureGenerator
from scan_search import mzML_Searcher

from logger_setup import setup_logger
setup_logger()
logger = logging.getLogger(__name__)


def check_sufficient_params(params, needed=None):
    """
    Helper function to check that necessary fields are provided for a given command. Used
    to print a more useful message than the default if a user does not give the required
    params. 

    Args:
        params (dict): Asari-X params dict
        needed (list, optional): list of required non-None fields in the params dict. Defaults to None.
    """
    if needed:
        for x in needed:
            logging.info(f"Checking params for {x}")
            assert x in params and params.get(x, None) is not None, f'{x} required for this operation!'

def process_params(params, args=None):
    """
    A shared pre-processing function for Asari-X. Using the default_parameters.py 
    allows flexibility and dynamic UI generation but is too unwieldly for passing around. 

    To remedy this, we will convert the verbose parameters into a concise key: value 
    representation and discard the information we do not need for computation. 

    Args:
        params (dict): Asari-X params dict
        args (namespace, optional): args namespace from the CLI. Defaults to None.
    """
    def __combine(params, args=None):
        """
        Combine takes a namespace (and potentially other inputs) and chains the 
        specified parameters in a fixed order to yield a combined params dictionary.

        The order is always:
            Parameters File > CLI/GUI Arguments > Defaults

        Thus if the default values for a field is X, and you pass Y on the CLI and pass Z 
        in the parameters file (specified by providing --parameters or -p) the field will 
        have value Z.

        Args:
            params (dict): Asari-X params dict
            args (namespace, optional): args namespace from the CLI. Defaults to None.

        Returns:
            (dict): updated Asari-X params dict 
        """
        if args is not None:
            for key, value in args.__dict__.items():
                params[key]['value'] = value

        for key, value in list(params.items()):
            if isinstance(value['value'], str) and value['value'].endswith(".json") and not value.get('skip_json', False):
                abs_path = os.path.abspath(value['value'])
                with open(abs_path) as json_fh:
                    data = json.load(json_fh)
                    params[key].update({
                        'value': data['data'],
                        'source': value['value'],
                        'types': [type(data['data'])]
                    })
                    if 'metadata' in data:
                        params[f"__{key}_metadata"] = data['metadata']

        if params.get('parameters', {}).get('value', None) is not None:
            with open(params['parameters']) as param_fh:
                for key, value in param_fh.items():
                    params[key]['value'] = value
        return params

    def __verify(params):
        """
        Using the syntax of default_parameters.py we can check the bounds on the provided
        parameters as well as their types. This is mainly to provide users feedback to 
        aid their own debugging and also provided some form of input sanity checking at 
        runtime. 

        Here we check that the types of the arguments are correct; however, future 
        iterations may expand upon this. Additionally, if the set of allowed
        values are provided, this will be enforced. 

        Args:
            params (dict): Asari-X params dict

        Raises:
            Warning: Invalid Type Encountered

        Returns:
            params: the Asari-X params dict
        """
        for k, v in [(k,v) for k,v in params.items() if not k.startswith('__')]:
            if not isinstance(v['value'], tuple(v['types'])):
                pass
                #print(type(v['value']))
                #raise Warning(f"WARNING: invalid type for {k}")
            if v.get("allowed", False):
                assert v['value'] in v["allowed"]
        return params

    def __reduce(params):
        """
        The full dictionary representation is no longer needed, we need only key value pairs
        for each field representing their ultimate values after configuration. We can discard
        the rest and keep the code simpler.

        Args:
            params (dict): Asari-X params dict (verbose)

        Returns:
            (dict): Concise Asari-X params dict
        """
        return {k: v.get('value', v) for k, v in params.items()}
    #print(json.dumps(__reduce(__verify(__combine(params, args))), indent=4))
    return __reduce(__verify(__combine(params, args)))

def main(params, dry_run=False):
    """
    This method uses the configuration provided in params to run the Asari-X job. 

    *All* UIs target this interface and this params configuration. Thus, this is 
    the entry point into Asari-X for all use cases, including, possibly the API. 

    The format of params can be found in a verbose format in default_parameters.py, 
    but by the time this function is called, you must first process the parameters 
    (see above function) to get into a simple key: value format. 

    The valid subcommands for Asari-X are defined by inspecting the functions defined
    under main. To get this list of functions, enable the dry_run argument. This is 
    complicated but not complex as it allows us to dynamically expand the UI, the 
    documentation, error messages, etc. as new functionality is added reducing 
    maintenance cost.

    Args:
        params (dict): Asari-X param dict
        dry_run (bool, optional): if true, return names of subcommands. Defaults to False.

    Returns:
        (set): names of valid subcommands if dry_run is True
        (bool): truth statement / exit code of requested command if dry_run was not False
    """
    # this runs Asari-X on a configured parameter payload
    # this function should be the target entry point for CLI / GUI / etc. 
    # there should be one function per allowed subcommand defined here.
    def build_signatures(params): 
        """
        This method takes a given set of compounds and reactions and expands them out to 
        a given reaction depth. The resulting products are organized based on the reactions
        that generated them and the compounds from which they are derived. The resulting
        signatures represent possible products given the input compounds and the reactions. 

        In other steps we can seach for and score these signatures to identify likely exposures
        to xenobiotic compounds. 

        Args:
            params (dict): the parameter dictionary, same payload for all methods
        """
        check_sufficient_params(params, ['compounds', 'reactions', 'signatures', 'reaction_depth'])
        SG = SignatureGenerator.from_compounds_reactions(params['compounds'], params['reactions'])
        SG.generate_signatures(reaction_depth=params['reaction_depth'])
        SG.save_signatures(signature_path=params['signatures'])

    def mzml_search(params): 
        """
        This method is one of two search functions in Asari-X. 

        The mzML_Searcher is configured using the params dictionary. Signatures must be provided, 
        either manually curated or generated automatically using the build_signatures command. 

        The input must be a path to either an individual mzML file or a directory with multiple 
        such files. 

        """
        check_sufficient_params(params, ['input', 'signatures'])
        if isinstance(params['signatures'], str) and params['signatures'].endswith('json'):
            params['signatures'] = json.load(open(params['signatures']))['data']
        XS = mzML_Searcher.from_params(params)
        XS.search()


        

    def ftable_search(params): 
        """
        placeholder
        """
        pass

    # Look at all callables in locals() and pick those whose qualname does not start with __"
    subcommand_func_objs = [x for x in locals().values() if callable(x) and not x.__name__.startswith("__")]
    if dry_run:
        logging.info(f"Dry Run of Asari-X main, this is normal")
        return {x.__name__ for x in subcommand_func_objs}
    logging.info(f"Starting Asari-X Run")
    logging.info(f"PARAMS for RUN:\n{json.dumps({k: v for k, v in params.items() if k not in {'compounds', 'reactions', 'signatures'}}, indent=4)}")
    assert params['run'] in {x.__name__ for x in subcommand_func_objs}
    return {x.__name__: x for x in subcommand_func_objs}[params["run"]](params)

def cli():
    """
    This function is called whenever Asari-X CLI is needed. This should only be ran from calling it 
    like a script, almost never programmatically. 
    """
    def __build_parser():
        """
        This function uses the specification in default_parameters.py to autopopulate the CLI

        Returns:
            parser (argparse.ArgumentParser): Asari-X ArgumentParser dynamically generated based on PARAMETERS
        """
        print(logo())
        parser = argparse.ArgumentParser(description="Asari-X, eXposome mining based on Asari")
        for key, parameter_spec in PARAMETERS.items():
            if parameter_spec.get("metavar", False):
                parser.add_argument(key, metavar=key, help=f"Valid Subcommands = {main({}, dry_run=True)}")
            else:
                if parameter_spec.get("short", None):
                    parser.add_argument(parameter_spec['short'], 
                                        f'--{key}', 
                                        default=parameter_spec['default'],
                                        type=parameter_spec['types'][0])
                else:
                    parser.add_argument(f'--{key}', 
                                        default=parameter_spec['default'],
                                        type=parameter_spec['types'][0])
        return parser

    parser = __build_parser()
    args = parser.parse_args()
    main(process_params(PARAMETERS, args))


if __name__ == '__main__':
    # when running Asari-X as a script, envoke the CLI.
    cli()