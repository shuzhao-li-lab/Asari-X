import argparse
import json
from default_parameters import PARAMETERS

import logging
from logger_setup import setup_logger
setup_logger()
logger = logging.getLogger(__name__)

LOGO = """ 
###############################
#                             #
#           ASARI-X           #
#                             #
###############################
"""

def process_params(params, args=None):
    def __combine(params, args=None):
        if args is not None:
            for key, value in args.__dict__.items():
                params[key]['value'] = value
        if params.get('parameters', False):
            with open(params['parameters']) as param_fh:
                for key, value in param_fh.items():
                    params[key]['value'] = value
        return params

    def __verify(params):
        for k, v in params.items():
            if not isinstance(v['value'], tuple(v['types'])):
                print(f"WARNING: invalid type for {k}")
            if v.get("allowed", False):
                assert v['value'] in v["allowed"]
        return params

    def __reduce(params):
        return {k: v['value'] for k, v in params.items()}
    
    return __reduce(__verify(__combine(params, args)))

def main(params, dry_run=False):
    # this runs Asari-X on a configured parameter payload
    # this function should be the target entry point for CLI / GUI / etc. 
    # there should be one function per allowed subcommand defined here.
    def build_signatures(params): 
        print(params)
        from signature_generator import SignatureGenerator
        SG = SignatureGenerator.from_compounds_reactions(params['compounds'], params['reactions'])
        signatures = SG.generate_signatures(reaction_depth=params['reaction_depth'])

    def mzml_search(params): 
        pass

    def ftable_search(params): 
        pass

    # Look at all callables in locals() and pick those whose qualname does not start with __"
    subcommand_func_objs = [x for x in locals().values() if callable(x) and not x.__name__.startswith("__")]
    if dry_run:
        return {x.__name__ for x in subcommand_func_objs}
    logging.info(f"PARAMS for RUN:\n{json.dumps(params, indent=4)}")
    assert params['run'] in {x.__name__ for x in subcommand_func_objs}
    return {x.__name__: x for x in subcommand_func_objs}[params["run"]](params)

def cli():
    """
    This function is called whenever Asari-X CLI is needed. This should only be ran from calling it 
    like a script, almost never programmatically. 

    Returns:
        _type_: _description_
    """
    def __build_parser():
        """
        This function uses the specification in default_parameters.py to autopopulate the CLI

        Returns:
            _type_: _description_
        """
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
        return parser.parse_args()
    return main(process_params(PARAMETERS, __build_parser()))

if __name__ == '__main__':
    cli()