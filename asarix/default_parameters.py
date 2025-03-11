"""
This defines the set of all parameters that can be required by Asari-X.

They are defined in this verbose format to allow for type and value checking
to be performed on the input values and to set sane defaults. Additionally, 
by organizing the parameters in this manner, we can automatically populate
the CLI and GUI, making the code less hard to maintain.

The format is key, value pairs where the key is the full name of the parameter. 
The value is a dictionary containing one or more of the following key value pairs:

default - required, this value will be used if a user does not provide one
types - required, the allowed types for the value
short - optional, the shorthand name for the option to be used in the CLI
allowed - optional, if provided only values in this list are permitted.
metavar - optional, if provided, the value is a metavar. Only used for run
help - optional, if provided, this is the help text used in the CLI for the 
    parameter's description
"""

PARAMETERS = {
    "run": {
        "default": None,
        "types": [str],
        "metavar": True,
    },
    "parameters": {
        "default": None,
        "types": [str, type(None)],
        "short": "-p",
        "help": "path to parameters.json"
    },
    "mode": {
        "default": "pos",
        "types": [str],
        "allowed": ["pos", "neg", "mixed"],
        "short": '-m',
    },
    "mz_tolerance_ppm": {
        "default": 10,
        "types": [int],
        "short": '-z',
    },
    "input": {
        "default": None,
        "types": [str, type(None)],
        "short": '-i'
    },
    "compounds": {
        "default": None,
        "types": [str, type(None)],
        "short": '-c'
    },
    "reactions": {
        "default": None,
        "types": [str, type(None)],
        "short": '-r'
    },
    "reaction_depth": {
        "default": 3,
        "types": [int, type(None)],
        "short": '-d'
    },
    "signatures": {
        "default": None,
        "types": [str, type(None)],
        "short": '-s',
        "skip_json": True
    },
    "snr_cutoff": {
        "default": 2.5,
        "types": [float, int],
    },
    "scan_cutoff": {
        "default": 0,
        "types": [int],
    }
}