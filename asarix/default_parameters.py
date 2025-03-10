# format is key: [default, type, short, long, ]

PARAMETERS = {
    "run": {
        "default": None,
        "types": [str],
        "metavar": True,
    },
    "mode": {
        "default": "pos",
        "types": [str],
        "allowed": ["pos", "neg", "mixed"],
        "short": '-m',
    },
    "input": {
        "default": None,
        "types": [str],
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
    }
}