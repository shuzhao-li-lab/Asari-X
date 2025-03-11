# Parmameters and Params

The payload to Asari-X is the parameters dictionary. All keys needed in that dictionary should be defined in parameters.py. This file is used to autogenerate the CLI and the GUI. The "run" field of the parameters file specifies which workflow or command to generate. These are defined within main.main, and the set of allowed commands autopopulated by dry_running the main method. This seems complicated but allows us to build a scalable and maintainable UI with minimal code changes. 

The global PARAMETERS is the default params configuration. When passed to a function, we call it 'params' to clarify that this instance of the parameters is not the global one. 

# Logging

Each component of Asari-X should import the logger from logger_setup.py. By doing this, all logs are unified and in one location. Avoid printing to terminal debug information when we can toggle it using logging instead. 