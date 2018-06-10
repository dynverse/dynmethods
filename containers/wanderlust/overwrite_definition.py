from yaml import load, dump

definition = load(open("/code/definition.yml", "r"))

definition["name"] = "wanderlust"
definition["parameters"]["branch"]["default"] = False

dump(definition, open("/code/definition.yml", "w"))
