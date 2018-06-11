from yaml import load, dump

definition = load(open("/code/definition.yml", "r"))
definition["name"] = "praga"
definition["output"]["outputs"] = ["dimred_projection"]

dump(definition, open("/code/definition.yml", "w"))
