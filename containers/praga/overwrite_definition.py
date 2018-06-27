from yaml import load, dump

definition = load(open("/code/definition.yml", "r"))
definition["name"] = "projected PAGA"
definition["short_name"] = "praga"
definition["output"]["outputs"] = ["dimred_projection", "dimred", "timings"]

dump(definition, open("/code/definition.yml", "w"))
