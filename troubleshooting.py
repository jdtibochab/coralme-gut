def open_bounds(model,met):
    rxnid = "TS_{}".format(met)
    if rxnid in model.reactions:
        rxn = model.get(rxnid)
    else:
        rxn = coralme.core.reaction.MEReaction(rxnid)
        model.add_reaction(rxn)
        rxn.add_metabolites({met:-1})
    rxn.bounds = (-1000,1000)

def close_bounds(model,met):
    rxnid = "TS_{}".format(met)
    rxn = model.get(rxnid)
    rxn.bounds = (0,0)
