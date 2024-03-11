import re
def filter_biomass(biomass_dct):
    exclude_lst = ['atp_c', 'dnarep_c', 'proteinsynth_c', 'rnatrans_c', 'ACP_c','PGP_c','h2o_c','gly_c']
    regex = "__L_c|d?[a,u,g,t,c][t,d]p_c|^[a-z]{1,2}\d?_c"
    dct = {}
    for k,v in biomass_dct.items():
        if v > 0:
            #print(1,k)
            continue
        if k in exclude_lst:
            #print(2,k)
            continue
        if regex and re.search(regex,k):
            #print(3,k)
            continue
        dct[k] = v
    return dct

import coralme
def correct_biomass(model,biomass_constituents):
    biomass_constituents = filter_biomass(biomass_constituents)
    model.global_info['flux_of_biomass_constituents'] = biomass_constituents
    # remove metabolites not in the model or without molecular weight
    biomass_constituents = { k:v for k,v in biomass_constituents.items()
                         if model.metabolites.has_id(k) and model.metabolites.get_by_id(k).formula_weight }

    problems = list(set(model.global_info.get('flux_of_biomass_constituents', {})).difference(biomass_constituents))
    if problems:
        print('The following biomass constituents are not in the ME-model or have no formula: {:s}.'.format(', '.join(problems)))
    model.remove_reactions(['biomass_constituent_demand'], remove_orphans = False)
    rxn = coralme.core.reaction.SummaryVariable('biomass_constituent_demand')
    model.add_reactions([rxn])
    rxn.add_metabolites({ k:-(abs(v)) for k,v in biomass_constituents.items() })
    rxn.lower_bound = model.mu # coralme.util.mu
    rxn.upper_bound = model.mu # coralme.util.mu
    constituent_mass = sum([model.metabolites.get_by_id(c).formula_weight / 1000. * abs(v) for c,v in biomass_constituents.items()])
    rxn.add_metabolites({model.metabolites.get_by_id('constituent_biomass'): constituent_mass})
