import os
from coralme.builder.main import MEBuilder
import anyconfig
import sys
import pandas as pd
def parse_diet(diet):
    d = {}
    for r,row in diet.iterrows():
        for i in r.split(","):
            d[i] = row["lb"]
    return d
    
def constrain_diet(model,dct):
    for r in model.reactions.query("EX_"):
        r.lower_bound = 0
    for r,lb in dct.items():
        if r not in model.reactions:
            print("{} not in model".format(r))
            continue
        print(r,lb)
        model.reactions.get_by_id(r).lower_bound = lb
def parser(args):
    org = args[1]
    param = {}
    if len(args)>2:
        for idx,a in enumerate(args):
            if '-' not in a: continue
            param[a] = args[idx+1]
    return org,param

def build_model(args):
    org,param = parser(args)
    config = anyconfig.load('input.json')
    dct = {
        'm-model-path' : './agora-models/AGORA_2_01_json_renamed/{:s}.json'.format(org),
        'genbank-path' : './agora-models/genbank-mixed-ncbi-agora-sources/{:s}.gb'.format(org),
        'ME-Model-ID' : "{}-ME".format(org),
        'out_directory' : "./me-models/{}/".format(org),
        'log_directory' : "./me-models/{}/".format(org),
        'df_gene_cplxs_mods_rxns' : "./me-models/{}/{}.xlsx".format(org,org)
    }
    config.update(dct)
    print(config)
    if '-g' not in param or param['-g'] == '1':
        if '-bbh' in param and param['-bbh'] == '0':
            config['run_bbh_blast'] = False
        else:
            config['run_bbh_blast'] = True
        builder = MEBuilder(*['parameters.json'], **config)
        builder.generate_files(overwrite = True)
        builder.save_builder_info()
    else:
        builder = MEBuilder(*['./me-models/{:s}/coralme-config.yaml'.format(org)])

    if '-b' not in param or param['-b'] == '1':   
        builder.build_me_model(overwrite = False)
    else:
        builder.me_model = builder.load(builder.configuration['out_directory']+'/MEModel-step2-{}-ME.pkl'.format(org))

    if '-ts' not in param or param['-ts'] == '1':
        diet = pd.read_csv("./diets/high_fiber_diet.txt",index_col=0,sep='\t',comment='#',header=None)
        diet.columns = ["lb"]
        diet_dct = parse_diet(diet)
        constrain_diet(builder.me_model,diet_dct)
        builder.troubleshoot(growth_key_and_value = { builder.me_model.mu : config.get("feasibility_mu",0.001) })

    if param.get("-return_builder",None) is not None:
        return builder