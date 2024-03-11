
import coralme
def load_model(org):
    return coralme.io.pickle.load_pickle_me_model("./me-models/{}/MEModel-step3-{}-ME-TS.pkl".format(org,org))

def get_bounds(solution,reactions,multiplier = 1e-3):
    dct = {}
    for reaction in reactions:
        flux = solution[reaction]
        bound = flux * multiplier
        dct[reaction] = bound
    return dct

def constrain_exchanges(model,dct):
    for r,v in dct.items():
        model.reactions.get_by_id(r).lower_bound = v
