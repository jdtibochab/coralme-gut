def get_exchange(reactions,flux_dict):
    l = []
    for r in reactions:
        if flux_dict[r.id] <= 0:
            continue
        met = next(i for i in r.reactants)
        if "C" not in met.elements:
            continue
        l.append(r)
    return l

def get_overflow(model,solution=None):
    if solution is None:
        solution = model.solution.to_frame()
    df = solution
    # Only exchange
    df = df[df.index.str.contains("^EX_|^TS_")]
    # Only secretion
    df = df[df>1e-6].fillna(0)
    df = df.loc[[i for i in df.index if hasattr(model.reactions,i)]]
    df0 = df.copy()[["fluxes"]]

    # To C-mole
    mets = [next(i for i in model.reactions.get_by_id(r).reactants) for r in df.index]
    carbons = [i.elements.get("C",0) for i in mets]
    df = (df.T * carbons).T
    # Only active
    df = df[df>1e-3].fillna(0.)
    df = df[df.any(axis=1)]
    # Not co2
    df = df[~df.index.str.contains("co2\(e\)")]

    # Only this column
    df = df[["fluxes"]]
    df.columns = ["C-fluxes"]
    return df0.join(df).dropna()

def get_normalized_overflow(model,solution=None):
    df = get_overflow(model,solution=solution)
    df = df[df.index.str.contains("EX_")]
    df = df/df.sum().max()
    return df
