from coralme.solver.solver import ME_NLP
import glob
import numpy
import pandas
import seaborn
import matplotlib.pyplot as plt
from tqdm import tqdm
import multiprocessing as mp
import os
from build import build_model
import coralme
from tqdm import tqdm
import cobra

def get_nlp(model):
    Sf, Se, lb, ub, b, c, cs, atoms, lambdas = model.construct_lp_problem(lambdify=True)
    me_nlp = ME_NLP(Sf, Se,b, c, lb, ub,  cs, atoms, lambdas)
    return me_nlp

def get_feasibility(me_nlp, basis=None):
    x_new,y_new,z_new,stat_new,hs_new = me_nlp.solvelp(0.001,basis,'quad')
    return (True,hs_new) if stat_new=="optimal" else (False,hs_new)

def knockout(index_dct,met,me_nlp,limit=lambda x:0):
    # met = "RNA_592010.4.peg.10"
    # rxns = [r.id for r in model.metabolites.get_by_id(met).reactions if "transcription" in r.id]
    # m_idx = model.metabolites.index(met)
    rxn = "transcription_TU_{}".format(met)
    r_idx = index_dct.get(rxn)
    me_nlp.xu[r_idx] = limit
    
    # for ri in r_idxs:
        # print(tmp.Sf[(m_idx,ri)])
        # me_nlp.Sf[(m_idx,ri)] = 0

def restore(index_dct,met,me_nlp):
    # met = "RNA_592010.4.peg.10"
    # rxns = [r.id for r in model.metabolites.get_by_id(met).reactions if "transcription" in r.id]
    # m_idx = model.metabolites.index(met)
    rxn = "transcription_TU_{}".format(met)
    r_idx = index_dct.get(rxn)
    me_nlp.xu[r_idx] = lambda x:1000
    
    # for ri in r_idxs:
        # print(tmp.Sf[(m_idx,ri)])
        # me_nlp.Sf[(m_idx,ri)] = 0
        
def printcoeff(model,met,me_nlp):
    # met = "RNA_592010.4.peg.10"
    # rxns = [r.id for r in model.metabolites.get_by_id(met).reactions if "transcription" in r.id]
    rxns = [model.get("transcription_TU_{}".format(met))]
    
    # m_idx = model.metabolites.index(met)
    r_idxs = [model.reactions.index(r) for r in rxns ]
    for ri in r_idxs:
        # print(tmp.Sf[(m_idx,ri)])
        print(me_nlp.Sf[(m_idx,ri)])

def get_targets(model,condition,NormalizedTCounts):
    model_genes = [i.id.split("RNA_")[1] for i in model.all_genes]
    return [g for g in model_genes if g not in NormalizedTCounts.index or NormalizedTCounts.loc[g][condition]<1e-16]

def get_killable(index_dct,me_nlp,basis,kill_genes,limit=lambda x:0,ListHandler=None,org="Model"):
    killable = []
    for idx,k in enumerate(kill_genes):
        # k = "RNA_" + k
        knockout(index_dct,k,me_nlp,limit=limit)
        f,new_basis = get_feasibility(me_nlp,basis=basis)
        if f:
            if ListHandler:ListHandler.print_and_log("{} feasible if {} knockout ({} out of {})".format(org,k,idx+1,len(kill_genes)))
            killable.append(k)
            basis = new_basis
            continue
        if ListHandler:ListHandler.print_and_log("{} not feasible if {} knockout ({} out of {})".format(org,k,idx+1,len(kill_genes)))
        restore(index_dct,k,me_nlp)
    if ListHandler:ListHandler.print_and_log("Done with {}".format(org))
    return killable

def optimize(index_dct,met_index_dct,me_nlp,max_mu = 2.8100561374051836, min_mu = 0., maxIter = 100, lambdify = True,
		tolerance = 1e-6, precision = 'quad', verbose = True):
    muopt, xopt, yopt, zopt, basis, stat = me_nlp.bisectmu(
				mumax = max_mu,
				mumin = min_mu,
				maxIter = maxIter,
				tolerance = tolerance,
				precision = precision,
				verbose = verbose)

    if stat == 'optimal':
        #f = sum([ rxn.objective_coefficient * xopt[idx] for idx, rxn in enumerate(self.reactions) ])
        x_primal = xopt[ 0:len(index_dct) ]   # The remainder are the slacks
        x_dict = { rxn : xopt[idx] for rxn,idx in index_dct.items() }
        #y = pi
        # J = [S; c]
        y_dict = { met : yopt[idx] for met,idx in met_index_dct.items() }
        z_dict = { rxn : zopt[idx] for rxn,idx in index_dct.items() }
        #y_dict['linear_objective'] = y[len(y)-1]

        #self.me.solution = Solution(f, x_primal, x_dict, y, y_dict, 'qminos', time_elapsed, status)
        return cobra.core.Solution(
            objective_value = muopt,
            status = stat,
            fluxes = x_dict, # x_primal is a numpy.array with only fluxes info
            reduced_costs = z_dict,
            shadow_prices = y_dict,
            )
    else:
        return "infeasible"