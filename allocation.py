from collections import defaultdict
import pandas
import coralme
def get_allocation(model):
    alloc = defaultdict(int)
    for r in model.reactions.query("^translation_"):
        gid = r.id.split("translation_")[1]
        #if gid == "dummy":
            #continue
        fl = model.solution.fluxes[r.id]
        for fu in model.metabolites.get_by_id("RNA_"+gid).functions:
            alloc[fu] += fl
    alloc = pandas.DataFrame.from_dict({"mass": alloc})
    return alloc[alloc["mass"] > 1e-16].sort_values("mass")

from coralme.builder.helper_functions import flux_based_reactions
def get_partition(m):
    df = flux_based_reactions(m._model,m.id)["met_flux"]
    df = df[df<0]
    df = df/df.sum()
    dct = df.to_dict()
    return {m._model.get(r):flux for r,flux in dct.items()}

from coralme.builder.helper_functions import get_next_from_type,substitute_value
def find_complexes(m, seen = dict(),multiplier=1):
    #print(1,type(m))
    #print(m,multiplier)
    if not m:
        return dict()
#     print(m.id)
    if m in seen:
        return dict()
    seen.update({m:multiplier})

    # Metabolite objects
    if isinstance(m,coralme.core.component.MEComponent):
        partition = get_partition(m)
        #print(2,partition)
        if isinstance(m,coralme.core.component.TranslatedGene):
            #print(1)
            cplxs = dict()
            for r,fraction in partition.items():
                cplxs.update(find_complexes(r, seen=seen, multiplier=multiplier*fraction))
            return cplxs
        if isinstance(m,coralme.core.component.TranscribedGene):
            #print(2)
            translated_protein = m.id.replace('RNA_','protein_')
            if translated_protein in m._model.metabolites:
                return find_complexes(m._model.metabolites.get_by_id(translated_protein), seen=seen, multiplier=multiplier)
            cplxs = dict()
            for r,fraction in partition.items():
                cplxs.update(find_complexes(r, seen=seen, multiplier=multiplier*fraction))
            return cplxs
        if isinstance(m,coralme.core.component.ProcessedProtein):
            #print(3)
            cplxs = dict()
            for r,fraction in partition.items():
                cplxs.update(find_complexes(r, seen=seen, multiplier=multiplier*fraction))
            return cplxs

        if isinstance(m,coralme.core.component.Complex) or isinstance(m,coralme.core.component.GenericComponent) or isinstance(m,coralme.core.component.GenerictRNA):
            #print(4)
            #print(partition)
            other_formations = {r:f for r,f in partition.items() if isinstance(r,coralme.core.reaction.ComplexFormation) or isinstance(r,coralme.core.reaction.GenericFormationReaction)  or isinstance(r,coralme.core.reaction.tRNAChargingReaction)}
            cplxs = {m:multiplier}
            if other_formations:
                for r,fraction in other_formations.items():
                    cplxs.update(find_complexes(r, seen=seen, multiplier=multiplier*fraction))
                del cplxs[m]
            return cplxs
    #         return set([m])

    # Reaction objects
    if isinstance(m,coralme.core.reaction.MEReaction):
        if isinstance(m,coralme.core.reaction.PostTranslationReaction):
            return find_complexes(get_next_from_type(m.metabolites,coralme.core.component.ProcessedProtein), seen=seen, multiplier=multiplier)
        if isinstance(m,coralme.core.reaction.ComplexFormation):
            return find_complexes(get_next_from_type(m.metabolites,coralme.core.component.Complex), seen=seen, multiplier=multiplier)
        if isinstance(m,coralme.core.reaction.GenericFormationReaction):
            return find_complexes(get_next_from_type(m.metabolites,coralme.core.component.GenericComponent), seen=seen, multiplier=multiplier)
        if isinstance(m,coralme.core.reaction.tRNAChargingReaction):
            return find_complexes(get_next_from_type(m.metabolites,coralme.core.component.GenerictRNA), seen=seen, multiplier=multiplier)
        if isinstance(m,coralme.core.reaction.MetabolicReaction):
            tmp1 = find_complexes(get_next_from_type(m.metabolites,coralme.core.component.Complex), seen=seen, multiplier=multiplier)
            tmp2 = find_complexes(get_next_from_type(m.metabolites,coralme.core.component.GenericComponent), seen=seen, multiplier=multiplier)
            tmp1.update(tmp2)
            return tmp1

    return dict()


def get_subsystem(r):
    if isinstance(r,coralme.core.reaction.MetabolicReaction) and hasattr(r,'subsystem'):
        if r.subsystem:
            return r.subsystem
        return "No_subsystem"
    if isinstance(r,coralme.core.reaction.TranslationReaction):
        return "Translation"
    elif isinstance(r,coralme.core.reaction.TranscriptionReaction):
        return "Transcription"
    elif isinstance(r,coralme.core.reaction.tRNAChargingReaction):
        return "tRNA-Charging"
    elif isinstance(r,coralme.core.reaction.PostTranslationReaction):
        return "Post-Translation"
    elif isinstance(r,coralme.core.reaction.SummaryVariable):
        return "Biomass"
