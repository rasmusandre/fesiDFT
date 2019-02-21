from ase.clease import CEBulk
from ase.clease import Concentration
from ase.clease import NewStructures
from ase.clease import BayesianCompressiveSensing
import numpy as np
import json
import time

def main():

    conc = Concentration(basis_elements=[['Fe','Si']])
    ceBulk = CEBulk(crystalstructure='bcc', a=2.876, db_name='FeSi_8atoms_12finished.db',
    max_cluster_size=4, size=[2,2,2],
    max_cluster_dia=[0.0, 10, 10, 7.5, 5.5], concentration=conc, cubic=True)
    #print(ceBulk.cluster_info_given_size(4))
    #ceBulk.reconfigure_settings()
    #ceBulk.view_clusters()
    struct_generator = NewStructures(ceBulk, struct_per_gen=10)
    #print(ceBulk.basis_functions)
    #reconfigure(ceBulk)
    evaluate(ceBulk)


def reconfigure(ceBulk):

    from ase.clease import CorrFunction
    cf = CorrFunction(ceBulk, parallel=True)
    cf.reconfigure_db_entries()


def evaluate(ceBulk):

    from ase.clease import Evaluate, BayesianCompressiveSensing
    compressive = BayesianCompressiveSensing(noise=0.05)#, variance_opt_start=0, lamb_opt_start=0)

    scond = [('converged','=',1), ('c1_0', '>', -0.05)]

    evaluator = Evaluate(ceBulk, fitting_scheme=compressive, parallel=False, alpha=1.29*1E-4,
    scoring_scheme="loocv_fast", max_cluster_dia=6, max_cluster_size=5, select_cond=scond)

    #evaluator.plot_CV()

    evaluator.plot_fit(interactive=True)#, show_hull=False)


    eci_names = evaluator.get_cluster_name_eci(return_type="dict")
    eci_fname = 'my_file_eci.json'
    with open(eci_fname,'w') as outfile:
        json.dump(eci_names, outfile, indent=2, separators=(",",":"))


def evaluate_GA(ceBulk):

    from ase.clease import GAFit, Evaluate

    scond = [('converged','=',1), ('c1_0', '>', -0.25)]

    #ga = GAFit(setting=ceBulk, alpha=1E-8, mutation_prob=0.01, num_individuals="auto",
    #fname="ga_fesi.csv", max_cluster_dia=6, cost_func='loocv', include_subclusters=False, select_cond=scond)

    #names = ga.run(min_change=0.001, gen_without_change=100)

    with open("ga_fesi_cluster_names.txt", 'r') as infile:
        lines = infile.readlines()
    names = [x.strip() for x in lines]

    evaluator = Evaluate(ceBulk, cluster_names=names, fitting_scheme="l1", parallel=False, alpha=1.3*1E-4,
    scoring_scheme="loocv_fast", max_cluster_dia=6, max_cluster_size=7, select_cond=scond)

    evaluator.plot_CV()

    evaluator.plot_fit(interactive=True, show_hull=False)
    eci_names = evaluator.get_cluster_name_eci(return_type="dict")
    eci_fname = 'my_file_eci_ga.json'
    with open(eci_fname,'w') as outfile:
        json.dump(eci_names, outfile, indent=2, separators=(",",":"))


main()
