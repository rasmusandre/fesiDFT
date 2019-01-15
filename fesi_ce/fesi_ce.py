from ase.clease import CEBulk
from ase.clease import Concentration
from ase.clease import NewStructures
from ase.clease import BayesianCompressiveSensing
import numpy as np
import json

def main():

    conc = Concentration(basis_elements=[['Fe','Si']])
    ceBulk = CEBulk(crystalstructure='bcc', a=2.876, db_name='FeSi_8atoms_12finished.db',
    max_cluster_size=4, size=[2,2,2],
    max_cluster_dia=[0.0, 10, 10, 7.5, 5.5], concentration=conc, cubic=True)
    #print(ceBulk.cluster_info_given_size(4))
    #ceBulk.reconfigure_settings()
    #ceBulk.view_clusters()
    struct_generator = NewStructures(ceBulk, struct_per_gen=10)

    #reconfigure(ceBulk)
    evaluate(ceBulk)
    #insert_experimental_fesi_structure(struct_generator)

def reconfigure(ceBulk):

    from ase.clease import CorrFunction
    cf = CorrFunction(ceBulk, parallel=True)
    cf.reconfigure_db_entries()


def evaluate(ceBulk):

    from ase.clease import GAFit, Evaluate
    ga = GAFit(setting=ceBulk, alpha=1E-8, mutation_prob=0.01, num_individuals="auto",
    change_prob=0.2, fname="ga_fesi.csv", max_cluster_dia=6)
    names = ga.run(min_change=0.001, gen_without_change=100)
    with open("ga_fesi_cluster_names.txt", 'r') as infile:
        lines = infile.readlines()
    names = [x.strip() for x in lines]
    compressive = BayesianCompressiveSensing(noise=0.1)
    evaluator = Evaluate(ceBulk, fitting_scheme="l2", parallel=False, alpha=1E-8,
    scoring_scheme="loocv_fast", cluster_names=names)
    #evaluator = Evaluate(ceBulk, fitting_scheme=compressive, parallel=False,
    #scoring_scheme="loocv_fast")
    #x = evaluator.cf_matrix
    #k = x.T.dot(x)
    #print(k.shape)
    #print(np.linalg.matrix_rank(x.T.dot(x)))
    evaluator.plot_fit(interactive=True)
    eci_names = evaluator.get_cluster_name_eci(return_type="dict")

    with open(eci_fname,'w') as outfile:
        json.dump(eci_name, outfile, indent=2, separators=(",",":"))

def insert_experimental_fesi_structure(struct_gen):

    from ase.calculators.singlepoint import SinglePointCalculator
    from ase.io import read

    init_structure = read('initial.xyz')
    final_structure = read('final.xyz')

    calc = SinglePointCalculator(final_structure, energy=-241.049)
    final_structure.set_calculator(calc)
    struct_gen.insert_structure(init_struct=init_structure, final_struct=final_structure, generate_template=True)


def create_xyz():

    from ase.db import connect
    from ase.io import write

    db = connect('FeSi_8atoms.db')

    initial = db.get_atoms(id=29)
    final = db.get_atoms(id=30)

    initial.write('initial.xyz')
    final.write('final.xyz')





main()
#create_xyz()
