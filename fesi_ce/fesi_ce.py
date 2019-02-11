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
    evaluate_GA(ceBulk)
    #insert_experimental_fesi_structure(struct_generator)

def reconfigure(ceBulk):

    from ase.clease import CorrFunction
    cf = CorrFunction(ceBulk, parallel=True)
    cf.reconfigure_db_entries()


def evaluate(ceBulk):

    from ase.clease import GAFit, Evaluate
    # ga = GAFit(setting=ceBulk, alpha=1E-8, mutation_prob=0.01, num_individuals="auto",
    # change_prob=0.2, fname="ga_fesi.csv", max_cluster_dia=6)
    # names = ga.run(min_change=0.001, gen_without_change=100)

    scond = [('converged','=',1), ('c1_0', '>', -0.4)]#, ('id','!=',15), ('id', '!=',16)]#, ('group','<',4), ('id','!=',15), ('id', '!=',16)] #(group, =, 0)



    with open("ga_fesi_cluster_names.txt", 'r') as infile:
        lines = infile.readlines()
    names = [x.strip() for x in lines]
    #compressive = BayesianCompressiveSensing(noise=0.1)
    # evaluator = Evaluate(ceBulk, fitting_scheme="l2", parallel=False, alpha=1E-8,
    # scoring_scheme="loocv_fast", cluster_names=names)
    evaluator = Evaluate(ceBulk, fitting_scheme="l1", parallel=False, alpha=1.3*1E-4,
    scoring_scheme="loocv_fast", max_cluster_dia=5, max_cluster_size=4, select_cond=scond)
    evaluator.plot_CV()
    #for c1_0 >-0.4, alpha=1.3E-4, dia=5, size=4
    #for group<3 alpha 0.2*1E-4
    #for l2, alpha=1E-2, max_dia=5, max_size=3
    #for group=0, alpha=1.3E-4, dia=6, size=4
    #for group=0,1, alpha=4.6E-4, dia=6, size=3
    #evaluator = Evaluate(ceBulk, fitting_scheme=compressive, parallel=False,
    #scoring_scheme="loocv_fast")
    #x = evaluator.cf_matrix
    #k = x.T.dot(x)
    #print(k.shape)
    #print(np.linalg.matrix_rank(x.T.dot(x)))
    evaluator.plot_fit(interactive=False)
    eci_names = evaluator.get_cluster_name_eci(return_type="dict")

    with open(eci_fname,'w') as outfile:
        json.dump(eci_name, outfile, indent=2, separators=(",",":"))

def convex_hull():
    from ase.clease import ConvexHull
    import matplotlib.pyplot as plt
    ch = ConvexHull('FeSi_8atoms_12finished.db')
    ch.plot()
    plt.show()

def evaluate_GA(ceBulk):

    from ase.clease import GAFit, Evaluate

    scond = [('converged','=',1)], ('c1_0', '>', -0.4)]

    ga = GAFit(setting=ceBulk, alpha=1E-8, mutation_prob=0.01, num_individuals="auto",
    fname="ga_fesi.csv", max_cluster_dia=6, include_subclusters=False, select_cond=scond)
    #ga = GAFit(setting=ceBulk, parallel=True)
    names = ga.run(min_change=0.001, gen_without_change=100)

    with open("ga_fesi_cluster_names.txt", 'r') as infile:
        lines = infile.readlines()
    names = [x.strip() for x in lines]

    #scond = [('converged','=',1), ('c1_0', '>', -0.4)]#, ('id','!=',15), ('id', '!=',16)]

    evaluator = Evaluate(ceBulk, cluster_names=names, fitting_scheme="l1", parallel=False, alpha=1.3*1E-4,
    scoring_scheme="loocv_fast", max_cluster_dia=5, max_cluster_size=4, select_cond=scond)

    #evaluator = Evaluate(ceBulk, fitting_scheme=compressive, parallel=False,
    #scoring_scheme="loocv_fast")

    evaluator.plot_fit(interactive=True)
    eci_names = evaluator.get_cluster_name_eci(return_type="dict")

    with open(eci_fname,'w') as outfile:
        json.dump(eci_name, outfile, indent=2, separators=(",",":"))

def insert_experimental_fesi_structure(struct_gen, struct_energy):

    from ase.calculators.singlepoint import SinglePointCalculator
    from ase.io import read

    init_structure = read('initial.xyz')
    final_structure = read('final.xyz')

    calc = SinglePointCalculator(final_structure, energy=struct_energy)
    final_structure.set_calculator(calc)
    struct_gen.insert_structure(init_struct=init_structure, final_struct=final_structure, generate_template=True)


def create_xyz(database_name, initial_id, final_id):

    from ase.db import connect
    from ase.io import write

    db = connect(database_name)

    initial = db.get_atoms(id=initial_id)
    final = db.get_atoms(id=final_id)

    initial.write('initial.xyz')
    final.write('final.xyz')

def insert_structures():

    database = 'FeSi_8atoms_12finished_cubic.db'
    #structure_ids = [(166,200,-15.324), (170,201,-30.649), (167,202,-33.294), (169,203,-30.178),
    #(174,204,-44.978), (175,205,-39.996), (177,206,-48.644), (173,207,-43.319), (172,208,-48.201), (168,209,-25.080)]
    structure_ids = [(169,203,-30.178)]
    #from ase.db import connect
    #from ase.io import write
    #db = connect('FeSi_8atoms_12finished.db')
    #del db[179]
    #del db[180]
    for data in structure_ids:
        create_xyz(database, data[0], data[1])


        conc = Concentration(basis_elements=[['Fe','Si']])
        ceBulk = CEBulk(crystalstructure='bcc', a=2.876, db_name='FeSi_8atoms_12finished.db',
        max_cluster_size=4, size=[2,2,2],
        max_cluster_dia=[0.0, 10, 10, 7.5, 5.5], concentration=conc, cubic=True)

        struct_generator = NewStructures(ceBulk, struct_per_gen=10)

        #insert_experimental_fesi_structure(struct_generator, data[2])

def create_db_organize():
    #create CE db from FeSi_8atoms and FeSi_16atoms
    from ase.db import connect
    from ase.io import write

    db_name = 'FeSi_16atoms.db'
    db = connect(db_name)
    start_id = 1
    end_id = 4

    conc = Concentration(basis_elements=[['Fe','Si']])
    ceBulk = CEBulk(crystalstructure='bcc', a=2.876, db_name='FeSi_test.db',
    max_cluster_size=4, size=[2,2,2],
    max_cluster_dia=[0.0, 10, 10, 7.5, 5.5], concentration=conc, cubic=True)
    ceBulk.reconfigure_settings()

    struct_generator = NewStructures(ceBulk, struct_per_gen=10)


    while start_id < end_id + 1:
        local_end_id = start_id + 1

        for obj in db.select(id=local_end_id):

            energy = obj['energy']
        print((start_id, local_end_id))
        create_xyz(db_name, start_id, local_end_id)
        insert_experimental_fesi_structure(struct_generator, energy)
        start_id += 2


def sym_check():

    from ase.utils.structure_comparator import SymmetryEquivalenceCheck

    from ase.db import connect

    db1 = connect('FeSi_8atoms.db')
    db2 = connect('FeSi_8atoms_12finished.db')

    init = db1.get_atoms(id=61)
    ids = []
    for k in db2.select():
        ids.append(k.id)

    for id in ids:

        fin_12 = db2.get_atoms(id=id)
        symcheck = SymmetryEquivalenceCheck()

        if (symcheck.compare(init,fin_12)):

            print(id)

#create_db_organize()
#insert_structures()
#create_xyz()
#sym_check()
main()
#convex_hull()
