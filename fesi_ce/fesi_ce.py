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

    from ase.clease import Evaluate

    scond = [('converged','=',1), ('c1_0', '>', -0.05)]#, ('id','!=',15), ('id', '!=',16)]#, ('group','<',4), ('id','!=',15), ('id', '!=',16)] #(group, =, 0)

    evaluator = Evaluate(ceBulk, fitting_scheme="l1", parallel=False, alpha=1.29*1E-4,
    scoring_scheme="loocv_fast", max_cluster_dia=6, max_cluster_size=5, select_cond=scond)

    evaluator.plot_CV()

    evaluator.plot_fit(interactive=True, show_hull=False)
    #evaluator.plot_ECI()
    eci_names = evaluator.get_cluster_name_eci(return_type="dict")
    eci_fname = 'my_file_eci.json'
    with open(eci_fname,'w') as outfile:
        json.dump(eci_names, outfile, indent=2, separators=(",",":"))

    #plot_eci(eci_fname, show_zero_one_body=False)
    #investigate_eci(eci_fname)
    #print(evaluator.cf_matrix.dot(evaluator.eci))

def ret_evaluator(c1_0):

    conc = Concentration(basis_elements=[['Fe','Si']])
    ceBulk = CEBulk(crystalstructure='bcc', a=2.876, db_name='FeSi_8atoms_12finished.db',
    max_cluster_size=4, size=[2,2,2],
    max_cluster_dia=[0.0, 10, 10, 7.5, 5.5], concentration=conc, cubic=True)
    struct_generator = NewStructures(ceBulk, struct_per_gen=10)

    from ase.clease import Evaluate

    #scond = [('converged','=',1)]#, ('c1_0', '>', c1_0)]#, ('id','!=',15), ('id', '!=',16)]#, ('group','<',4), ('id','!=',15), ('id', '!=',16)] #(group, =, 0)

    evaluator = Evaluate(ceBulk, fitting_scheme="l1", parallel=False, alpha=1.29*1E-4,
    scoring_scheme="loocv_fast", max_cluster_dia=6, max_cluster_size=5)#, select_cond=scond)
    #evaluator.plot_fit()
    eci_names = evaluator.get_cluster_name_eci(return_type="dict")
    eci_fname = 'my_file_eci.json'
    with open(eci_fname,'w') as outfile:
        json.dump(eci_names, outfile, indent=2, separators=(",",":"))

    return evaluator

def ret_evaluator_bcs(c1_0):

    conc = Concentration(basis_elements=[['Fe','Si']])
    ceBulk = CEBulk(crystalstructure='bcc', a=2.876, db_name='FeSi_8atoms_12finished.db',
    max_cluster_size=4, size=[2,2,2],
    max_cluster_dia=[0.0, 10, 10, 7.5, 5.5], concentration=conc, cubic=True)
    struct_generator = NewStructures(ceBulk, struct_per_gen=10)

    from ase.clease import Evaluate

    from ase.clease import BayesianCompressiveSensing
    compressive = BayesianCompressiveSensing(noise=0.05)

    scond = [('converged','=',1), ('c1_0', '>', c1_0)]#, ('id','!=',15), ('id', '!=',16)]#, ('group','<',4), ('id','!=',15), ('id', '!=',16)] #(group, =, 0)

    evaluator = Evaluate(ceBulk, fitting_scheme=compressive, parallel=False, alpha=1.29*1E-4,
    scoring_scheme="loocv_fast", max_cluster_dia=6, max_cluster_size=5, select_cond=scond)#, select_cond=scond)
    #evaluator.plot_fit()
    eci_names = evaluator.get_cluster_name_eci(return_type="dict")
    eci_fname = 'my_file_eci.json'
    with open(eci_fname,'w') as outfile:
        json.dump(eci_names, outfile, indent=2, separators=(",",":"))

    return evaluator

def evaluate_GA(ceBulk):

    from ase.clease import GAFit, Evaluate

    scond = [('converged','=',1), ('c1_0', '>', -0.25)]#, ('id', '!=', '15'), ('id', '!=', '16')]#, ('group', '<=', 2)]

    #ga = GAFit(setting=ceBulk, alpha=1E-8, mutation_prob=0.01, num_individuals="auto",
    #fname="ga_fesi.csv", max_cluster_dia=6, cost_func='loocv', include_subclusters=False, select_cond=scond)

    #names = ga.run(min_change=0.001, gen_without_change=100)

    with open("ga_fesi_cluster_names.txt", 'r') as infile:
        lines = infile.readlines()
    names = [x.strip() for x in lines]

    #scond = [('converged','=',1), ('c1_0', '>', -0.4)]#, ('id','!=',15), ('id', '!=',16)]

    evaluator = Evaluate(ceBulk, cluster_names=names, fitting_scheme="l1", parallel=False, alpha=1.3*1E-4,
    scoring_scheme="loocv_fast", max_cluster_dia=6, max_cluster_size=7, select_cond=scond)

    #alp = evaluator.plot_CV()

    evaluator.plot_fit(interactive=True, show_hull=False)
    eci_names = evaluator.get_cluster_name_eci(return_type="dict")
    eci_fname = 'my_file_eci_ga.json'
    with open(eci_fname,'w') as outfile:
        json.dump(eci_names, outfile, indent=2, separators=(",",":"))

    #investigate_eci(eci_fname)
    #plot_eci(eci_fname, show_zero_one_body=False)

def plot_eci(name_of_json, show_zero_one_body=True):

    with open(name_of_json) as f:
        data = json.load(f)


    import matplotlib.pyplot as plt
    plt.figure(0)
    for key in data:
        if (show_zero_one_body):
            if (len(key) > 5):
                plt.bar(key[:9],data[key])
            else:
                plt.bar(key,data[key])
        else:
            if (len(key) > 5) and key != 'c0' and key != 'c1_0':
                plt.bar(key[:10],data[key])
            if (len(key) <= 5) and key != 'c0' and key != 'c1_0':
                plt.bar(key,data[key])

    plt.xticks(rotation=90)
    plt.ylabel('ECI [eV/atom]')
    plt.show()

def investigate_eci(name_of_json):

    negative_corr_val = 'Fe'
    positive_corr_val = 'Si'

    nn_configurations = []

    with open(name_of_json) as f:
        data = json.load(f)

    for z, key in enumerate(data):

        name = key
        name.split('_')
        number_of_atoms_in_cluster = int(name[1])
        num_conf = []
        if data[key] > 0:


            print(key)

            if number_of_atoms_in_cluster%2 == 0:

                for i in range(1,number_of_atoms_in_cluster,2):

                    num_conf.append(i)
            else:

                for i in range(1, number_of_atoms_in_cluster+1, 2):

                    num_conf.append(i)
        if data[key] < 0:

            print(key)

            if number_of_atoms_in_cluster%2 == 0:

                for i in range(0,number_of_atoms_in_cluster+1,2):

                    num_conf.append(i)

            else:

                for i in range(0, number_of_atoms_in_cluster, 2):

                    num_conf.append(i)

        nn_configurations.append({'key': key, 'energy': data[key], 'number of atoms': number_of_atoms_in_cluster, 'data': num_conf})

    for c in nn_configurations:

        for num_neg in c['data']:

            c[str(num_neg)] = negative_corr_val+str(num_neg) + positive_corr_val + str(c['number of atoms'] - num_neg)

    for c in nn_configurations:
        print(c)



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

    database = 'FeSi_27atoms_final_v3.db'
    #structure_ids = [(166,200,-15.324), (170,201,-30.649), (167,202,-33.294), (169,203,-30.178),
    #(174,204,-44.978), (175,205,-39.996), (177,206,-48.644), (173,207,-43.319), (172,208,-48.201), (168,209,-25.080)]
    #structure_ids = [(169,203,-30.178)]
    #structure_ids = [(187,216,-64.360), (181,217,-62.878), (194,218,-66.587), (195,219,-66.965), (197,220,-58.033)]
    #structure_ids = [(17,18,-205.707), (19,20,-177.696), (21,22,-172.908), (23,24,-173.178), (25,26,-238.484)]
    #structure_ids = [(188,221,-54.473), (190,222,-55.886), (189,223,-59.582), (193,224,-50.170), (199,225,-61.299)]
    #structure_ids = [(27,28,-168.263),(29,30,-197.639),(31,32,-235.868),(33,34,-229.364),(35,36,-238.484)]
    #structure_ids = [(191, 226, -59.257), (192,227, -49.689 ), (229,198,-61.304)]
    structure_ids = [(1,2,-223.448), (3,4,-225.704), (5,6,-223.011), (7,8,-203.672), (9,10,-216.070)]
    #cubic: from 215-220
    #27 new: 18-26
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

        insert_experimental_fesi_structure(struct_generator, data[2])

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

def check_databases():
    from ase.db import connect
    db1 = connect('FeSi_8atoms.db')
    db2 = connect('FeSi_16atoms.db')
    db3 = connect('FeSi_8atoms_12finished_v2.db')
    db4 = connect('FeSi_8atoms_12finished_cubic.db')
    db5 = connect('FeSi_27atoms_final.db')
    db_ce = connect('FeSi_8atoms_12finished.db')

    num_ids = 0

    for obj in db_ce.select(struct_type='final'):

        energy = obj['energy']
        volume = obj['volume']
        if not compare_dbs(db1, energy, volume) and not compare_dbs(db2, energy, volume) and not compare_dbs(db3, energy, volume) and not compare_dbs(db4, energy, volume) and not compare_dbs(db5, energy, volume):

            print('Did not find energy or volume for id: ' + str(obj['id']))
            num_ids += 1

    print(num_ids)


def compare_dbs(db, energy, volume):

    found_energy = False

    for obj in db.select(is_initial=0):

        if round(energy,3) == round(obj['energy'],3):
            #print(round(energy,3), round(obj['energy'],3))
            if (volume == obj['volume']):

                found_energy = True

    return found_energy

if __name__ == '__main__':
    #check_databases()
    #create_db_organize()
    #insert_structures()
    #create_xyz()
    #sym_check()
    main()
    #convex_hull()
