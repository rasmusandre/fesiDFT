from ase.build import bulk
from ase import Atom, Atoms
from ase.visualize import view
from ase.db import connect
from ase.utils.structure_comparator import SymmetryEquivalenceCheck
from itertools import product
from random import randint

def prepare_structures():

    atoms = bulk('Fe', 'bcc', a = 2.876)*(3,3,3)
    symbols = ['Fe', 'Si']
    #SET MAGMOM!
    db = connect('FeSi_8atoms_initial.db')
    for s in product(symbols, repeat = 8):
        for atom in atoms:
            atom.symbol = s[atom.index]
        if not exists_in_db(atoms):
            db.write(atoms, struct_type = 'initial')

def prepare_random_structure():

    num_fe = 0
    num_si = 10

    while num_si/(num_si + num_fe) > 0.4:
        num_fe = randint(0,27)
        num_si = 27 - num_fe

    return [num_fe, num_si]

def prepare_random_structures():


    db = connect('FeSi_27atoms_initial_v3.db')
    added_to_db = 0

    while added_to_db < 1000:
        atoms = bulk('Fe', 'bcc', a = 2.876)*(3,3,3)
        atom_config = prepare_random_structure()

        num_si = 0
        while num_si < atom_config[1]:
            for atom in atoms:
                if atom.symbol != 'Si' and randint(0,10)/10 <= atom_config[1]/27:
                    atom.symbol = 'Si'
                    num_si += 1


        if not exists_in_db(atoms, 'FeSi_27_atoms_initial_v3.db') and not exists_in_db(atoms, 'FeSi_27_atoms_initial_v2.db' ) and not exists_in_db(atoms, 'FeSi_27atoms_initial.db') and not exists_in_db(atoms, 'FeSi_16_atoms_initial.db') and not exists_in_db(atoms, 'FeSi_8atoms_initial.db'):
            db.write(atoms, struct_type = 'initial')
            added_to_db += 1

def exists_in_db(a, db_name):

    db = connect(db_name)
    atoms = []
    print('Checking database ' + db_name)
    for row in db.select(struct_type = 'initial'):

        atoms.append(row.toatoms())

    symcheck = SymmetryEquivalenceCheck()

    return symcheck.compare(a, atoms)

prepare_random_structures()
