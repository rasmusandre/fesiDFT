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

def prepare_random_structure(num_atoms):

    num_fe = 0
    num_si = 10

    while num_si/(num_si + num_fe) > 0.75 and num_si != num_atoms and num_fe != num_atoms:
        num_fe = randint(1,num_atoms - 1)
        num_si = num_atoms - num_fe

    return [num_fe, num_si]

def prepare_random_structures():


    db = connect('FeSi_8atoms_12finished.db')
    added_to_db = 0
    trial = 0
    i = 1
    j = 3
    k = 3
    atom_conf = (i,j,k)
    num_atoms = i*j*k
    while added_to_db < 10 and trial < 100:
        trial += 1
        atoms = bulk('Fe', 'bcc', a = 2.876)*atom_conf
        atom_config = prepare_random_structure(num_atoms)

        num_si = 0
        while num_si < atom_config[1]:
            for atom in atoms:
                if atom.symbol != 'Si' and randint(0,10)/10 <= atom_config[1]/num_atoms:
                    atom.symbol = 'Si'
                    num_si += 1
                if num_si == atom_config[1]:
                    break
            print(num_si)

        if not exists_in_db(atoms, 'FeSi_8atoms_12finished.db' ) and not exists_in_db(atoms, 'FeSi_27atoms_initial.db') and not exists_in_db(atoms, 'FeSi_16_atoms_initial.db') and not exists_in_db(atoms, 'FeSi_8atoms_initial.db'):
            db.write(atoms, struct_type = 'initial')
            added_to_db += 1

def exists_in_db(a, db_name):

    db = connect(db_name)
    atoms = []
    #print('Checking database ' + db_name)
    for row in db.select(struct_type = 'initial'):

        atoms.append(row.toatoms())

    symcheck = SymmetryEquivalenceCheck()

    if symcheck.compare(a, atoms):
        print(db_name)

    return symcheck.compare(a, atoms)

prepare_random_structures()
