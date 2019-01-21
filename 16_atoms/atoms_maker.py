from ase.build import bulk
from ase import Atom, Atoms
from ase.visualize import view
from ase.db import connect
from ase.utils.structure_comparator import SymmetryEquivalenceCheck
from itertools import product
from random import randint

def prepare_structures():

    atoms = bulk('Fe', 'bcc', a = 2.876)*(2,2,4)
    symbols = ['Fe', 'Si']
    #SET MAGMOM!
    db = connect('FeSi_8atoms_initial.db')
    for s in product(symbols, repeat = 8):
        for atom in atoms:
            atom.symbol = s[atom.index]
        if not exists_in_db(atoms):
            db.write(atoms, struct_type = 'initial')

def prepare_random_structure():

    num_fe = randint(0,16)
    num_si = 16 - num_fe

    return [num_fe, num_si]

def prepare_random_structures():


    db = connect('FeSi_16atoms_initial.db')
    added_to_db = 0

    while added_to_db < 10:
        atoms = bulk('Fe', 'bcc', a = 2.876)*(2,2,4)
        atom_config = prepare_random_structure()
        num_si = 0
        while num_si < atom_config[1]:
            for atom in atoms:
                if atom.symbol != 'Si' and randint(0,10)/10 <= atom_config[1]/16:
                    atom.symbol = 'Si'
                    num_si += 1


        if not exists_in_db(atoms) and num_si != 0 and num_si != 16:
            db.write(atoms, struct_type = 'initial')
            added_to_db += 1

    added_to_db = 10
    while added_to_db < 10:

        atoms = bulk('Fe', 'bcc', a = 2.876)*(2,2,4)
        max_si = randint(13,14)
        print(max_si)
        num_si = 0
        while num_si < max_si:

            for atom in atoms:
                if atom.symbol != 'Si' and randint(0,10)/10 <= 0.5 and num_si < max_si:
                    atom.symbol = 'Si'
                    num_si += 1

        if not exists_in_db(atoms):
            db.write(atoms, struct_type = 'initial')
            added_to_db += 1

def exists_in_db(a):

    db = connect('FeSi_16atoms_initial.db')
    atoms = []

    for row in db.select(struct_type = 'initial'):

        atoms.append(row.toatoms())

    symcheck = SymmetryEquivalenceCheck()

    return symcheck.compare(a, atoms)

prepare_random_structures()
