from ase.db import connect
import matplotlib.pyplot as plt


def print_energies(database_name):

    db = connect(database_name)
    bulk_fe = 0
    bulk_si = 0
    form_energies = []
    concentrations = []

    for obj in db.select(is_initial = 0):

        if not 'Si' in obj['formula']:
            bulk_fe = obj['energy']
        if not 'Fe' in obj['formula']:
            bulk_si = obj['energy']


    for obj in db.select(is_initial = 0):

        if not 'Si' in obj['formula']:
            conc = 0
        if not 'Fe' in obj['formula']:
            conc = 1

        if 'Si' in obj['formula'] and 'Fe' in obj['formula']:

            formula = obj['formula']
            fe, si = formula.split('Si')
            fe = fe.split('Fe')[1]

            if not fe:
                fe = 1
            if not si:
                si = 1

            conc = int(si)/(int(si) + int(fe))

        form_energies.append(obj['energy'] - conc*bulk_si - (1-conc)*bulk_fe)
        concentrations.append(conc)

    
    plt.figure(0)
    plt.plot(concentrations, form_energies, '*')
    plt.ylabel('Formation energy [eV]')
    plt.xlabel('X(Si)')
    plt.grid(True)
    plt.show()

print_energies('FeSi_8atoms_12finished.db')
