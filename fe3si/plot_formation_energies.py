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

        if obj['formula'] == 'Fe12Si4':
            #Grab fe3si data and renormalize bulk energies
            fe3si_en = (obj['energy'] - conc*bulk_si*(16/8) - (1-conc)*bulk_fe*(16/8))/16
            fe3si_co = conc
        else:
            form_energies.append((obj['energy'] - conc*bulk_si - (1-conc)*bulk_fe)/8)
            concentrations.append(conc)




    plt.figure(0)
    plt.plot(concentrations, form_energies, 'b*')
    plt.plot(fe3si_co, fe3si_en, 'r*')
    plt.legend(['8 atom structures', 'Fe3Si structure (Fe12Si4)'])

    #minima_c = [0, fe3si_co, 0.5, 0.625, 0.75, 0.875, 1]
    #minima_e = [0, fe3si_en, -0.7155, -0.2909, -0.0723, 0.0402, 0]

    minima_c = [0, fe3si_co, 0.5, 1]
    minima_e = [0, fe3si_en, -0.7155, 0]

    plt.plot(minima_c, minima_e, 'y')

    print(form_energies)
    print(concentrations)



    plt.ylabel('Formation energy [eV/atom]')
    plt.xlabel('X(Si)')
    plt.grid(True)
    plt.show()

print_energies('FeSi_8atoms_withfe3si.db')
