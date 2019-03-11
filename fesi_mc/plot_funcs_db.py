import dataset
import time
import numpy as np

def plot_energies():
    import matplotlib.pyplot as plt

    db_name = "sqlite:///fesi_simannealing_888.db"

    db_set = dataset.connect(db_name)

    table = 'simulations'

    """Find the DFT energies of pure Fe and Si"""
    from ase.db import connect

    db = connect('FeSi_8atoms_12finished.db')

    for obj in db.select(formula='Fe8',struct_type='final'):
        fe_ene = obj['energy']/8
    for obj in db.select(formula='Si8',struct_type='final'):
        si_ene = obj['energy']/8

    plt.figure(0)

    temps = [4000, 3000, 2500, 2000, 1500, 1000, 750, 500, 250, 200, 150, 100]

    for t in temps:
        concs = []
        energies = []
        for obj in db_set[table]:
            if obj['temperature'] == t:
                concs.append(obj['Si_conc'])
                energies.append(obj['energy']/(1024)-obj['Fe_conc']*fe_ene - obj['Si_conc']*si_ene)

        plt.plot(concs,energies, '-x')

    plt.legend(temps)
    plt.xlabel('Si conc')
    plt.ylabel('Formation energy')
    #plt.colorbar()
    plt.grid(True)

    plt.show()

def plot_heat_capacity():
    import matplotlib.pyplot as plt

    db_name = "sqlite:///fesi_simannealing_888.db"

    db_set = dataset.connect(db_name)

    table = 'simulations'

    plt.figure(0)

    temps = [4000, 3000, 2500, 2000, 1500, 1000, 750, 500, 250, 200, 150, 100]
    plot_temps = temps
    plot_concs = [0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65]
    plot_hc = []

    for t in temps:
        concs = []
        heat_cap = []
        for obj in db_set[table]:
            if obj['temperature'] == t:

                if obj['Fe_conc'] > 0.64:
                    plot_hc.append(obj['heat_capacity'])

                concs.append(obj['Si_conc'])
                heat_cap.append(obj['heat_capacity'])

        plt.plot(concs,heat_cap, '-x')


    plt.legend(temps)
    plt.xlabel('Si conc')
    plt.ylabel('Heat capacity')

    plot_concs = [1-i for i in plot_concs]

    #plt.figure(1)
    # z = np.array(plot_hc).reshape((len(plot_concs), len(plot_temps)))
    # X, Y = np.meshgrid(plot_temps, plot_concs)
    # cp = plt.contour(X,Y,z, 100, cmap='RdGy')
    # plt.colorbar(cp)

    plt.show()

def plot_heat_capacity_vs_t():
    import matplotlib.pyplot as plt

    db_name = "sqlite:///fesi_simannealing_888.db"

    db_set = dataset.connect(db_name)

    table = 'simulations'

    plt.figure(0)

    temps = [4000, 3000, 2500, 2000, 1500, 1000, 750, 500, 250, 200, 150, 100]
    concs = [0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65]
    plot_temps = temps
    plot_concs = []


    for c in concs:
        concs = []

        for obj in db_set[table]:

            if obj['Si_conc'] not in plot_concs:

                temp = []
                heat_cap = []
                plot_hc = []

                plot_concs.append(obj['Si_conc'])
                loc_conc = obj['Si_conc']

                for obj2 in db_set[table]:
                    if obj2['Si_conc'] == loc_conc:
                        plot_hc.append(obj2['heat_capacity'])
                        temp.append(obj2['temperature'])



                plt.plot(temp,plot_hc, '-x')


    plt.legend([round(i,2) for i in plot_concs])
    plt.xlabel('Temperature')
    plt.ylabel('Heat capacity')

    plt.show()

if __name__ == '__main__':

    #plot_energies()
    #plot_heat_capacity()
    plot_heat_capacity_vs_t()
