import dataset
import time

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
                energies.append(obj['energy']/(1024)-obj['Fe_conc']*fe_ene-obj['Si_conc']*si_ene)

        plt.plot(concs,energies, '-x')

    plt.legend(temps)
    plt.xlabel('Si conc')
    plt.ylabel('Formation energy')

    plt.show()

def run_job():
    from fesi_mc import simulated_ann

    concs = [0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5]

    size = [8,8,8]
    temps = [4000, 3000, 2500, 2000, 1500, 1000, 750, 500, 250, 200, 150, 100]
    db_name = "sqlite:///fesi_simannealing_888.db"

    for conc in concs:
        print('---------------------------------------------------')
        print('The concentration is now {}'.format(conc))
        print('---------------------------------------------------')
        start = time.time()
        simulated_ann(conc, size, temps, db_name)
        end = time.time()
        time_delta = end - start
        print('The time used for this concentration was {}'.format(time_delta))


if __name__ == '__main__':

    plot_energies()
    #run_job()
