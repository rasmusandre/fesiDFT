from ase.build import bulk
from ase import Atom, Atoms
from ase.visualize import view
from ase.db import connect
import matplotlib.pyplot as plt

def get_energies():

    db = connect('FeSi_8atoms_12finished.db')
    energies = []
    cu_conc = []
    fe_bulk = 0
    cu_bulk = 0

    for obj in db.select(is_initial = 0):

        if obj['formula'] == 'Fe8':
            fe_bulk = obj['energy']
        if obj['formula'] == 'Si8':
            cu_bulk = obj['energy']


    for obj in db.select(is_initial = 0):

        formula = obj['formula']
        if formula == 'Fe8':
            numcu = 0
        if formula == 'Si8':
            numcu = 1
        elif len(formula) > 3:
            formula.split('Fe')
            print(formula)

        cu_conc.append(numcu/obj['natoms'])

        energy = obj['energy'] - cu_bulk*numcu/obj['natoms'] - fe_bulk*(1-numcu/obj['natoms'])

        energies.append(energy)

    plt.figure(0)
    plt.plot(energies, cu_conc)
    plt.grid(True)
    plt.show()

def get_ids():

    db = connect('FeSi_8atoms_12finished.db')

    ids = []
    final_ids = []

    for obj in db.select():

        if 'energy' in obj:
            ids.append(obj['id'])
            print(obj['formula'])

    return ids

def plot_lattice_const():

    db = connect('FeSi_8atoms_12finished.db')

    lattice_const = []
    mag_mom = []
    si_conc = []
    ids = get_ids()
    plt.figure(2)

    for id in ids:

        bulk = db.get_atoms(id=id)

        symbols = bulk.get_chemical_symbols()
        num_fe = 0
        num_si = 0

        for sym in symbols:

            if sym == 'Fe':
                num_fe += 1
            else:
                num_si += 1

        si_conc.append(num_si/(num_si+num_fe))

        #print(bulk.get_cell())
        #print('----------')
        #for obj in db.select(id=id):
            #mag_val = obj['magmom']/(num_si+num_fe)
        #mag_val = bulk.get_magnetic_moment()
        #mag_mom.append(mag_val)
        num_atoms = num_si+num_fe
        if num_atoms==27:
            val = (1/3)*(2*bulk.get_volume())**(1/3)
            plt.plot(num_si/num_atoms, val, 'ro')
        if num_atoms==16:
            val = (1/2)*(2*bulk.get_volume())**(1/3)
            plt.plot(num_si/num_atoms, val, 'bo')
        if num_atoms==8:
            val = (2/5)*(2*bulk.get_volume())**(1/3)
            plt.plot(num_si/num_atoms, val, 'go')
        
        lattice_const.append(val)

    si_conc2 = si_conc
    #si_conc2, mag_mom = (list(t) for t in zip(*sorted(zip(si_conc2, mag_mom))))
    si_conc, lattice_const = (list(t) for t in zip(*sorted(zip(si_conc, lattice_const))))
    plt.figure(0)
    plt.plot(si_conc, lattice_const,'o')
    #plt.figure(1)
    #plt.plot(si_conc2, mag_mom, 'o')
    #plt.plot([0,1],[mag_mom[0],mag_mom[-1]],'r')
    plt.xlabel('X(Si)')
    plt.ylabel('Magnetic moment per atom')
    plt.show()




def find_atoms():

    db = connect('fecu_8atoms_first.db')

    atom1 = db.get_atoms(id=24)
    view(atom1)




def brute_force():
    efe = -72.032
    ecu = -28.897
    magmoms = [17.267, 16.077, 0, 10.628, 13.002, 14.416, 2.961, 10.227, 5.491, 10.416, 8.913, 14.793, 7.816, 10.349, 5.346, 12.546]
    energies = [-72.032, -66.090, -28.897, -48.475, -54.355, -60.329, -33.341, -49.211, -38.469, -48.967, -43.276, -60.259, -43.730, -49.158, -38.500, -54.701]
    conc = [0, 1/8, 1, 1/2, 3/8, 2/8, 7/8, 1/2, 6/8, 1/2, 5/8, 2/8, 5/8, 1/2, 6/8, 3/8]
    conc2 = conc
    conc2, magmoms = (list(t) for t in zip(*sorted(zip(conc2, magmoms))))

    conc, energies = (list(t) for t in zip(*sorted(zip(conc, energies))))

    energies = [energies[i] - conc[i]*ecu - (1-conc[i])*efe for i in range(0,len(energies))]

    plt.figure(0)
    plt.plot(conc, energies, 'bo')
    plt.xlabel('X(Cu)')
    plt.ylabel('Energy of formation [eV]')
    plt.grid(True)

    plt.figure(1)
    plt.plot(conc, magmoms, 'r*')
    plt.plot([0,1], [17.267,0], 'b')
    plt.xlabel('X(Cu)')
    plt.ylabel('Total magnetic moment [$\mu_B$]')
    plt.legend(['Calculated magnetic moment', 'Zero induced magnetic moment'],fancybox=True, framealpha=1,shadow=True,prop={'size': 10})
    plt.grid(True)
    plt.show()

plot_lattice_const()
k = 2
mat = bulk('Fe','bcc',a=2.86)*(2,2,4)
#view(mat)
print(mat.get_volume())
lc = (2*mat.get_volume())**(1/3)
print(lc)











#get_energies()
