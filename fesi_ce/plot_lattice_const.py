from ase.build import bulk
from ase import Atom, Atoms
from ase.visualize import view
from ase.db import connect
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import ConvexHull

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
            #print(obj['formula'])

    return ids

def plot_lattice_const():

    db = connect('FeSi_8atoms.db')

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
        try:

            mag_val = bulk.get_magnetic_moment()

            print(mag_val)
        except:
            raise ValueError('Could not find magmom for id ' + str(id))
        num_atoms = num_si+num_fe
        mag_mom.append(mag_val/num_atoms)
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
    si_conc2, mag_mom = (list(t) for t in zip(*sorted(zip(si_conc2, mag_mom))))
    si_conc, lattice_const = (list(t) for t in zip(*sorted(zip(si_conc, lattice_const))))
    plt.figure(0)
    plt.plot(si_conc, lattice_const,'o')
    plt.figure(1)
    plt.plot(si_conc2, mag_mom, 'o')
    plt.plot([0,1],[mag_mom[0],mag_mom[-1]],'r')
    plt.xlabel('X(Si)')
    plt.ylabel('Magnetic moment per atom')
    plt.show()


def get_lattice_constants():

    db = connect('FeSi_8atoms_12finished.db')
    ids = get_ids()
    si_conc = []
    lattice_constants = []
    min_distances = []

    for num, id in enumerate(ids, 1):

        bulk = db.get_atoms(id=id)
        symbols = bulk.get_chemical_symbols()
        num_fe = 0
        num_si = 0

        for sym in symbols:
            if sym == 'Fe':
                num_fe += 1
            else:
                num_si += 1

        num_atoms = num_fe + num_si

        unit_cell_vectors = bulk.get_cell()
        vol = np.dot(unit_cell_vectors[0], np.cross(unit_cell_vectors[1],unit_cell_vectors[2]))
        lattice_param = (2*vol)**(1/3)
        if num_atoms == 16:
            lattice_constants.append(lattice_param)
            si_conc.append(num_si/num_atoms)

        positions = bulk.get_positions()

        distances = []

        for i in range(len(positions)-1):
            for j in range(i+1, len(positions)):
                distances.append(get_distance(positions[i],positions[j]))

        min_distances.append(min(distances))

    plt.figure(0)
    plt.plot(si_conc, lattice_constants, '*')
    plt.show()

def get_distance(pt1, pt2):

    return ((pt1[0]-pt2[0])**2 + (pt1[1]-pt2[1])**2 + (pt1[2]-pt2[2])**2)**(1/3)


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



def plot_formation_energy():

    db = connect('FeSi_8atoms_12finished.db')
    ids = get_ids()
    si_conc = []
    energies = []

    for obj in db.select(formula='Fe8',struct_type='final'):
        fe_ene = obj['energy']/8
    for obj in db.select(formula='Si8',struct_type='final'):
        si_ene = obj['energy']/8

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

        num_atoms = num_si + num_fe
        si_conc.append(num_si/num_atoms)

        for obj in db.select(id=id):

            energies.append(obj['energy']/(num_atoms) - (num_si/num_atoms)*si_ene - (num_fe/num_atoms)*fe_ene)


    si_conc, energies = (list(t) for t in zip(*sorted(zip(si_conc, energies))))

    ch_conc = []
    ch_en = []

    ch_conc.append(si_conc[0])
    ch_en.append(energies[0])

    for i in range(1,len(si_conc)):

        if ch_conc[-1] != si_conc[i]:

            ch_conc.append(si_conc[i])
            ch_en.append(energies[i])

        if ch_conc[-1] == si_conc[i] and energies[i] < ch_en[-1]:

            ch_en[-1] = energies[i]

    ch_conc, ch_en = (list(t) for t in zip(*sorted(zip(ch_conc, ch_en))))


    conv_hull_x = []
    conv_hull_y = []

    conv_hull_x.append(ch_conc[0])
    conv_hull_y.append(ch_en[0])
    data = np.column_stack((ch_conc,ch_en))
    hull = ConvexHull(data)


    #plt.plot(data[hull.vertices[0],0], data[hull.vertices[0],1], 'ro')

    x = list(data[hull.vertices,0])

    y = list(data[hull.vertices,1])

    x,y = (list(t) for t in zip(*sorted(zip(x, y))))

    plt.figure(0)
    #plt.plot(data[hull.vertices,0], data[hull.vertices,1], 'r--', lw=2)
    #print(x)
    plt.plot(x, y, 'r--', lw=2)

    plt.plot(si_conc, energies,'bo')
    plt.xlabel('X(Si)')
    plt.ylabel('Formation energy [eV/atom]')
    plt.title('Formation energy and convex hull for ' + str(len(ids)) + ' structures')
    #plt.fill_between([0.5375,1],[-0.75,-0.75], facecolor='grey', alpha=0.5)
    plt.fill_between([0.5375,1.01],[-0.75,-0.75], [0.075, 0.075], facecolor='grey', alpha=0.3)
    plt.show()



get_lattice_constants()
#plot_formation_energy()
#get_energies()
