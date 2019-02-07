from ase.db import connect
from ase.optimize import BFGS
from ase.constraints import UnitCellFilter
from ase.io.trajectory import Trajectory
from ase.build import niggli_reduce
from gpaw import GPAW, PW, FermiDirac, restart
import sys
import os
import time

#This script is used to relax 8 atom structures of FeSi

class Gpw_save(object):

    def __init__(self, calc, txtname):

        self.calc = calc
        self.txtname = txtname

    def __call__(self):

        self.calc.write(self.txtname)

def save_atoms(my_atoms, E_c, Nbands, Kpts, Fermi_dirac, Magnetic_moment, Si_concentration, Initial_id, Is_initial, database):

    db = connect(database)
    uid = db.write(my_atoms, energy_cutoff = E_c,
             nbands = Nbands, k_points = Kpts,
             smearing_factor = Fermi_dirac,
             initial_magnetic_moment = Magnetic_moment,
             si_concentration = Si_concentration,
             initial_id = Initial_id,
             is_initial = Is_initial,
             struct_type = "final")

    db.update(Initial_id, converged=1, final_struct_id=uid)

if __name__ == "__main__":


    #Main arguments:
    e_cut = 600
    nbands = -100
    xc = 'PBE'
    density_number = int(sys.argv[3])
    #density_number = 5.4, 1.37
    k_pts = {'density': density_number, 'even': True}
    smear = 0.1
    system_id = int((sys.argv[1]))
    system_name = 'fesi_id' + str(system_id)
    final_db = str(sys.argv[2])

    print (time.strftime("%d/%m/%Y"))
    print(time.strftime("%H:%M:%S"))
    print('The id for this script is: ' + str(system_id))


    calc = GPAW(mode=PW(e_cut), nbands = nbands,
                xc='PBE', spinpol=True, kpts=k_pts,
                occupations=FermiDirac(smear), txt=system_name + '.out')

    if os.path.isfile('./' + system_name + '.gpw'):
        #Recover past done calculations if available
        #bulk_mat, calc = restart(system_name + '.gpw', txt = None)
        #bulk_mat.set_calculator(calc)
        traj = Trajectory(system_name + '.traj')
        bulk_mat = traj[-1]
        bulk_mat.set_calculator(calc)


    else:
        #Initialize new calculations
        db = connect(final_db)
        bulk_mat = db.get_atoms(id = system_id)

        for atom in bulk_mat:
            #initialize magnetic moments
            if atom.symbol == 'Fe':
                atom.magmom = 2.2
            else:
                atom.magmom = 0.0


        bulk_mat.set_calculator(calc)


    number_of_atoms = bulk_mat.get_number_of_atoms()
    number_of_si = 0
    initial_magmom = 0

    for atom in bulk_mat:
        if atom.symbol == 'Fe':
            initial_magmom += 2.2
        else:
            number_of_si += 1

    si_concentration = number_of_si/number_of_atoms

    saver = Gpw_save(calc, system_name + '.gpw')
    #Check if there is an existing .traj, if so append
    if os.path.isfile('./' + system_name + '.traj'):
        traj = Trajectory(system_name + '.traj', 'a', bulk_mat)
    else:
        traj = Trajectory(system_name + '.traj', 'w', bulk_mat)
    ucf = UnitCellFilter(bulk_mat)
    relaxer = BFGS(ucf, logfile = system_name + '.txt')
    relaxer.attach(traj)
    relaxer.attach(saver)
    relaxer.run(fmax = 0.025)
    bulk_mat.get_potential_energy()

    #Save the final state of the calculations
    if density_number == 5.4:

    	calc.write(system_name + '_final.gpw')
        save_atoms(bulk_mat, e_cut, nbands, density_number, smear, initial_magmom, si_concentration, system_id, 0, final_db)
