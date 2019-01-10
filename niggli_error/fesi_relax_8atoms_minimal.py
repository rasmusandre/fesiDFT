from ase.optimize import BFGS
from ase.constraints import UnitCellFilter
from ase.io import read
from gpaw import GPAW, PW, FermiDirac
import sys
import os

my_atoms = read('fesi_system.xyz')

calc = GPAW(mode=PW(600), nbands = -100,
            xc='PBE', spinpol=True, kpts={'density': 1.37, 'even': True},
            occupations=FermiDirac(0.1))

my_atoms.set_calculator(calc)
ucf = UnitCellFilter(bulk_mat)
relaxer = BFGS(ucf, logfile = 'relaxer_log.txt')
relaxer.run(fmax = 0.025)
bulk_mat.get_potential_energy()
