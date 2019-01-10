from ase.io import read
from gpaw import GPAW, PW
from ase.build import niggli_reduce

atoms = read('fesi_system.xyz')
#niggli_reduce(atoms)
calc = GPAW(mode=PW(300),
            kpts=[2,2,2],
            symmetry={'do_not_symmetrize_the_density': True})
atoms.calc = calc
atoms.get_potential_energy()
