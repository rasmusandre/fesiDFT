from cemc.mcmc import Montecarlo
from cemc import get_atoms_with_ce_calc
from cemc.mcmc import LowestEnergyStructure
from cemc.mcmc import Snapshot
from cemc.mcmc import SiteOrderParameter
import json
import dataset

def get_atoms(size):

    from ase.clease import Concentration, CEBulk

    conc = Concentration(basis_elements=[['Fe','Si']])

    args = dict(crystalstructure='bcc', a=2.876, db_name='FeSi_mc.db',
    max_cluster_size=4, size=[2,2,2],
    max_cluster_dia=[0.0, 10, 10, 7.5, 5.5], concentration=conc, cubic=True)

    ceBulk = CEBulk(**args)

    with open('my_file_eci_4nn.json') as infile:

        eci = json.load(infile)

    atoms = get_atoms_with_ce_calc(ceBulk, args, eci=eci, size=size, db_name='FeSi_bigcell.db')
    #if 4,4,4 success go to 8,8,8

    return atoms

def simulated_ann(iron_conc, size, temps, db_name):

    conc = {'Fe': iron_conc, 'Si':1 - iron_conc}

    atoms = get_atoms(size)

    atoms.get_calculator().set_composition(conc)

    low_en_struct = LowestEnergyStructure(atoms.get_calculator(), None)

    db = dataset.connect(db_name)
    tbl = db['simulations']

    for t in temps:

        print(t)

        mc = Montecarlo(atoms, t)
        low_en_struct.mc_obj = mc
        mc.attach(low_en_struct, interval=1)
        mc.runMC(mode='fixed', equil=True, steps=1000*len(atoms))

        thermo = mc.get_thermodynamic()
        thermo['size'] = size[0]

        tbl.insert(thermo)

    from ase.io import write
    write('groundstate.xyz', low_en_struct.atoms)

def dissolution(temps, size):

    conc = {'Fe': iron_conc, 'Si':1 - iron_conc}

    atoms = get_atoms(size)

    atoms.get_calculator().set_composition(conc)

    #Sett in fe3si i atoms, 10x10x10 fe3si, insert -> for hvert atom i fe3si set inn,
    #calc.set_symbols

    snap = Snapshot('my_snap.traj', atoms)

    db = dataset.connect(db_name)
    tbl = db['simulations']
    mc = Montecarlo(atoms, t)
    mc.attach(snap, interval = 10**4)
    site_order_param = SiteOrderParameter(atoms)
    mc.attach(site_order_param)


    for t in temps:

        mc.T = t

        mc.runMC(mode='fixed', equil=True, steps=1000*len(atoms))

        thermo = mc.get_thermodynamic()
        thermo['size'] = size[0]

        tbl.insert(thermo)
        mc.reset()

    from ase.io import write
    write('groundstate.xyz', low_en_struct.atoms)


if __name__ == '__main__':

    iron_conc = 0.75
    size = [30,30,30]
    temps = [100,250,300,400,500,750,1000,1250,1500,2000,3000,4000,5000]#[5000, 4000, 3000, 2500, 2000, 1500, 1000, 750, 500, 250]
    db_name = "sqlite:///fesi_simannealing.db"

    #simulated_ann(iron_conc, size, temps)
    #[1,0.95,0.9,0.85....]
