from openmm import *
from openmm.app import *
from openmm.unit import *

ff = ForceField('amber14-all.xml', 'implicit/obc2.xml')

pdb = PDBFile('./1gcn_capped.pdb')
mod = Modeller(pdb.topology, pdb.positions)
variants = mod.addHydrogens()

system = ff.createSystem(mod.topology)

#optional: fix non-cap atoms:
for residue in mod.topology.residues():
    if residue.name not in ['ACE','NME']:
        for atom in residue.atoms():
            system.setParticleMass(atom.index, 0)
            
integrator = LangevinIntegrator(250*kelvin,1/picosecond,1*femtosecond)
sim = Simulation(mod.topology, system, integrator, openmm.Platform.getPlatformByName('CPU'))
sim.context.setPositions(mod.positions)
sim.minimizeEnergy()

#optional:
#sim.step(10_000)

PDBFile.writeFile(
    mod.topology, 
    sim.context.getState(getPositions=True).getPositions(),
    open('./1gcn_minimized.pdb', 'w'),
    )
