#!/usr/bin/env python3
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os

try:
  platform = Platform.getPlatformByName("CUDA")
except Exception:
  platform = Platform.getPlatformByName("OpenCL")

# already folded protein
#protein_pdb = "proteins/villin/1vii.pdb"
#pdb = PDBFile(protein_pdb)

# unfolded protein
protein_fasta = "proteins/villin/1vii.fasta"
fasta = open(protein_fasta).read().split("\n")[1]
print("folding %s" % fasta)
from lib import write_unfolded
write_unfolded(fasta, "/tmp/unfolded.pdb")
pdb = PDBFile("/tmp/unfolded.pdb")

forcefield = ForceField('amber99sb.xml', 'tip3p.xml')

modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield)
print(modeller.topology)

system = forcefield.createSystem(modeller.topology, nonbondedMethod=NoCutoff, nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(modeller.topology, system, integrator, platform)
simulation.context.setPositions(modeller.positions)
simulation.minimizeEnergy()

#steps = 100000
# 2500000000

# still off by factor of 100!
steps = 25000000

steps_write = steps//1000
print("writing every %d steps" % steps_write)
simulation.reporters.append(PDBReporter('/tmp/output.pdb', steps_write))
simulation.reporters.append(StateDataReporter(stdout, steps_write, step=True, potentialEnergy=True, temperature=True))
simulation.step(steps)

