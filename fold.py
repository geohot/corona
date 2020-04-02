#!/usr/bin/env python3
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os

#platform = Platform.getPlatformByName("CUDA")
platform = Platform.getPlatformByName("OpenCL")

protein_pdb = "proteins/villin/1vii.pdb"
pdb = PDBFile(protein_pdb)
forcefield = ForceField('amber99sb.xml', 'tip3p.xml')

modeller = Modeller(pdb.topology, pdb.positions)
print(modeller.topology)

system = forcefield.createSystem(modeller.topology, nonbondedMethod=NoCutoff, nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(modeller.topology, system, integrator, platform)
simulation.context.setPositions(modeller.positions)
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter('/tmp/output.pdb', 25))
simulation.reporters.append(StateDataReporter(stdout, 100, step=True, potentialEnergy=True, temperature=True))
simulation.step(5000)

