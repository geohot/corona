#!/usr/bin/env python3
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os
import sys
import argparse

parser = argparse.ArgumentParser(description='Fold some proteins.')
parser.add_argument('--scratch', action='store_true')
parser.add_argument('--steps', type=int, default=100000, help="2500000000 should fold the protein")
parser.add_argument('--writes', type=int, default=1000, help="default is 1000")
parser.add_argument('--out', type=str, default="/tmp/output.pdb")
args = parser.parse_args(sys.argv[1:])

try:
  platform = Platform.getPlatformByName("CUDA")
except Exception:
  platform = Platform.getPlatformByName("OpenCL")

if args.scratch:
  # unfolded protein
  protein_fasta = "proteins/villin/1vii.fasta"
  fasta = open(protein_fasta).read().split("\n")[1]
  print("folding %s" % fasta)
  from lib import write_unfolded
  write_unfolded(fasta, "/tmp/unfolded.pdb")
  pdb = PDBFile("/tmp/unfolded.pdb")
else:
  # already folded protein
  protein_pdb = "proteins/villin/1vii.pdb"
  pdb = PDBFile(protein_pdb)

forcefield = ForceField('amber99sb.xml', 'tip3p.xml')

modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield)
print(modeller.topology)

system = forcefield.createSystem(modeller.topology, nonbondedMethod=NoCutoff, nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(modeller.topology, system, integrator, platform)
simulation.context.setPositions(modeller.positions)
simulation.minimizeEnergy()

steps = args.steps
steps_write = steps//args.writes
print("writing every %d steps" % steps_write)
simulation.reporters.append(PDBReporter(args.out, steps_write))
simulation.reporters.append(StateDataReporter(stdout, steps_write, step=True, potentialEnergy=True, temperature=True))
simulation.step(steps)

