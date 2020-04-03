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
parser.add_argument('--temp', type=int, default=300)
parser.add_argument('--steps', type=int, default=100000, help="2500000000 should fold the protein")
parser.add_argument('--writes', type=int, default=1000, help="default is 1000")
parser.add_argument('--out', type=str, default="/tmp/output.pdb")
parser.add_argument('--pdb', type=str, default="proteins/villin/1vii.pdb")
parser.add_argument('--fasta', type=str, default=None)
args = parser.parse_args(sys.argv[1:])

try:
  platform = Platform.getPlatformByName("CUDA")
except Exception:
  platform = Platform.getPlatformByName("OpenCL")

if args.scratch:
  # unfolded protein
  if args.fasta is not None:
    fasta = args.fasta
  else:
    protein_fasta = "proteins/villin/1vii.fasta"
    fasta = open(protein_fasta).read().split("\n")[1]
  print("folding %s" % fasta)
  from lib import write_unfolded
  write_unfolded(fasta, "/tmp/unfolded.pdb")
  pdb = PDBFile("/tmp/unfolded.pdb")
else:
  # already folded protein
  pdb = PDBFile(args.pdb)

#forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
forcefield = ForceField('amber03.xml', 'amber03_obc.xml')

modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield)
print(modeller.topology)

system = forcefield.createSystem(modeller.topology,
  implicitSolvent=OBC2,   # matches https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2980750/#bib39
  nonbondedMethod=NoCutoff, nonbondedCutoff=1*nanometer,
  constraints=HBonds)
integrator = LangevinIntegrator(args.temp*kelvin, 1/picosecond, 2*femtoseconds)
simulation = Simulation(modeller.topology, system, integrator, platform)
simulation.context.setPositions(modeller.positions)
simulation.minimizeEnergy()

steps = args.steps
steps_write = max(1, steps//args.writes)
print("writing every %d steps" % steps_write)
simulation.reporters.append(PDBReporter(args.out, steps_write))
simulation.reporters.append(StateDataReporter(stdout, steps_write, step=True, potentialEnergy=True, temperature=True))
simulation.step(steps)

