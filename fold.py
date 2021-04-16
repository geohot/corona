#!/usr/bin/env python3
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os
import sys
import argparse

from parser import Parser
from unfold import Unfold
from simulate import Simulate


args = Parser(argparse).parse()

try:
    platform = Platform.getPlatformByName("CUDA")
except Exception:
    platform = Platform.getPlatformByName("OpenCL")

pdb = Unfold(args).unfold

#forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
forcefield = ForceField('amber03.xml', 'amber03_obc.xml')

modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield)
print(modeller.topology)

Simulate(forcefield, modeller, args).simulate()
