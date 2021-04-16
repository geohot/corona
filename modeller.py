import sys


class ModellerF:
    def __init__(self, forcefield):
        self.forcefield = forcefield

    def model(self):
       #forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
        modeller = Modeller(pdb.topology, pdb.positions)
        modeller.addHydrogens(self.forcefield)
        print(modeller.topology)
