import sys


class Simulate:
    def __init__(self, forcefield, modeller, args):
        self.forcefield = forcefield
        self.modeller = modeller
        self.args = args

    def simulate(self):
        system = self.forcefield.createSystem(self.modeller.topology,
                                              # matches https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2980750/#bib39
                                              implicitSolvent=OBC2,
                                              nonbondedMethod=NoCutoff, nonbondedCutoff=1*nanometer,
                                              constraints=HBonds)
        integrator = LangevinIntegrator(
            self.args.temp*kelvin, 1/picosecond, 2*femtoseconds)
        simulation = Simulation(
            self.modeller.topology, system, integrator, platform)
        simulation.context.setPositions(self.modeller.positions)
        simulation.minimizeEnergy()

        steps = self.args.steps
        steps_write = max(1, steps//self.args.writes)
        print("writing every %d steps" % steps_write)
        simulation.reporters.append(PDBReporter(args.out, steps_write))
        simulation.reporters.append(StateDataReporter(
            stdout, steps_write, step=True, potentialEnergy=True, temperature=True))
        simulation.step(steps)
