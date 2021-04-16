# A genome sequence is the complete list of the nucleotides
# (A, C, G, and T for DNA genomes) that make up all the
# chromosomes of an individual or a species

# TODO: make nucleotides have the apropriate property...
# sepcify the nucleotides is a string containng only A,C,G,T

# nucleotides -- are all sequence of ACGT valid genome ?

class Genome:
    def __init__(self, nucleotides):
        self.nucleotides = nucleotides

    def get_nucleotides(self):
        return self.nucleotides


class GenomeBuilder():
    def __init__(self, nucleotides):
        self.genome = Genome(nucleotides)

    def build(self):
        return self.genome
