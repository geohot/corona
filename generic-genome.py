class Genome:
    def __init__(self, chromosomes, genes, dnscode):
        self.chromosomes = chromosomes

    def get_nucleotideSequence(self):
        return self.chromosomes


class GenomeBuilder():
    def add_chromosomes(self, chromosomes):
        self.chromosomes = chromosomes

    def add_genes(self, genes):
        self.genes = genes

    def add_dnscode(self, dnscode):
        self.dnscode = dnscode

    def build(self):
        return Genome(self.chromosomes, self.genes, self.dnscode)
