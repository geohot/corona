class Genome:
    def __init__(self, nucleotideSequence):
        self.nucleotideSequence = nucleotideSequence

    def get_nucleotideSequence(self):
        return self.nucleotideSequence


class GenomeBuilder():
    def __init__(self, nucleotideSequence):
        self.genome = Genome(nucleotideSequence)

    def build(self):
        if(self.validate()):
            return self.genome
        else:
            return None

    def validate(self):
        my_set = set('ACGTU')
        return set(self.genome.get_nucleotideSequence).issubset(my_set)
