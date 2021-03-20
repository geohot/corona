from lib import cc, translate


class gene:
    def __init__(self, id, sequence_range):
        self.id = id
        self.sequence_range = sequence_range

class rna_sequence:
    def __init__(self, id, sequence_range, is_protein):
        self.id = id
        self.sequence_range = sequence_range
        self.is_protein = is_protein
    
    def translate_rna(self):
        return translate(self.sequence_range, self.is_protein)
