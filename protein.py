from lib import cc, translate
from fold import *

class protein:
    def __init__(self, id, aminoacid_sequence):
        self.id = id
        self.aminoacid_sequence = aminoacid_sequence
        
    
    def structure(self):
        return fold (self.sequence_range, self.is_protein)
