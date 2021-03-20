#!/usr/bin/env python3
from lib import cc as virus
from vaccine.load import dat as vaccine
from corona import *

# virus = virus.replace("T", "U")
vaccine = vaccine.replace("Ψ", "U")

"""
for i in range(len(virus)-len(vaccine)):
  mm = virus[i:i+len(vaccine)]
  mr = sum([c1 == c2 for c1,c2 in zip(mm, vaccine)]) / len(vaccine)
  if mr > 0.5:
    print(i, mr)
exit(0)
"""

# vaccine starts at 21508 with a 67% match
# spike protein starts at 21562
# vvirus = virus[21508:21508+len(vaccine)]

# changed the RNA sequence of the corona virus where we have the overlap with the vaccine 
# so that it is an rna_sequence object from corona.py so that this can have the full access to the method of the rna_sequence object 

vaccine_coverage = rna_sequence('vaccine_coverage', cc[21508:21508+len(vaccine)], True)
vvirus = vaccine_coverage.rna_sequence_range.replace("T", "U")


print(vvirus)
print(vaccine)

