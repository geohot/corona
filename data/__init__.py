import os
import pathlib

with open(os.path.join(pathlib.Path(__file__).parent.absolute(), "nl63.fasta")) as f:
  nl63 = ''.join(f.read().split("\n")[1:]).upper()
  
