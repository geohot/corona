import os
import pathlib

BASEDIR = pathlib.Path(__file__).parent.absolute()

with open(os.path.join(BASEDIR, "nl63.fasta")) as f:
  nl63 = ''.join(f.read().split("\n")[1:]).upper()

def get_amber99sb():
  from xml.dom import minidom
  amber99sb = minidom.parse(os.path.join(BASEDIR, "amber99sb.xml"))
  return amber99sb

