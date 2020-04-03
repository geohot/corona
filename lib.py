import random

# Asn or Asp / B  AAU, AAC; GAU, GAC
# Gln or Glu / Z  CAA, CAG; GAA, GAG
# START AUG
tt = """Ala / A GCU, GCC, GCA, GCG
Ile / I AUU, AUC, AUA
Arg / R CGU, CGC, CGA, CGG; AGA, AGG
Leu / L CUU, CUC, CUA, CUG; UUA, UUG
Asn / N AAU, AAC
Lys / K AAA, AAG
Asp / D GAU, GAC
Met / M AUG
Phe / F UUU, UUC
Cys / C UGU, UGC
Pro / P CCU, CCC, CCA, CCG
Gln / Q CAA, CAG
Ser / S UCU, UCC, UCA, UCG; AGU, AGC
Glu / E GAA, GAG
Thr / T ACU, ACC, ACA, ACG
Trp / W UGG
Gly / G GGU, GGC, GGA, GGG
Tyr / Y UAU, UAC
His / H CAU, CAC
Val / V GUU, GUC, GUA, GUG
STOP    UAA, UGA, UAG
""".strip()
dec = {}
for t in tt.split("\n"):
  k = t[:len("Val / V")].strip()
  v = t[len("Val / V "):]
  if '/' in k:
    k = k.split("/")[-1].strip()
  k = k.replace("STOP", "*")
  v = v.replace(",", "").replace(";", "").lower().replace("u", "t").split(" ")
  for vv in v:
    if vv in dec:
      print("dup", vv)
    dec[vv.strip()] = k

def translate(x, protein=False):
  x = x.lower()
  aa = []
  for i in range(0, len(x)-2, 3):
    aa.append(dec[x[i:i+3]])
  aa = ''.join(aa)
  if protein:
    if aa[0] != "M" or aa[-1] != "*":
      print("BAD PROTEIN")
      print(aa)
      return None
    aa = aa[:-1]
  return aa

ltl = 'Asp D Glu E Arg R Lys K His H Asn N Gln Q Ser S Thr T Tyr Y Ala A Gly G Val V Leu L Ile I Pro P Phe F Met M Trp W Cys C'
ltl = ltl.split(" ")
ltl = dict(zip(ltl[1::2], ltl[0::2]))

def get_atoms():
  from data import get_amber99sb
  amber99sb = get_amber99sb()
  residues = amber99sb.getElementsByTagName("Residue")
  atoms = {}
  for r in residues:
    name = r.attributes['name'].value
    atoms[name] = [x.attributes['name'].value for x in r.getElementsByTagName("Atom")]
  return atoms

def write_unfolded(fasta, fn):
  atoms = get_atoms()
  atom_num = 1
  res_num = 1
  ss = []
  random.seed(1337)
  for i, aa in enumerate(fasta):
    tl = ltl[aa].upper()
    for a in atoms[tl] + ([] if i != len(fasta)-1 else ["OXT"]):
      if len(a) < 4:
        pa = " " + a
      else:
        pa = a
      gr = lambda: 1.0*(random.random()-0.5)
      x,y,z = gr(), gr(), gr()
      x += res_num*5
      s = "ATOM %6d %-4s %3s A %3d    %8.3f%8.3f%8.3f  1.00  1.00           %s" % \
        (atom_num, pa, tl, res_num, x, y, z, a[0:1])
      ss.append(s)
      atom_num += 1
    res_num += 1

  with open(fn, "w") as f:
    f.write('\n'.join(ss))
  
def invert(dd):
  dd = dd.upper()
  def _invert(x):
    if x == 'A':
      return 'T'
    elif x == 'T':
      return 'A'
    elif x == 'C':
      return 'G'
    elif x == 'G':
      return 'C'
  return (''.join([_invert(x) for x in dd]))[::-1]

import pathlib
import os
import json
with open(os.path.join(pathlib.Path(__file__).parent.absolute(), "data", "allseq.json")) as f:
  allseq = json.load(f)
cc = allseq['MN908947']

