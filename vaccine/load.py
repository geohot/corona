#!/usr/bin/env python3

import pathlib
import os

# code from world health organization
# https://mednet-communities.net/inn/db/media/docs/11889.doc

"""
1-2        -- cap A modified 5’-cap1 structure (m7G+m3'-5'-ppp-5'-Am)
3-54       -- 5’-UTR  5´-untranslated region derived from human alpha-globin RNA with an optimized Kozak sequence
55-102     -- sig S glycoprotein signal peptide (extended leader sequence), which guides translocation of the nascent polypeptide chain into the endoplasmic reticulum.
103-3879   -- S protein_mut Codon-optimized sequence encoding full-length SARS-CoV-2 spike (S) glycoprotein containing mutations K986P and V987P to ensure the S glycoprotein remains in an antigenically optimal pre-fusion conformation; stop codons: 3874-3879 (underlined)
3880-4174  -- 3’-UTR  The 3´ untranslated region comprises two sequence elements derived from the amino-terminal enhancer of split (AES) mRNA and the mitochondrial encoded 12S ribosomal RNA to confer RNA stability and high total protein expression.
4175-4284  -- poly(A) A 110-nucleotide poly(A)-tail consisting of a stretch of 30 adenosine residues, followed by a 10-nucleotide linker sequence and another 70 adenosine residues.
"""

with open(os.path.join(pathlib.Path(__file__).parent.absolute(), "code")) as f:
  dat = f.read().strip().replace("\n", "").replace(" ", "")

if __name__ == "__main__":
  print(len(dat))
  print(dat)



