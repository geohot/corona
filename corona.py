from lib import cc, translate
# entire diff: https://www.ncbi.nlm.nih.gov/projects/msaviewer/?rid=7FYNU14F01R&coloring=
# protein alignments: http://virological.org/t/alignment-of-58-sarbecovirus-genomes-for-conservation-analysis-of-sars-cov-2/430

# https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2
# https://www.ncbi.nlm.nih.gov/nuccore/MN908947
# https://zhanglab.ccmb.med.umich.edu/C-I-TASSER/2019-nCov/

# whole thing has a "lipid bilayer envelope", with S E M sticking out
# the orf proteins are "non structural" and form a "replicase-transcriptase complex"

# copy machine -- https://www.uniprot.org/uniprot/Q0ZJN1
# zhanglab breaks this down into many more proteins

corona = {}

# begin: 266 base pairs "untranslated"
# https://en.wikipedia.org/wiki/Five_prime_untranslated_region

# reference genome parsed: https://www.ncbi.nlm.nih.gov/genbank/sars-cov-2-seqs/

# 1ab = replicase polyprotein, https://www.ncbi.nlm.nih.gov/protein/YP_009724389.1?report=graph
# (same one for SARS v1 https://www.ncbi.nlm.nih.gov/protein/NP_828849.2?report=graph)
# 1-10 are in orf1a
#    1 = Host translation inhibitor nsp1, leader protein?
#    2 = ???
#    3 = Papain-like proteinase (includes PLpro)
#        see diff https://www.ncbi.nlm.nih.gov/projects/msaviewer/?rid=7FXGTZFN016&coloring=cons
#        Ac + ADRP + SUD + PLpro(1541-1855) + TM + Y domain (from SARS)
#    4 = nsp4B_TM; contains transmembrane domain 2 (TM2); produced by both pp1a and pp1ab
#    5 = 3CLpro, Proteinase 3CL-PRO
#    6 = putative transmembrane domain
#    7 = hexadecamer with nsp8
#    8 = primase, makes primers for RdRp
#    9 = ssRNA-binding protein; produced by both pp1a and pp1ab
#   10 = nsp10_CysHis; formerly known as growth-factor-like protein (GFL)
#   11 = only 13 aa
# 12 is mostly and 13-16 in orf1b
#   12 = https://en.wikipedia.org/wiki/RNA-dependent_RNA_polymerase
#   13 = Helicase (Hel).
#   14 = Guanine-N7 methyltransferase (ExoN) or maybe 3'-to-5' exonuclease
#   15 = Uridylate-specific endoribonuclease (NendoU), endoRNAse
#   16 = 2'-O-methyltransferase (2'-O-MT)
#        https://en.wikipedia.org/wiki/MRNA_(nucleoside-2%27-O-)-methyltransferase

# in front "the untranslated leader sequence that ends with the Transcription Regulation Sequence"
corona['untranslated_region'] = cc[0:265]

corona['orf1a'] = translate(cc[266-1:13483], True)

# cc[266-1+4398*3:13468] = 'TTT_TTA_AAC' aka 'X_XXY_YYZ'
# https://en.wikipedia.org/wiki/Ribosomal_frameshift
# Programmed −1 Ribosomal Frameshifting
# TODO: add this to the translate function with automatic detection
corona['orf1b'] = translate(cc[13468-1:21555], False).strip("*")  # chop off the stop, note this doesn't have a start

# exploit vector, this attaches to ACE2. also called "surface glycoprotein"
# https://www.ncbi.nlm.nih.gov/Structure/pdb/6VYB -- open state
# https://www.ncbi.nlm.nih.gov/Structure/pdb/6VXX -- closed state
# 1273 amino acids
#   S1  = 14-685
#   S2  = 686-1273
#   S2' = 816-1273
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2750777/
corona['spike_glycoprotein'] = translate(cc[21563-1:25384], True)

# Forms homotetrameric potassium sensitive ion channels (viroporin) and may modulate virus release.
corona['orf3a'] = translate(cc[25393-1:26220], True)

# these two things stick out, used in assembly aka they package the virus
corona['envelope_protein'] = translate(cc[26245-1:26472], True)  # also known as small membrane
corona['membrane_glycoprotein'] = translate(cc[26523-1:27191], True)

corona['orf6'] = translate(cc[27202-1:27387], True)

corona['orf7a'] = translate(cc[27394-1:27759], True)
corona['orf7b'] = translate(cc[27756-1:27887], True)  # is this one real?

corona['orf8'] = translate(cc[27894-1:28259], True)

# https://en.wikipedia.org/wiki/Capsid
# Packages the positive strand viral genome RNA into a helical ribonucleocapsid
# Includes the "internal" protein (from Coronavirus Pathogenesis)
# https://www.sciencedirect.com/topics/veterinary-science-and-veterinary-medicine/human-coronavirus-oc43
corona['nucleocapsid_phosphoprotein'] = translate(cc[28274-1:29533], True)

# might be called the internal protein (Coronavirus Pathogenesis)
corona['orf10'] = translate(cc[29558-1:29674], True)





# Component
# interface for objects in the composition of 
# leaf and composite nodes

class Component:
    # size of sequence
    def display(self):
        pass

# RNA sequence 
class RNALeaf(Component):
    def __init__(self, name, rna_sequence):
        self.__rna_sequence = rna_sequence
        self.__name = name
    def display(self):
        print("RNA sequence: " + self.__name + " " + str(self.__rna_sequence))


# Amino acid sequence        
class AminoAcidLeaf(Component):
    def __init__(self, name, amino_acid_sequence):
        self.__amino_acid_sequence = amino_acid_sequence
        self.__name = name

    def display(self):
        print("Amino Acid sequence: " + self.__name + " " + str(self.__amino_acid_sequence))

class Composite(Component): 
    def __init__(self):
        Component.__init__(self)
        self._children = list()

    def display(self):
        for child in self._children:
            child.operation()

    def add(self, component):
        self._children.append(component)

    def remove(self, component):
        self._children.remove(component)


def main ():
    # Creating a composite consisting of covid known RNA segments and Amino Acid 
    # segments. A different strain of corona virus with a mutated RNA and its
    # corresponding amino acid sequence can be created by constructing a new
    # RNA and amino acid leaf with the mutation and adding (replacing) the 
    # original sequence.
    rna_one = RNALeaf('untranslated_region', cc[0:265])
    rna_two = RNALeaf('orf1a',cc[266-1:13483])
    rna_three = RNALeaf('orf1b',cc[13468-1:21555])
    rna_four = RNALeaf('spike_glycoprotein',cc[21563-1:25384])
    rna_five = RNALeaf('orf3a',cc[25393-1:26220])
    rna_six = RNALeaf('envelope_protein',cc[26245-1:26472])
    rna_seven = RNALeaf('membrane_glycoprotein',cc[26523-1:27191])
    rna_eight = RNALeaf('orf6',cc[27202-1:27387])
    rna_nine = RNALeaf('orf7a',cc[27394-1:27759])
    rna_ten = RNALeaf('orf7b',cc[27756-1:27887])
    rna_eleven = RNALeaf('orf8',cc[27894-1:28259])
    rna_twelve = RNALeaf('nucleocapsid_phosphoprotein',cc[28274-1:29533])
    rna_thirteen = RNALeaf('orf10',cc[29558-1:29674])


  
    aa_one = AminoAcidLeaf('orf1a',translate(cc[266-1:13483], True))
    aa_two = AminoAcidLeaf('orf1b',translate(cc[13468-1:21555], False).strip("*"))
    aa_three = AminoAcidLeaf('spike_glycoprotein',translate(cc[21563-1:25384], True))
    aa_four = AminoAcidLeaf('orf3a',translate(cc[25393-1:26220], True))
    aa_five = AminoAcidLeaf('envelope_protein',translate(cc[26245-1:26472], True))
    aa_six = AminoAcidLeaf('membrane_glycoprotein',translate(cc[26523-1:27191], True))
    aa_seven = AminoAcidLeaf('orf6',translate(cc[27202-1:27387], True))
    aa_eight = AminoAcidLeaf('orf7a',translate(cc[27394-1:27759], True))
    aa_nine = AminoAcidLeaf('orf7b',translate(cc[27756-1:27887], True))
    aa_ten = AminoAcidLeaf('orf8',translate(cc[27894-1:28259], True))
    aa_eleven = AminoAcidLeaf('nucleocapsid_phosphoprotein',translate(cc[28274-1:29533], True))
    aa_twelve = AminoAcidLeaf('orf10',translate(cc[29558-1:29674], True))


    corona_composite = Composite()
    corona_composite.add(rna_one)
    corona_composite.add(rna_two)
    corona_composite.add(rna_three)
    corona_composite.add(rna_four)
    corona_composite.add(rna_five)
    corona_composite.add(rna_six)
    corona_composite.add(rna_seven)
    corona_composite.add(rna_eight)
    corona_composite.add(rna_nine)
    corona_composite.add(rna_ten)
    corona_composite.add(rna_eleven)
    corona_composite.add(rna_twelve)
    corona_composite.add(rna_thirteen)
    corona_composite.add(aa_one)
    corona_composite.add(aa_two)
    corona_composite.add(aa_three)
    corona_composite.add(aa_four)
    corona_composite.add(aa_five)
    corona_composite.add(aa_six)
    corona_composite.add(aa_seven)
    corona_composite.add(aa_eight)
    corona_composite.add(aa_nine)
    corona_composite.add(aa_ten)
    corona_composite.add(aa_eleven)
    corona_composite.add(aa_twelve)
