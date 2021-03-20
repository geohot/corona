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
# corona['spike_glycoprotein'] = translate(cc[21563-1:25384], True)

# # Forms homotetrameric potassium sensitive ion channels (viroporin) and may modulate virus release.
# corona['orf3a'] = translate(cc[25393-1:26220], True)

# # these two things stick out, used in assembly aka they package the virus
# corona['envelope_protein'] = translate(cc[26245-1:26472], True)  # also known as small membrane
# corona['membrane_glycoprotein'] = translate(cc[26523-1:27191], True)

# corona['orf6'] = translate(cc[27202-1:27387], True)

# corona['orf7a'] = translate(cc[27394-1:27759], True)
# corona['orf7b'] = translate(cc[27756-1:27887], True)  # is this one real?

# corona['orf8'] = translate(cc[27894-1:28259], True)

# # https://en.wikipedia.org/wiki/Capsid
# # Packages the positive strand viral genome RNA into a helical ribonucleocapsid
# # Includes the "internal" protein (from Coronavirus Pathogenesis)
# # https://www.sciencedirect.com/topics/veterinary-science-and-veterinary-medicine/human-coronavirus-oc43
# corona['nucleocapsid_phosphoprotein'] = translate(cc[28274-1:29533], True)

# # might be called the internal protein (Coronavirus Pathogenesis)
# corona['orf10'] = translate(cc[29558-1:29674], True)


class rna_sequence:
    def __init__(self, id, rna_sequence_range, is_protein):
        self.id = id
        self.rna_sequence_range = rna_sequence_range
        self.is_protein = is_protein
    #RNA sequence can be translated into amino acid sequence. So, each rna_sequence has its corresponding amino acid sequence.
    def amino_acid_sequence(self):
        return translate(self.rna_sequence_range, self.is_protein)



#RNA sequences and their gene 'names' 
untranslated_region = rna_sequence('untranslated_region', cc[0:265], False)
orf1a = rna_sequence('orf1a', cc[266-1:13483], True)
orf1b = rna_sequence('orf1b', cc[13468-1:21555], False)
spike_glycoprotein = rna_sequence('spike_glycoprotein', cc[21563-1:25384], True)
orf3a = rna_sequence ('orf3a', cc[25393-1:26220], True)
envelope_protein = rna_sequence ('envelope_protein', cc[26245-1:26472], True)
membrane_glycoprotein = rna_sequence('membrane_glycoprotein', cc[26523-1:27191], True) 
orf6 = rna_sequence('orf6', cc[27202-1:27387], True)
orf7a = rna_sequence ('orf7a', cc[27394-1:27759], True)
orf7b = rna_sequence ('orf7b', cc[27756-1:27887], True)
orf8 = rna_sequence ('orf8', cc[27894-1:28259], True)
nucleocapsid_phosphoprotein = rna_sequence ('nucleocapsid_phosphoprotein', cc[28274-1:29533], True)
orf10 = rna_sequence ('orf10', cc[29558-1:29674], True)

#the corona array contains: an RNA sequence of the untranslated region and an amino acid sequence of rest of the genes
corona['untranslated_region'] = untranslated_region.rna_sequence_range
corona['orf1a'] = orf1a.amino_acid_sequence
corona['orf1b'] = orf1b.amino_acid_sequence.strip("*")  # chop off the stop, note this doesn't have a start
corona['spike_glycoprotein'] = spike_glycoprotein.amino_acid_sequence
corona['orf3a'] = orf3a.amino_acid_sequence
corona['envelope_protein'] = envelope_protein.amino_acid_sequence
corona['membrane_glycoprotein'] = membrane_glycoprotein.amino_acid_sequence
corona['orf6'] = orf6.amino_acid_sequence
corona['orf7a'] = orf7a.amino_acid_sequence
corona['orf7b'] = orf7b.amino_acid_sequence
corona['orf8'] = orf8.amino_acid_sequence
corona['nucleocapsid_phosphoprotein'] = nucleocapsid_phosphoprotein.amino_acid_sequence
corona['orf10'] = orf10.amino_acid_sequence





