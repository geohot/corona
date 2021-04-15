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




# The director which controls the construction process
class Director:
    __builder = None

    def setBuilder(self, builder):
        self.__builder = builder

    def getOrganism(self):
        organism = Organism()

        # name of organism 
        name = self.__builder.get_name()
        organism.set_name(name)

        # DNA sequence
        dna_seq = self.__builder.get_dna_seq()
        organism.set_dna_seq(dna_seq)

        # RNA sequence 
        rna_seq = self.__builder.get_rna_seq()
        organism.set_rna_seq(rna_seq)

        # Known Gene list
        genes = self.__builder.getGenes()
        organism.setGenes(genes)
        
        return organism

# The whole organism
class Organism: 
    def __init__(self):
        self.__name = None
        self.__dna_seq = None
        self.__rna_seq = None
        self.__genes = None
    
    def set_name(self, name):
        self.__name = name

    def set_dna_seq(self,dna_seq):
        self.__dna_se = dna_seq
    
    def set_rna_seq(self,rna_seq):
        self.__rna_seq = rna_seq
    
    def setGenes(self, genes):
        self.__genes = genes
    
    def specification(self):
        print "organism name: %s" % self.__name.organismName
        print "number of genes: %s" %self.__genes.len(geneNameSeq_tuple)
        
#    abstract interface builder
class Builder:
    def get_name(self): pass
    def get_dna_seq(self): pass
    def get_rna_seq(self): pass
    def getGenes(self): pass



# Organism parts
class Name:
    organismName = None
class DNA_seq:
    sequence = None
class RNA_seq:
    sequence = None
class Genes:
    geneNameList = list()
    geneSeqList = list()
    geneNameSeq_tuple = list(zip(geneNameList, geneSeqList))


class coronaVirusBuilder(Builder):
    # Concrete Builder implementation.
    # Builds parts for an RNA only organism

    def getName(self):
        organismName = Name()
        organismName.name = "COVID-19: strain 1"


    def get_dna_seq(self):
        dna_seq = DNA_seq()
        dna_seq.sequence = None
        return dna_seq

    def get_rna_seq(self):
        rna_seq = RNA_seq()
        rna_seq.sequence = cc[0:29674]
        return rna_seq
    
    def get_genes(self):
        genes = Genes()
        genes.geneNameList.extend('untranslated_region','orf1a','orf1b', 'spike_glycoprotein',
        'orf3a','envelope_protein','membrane_glycoprotein', 'orf6', 'orf7a', 'orf7b', 'orf8',
        'nucleocapsid_phosphoprotein', 'orf10')
        genes.geneSeqList.extend(cc[0:265], cc[266-1:13483], cc[13468-1:21555], cc[21563-1:25384],
        cc[25393-1:26220], cc[26245-1:26472], cc[26523-1:27191], cc[27202-1:27387], cc[27394-1:27759],
        cc[27756-1:27887], cc[27894-1:28259], cc[28274-1:29533], cc[29558-1:29674])
        return genes

# construct a corona virus object containing: the name of the organism,
# the entire RNA sequence, and a tuple-list of the names of genes and corresponding RNA sequence 
def main():
    covidBuilder = coronaVirusBuilder()

    director = Director()

    # Build covid19
    director.setBuilder(covidBuilder)
    covid = director.getOrganism()
    covid.specification()
