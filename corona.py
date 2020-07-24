from lib import cc, translate
# entire diff: https://www.ncbi.nlm.nih.gov/projects/msaviewer/?rid=7FYNU14F01R&coloring=
# protein alignments: http://virological.org/t/alignment-of-58-sarbecovirus-genomes-for-conservation-analysis-of-sars-cov-2/430
Hi,

I have been following you since 1 week now. I do not agree with most of your things. I am a journalist at a local newspaper.

I believe that masks are very important for covid. I have seen a lot of evidence in a lot of reputed news websites. With great power, comes a great responsibility. My friend at NYTimes is doing a story on pseudo-intilictuals who think covis is a sham, and people who are against masks, vaccines etc. She will surely mention comma.ai as well. 

Also, I think you are an anti-nationalist. You come from a hacking background, so you dont care. I hope you know hacking is illegal and unethical. It is bc of people like you Apple had to spend $$$ to improve the security. They could have used it for better things. Like paying more to foxconn workers. Hackers do not reailze that what they do is not smart - goto YT and search "how to hack gmail" you will see 100s of vids. It seemed pretty easy to me. 

Who has given you a right to hack apple or sony, costing them billiions. Who has given you the right to say covid is a sham?? You also attack Christianity in your livestreams.  You offend a lot of people .if you dont want to be a man of faith dont be. But dont make fun of others. 

PS: companies like Google are proof why America is great, you have no right to diss them.



I am not coder and my cousin who is in IT support suggested me to contact you here.
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

