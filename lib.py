# Asn or Asp / B	AAU, AAC; GAU, GAC
# Gln or Glu / Z	CAA, CAG; GAA, GAG
# START	AUG
tt = """Ala / A	GCU, GCC, GCA, GCG
Ile / I	AUU, AUC, AUA
Arg / R	CGU, CGC, CGA, CGG; AGA, AGG
Leu / L	CUU, CUC, CUA, CUG; UUA, UUG
Asn / N	AAU, AAC
Lys / K	AAA, AAG
Asp / D	GAU, GAC
Met / M	AUG
Phe / F	UUU, UUC
Cys / C	UGU, UGC
Pro / P	CCU, CCC, CCA, CCG
Gln / Q	CAA, CAG
Ser / S	UCU, UCC, UCA, UCG; AGU, AGC
Glu / E	GAA, GAG
Thr / T	ACU, ACC, ACA, ACG
Trp / W	UGG
Gly / G	GGU, GGC, GGA, GGG
Tyr / Y	UAU, UAC
His / H	CAU, CAC
Val / V	GUU, GUC, GUA, GUG
STOP	UAA, UGA, UAG
""".strip()
dec = {}
for t in tt.split("\n"):
  k,v = t.split("\x09")
  if '/' in k:
    k = k.split("/")[-1].strip()
  k = k.replace("STOP", "*")
  v = v.replace(",", "").replace(";", "").lower().replace("u", "t").split(" ")
  for vv in v:
    if vv in dec:
      print("dup", vv)
    dec[vv] = k

def translate(x, protein=False):
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

