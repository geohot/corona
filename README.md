## Reverse engineering the Coronavirus

See [`corona.py`](corona.py)

### Open questions
- How is orf1ab cleaved into polypeptides? Can we predict this from the sequence?
- How do the researchers know (guess?) where orf1ab cleaves?
  - nsp3 and nsp5 do it -- https://www.pnas.org/content/pnas/103/15/5717.full.pdf
- Which protein is the immune system responding to?
  - "spike" and "nucleocapsid" -- http://www.cmi.ustc.edu.cn/1/3/193.pdf
- Find the "furin cleavage site" in the "spike glycoprotein"
  - It might be at the "PRRA" -- https://www.sciencedirect.com/science/article/pii/S0166354220300528
  - Use ProP or PiTou to predict? -- https://en.wikipedia.org/wiki/Furin

### Work to be done
- Automatic extraction of genes from different coronaviruses
- Good multisequence compare tool
- Molecular dynamics?

### Homemade Vaccine
- https://siasky.net/bACLKGmcmX4NCp47WwOOJf0lU666VLeT5HRWpWVtqZPjEA
- Based on injecting DNA (plasmid) that expresses the spike protein


