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

### Bio vs IDA
- Each letter in genome is a "byte at address x"
- Translation = Disassembly, it's a 3 byte wide instruction set with arbitrary "reading frames"
- Protein ~ Function. polyprotein = Function with multiple pieces
- Proteins appear to have "basic blocks" = "Secondary Structure"
- Gene = library (bacteria are static linked, viruses are dynamically linked)

### Analysis
- There is no equivalent to execution, we are reverse engineering a CAD format
- Static analysis, looking at the DNA, protein structure prediction, FLIRT signatures, etc...
- Simulation doesn't seem to work yet
- Tons of in system dynamic analysis, but the tools are crap
- Runs more like FPGA code, all at once, no serial execution (what are the FPGA re tools?)

