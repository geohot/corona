## Reverse engineering the Coronavirus

See [`corona.py`](corona.py)

### Open questions
- How is orf1ab cleaved into polypeptides? Can we predict this from the sequence?
- How do the researchers know (guess?) where orf1ab cleaves?
  - nsp3 and nsp5 do it -- https://www.pnas.org/content/pnas/103/15/5717.full.pdf
- Which protein is the immune system responding to?
  - "spike" and "nucleocapsid" -- http://www.cmi.ustc.edu.cn/1/3/193.pdf
  - Are some people already immune from exposure to other coronavirus?
- Find the "furin cleavage site" in the "spike glycoprotein"
  - It might be at the "PRRA" -- https://www.sciencedirect.com/science/article/pii/S0166354220300528
  - Use ProP or PiTou to predict? -- https://en.wikipedia.org/wiki/Furin
- How similar are the other coronaviruses? (causes colds, not either SARS or MERS)
  - alpha
    - https://en.wikipedia.org/wiki/Human_coronavirus_229E (simpler, though targets APN)
    - https://en.wikipedia.org/wiki/Human_coronavirus_NL63 (targets ACE2!)
      - https://www.ncbi.nlm.nih.gov/nuccore/MG772808
  - beta
    - https://en.wikipedia.org/wiki/Human_coronavirus_OC43 (targets Neu5Ac)
      - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2095096/pdf/JIDMM17330.pdf
      - Specifically, how similar is the N protein OC43, SARS v1, and SARS v2?
    - https://en.wikipedia.org/wiki/Human_coronavirus_HKU1 (targets Neu5Ac)
    - MERS-CoV
    - SARS-CoV
    - SARS-CoV-2
- What adds the phosphate group to the N protein? Kinase?

### Work to be done
- Automatic extraction of genes from different coronaviruses
- Good multisequence compare tool
- Molecular dynamics?
- Secondary Structure prediction on orf1a?

### Homemade Vaccine
- https://siasky.net/bACLKGmcmX4NCp47WwOOJf0lU666VLeT5HRWpWVtqZPjEA
- Based on injecting DNA (plasmid) that expresses the spike protein

### Bio vs IDA
- Each letter in genome is a "byte at address x"
- Translation = Disassembly, it's a 3 byte wide instruction set with arbitrary "reading frames"
- Protein ~ Function. polyprotein = Function with multiple pieces
- Proteins appear to have "basic blocks" = "Secondary Structure"
  - 80% accuracy in prediction: https://en.wikipedia.org/wiki/Protein_structure_prediction#Secondary_structure
- "Tertiary Structure" forms a function
  - This seems like the hard one to predict: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0205819
- "Quaternary Structure" is like compiled function with inlining
  - https://en.wikipedia.org/wiki/Protein%E2%80%93protein_interaction_prediction
- Gene ~ library (bacteria are static linked, viruses are dynamically linked)
- Transcription = loading off disk

### Analysis
- There is no equivalent to execution, we are reverse engineering a CAD format
- Static analysis, looking at the DNA, protein structure prediction, FLIRT signatures, etc...
- Simulation doesn't seem to work yet
- Tons of in system dynamic analysis, but the tools are crap
- Runs more like FPGA code, all at once, no serial execution (what are the FPGA re tools?)

### Tests (how they work)
- All based on https://en.wikipedia.org/wiki/Reverse_transcription_polymerase_chain_reaction
- USA -- https://www.fda.gov/media/134922/download
  - selected from regions of the virus nucleocapsid (N) gene
  - 28286---28308--28332---28358
  - 29163---29187--29210---29230
  - https://biosearchtech.a.bigcontent.io/v1/static/coa_KIT-NCOV-PP1-1000_Lot-No-143503
- South Korea -- http://www.kogene.co.kr/eng/about_us/news/listbody.php?h_gcode=board&h_code=7&po_no=288
  - E gene detection (same for all coronavirus)
  - specific RdRp detection

### Homemade test?
- Isolation of viral RNA (no matter what)
  - https://www.qiagen.com/us/products/diagnostics-and-clinical-research/sample-processing/qiaamp-viral-rna-mini-kit/#orderinginformation
- Primers and probes (to detect SARS-CoV-2)
  - https://www.biosearchtech.com/products/pcr-kits-and-reagents/pathogen-detection/2019-ncov-cdc-probe-and-primer-kit-for-sars-cov-2
  - Wouldn't need if using a nanopore sequencer (nanopore MinION)
- RT-qPCR Master Mix (to PCR)
  - https://www.thermofisher.com/order/catalog/product/A15300#/A15300
  - Probably wouldn't need if using a nanopore sequencer
- All in one?
  - https://www.chaibio.com/coronavirus
  - Open qPCR, understand https://www.chaibio.com/openqpcr
  - FAM and HEX fluorophores?

### Treatments 
- Hydroxychloroquine + Zinc
  - Zinc blocks RdRp
    - https://jvi.asm.org/content/91/21/e00754-17 -- how similar is Hep E RdRp?
    - https://www.ncbi.nlm.nih.gov/pubmed/21079686
  - Chloroquine Is a Zinc Ionophore (allows zinc into the cell)
    - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4182877/
- Novel RdRp inhibitors
  - Favipiravir (prodrug for favipiravir-RTP)
- Adenosine Analog
  - Remdesivir (prodrug for GS-441524)
  - Galidesivir

### Resources
- corona
  - Chapter 4 - Coronavirus Pathogenesis -- https://www.sciencedirect.com/science/article/pii/B9780123858856000092
  - https://www.futuremedicine.com/doi/pdf/10.2217/fvl-2018-0008
  - https://www.sciencedirect.com/science/article/pii/S2095177920302045
- textbooks
  - Molecular Biology of the Cell
- classes 
  - better tests - https://ocw.mit.edu/courses/biology/7-012-introduction-to-biology-fall-2004/index.htm
  - suspected better lectures - https://ocw.mit.edu/courses/biology/7-014-introductory-biology-spring-2005/index.htm
- Analysis of H1N1 - https://www.bunniestudios.com/blog/?p=353

