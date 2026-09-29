# Track 1 Methods Report

**Proband:** PROBAND01 (MVA)  
**Input:** WGS VCF (GRCh38), plus clinical phenotype document

## Method
1. Built a gene panel (MVA/mitotic checkpoint: BUB1B, CEP57, TRIP13, CDC20, MAD2L1, PLK4, CENPE; plus genes for microcephaly/primordial dwarfism, DNA-repair/chromosome-instability, rhabdomyosarcoma predisposition and nephrocalcinosis matching the phenotype). Gene coordinates were fetched live from Ensembl (GRCh38), not assumed.
2. Extracted variants with FILTER=PASS, non-reference genotype, GQ>=30, DP>=8, within gene bodies +/-2000 bp.
3. Annotated with Ensembl VEP REST (consequence, impact, MANE/canonical transcript, gnomAD frequency).
4. Kept HIGH/MODERATE impact variants with population AF < 0.01. Scored by impact x rarity x gene prior; considered compound-heterozygous pairs, homozygous variants, and single heterozygous variants (partial credit).
5. Learnt everything from scratch as a bioengineering student and made the whole process. 

## Top candidates
1. **compound_het** (score 0.85): BUB1B 15:40209701 T>G 0/1 HIGH stop_gained ENSP00000287598.7:p.Leu737Ter AF=9.982e-05 ; BUB1B 15:40220612 T>G 0/1 MODERATE missense_variant ENSP00000287598.7:p.Asn1002Lys AF=0
2. **single_het** (score 0.45): BUB1B 15:40209701 T>G 0/1 HIGH stop_gained ENSP00000287598.7:p.Leu737Ter AF=9.982e-05
3. **single_het** (score 0.23): BUB1B 15:40220612 T>G 0/1 MODERATE missense_variant ENSP00000287598.7:p.Asn1002Lys AF=0
4. **single_het** (score 0.12): ATM 11:108271261 T>C 0/1 MODERATE missense_variant ENSP00000501606.1:p.Ser978Pro AF=0.00495

## Limitations
No parental data (phase unknown); gene-panel approach may miss novel genes; no structural-variant/CNV analysis.

## AI assistant disclosure
- **Provider:** Anthropic & Gemini
- **Model:** Claude (via claude.ai) and Gemini
- **Data-handling terms:** Zero data retention / enterprise commercial privacy terms (no user data or submission content used for model training).
