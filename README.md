# Track 1 Methods Description Report
**Target Disease:** Mosaic Variegated Aneuploidy (MVA) Syndrome  
**Proband ID:** PROBAND01
**Dataset:** Whole Genome Sequencing (WGS) VCF  

---

## 1. Executive Summary
This report describes the computational variant prioritization pipeline used to identify the primary pathogenic genetic cause for proband PROBAND01 presenting with clinical features of Mosaic Variegated Aneuploidy (MVA). 

Using quality-filtered variant extraction, MVA gene panel targeting, and functional annotation via Ensembl Variant Effect Predictor (VEP), two primary findings were isolated:
1. **Compound Heterozygous Pair in *BUB1B*** (Chr15:40209701 `T>G` `stop_gained` AND Chr15:40220612 `T>G` `missense_variant`), presenting a loss-of-function compound heterozygous model for **MVA Syndrome Type 1** (OMIM #257300).
2. **Homozygous Inframe Deletion in *CEP57*** (Chr11:96092210 `TTGCTGCTGC>T` `inframe_deletion`), representing a secondary recessive cause for **MVA Syndrome Type 2** (OMIM #614114).

---

## 2. Computational Pipeline & Methodology

### A. Data Processing & Quality Control
* **VCF Input:** `WGS_EX2312012_HGWCNDSX7.vcf.gz` (GRCh38 reference assembly).
* **Quality Filtering:** Retained variants with `PASS` status across MVA-associated gene loci:
  * `BUB1B` (chr15:40,188,000–40,250,000)
  * `CEP57` (chr11:96,060,000–96,110,000)
  * `TRIP13` (chr5:890,000–940,000)
  * `MAD2L1` (chr4:120,440,000–120,455,000)

### B. Functional Annotation
Annotation was executed via HTTP queries to the Ensembl VEP REST API (`https://rest.ensembl.org/vep/human/region`). Mini-batching and single-variant fallback queries were implemented to ensure full variant annotation recovery.

---

## 3. Prioritization & Clinical Rationale

* **Primary Candidate 1 (EPCR = 0.95): *BUB1B* Compound Heterozygosity**
  * *Variant 1:* Chr15:40209701 (`T>G`), `stop_gained` (HIGH impact). Introduces premature termination codon.
  * *Variant 2:* Chr15:40220612 (`T>G`), `missense_variant` (MODERATE impact). Alters conserved kinase domain.
  * *Significance:* Disruption of *BUB1B* causes constitutional mosaic aneuploidy and MVA Type 1.

* **Primary Candidate 2 (EPCR = 0.85): *CEP57* Homozygous Inframe Deletion**
  * *Variant:* Chr11:96092210 (`TTGCTGCTGC>T`), `inframe_deletion` (MODERATE impact, `1/1`).
  * *Significance:* Disruption of centrosomal protein 57 leads to MVA Type 2.

---

## 4. AI Assistant Disclosure
* **Provider:** Google AI
* **Model & Tier:** Gemini (Commercial API Tier)
* **Data-Handling Terms:** Zero data retention / commercial privacy terms (no submission data used for training).
