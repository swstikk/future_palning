# 🟡 TIER 3: ADVANCED & INTEGRATIVE
## Aging Biology → Immunology → Genomics Tech → Biostatistics → Computational Prerequisites

> **PRIORITY:** IMPORTANT — Learn within Month 6-12
> **PREREQUISITE:** Tier 1 + Tier 2 checklists complete, Tier 1B in progress
> **FOCUS:** This tier connects biology to your computational tools

---

## 3.1 THE 12 HALLMARKS OF AGING — THE COMPLETE MAP
**Why:** This IS the aging field. Every intervention targets one or more hallmarks.

```
López-Otín et al. 2023 — "Hallmarks of Aging: An Expanding Universe" (Cell)
→ READ THIS PAPER. It's free on PubMed. It's the bible of aging research.

THE 12 HALLMARKS:

═══════════════════════════════════════════
PRIMARY HALLMARKS (damage accumulation):
═══════════════════════════════════════════

1. GENOMIC INSTABILITY:
├── DNA damage accumulates: ~10,000 lesions/cell/day → most repaired, but not all
├── Somatic mutation burden increases linearly with age
├── Transposon reactivation: heterochromatin erosion → LINE1 elements activate
├── Clonal hematopoiesis: blood cells with somatic mutations expand with age
│   └── CHIP (Clonal Hematopoiesis of Indeterminate Potential) → cardiovascular risk!
└── Intervention: boost DNA repair (SIRT6, WRN helicase), silence transposons

2. TELOMERE ATTRITION:
├── Telomeres: TTAGGG repeats at chromosome ends → shorten each division (~50-200bp)
├── Shelterin complex: protects telomere ends (TRF1, TRF2, POT1, TIN2, TPP1, RAP1)
├── When too short → DNA damage response → senescence or apoptosis
├── Telomerase: TERT + TERC → extends telomeres (active in stem cells, not somatic)
├── Measurement: qPCR, Flow-FISH, TRF by Southern blot, nanopore
└── Intervention: TERT activation (AAV-TERT in mice extended lifespan!)

3. EPIGENETIC ALTERATIONS:
├── DNA methylation drift → aging clocks (covered in Tier 2.4)
├── Histone mark redistribution → gene derepression
├── Heterochromatin loss → transposon activation → genomic instability (connects #1!)
└── Intervention: epigenetic reprogramming (OSKM, dCas9-TET1)

4. LOSS OF PROTEOSTASIS:
├── Chaperone decline (HSP70/90 expression drops)
├── Proteasome activity drops ~50% with age
├── Autophagy-lysosome pathway declines
├── Result: misfolded protein AGGREGATION → Alzheimer's (Aβ), Parkinson's (α-synuclein)
└── Intervention: autophagy inducers (rapamycin, spermidine), proteasome activators

═══════════════════════════════════════════
ANTAGONISTIC HALLMARKS (responses that become harmful):
═══════════════════════════════════════════

5. DEREGULATED NUTRIENT SENSING:
├── The 4 pathways: mTOR↑, AMPK↓, Sirtuins↓, Insulin/IGF-1↑ with age
├── All covered in detail in Tier 2.3
└── Intervention: caloric restriction, rapamycin, metformin, NMN/NR

6. MITOCHONDRIAL DYSFUNCTION:
├── ETC efficiency drops → more ROS, less ATP
├── mtDNA mutations accumulate (no histones, near ROS source, poor repair)
├── Mitochondrial dynamics shift: less fusion, more fission → fragmented mito
├── NAD⁺ decline → ETC Complex I underperforming
└── Intervention: mitophagy enhancers, DdCBE for mtDNA repair, NMN/NR

7. CELLULAR SENESCENCE:
├── All covered in Tier 1.6 (SASP, markers, senolytics, senomorphics)
├── Baker 2016: clearing p16+ cells → lifespan extension in mice
└── Intervention: senolytics (D+Q), senomorphics (rapamycin), CAR-T senolytics

═══════════════════════════════════════════
INTEGRATIVE HALLMARKS (consequences):
═══════════════════════════════════════════

8. STEM CELL EXHAUSTION:
├── Stem cells decline in number AND function with age
├── HSCs (hematopoietic): bias toward myeloid lineage (fewer lymphoid → immunosenescence)
├── Muscle satellite cells: fewer → sarcopenia (muscle wasting)
├── Neural stem cells: decline → less neurogenesis → cognitive decline
├── Intestinal stem cells: Lgr5+ cells → decline → gut barrier dysfunction
└── Intervention: ex vivo expansion + CRISPR editing → autologous transplant

9. ALTERED INTERCELLULAR COMMUNICATION:
├── SASP: senescent cells poison neighboring cells (paracrine senescence!)
├── Inflammaging: chronic low-grade systemic inflammation
├── Extracellular vesicles: change cargo with age
├── Blood factors: young blood factors (GDF11? — debated) vs old blood (CCL11, β2M)
│   └── Heterochronic parabiosis: connect young+old mouse → old mouse rejuvenates!
└── Intervention: plasma dilution, young plasma factors, anti-SASP therapy

10. DISABLED MACROAUTOPHAGY:
├── Autophagy: cell self-eating → degrades damaged organelles + proteins
├── Process: cargo → autophagosome → fusion with lysosome → degradation
├── Key regulators: ULK1, Beclin-1, ATG proteins, LC3 (autophagosome marker)
├── mTOR inhibits ULK1 → autophagy blocked
├── AMPK activates ULK1 → autophagy promoted
└── Intervention: autophagy inducers → rapamycin, spermidine, fasting

11. CHRONIC INFLAMMATION (Inflammaging):
├── NF-κB chronically active (covered in Tier 2.3)
├── SASP factors: IL-6, TNF-α, IL-1β → systemic inflammation
├── NLRP3 inflammasome: activated by cellular debris → IL-1β maturation
├── cGAS-STING pathway: cytoplasmic DNA (from damaged nuclei/mito) → type I IFN → inflammation
│   └── This is a MAJOR new target for aging intervention!
└── Intervention: anti-inflammatory (targeted NF-κB, NLRP3 inhibitors, STING inhibitors)

12. DYSBIOSIS:
├── Gut microbiome changes (covered in Tier 2.5)
├── Leaky gut → LPS enters bloodstream → systemic inflammation
├── Microbiome-gut-brain axis: affects cognition in aging
└── Intervention: FMT from young donors, probiotics, dietary fiber

HALLMARK INTERCONNECTIONS (this is key — they're not independent!):
├── Genomic instability (#1) → triggers DDR → senescence (#7)
├── Senescence (#7) → SASP → inflammation (#11) → altered communication (#9)
├── Epigenetic alterations (#3) → gene derepression → transposon activation → genomic instability (#1)
├── Mitochondrial dysfunction (#6) → ROS → DNA damage (#1) + proteotoxic stress (#4)
├── Nutrient sensing (#5) → mTOR↑ → autophagy↓ (#10) → proteostasis loss (#4)
└── Everything converges → inflammaging (#11) → tissue dysfunction → disease!
```

---

## 3.2 IMMUNOLOGY — THE AGING IMMUNE SYSTEM
**Why:** Immune rejuvenation is a major strategy. CAR-T senolytics require immunology.

```
═══════════════════════════════════════════
A. INNATE vs ADAPTIVE IMMUNITY:
═══════════════════════════════════════════

Innate (fast, non-specific, no memory):
├── Barriers: skin, mucus, stomach acid, antimicrobial peptides
├── Cells:
│   ├── Neutrophils: first responders, phagocytosis, NETs (neutrophil extracellular traps)
│   ├── Macrophages: phagocytosis, antigen presentation, M1 (pro-inflammatory) vs M2 (repair)
│   ├── Dendritic cells: antigen presentation to T cells (bridge innate↔adaptive)
│   ├── NK cells: kill stressed/infected cells without prior sensitization
│   │   └── Recognize "missing self" (cells that downregulate MHC-I → cancer/virus)
│   └── Mast cells: histamine release, allergy
├── Pattern Recognition Receptors (PRRs):
│   ├── TLRs (Toll-like receptors): detect pathogen patterns (LPS → TLR4)
│   ├── cGAS-STING: detect cytoplasmic DNA → type I interferon
│   │   └── AGING: damaged nuclear/mitochondrial DNA leaks → chronic STING activation!
│   └── NLRP3 inflammasome: detect damage signals → cleave IL-1β → inflammation
│       └── AGING: constitutively active → inflammaging

Adaptive (slow, specific, has MEMORY):
├── T cells (mature in Thymus):
│   ├── CD4+ Helper: coordinate response, activate B cells and macrophages
│   │   └── Th1 (antiviral), Th2 (antiparasitic), Th17 (autoimmune), Treg (suppress)
│   ├── CD8+ Cytotoxic (CTL): kill infected/cancer cells via perforin/granzymes
│   ├── Tregs (CD4+CD25+FoxP3+): suppress immune response → prevent autoimmunity
│   ├── Memory T cells: long-lived, fast response on re-exposure
│   └── T cell receptor (TCR): recognizes antigen presented on MHC
├── B cells (mature in Bone marrow):
│   ├── Produce antibodies (immunoglobulins): IgM, IgG, IgA, IgE, IgD
│   ├── Antigen presentation to T cells
│   ├── Somatic hypermutation: V(D)J recombination → diverse antibody repertoire
│   └── Memory B cells + plasma cells → long-term antibody production

MHC (Major Histocompatibility Complex):
├── MHC-I: on ALL nucleated cells → presents intracellular peptides to CD8+ T cells
│   └── If cell infected → viral peptides on MHC-I → CTL kills cell
├── MHC-II: on APCs (dendritic, macrophage, B cell) → presents extracellular peptides to CD4+
└── HLA (Human Leukocyte Antigen) = human MHC → crucial for transplant compatibility

═══════════════════════════════════════════
B. IMMUNOSENESCENCE (aging of immune system):
═══════════════════════════════════════════
├── Thymic involution: thymus shrinks dramatically after puberty
│   └── By age 50: >90% replaced by fat → very few new naive T cells produced
├── Naive T cell pool → depleted → less ability to respond to NEW threats
├── Memory T cell accumulation: some become SENESCENT (express p16, SASP-like)
│   └── CMV-driven: cytomegalovirus causes massive clonal expansion of CD8+ T cells
├── B cell defects: fewer naive B cells, reduced antibody diversification
├── NK cell: number maintained but FUNCTION declines (less cytotoxicity)
├── Innate immune cells: macrophages become more inflammatory (M1-skewed)
├── Vaccine responses decline: flu vaccine only ~30% effective in elderly (vs ~70% young)
└── NET result: IMMUNOSENESCENCE + chronic inflammation = "inflamm-aging"

═══════════════════════════════════════════
C. IMMUNE INTERVENTIONS FOR AGING:
═══════════════════════════════════════════
├── Thymus rejuvenation: FOXN1 overexpression → regrow thymus (proof of concept in mice!)
├── CAR-T senolytics (Amor et al., 2024):
│   ├── Engineer T cells with CAR targeting uPAR (senescent cell surface marker)
│   ├── CAR-T cells selectively kill senescent cells → reverse aging phenotypes!
│   └── This is immunology + genetic engineering + aging biology COMBINED
├── NK cell rejuvenation: expand/activate NK cells ex vivo → reinfuse
├── Immune checkpoint: PD-1/PD-L1 blockade (cancer) → might clear senescent cells too
├── Plasma exchange/dilution: Irina Conboy's work → dilute old blood factors → rejuvenate
└── Rapamycin: low-dose → improves vaccine response in elderly (Mannick 2014)

🔗 BIOINFORMATICS: 
├── scRNA-seq of immune cells → identify senescent T cell populations
├── CITE-seq: surface marker + transcriptome simultaneously
├── TCR/BCR repertoire sequencing: measure immune diversity (declines with age)
└── Module 20 in DEEP SYLLABUS: full immunosenescence analysis
```

---

## 3.3 GENOMICS TECHNOLOGIES — THE DATA YOU'LL ANALYZE
**Why:** You must understand HOW data is generated to analyze it properly.

```
═══════════════════════════════════════════
A. DNA SEQUENCING TECHNOLOGIES:
═══════════════════════════════════════════

Sanger Sequencing (1977):
├── Chain termination method: ddNTPs → fluorescent termination
├── Read length: up to ~1000bp
├── Throughput: 1 sequence at a time → slow
├── Use now: verification of plasmid clones, single-gene diagnostics
└── Gold standard for accuracy → still used to validate NGS findings

Illumina (Short-read NGS) — WORKHORSE:
├── Sequencing by synthesis (SBS):
│   ├── Library prep: fragment DNA → add adapters → amplify on flow cell
│   ├── Bridge amplification: create clusters of identical fragments
│   ├── Fluorescent reversible terminators: one base added per cycle
│   └── Image after each cycle → identify base by color
├── Read length: 1×75bp, 2×150bp, 2×300bp (paired-end most common)
├── Output: millions to BILLIONS of reads per run
├── Cost: ~$200-600 for 30× human genome (2025 prices)
├── Platforms: MiSeq (small), NextSeq (medium), NovaSeq (high throughput)
├── Applications: WGS, WES, RNA-seq, ChIP-seq, ATAC-seq, WGBS, scRNA-seq
└── Limitation: short reads → repetitive regions, structural variants challenging

PacBio (Long-read):
├── SMRT (Single Molecule Real-Time) sequencing
├── Read length: 10-30kb (HiFi mode: 15-20kb at >99.9% accuracy!)
├── Sees through: repeats, structural variants, haplotype phasing
├── CCS (Circular Consensus Sequencing): read same molecule multiple times → HiFi reads
└── Applications: de novo assembly, full-length transcript (Iso-Seq), methylation detection

Oxford Nanopore:
├── Real-time: ssDNA threaded through protein pore → current changes = base identity
├── Read length: 10kb-1Mb+ (longest reads of ANY technology!)
├── Portable: MinION = USB-stick-sized sequencer → field deployable!
├── Detects base modifications DIRECTLY: 5mC, 6mA (no bisulfite needed!)
├── Platforms: MinION (portable), GridION (5 flow cells), PromethION (high throughput)
└── Applications: structural variants, metagenomics, real-time surveillance, methylation

═══════════════════════════════════════════
B. FUNCTIONAL GENOMICS ASSAYS:
═══════════════════════════════════════════

Transcriptomics:
├── RNA-seq (bulk): measure gene expression across ALL cells → averages
│   └── Pipeline: FASTQ → QC (FastQC) → Trim → Align (STAR) → Count (featureCounts) → DE (DESeq2)
├── scRNA-seq (single-cell): expression in INDIVIDUAL cells → cell types!
│   ├── 10x Genomics Chromium: droplet-based, ~10,000 cells/run
│   ├── Smart-seq2/3: plate-based, full-length transcripts, fewer cells
│   └── Pipeline: CellRanger → Seurat (R) or Scanpy (Python)
├── Spatial transcriptomics: gene expression WITH tissue location
│   ├── Visium (10x): spots on tissue section → ~55µm resolution
│   ├── MERFISH: fluorescence in situ → subcellular resolution
│   └── Slide-seq: bead-based → near single-cell spatial
└── Ribo-seq: measure actively TRANSLATED mRNAs (ribosome profiling)

Epigenomics:
├── WGBS: whole genome bisulfite sequencing → single-CpG methylation (gold standard)
├── RRBS: reduced representation → enriched for CpG islands (cheaper)
├── Infinium 450K/EPIC array: probe-based → 450K/850K CpGs (clock data!)
├── ChIP-seq: protein-DNA binding (TFs, histone marks)
│   └── Antibody → immunoprecipitate → sequence bound DNA
├── CUT&RUN / CUT&TAG: lower input alternatives (100 cells vs millions)
├── ATAC-seq: open chromatin → Tn5 transposase inserts adapters into accessible regions
└── Hi-C / Micro-C: 3D genome organization (crosslink → ligate → sequence)

Proteomics & Metabolomics:
├── Mass spectrometry (MS): identify + quantify proteins/metabolites
│   ├── LC-MS/MS: liquid chromatography + tandem MS → standard proteomics
│   ├── TMT (tandem mass tags): multiplex → compare 16+ samples simultaneously
│   └── Metabolomics: untargeted (discover) vs targeted (quantify known metabolites)
├── Protein arrays: antibody-based, high throughput
└── Proximity ligation assay (PLA): detect protein-protein interactions in situ

═══════════════════════════════════════════
C. GENOME ASSEMBLY — HOW GENOMES ARE BUILT:
═══════════════════════════════════════════
⚠️ GAP FIX: Was missing from original roadmap.

├── Reference-based: align reads to EXISTING reference genome
│   └── Most common for human/mouse (reference genomes available)
│
├── De novo assembly: NO reference available → build from scratch!
│   ├── de Bruijn graph: reads → k-mers → graph → traverse paths → contigs
│   │   └── Assemblers: SPAdes (bacterial), MEGAHIT (metagenomics)
│   ├── OLC (Overlap-Layout-Consensus): for long reads
│   │   └── Assemblers: Canu, Flye (long-read de novo)
│   ├── Steps: reads → contigs → scaffolds → chromosome-level assembly
│   │   └── Hi-C data helps scaffold contigs into chromosomes!
│   └── Quality metrics:
│       ├── N50: 50% of assembly is in contigs ≥ this length (higher = better)
│       ├── BUSCO: completeness check → universally conserved genes present?
│       └── Misassemblies: QUAST → detect structural errors
│
└── WHY for aging: 
    ├── Non-model aging organisms: naked mole rat, bowhead whale, Hydra → need assembly
    ├── Metagenomics: assemble unknown microbiome species
    └── Understanding assembly → better understanding of alignment/mapping
```

---

## 3.4 BIOSTATISTICS FOR BIOLOGY — YOUR DATA WILL LIE WITHOUT THIS
**Why:** Every -omics experiment generates thousands of statistical tests. Wrong stats = wrong conclusions.

```
⚠️ GAP FIX: Biology-specific statistics were missing.

═══════════════════════════════════════════
A. MULTIPLE TESTING CORRECTION — THE BIGGEST TRAP:
═══════════════════════════════════════════
├── Problem: RNA-seq tests ~20,000 genes → 20,000 p-values
│   └── At α=0.05: expect 1,000 FALSE positives by chance alone!
│
├── Bonferroni correction: p_adjusted = p × n_tests
│   └── Too conservative: almost nothing survives for large n
│
├── False Discovery Rate (FDR) — Benjamini-Hochberg:
│   ├── Sort p-values smallest→largest
│   ├── Adjusted p_i = p_i × (n/rank_i)
│   ├── FDR 5% = expect ≤5% of reported results are false positives
│   └── THIS IS THE STANDARD for -omics data!
│
├── In practice:
│   ├── padj < 0.05 (after BH correction) = significant
│   ├── Report BOTH raw p-value and adjusted p-value
│   └── Tools: p.adjust() in R, statsmodels.multipletests() in Python
│
└── NEVER publish -omics results without multiple testing correction!

═══════════════════════════════════════════
B. SURVIVAL ANALYSIS — HOW AGING STUDIES ARE EVALUATED:
═══════════════════════════════════════════
├── Kaplan-Meier curves: plot survival probability over time
│   ├── Y-axis: fraction surviving, X-axis: time
│   ├── Steps (drops) at each death event
│   ├── Censoring: animal still alive at end of study → line continues
│   └── Compare curves: log-rank test (p-value for difference)
│
├── Median survival: time at which 50% have died → primary aging endpoint!
│
├── Maximum lifespan: usually measured as 90th percentile or last survivor
│
├── Cox Proportional Hazards regression:
│   ├── Hazard ratio (HR): risk of death in treatment vs control
│   ├── HR < 1: treatment REDUCES death risk → GOOD
│   ├── HR > 1: treatment INCREASES death risk → BAD
│   └── Can include covariates (age, sex, genotype → adjust for confounders)
│
├── Tools: lifelines (Python), survival (R)
└── EVERY aging intervention paper uses Kaplan-Meier + log-rank test

═══════════════════════════════════════════
C. BATCH EFFECTS — THE SILENT KILLER OF GENOMICS:
═══════════════════════════════════════════
├── What: technical variation between experiment runs (different days, reagent lots, operators)
├── Problem: batch effect can be LARGER than biological signal!
│   └── Example: PCA shows samples cluster by processing date, NOT by treatment!
├── Detection: PCA → color by batch → do batches separate? → PROBLEM
├── Correction:
│   ├── ComBat (sva package in R): empirical Bayes batch correction
│   ├── Harmony: integration method for scRNA-seq batches
│   └── limma::removeBatchEffect(): simple removal for visualization
├── Prevention: randomize samples across batches, include controls in each batch
└── ALWAYS check for batch effects BEFORE differential expression analysis!

═══════════════════════════════════════════
D. EFFECT SIZES & POWER:
═══════════════════════════════════════════
├── P-value ≠ biological importance!
│   └── Very large dataset → tiny meaningless difference can be "significant"
├── Effect size:
│   ├── log2 Fold Change (log2FC): common for gene expression
│   │   └── |log2FC| > 1 (2-fold change) = biologically meaningful threshold
│   ├── Cohen's d: standardized mean difference
│   └── Odds ratio / Hazard ratio: for categorical/survival outcomes
├── Volcano plot: log2FC (x) vs -log10(p) (y) → see significance AND effect size
│   └── Upper-right/upper-left quadrants = significant AND large effect → TRUE hits
├── Power analysis: how many samples do you need?
│   └── Too few → miss real effects (underpowered)
│   └── Too many → waste resources
│   └── Tools: G*Power, statsmodels in Python
└── Always report: adjusted p-value + effect size + confidence interval
```

---

## 3.5 COMPUTATIONAL PREREQUISITES — THE TOOLS YOU NEED
**Why:** You can't do bioinformatics without these. Learn alongside the biology.

```
⚠️ GAP FIX: Linux/CLI and R were missing from the biology roadmap.

═══════════════════════════════════════════
A. LINUX COMMAND LINE:
═══════════════════════════════════════════
├── WHY: 90%+ of bioinformatics tools run on Linux CLI
│   ├── GATK, BWA, STAR, samtools, bcftools, Picard → ALL CLI
│   ├── HPC clusters → Linux servers → no GUI
│   └── Docker/Singularity containers → Linux
│
├── ESSENTIAL COMMANDS:
│   ├── Navigation: cd, ls, pwd, mkdir, rmdir
│   ├── Files: cp, mv, rm, touch, cat, head, tail, less, wc
│   ├── Text processing: grep, awk, sed, cut, sort, uniq, tr
│   ├── Pipes: command1 | command2 (chain operations)
│   ├── Redirect: > (overwrite), >> (append), 2> (error redirect)
│   ├── File permissions: chmod, chown
│   ├── Process: ps, top, kill, nohup, screen, tmux
│   └── Remote: ssh, scp, rsync
│
├── Package managers:
│   ├── apt (Ubuntu/Debian): sudo apt install samtools
│   ├── conda: conda install -c bioconda star
│   └── pip: pip install scanpy
│
├── WHERE TO LEARN:
│   ├── Linux Journey (linuxjourney.com) → free, beginner-friendly
│   ├── Software Carpentry shell lesson → bioinformatics-focused
│   └── WSL (Windows Subsystem for Linux) → run Linux on your Windows PC!
│
└── You're on Windows → Install WSL2 → Ubuntu → then all bio tools work!

═══════════════════════════════════════════
B. R PROGRAMMING & BIOCONDUCTOR:
═══════════════════════════════════════════
├── WHY R (in addition to Python)?
│   ├── DESeq2 (differential expression) → R only
│   ├── Seurat (single-cell analysis) → R (Scanpy is Python alternative)
│   ├── limma, edgeR → statistical testing → R only
│   ├── ggplot2 → publication-quality figures → R standard
│   ├── Bioconductor: 2,200+ packages for biological data analysis
│   └── Many tutorials and papers use R → you MUST be able to read R code
│
├── WHAT TO LEARN:
│   ├── Base R: data types, vectors, data frames, functions, loops, apply family
│   ├── tidyverse: dplyr (data manipulation), tidyr (reshaping), ggplot2 (plotting)
│   ├── Bioconductor basics: install packages, use SummarizedExperiment object
│   ├── DESeq2: full RNA-seq differential expression pipeline
│   ├── GenomicRanges: manipulate genomic intervals (BED, GFF)
│   └── ComplexHeatmap: heatmaps for -omics data
│
├── WHERE TO LEARN:
│   ├── Bioconductor vignettes → best documentation for each package
│   ├── R for Data Science (r4ds.hadley.nz) → free online book
│   └── Orchestrating Single-Cell Analysis with Bioconductor (OSCA) → free online book
│
└── Priority: Learn R AFTER Python basics, focus on DESeq2 + ggplot2 first

═══════════════════════════════════════════
C. DATABASES YOU'LL USE DAILY:
═══════════════════════════════════════════
├── SEQUENCE:
│   ├── NCBI GenBank: all publicly available sequences
│   ├── UniProt: protein sequences + annotation (SwissProt = curated)
│   ├── Ensembl: genome browser + gene models + comparative genomics
│   └── UCSC Genome Browser: visualization + tracks + tools
│
├── STRUCTURE:
│   ├── PDB (Protein Data Bank): experimentally determined 3D structures
│   ├── AlphaFold DB: predicted structures for ~200M proteins
│   └── Swiss-Model: homology-based structure prediction
│
├── PATHWAYS & FUNCTIONS:
│   ├── KEGG: metabolic/signaling pathways
│   ├── GO (Gene Ontology): functional annotation categories
│   ├── Reactome: pathway database (more curated than KEGG)
│   └── STRING: protein-protein interaction networks
│
├── VARIANTS & DISEASE:
│   ├── ClinVar: clinical significance of variants
│   ├── gnomAD: population allele frequencies
│   ├── COSMIC: somatic mutations in cancer
│   └── OMIM: Mendelian disease genetics
│
├── AGING-SPECIFIC:
│   ├── GenAge: genes associated with aging
│   ├── AnAge: animal aging and longevity database
│   ├── HAGR (Human Ageing Genomic Resources)
│   └── SeneQuest: senescence gene database
│
└── EPIGENETICS:
    ├── ENCODE: encyclopedia of regulatory elements
    ├── Roadmap Epigenomics: tissue-specific epigenetic maps
    └── ChIP-Atlas: comprehensive ChIP-seq database
```

---

# ✅ TIER 3 COMPLETION CHECKLIST

```
═══════════════════════════════════════════════════════════════════════
 SECTION 3.1 — 12 HALLMARKS OF AGING
═══════════════════════════════════════════════════════════════════════
 □ Can name all 12 hallmarks and classify as primary/antagonistic/integrative
 □ Can explain: genomic instability → somatic mutations accumulate with age
 □ Can explain telomere biology: TTAGGG, shelterin, telomerase (TERT + TERC)
 □ Know: epigenetic alterations → methylation drift → aging clocks
 □ Can explain proteostasis: chaperones + proteasome + autophagy all decline
 □ Know the 4 nutrient sensing pathways and how they're deregulated (mTOR↑ AMPK↓ Sirtuins↓ IIS↑)
 □ Can explain mitochondrial dysfunction: ROS↑, ATP↓, mtDNA mutations, NAD⁺↓
 □ Know cellular senescence markers: p16, SA-β-gal, SASP, γH2AX
 □ Can explain stem cell exhaustion in at least 2 tissues (blood, muscle)
 □ Know: SASP → paracrine senescence → altered intercellular communication
 □ Can explain cGAS-STING pathway and its role in inflammaging
 □ Know: all hallmarks are INTERCONNECTED (can draw at least 5 connections)
 □ Have READ the López-Otín 2023 paper (even partially)

═══════════════════════════════════════════════════════════════════════
 SECTION 3.2 — IMMUNOLOGY
═══════════════════════════════════════════════════════════════════════
 □ Can explain innate vs adaptive immunity with cell types
 □ Know: CD4+ (helper), CD8+ (cytotoxic), Tregs, NK cells and their roles
 □ Can explain MHC-I (all cells → CD8+) vs MHC-II (APCs → CD4+)
 □ Know: thymic involution → fewer naive T cells → immunosenescence
 □ Can explain the concept of inflammaging (immunosenescence + chronic inflammation)
 □ Know: cGAS-STING and NLRP3 inflammasome roles in aging
 □ Can explain CAR-T senolytics (Amor 2024): T cells targeting uPAR on senescent cells
 □ Know: rapamycin improves vaccine response in elderly

═══════════════════════════════════════════════════════════════════════
 SECTION 3.3 — GENOMICS TECHNOLOGIES
═══════════════════════════════════════════════════════════════════════
 □ Can compare: Sanger vs Illumina vs PacBio vs Nanopore (read length, accuracy, cost)
 □ Know Illumina mechanism: bridge amplification, SBS, reversible terminators
 □ Know: PacBio HiFi = long + accurate, Nanopore = longest reads + portable
 □ Can list functional assays: RNA-seq, scRNA-seq, ChIP-seq, ATAC-seq, Hi-C, WGBS
 □ Know the basic RNA-seq pipeline: FASTQ → QC → Trim → Align → Count → DE
 □ Can explain de novo assembly: reads → k-mers → de Bruijn graph → contigs → scaffolds
 □ Know assembly quality metrics: N50, BUSCO
 □ Know: spatial transcriptomics adds tissue location to gene expression data

═══════════════════════════════════════════════════════════════════════
 SECTION 3.4 — BIOSTATISTICS
═══════════════════════════════════════════════════════════════════════
 □ Can explain the multiple testing problem and why FDR correction is needed
 □ Know Benjamini-Hochberg procedure and can apply it conceptually
 □ Can interpret a Kaplan-Meier survival curve and explain log-rank test
 □ Know: hazard ratio < 1 = treatment reduces death risk (good for aging intervention)
 □ Can explain batch effects and why PCA is used to detect them
 □ Know ComBat and Harmony as batch correction methods
 □ Can read a volcano plot: log2FC (x) vs -log10(padj) (y)
 □ Know: always report padj + effect size + confidence interval, not just p-value

═══════════════════════════════════════════════════════════════════════
 SECTION 3.5 — COMPUTATIONAL PREREQUISITES
═══════════════════════════════════════════════════════════════════════
 □ Can navigate Linux CLI: cd, ls, grep, awk, head, tail, pipes (|), redirect (>)
 □ Have WSL2 installed on Windows with Ubuntu
 □ Know: conda for managing bioinformatics packages
 □ Can write basic R code: vectors, data frames, ggplot2 basics
 □ Know what Bioconductor is and can install packages from it
 □ Can name at least 3 databases for: sequences (GenBank, UniProt, Ensembl),
   variants (ClinVar, gnomAD), pathways (KEGG, GO), aging (GenAge, HAGR)
 □ Know: DESeq2 (R) for differential expression, Scanpy (Python) for single-cell

═══════════════════════════════════════════════════════════════════════
 TOTAL: 50 items. All checked = TIER 3 COMPLETE → YOU ARE RESEARCH READY!
═══════════════════════════════════════════════════════════════════════
```
